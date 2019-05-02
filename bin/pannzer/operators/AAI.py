import sys
from myoperator import BlockOperator
import DictServer

# Clone a rowwise operator. Class name and file name (operators/classname.py) must be identical.
class AAI(BlockOperator) :
    """
Output_1 = SANS result plus lineage, genus, kingdom
Output_2 = taxonomic summary of SANS hits per species: weighted AUC, count, genus, kingdom, coverage, average, median, AAI histogram

Filter Output_2 on coverage range, average range, max number of hits to generate Plotly barcharts.
Filter Output_2 on minimum AAI to generate Plotly scatterplot.
Exclude query-species (top .tax species) from Output_1 and select best hit per query to generate Krona taxonomic profile piechart.

"""
    def __init__(self,glob):
        self.glob=glob
        [self.data,self.summary_data]=glob.use_sheets(['data','summary'])
        [self.qpid_col,self.spid_col,self.isquery_col,self.pide_col,self.bits_col,self.lineage_col,self.species_col,self.genus_col,self.kingdom_col]=\
                self.data.use_columns(['qpid','spid','isquery','pide','bits','lineage','species','taxon','kingdom'])
        self.data.use_columns(['nid','qcov','scov','lali','desc','qseq','vote','genename'])
        try:
                self.kt_format=glob.kt_format
        except:
                self.kt_format=False
        if self.kt_format:
                self.summary_data.create_columns(['#wsum','lineage','species'])
        else:
                self.summary_data.create_columns(['wsum','count','target_count','multiplicity','kingdom','genus','species','matched_fraction','average_identity','median','pide_bins','lineage'])
        # do sums over input stream
        self.counts={} # key = species, counts iindex = pide/100
        self.total={} # key = species, number of hits
        self.wsum={} # key = species, pide-sum
        self.nprot=0
        # global parameters used: glob.param['PANZ_MAXHITS'], glob.param[input_QUERYSPECIES'], glob.kt_format (no default)
        try:
                self.kt_format=glob.kt_format
        except:
                self.kt_format=False
        self.MAXHITS=self.glob.param['PANZ_MAXHITS']
        self.topspecies=self.glob.param['input_QUERYSPECIES']
        # count many-to-one matches
        self.pidcount={} # key = spid, count queries
        self.pid2species={} # key = spid, map to species
        # remote access
        self.REMOTE=self.glob.param['CONN_REMOTE']
        self.DICTHOST=self.glob.param['CONN_HOSTNAME']

    def process(self,block):
        seen={}
        # to exclude query species if defined
        try:
                seen[self.topspecies]=1
        except:
                pass
        self.data.sort_block(self.bits_col,reverse=True)
        if len(block)>0: self.nprot+=1
        nhits=0
        for row in block:
                if row[self.isquery_col]=="1": continue # exclude query
                # only count best hit per species
                species=row[self.species_col].rstrip()
                #print >> sys.stderr, nhits,(species in seen),species
                if species in seen: continue
                seen[species]=1
                nhits+=1
                if nhits>self.MAXHITS: break
                # init counts of new pid
                pid=row[self.spid_col]
                if not pid in self.pid2species:
                        self.pid2species[pid]=species
                        self.pidcount[pid]=0
                self.pidcount[pid]+=1
                # init counts of new species
                if not species in self.counts:
                        self.counts[species]=[0]*101
                        self.total[species]=0.0
                        self.wsum[species]=0.0
                pide=0.005+float(row[self.pide_col])
                x=int(pide*100)
                self.counts[species][x]+=1
                self.total[species]+=1.0
                self.wsum[species]+=pide

    def finalise(self):
        # get lineage for species in data
        msg=''
        for species in self.counts: msg+="LINEAGE"+"\t"+species.upper()+"\n"
        tmp=DictServer.DICTquery(msg,sys.version_info[0],HOSTNAME=self.DICTHOST,REMOTE=self.REMOTE).split("\n")
        # load values to dictionaries (keys occur in chunk)
        lineage={}
        for row in tmp:
                row.encode('utf-8')
                if not row: continue
                (table,key,value)=row.split("\t")
                lineage[key]=value
        # short form output if kt_format
        if self.kt_format:
                for species in self.counts:
                        ucspecies=species.upper()
                        if not ucspecies in lineage: continue
                        x=str(self.wsum[species])
                        y="\t".join(lineage[ucspecies].split('; '))
                        datarow=[x,y,species]
                        self.summary_data.append_row(datarow)
                self.summary_data.sort_block(0,reverse=True)
                self.summary_data.output(result=True)
                return
        # how many target proteins were matched?
        speciescount={}
        for pid in self.pidcount:
                species=self.pid2species[pid]
                if not species in speciescount: speciescount[species]=0
                speciescount[species]+=1
        # long form output default
        for species in self.counts:
                w=0
                x=self.counts[species]
                cumul=[0]*101
                for i in range(0,101):
                        n=x[i]
                        w+=i*n
                        if i>0: cumul[i]=cumul[i-1]+n
                target=self.total[species]/2
                median=0
                while cumul[median]<target:
                        if median==100: break
                        median+=1
                ucspecies=species.upper()
                if ucspecies in lineage:
                        line=lineage[ucspecies]
                        genus,kingdom=self.genus_kingdom(lineage[ucspecies])
                else:
                        line=''
                        genus=''
                        kingdom=''
                tot=self.total[species]
                ntarget=speciescount[species]
                multiplicity=tot/ntarget
                datarow=[str(float(w)/100),str(tot),str(ntarget),str(multiplicity),kingdom,genus,species,str(tot/self.nprot),str(float(w)/100/tot),str(float(median)/100),",".join(str(a) for a in x),line]
                self.summary_data.append_row(datarow)
        self.summary_data.sort_block(0,reverse=True)
        self.summary_data.output(result=True)

    def genus_kingdom(self,lineage):
        tmp=lineage.split("; ")
        taxon_depth=3
        lentmp=len(tmp)
        if lentmp<1: return('','')
        kingdom=tmp[0]
        if kingdom == "Eukaryota":
          try:
                if tmp[1] == "Viridiplantae":
                        taxon_depth=6
                        if tmp[taxon_depth] == "Magnoliophyta": taxon_depth=10
                elif tmp[6] == "Mammalia":
                        taxon_depth=9
                elif tmp[12] == "Aves":
                        taxon_depth=14
          except:
                pass
        else:
                taxon_depth=lentmp-1
        if lentmp > taxon_depth:
                genus=tmp[taxon_depth].replace("Candidatus ","")
        else:
                genus=tmp[-1]
        return(genus,kingdom)

