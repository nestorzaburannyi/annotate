# TextOperator can be used to parse text files and output tab

from myoperator import BlockOperator
import sys,re

class taxinspector(BlockOperator):
    """
Input data stream is one proteome. Count histogram of identities of best hit from each query protein to best hit in any target species.

Input: qpid, pide, species
Output is written directly to STDOUT: totalcount species genus kingdom count(pide, pide=0.00 to 1.00 in steps of 0.01)

genus and kingdom are used to group data in plotly scatterplot
    """
    def __init__(self,glob):
        # set parameter
        self.glob=glob
        [self.data,self.summary_data]=glob.use_sheets(['data','summary'])
        [self.qpid_col,self.isquery_col,self.pide_col,self.species_col,self.lineage_col,self.genus_col,self.kingdom_col]=\
                self.data.use_columns(['qpid','isquery','pide','species','lineage','taxon','kingdom'])
        self.summary_data.use_columns(['wsum','count','kingdom','genus','species','matched_fraction','average_identity','median','pide_bins','lineage'])
        # do sums over input stream
        self.counts={} # key = species, counts iindex = pide/100
        self.nprot=0
        self.kingdom={} # key = species
        self.genus={} # key = species
        self.lineage={}
        self.nprot=0
        # add lineage, genus, kingdom
        [self.a,self.b]=glob.use_operators(['lineage','genus'])

    def process(self, block):
        # only count best hit per species
        seen={}
        self.data.sort_block(self.pide_col,reverse=True)
        if len(block)>0: self.nprot+=1
        for row in block:
                if row[self.isquery_col]=="1": continue # exclude query
                species=row[self.species_col]
                if species in seen: continue
                seen[species]=1
                # init counts
                if not species in self.counts:
                        self.counts[species]=[0]*101
                        self.a.process(row) # add lineage
                        self.b.process(row) # add genus, kingdom
                        self.kingdom[species]=row[self.kingdom_col]
                        self.genus[species]=row[self.genus_col]
                        self.lineage[species]=row[self.lineage_col]
                pide=row[self.pide_col]
                x=int(float(pide)*100)
                self.counts[species][x]+=1

    def finalise(self):
        for species in self.counts:
                total=0
                w=0
                x=self.counts[species]
                cumul=[0]*101
                for i in range(0,101):
                        n=x[i]
                        total+=n
                        w+=i*n
                        if i>0: cumul[i]=cumul[i-1]+n
                target=total/2
                median=0
                while cumul[median]<target:
                        if median==100: break
                        median+=1
                #genus=species.split()[0]
                genus=self.genus[species]
                kingdom=self.kingdom[species]
                lineage=self.lineage[species]
                datarow=[str(float(w)/100),str(total),kingdom,genus,species,str(float(total)/self.nprot),str(float(w)/100/total),str(float(median)/100),",".join(str(a) for a in x),lineage]
                self.summary_data.append_row(datarow)
        self.summary_data.sort_block(0,reverse=True)
        self.summary_data.output(result=True)

