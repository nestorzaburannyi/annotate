# cut -f 2,5,7 goa_uniprot_noIEA.gaf

from myoperator import BlockOperator
from Read_and_Print import read_dict_data
from PannzerFunctions import Propagate
import sys,re,math

class gaf2tab(BlockOperator):
    """Extract subset of rows from gaf file, write (acc, list of goids)"""
    def __init__(self,glob):
        # set parameters
        self.glob=glob
        # propagate assigned GOids to parents
        obofile=glob.param['eval_OBOTAB']
        self.mapped_to=read_dict_data(obofile,'goid','mapped_to',header=True,UseColNames=True) # goid, mapped_to
        self.directparents=read_dict_data(obofile,'goid','directParentSet',header=True,UseColNames=True) # goid,parentset
        self.coparents=read_dict_data(obofile,'goid','CoParentSets',header=True,UseColNames=True) # goid, coparentset
        # initialise counts of v, P(v)
        self.count={}
        pseudocount=1
        for goid in self.directparents.keys():
                self.count[self.mapped_to[goid]]=pseudocount
                self.count[self.directparents[goid]]=pseudocount
        # write blockwise summmary
        [self.data,self.ic_data,self.godict_data]=glob.use_sheets(["data","godict","ic"])
        # column indices in GAF file
        [self.acc_col1,self.goid_col1]=self.data.use_columns(["acc","goid"])
        self.ic_data.use_columns(["goid","goidcount","parentsetcount","IC"])
        self.godict_data.use_columns(["acc","propagated"])

    def process(self,block):
        try:
                acc=block[0][self.acc_col1]
        except:
                return
        if acc[0]=='!': return # gaf comment line
        goidlist=[]
        directparents={}
        for row in block:
                goid=row[self.goid_col1]
                goid=re.sub(r'GO:','',goid)
                # direct assignment + parents
                directparents[goid]=1
        # propagate to all parents, write to godict sheet
        tmp=directparents.keys()
        propagated=Propagate(tmp,self.directparents)
        propagated.sort()
        self.godict_data.append_row([acc,','.join(propagated)])
        # accumulate counts(v), counts(P(v)) for IC calculation
        done={}
        for g in propagated:
                goid=self.mapped_to[g] # add alt_id counts to representative
                if goid in done: continue
                self.count[goid]+=1 # go-term
                done[goid]=1
                for x in self.coparents[goid].split(';'): # stringified coparent sets
                        if x in done: continue
                        ok=True
                        for y in x.split(','):
                                if not y in propagated:
                                        ok=False
                                        break
                        if ok: self.count[x]+=1
                        done[x]=1 # only count coparentset once!

    def finalise(self):
        # IC(v) = -log2 ( count(v)/count(Pv)
        for goid in self.directparents.keys():
                countrep=self.mapped_to[goid]
                goidcount=self.count[countrep]
                parentsetcount=self.count[self.directparents[countrep]]
                ic=-math.log(float(goidcount)/parentsetcount,2)
                if ic<0: ic=0.0 # root nodes
                self.ic_data.append_row([goid,str(goidcount),str(parentsetcount),str(ic)])
        # output ic table
        self.ic_data.output(result=True)

