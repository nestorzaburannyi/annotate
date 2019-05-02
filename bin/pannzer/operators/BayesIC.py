# cut -f 2,5,7 goa_uniprot_noIEA.gaf

from myoperator import RowOperator
from Read_and_Print import read_dict_data
from PannzerFunctions import Propagate
import sys,re,math

class BayesIC(RowOperator):
    """
Input: count parentset
Data: obo.tab
Output: goid children parentset ontology desc goidcount parentsetcount ic
    """
    def __init__(self,glob):
        # set parameters
        self.glob=glob
        # propagate assigned GOids to parents
        obotabfile=glob.param['eval_OBOTAB']
        self.mapped_to=read_dict_data(obotabfile,'goid','mapped_to',header=True,UseColNames=True) # goid, mapped_to
        self.children=read_dict_data(obotabfile,'goid','children',header=True,UseColNames=True)
        self.directparents=read_dict_data(obotabfile,'goid','directParentSet',header=True,UseColNames=True) # goid,parentset
        self.coparents=read_dict_data(obotabfile,'goid','CoParentSets',header=True,UseColNames=True) # goid, coparentset
        self.ontology=read_dict_data(obotabfile,'goid','ontology',header=True,UseColNames=True)
        self.desc=read_dict_data(obotabfile,'goid','desc',header=True,UseColNames=True)
        # initialise counts of v, P(v)
        self.count={}
        pseudocount=1
        for goid in self.directparents.keys():
                self.count[self.mapped_to[goid]]=pseudocount
                self.count[self.directparents[goid]]=pseudocount
        # write blockwise summmary
        [self.data,self.ic_data]=glob.use_sheets(["data","ic"])
        # column indices in GAF file
        [self.acc_col1,self.goid_col1]=self.data.use_columns(["qpid","propagated"])
        self.ic_data.use_columns(["goid","IC","ontology","children","parentset","goidcount","parentsetcount","desc"])

    def process(self,row):
        cnt=int(row[0])
        propagated=row[1].split(',')
        # accumulate counts(v), counts(P(v)) for IC calculation
        done={}
        for g in propagated:
                if not g in self.mapped_to: continue # obsolete term
                goid=self.mapped_to[g] # add alt_id counts to representative
                if goid in done: continue
                self.count[goid]+=cnt # go-term
                done[goid]=1
                for x in self.coparents[goid].split(';'): # stringified coparent sets
                        if x in done: continue
                        ok=True
                        for y in x.split(','):
                                if not y in propagated:
                                        ok=False
                                        break
                        if ok: self.count[x]+=cnt
                        done[x]=1 # only count coparentset once!

    def finalise(self):
        # IC(v) = -log2 ( count(v)/count(Pv)
        for goid in self.directparents.keys():
                countrep=self.mapped_to[goid]
                goidcount=self.count[countrep]
                parentset=self.directparents[goid]
                parentsetcount=self.count[parentset]
                ic=-math.log(float(goidcount)/parentsetcount,2)
                if ic<0: ic=0.0 # root nodes
                self.ic_data.append_row([goid,str(ic),self.ontology[goid],self.children[goid],parentset,str(goidcount),str(parentsetcount),self.desc[goid]])
        # output ic table
        self.ic_data.output(result=True)

