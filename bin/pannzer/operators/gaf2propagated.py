# cut -f 2,5,7 goa_uniprot_noIEA.gaf

from myoperator import BlockOperator
from Read_and_Print import read_dict_data
from PannzerFunctions import Propagate
import sys,re,math

class gaf2propagated(BlockOperator):
    """
Input: qpid, goid (multiline block)
Output: qpid, propagated_goids (one line per qpid)
    """
    def __init__(self,glob):
        # set parameters
        self.glob=glob
        # propagate assigned GOids to parents
        obotabfile=glob.param['eval_OBOTAB']
        self.mapped_to=read_dict_data(obotabfile,'goid','mapped_to',header=True,UseColNames=True) # goid, mapped_to
        self.directparents=read_dict_data(obotabfile,'goid','directParentSet',header=True,UseColNames=True) # goid,parentset
        # write blockwise summmary
        [self.data,self.godict_data]=glob.use_sheets(["data","godict"])
        # column indices in GAF file
        [self.acc_col1,self.goid_col1]=self.data.use_columns(["qpid","goid"])
        self.godict_data.use_columns(["qpid","propagated"])

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
