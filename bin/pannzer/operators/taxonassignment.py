from myoperator import BlockOperator
import operator

# Clone a blockwise operator. Class name and file name (operators/classname.py) must be identical.
class taxonassignment(BlockOperator) :
    """Outputs most frequent taxon per block."""
    def __init__(self,glob, maxrank=1):
        # parameters used by self.process()
        self.glob=glob
        self.maxrank=maxrank
        # define private object handle to spreadsheets (glob maps handles by name)
        [self.data,self.summary_data]=glob.use_sheets(["data","summary"])
        # define nicknames for column indices
        [self.qpid_col1,self.isquery_col,self.taxon_col1]=self.data.use_columns(["qpid","isquery","taxon"])
        self.summary_data.use_columns(["qpid","taxon","count","total","frequency","rank"])
        # this is a composite operator
        [self.a,self.b]=glob.use_operators(['lineage','taxon'])

    def process(self, block):
        # call preprocessing operators
        for row in block:
                self.a.process(row) # add lineage
                self.b.process(row) # add taxon
        # get query
        qpid='unknown'
        for row in block:
                if row[self.isquery_col]=="1": qpid=row[self.qpid_col1] # string comparison
        # compute sum statistics
        n=0
        cnt={}
        for row in block:
                n+=1
                taxon=row[self.taxon_col1]
                if not taxon in cnt: cnt[taxon]=0
                cnt[taxon]+=1
        # reverse sort by counts
        sorted_cnt = sorted(cnt.items(), key=operator.itemgetter(1), reverse=True)
        # save top ranks in summary table
        rank=0
        for (taxon,cnt) in sorted_cnt:
                rank+=1
                freq=0
                if n>0: freq=cnt/float(n)
                # assume fixed column order
                datarow=[qpid,taxon,str(cnt),str(n),str(freq),str(rank)]
                self.summary_data.append_row(datarow)
                if rank>=self.maxrank: break

