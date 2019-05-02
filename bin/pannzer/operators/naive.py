import sys
from myoperator import BlockOperator

# Clone a rowwise operator. Class name and file name (operators/classname.py) must be identical.
class naive(BlockOperator) :
    """Naive predictor as in CAFA: p(goid) is its base frequency in GOA. Input: qpid list"""
    def __init__(self,glob):
        self.glob=glob
        # define nicknames for column indices
        [self.data,self.prediction_data]=self.glob.use_sheets(["data","prediction"])
        [self.qpid_col]=self.data.use_columns(["qpid"])
        self.prediction_data.use_columns(["qpid","goid","frequency","ontology","description"])
        # use online dictionary. Object handles in glob are hardcoded
        self.glob.use_online_dictionaries(['GOIDELIC']) # downloads copy in Runner
        self.pred=None

    def process(self,block):
        try:
                row=block[0]
        except:
                return # empty block
        if not self.pred: self.pred=self.naive() # load base frequencies
        qpid=row[self.qpid_col] # input protein identifier
        # multiply predictions for each qpid
        for x in self.pred:
                datarow=[qpid] + x
                self.prediction_data.append_row(datarow)

    def naive(self):
        "base frequencies of GO terms, truncated to 2 decimals"
        pred=[] # tuples (goid,freq)
        for goid in self.glob.GOcounts.keys():
                ontology=self.glob.ontology[goid]
                rootcount=self.glob.rootcount[ontology]
                if rootcount<1: continue
                count=self.glob.GOcounts[goid]
                freq=float(count)/rootcount
                if freq >= 0.01: pred.append( [goid, "%.2f" %freq, ontology, self.glob.godesc[goid] ] )
        return(pred)
