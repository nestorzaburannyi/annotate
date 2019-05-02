from myoperator import RowOperator
import sys

# Clone a rowwise operator. Class name and file name (operators/classname.py) must be identical.
class taxon(RowOperator) :
    """Extract taxon from lineage to taxon_depth."""
    def __init__(self,glob):
        # define nicknames for column indices
        [self.lineage_col, self.taxon_col]=glob.use_sheet("data").use_columns(["lineage","taxon"])
        self.taxon_depth=int(glob.param['TAXON_DEPTH'])

    def process(self, row):
        x=row[self.lineage_col]
        #print >> sys.stderr, '#taxon:',x,self.lineage_col,row
        tmp=x.split(";")
        if len(tmp) > self.taxon_depth:
                row[self.taxon_col]=";".join(tmp[0:self.taxon_depth])
        else:
                row[self.taxon_col]=";".join(tmp)
