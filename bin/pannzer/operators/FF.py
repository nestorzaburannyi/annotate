from myoperator import RowOperator
from PannzerFunctions import FormFactor
import sys

class FF(RowOperator):
        """
        Adds FormFactor column to block.

        Creates data column 'FF'
        Input: data column 'desc'
        """
        def __init__(self,glob):
                sys.stderr.write("# Init FF\n")
                [self.ff_col,self.desc_col]=glob.use_sheet("data").use_columns(["FF","desc"])
                self.f=FormFactor()

        def process(self,row):
                desc=row[self.desc_col].strip()
                ff=str(self.f.formfactor(desc))
                row[self.ff_col]=ff
                row[self.desc_col]=desc

