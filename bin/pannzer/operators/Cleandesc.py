from myoperator import RowOperator
import sys
from PannzerFunctions import Cleaner

class Cleandesc(RowOperator):
        """
        Generate cleaned descriptions.

        Creates data column 'cleandesc'
        Inputs: data column 'desc'
        """
        def __init__(self, glob):
                sys.stderr.write("# Init Cleandesc\n")
                [self.desc_col,self.cleandesc_col]=glob.use_sheet("data").use_columns(["desc","cleandesc"])
                # nicknmae parameter
                self.remove_abbr=glob.param['PANZ_REMOVE_ABBR']

        def process(self,row):
                x=Cleaner(row[self.desc_col], remove_abbr=self.remove_abbr).upper()
                row[self.cleandesc_col]=x
