from myoperator import RowOperator
import sys

class DE_status(RowOperator):
        """Sets status = True if sbjct passes alignment and form-factor criteria.
        Creates data column 'DE_status'
        Input: data columns 'FF','status'
        """
        def __init__(self,glob):
                sys.stderr.write("# Init DE_status\n")
                self.blockwise=False
                [self.status_col,self.ff_col,self.DE_status_col]=glob.use_sheet("data").use_columns(['status','FF','DE_status'])
                self.FFCUTOFF=glob.param['PANZ_FFCUTOFF']

        def process(self,row):
                if row[self.status_col]=="True" and float(row[self.ff_col])>=self.FFCUTOFF:
                        row[self.DE_status_col]="True"
                else:
                        row[self.DE_status_col]="False"

