from myoperator import RowOperator
import sys,math

class RM1(RowOperator):
        """
        Calculate Regression Model 1 score.

        RM1 = 0.2 + log_10(pide) + 0.16 * qcov * scov + 0.40 * scov * taxdist

        RM1 is set to zero if status is "False".

        Creates data column 'RM1'.
        Inputs: data columns 'qcov', 'scov', 'pide','taxdist','status'
        """
        def __init__(self, glob):
                sys.stderr.write("# Init RM1\n")
                self.glob = glob
                [self.pide_col,self.qcov_col,self.scov_col,self.taxdist_col,self.status_col,self.rm1_col]=glob.use_sheet("data").use_columns(['qcov', 'scov', 'pide','taxdist','DE_status','RM1'])

        def process(self,row):
                if row[self.status_col]=="False":
                        row[self.rm1_col]="0.0"
                        return
                percent_identity=float(row[self.pide_col])*100.0
                qcov=float(row[self.qcov_col])
                scov=float(row[self.scov_col])
                taxdist=float(row[self.taxdist_col])
                x=0.21*math.log(percent_identity)/math.log(10.0)+0.26*qcov*scov+0.40*scov*taxdist
                row[self.rm1_col]=str(x) # data is string

