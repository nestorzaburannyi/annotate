from __future__ import print_function
from myoperator import BlockOperator
import sys

class Filter(BlockOperator):
        """
        set status flag = true if row passes filter criteria: qcov, scov, pide.

        set status = false after glob.param.maxhits informative hits.

        creates data column 'status'
        input: data columns 'qcov','scov','pide','desc'
        """
        def __init__(self, glob,verbose=False):
                sys.stderr.write("# init filter\n")
                [self.status_col,self.qcov_col,self.scov_col,self.pide_col,self.lali_col]=glob.use_sheet("data").use_columns(['status','qcov','scov','pide','lali'])
                # define nicknames for parameters
                self.MAXHITS=glob.param['PANZ_MAXHITS']
                self.FILTER_PERMISSIVE=glob.param['PANZ_FILTER_PERMISSIVE']
                self.QCOVCUTOFF=glob.param['PANZ_QCOVCUTOFF']
                self.SCOVCUTOFF=glob.param['PANZ_SCOVCUTOFF']
                self.MINPIDECUTOFF=glob.param['PANZ_MINPIDECUTOFF']
                self.MAXPIDECUTOFF=glob.param['PANZ_MAXPIDECUTOFF']
                self.MINLALI=glob.param['PANZ_MINLALI']
                if verbose: print("Filter parameters",self.MAXHITS,self.FILTER_PERMISSIVE,self.QCOVCUTOFF,self.SCOVCUTOFF,self.MINPIDECUTOFF,self.MAXPIDECUTOFF,self.MINLALI, file=sys.stderr)

        def process(self,block):
                # accept maxhits informative
                naccepted=0
                for row in block:
                        if naccepted >= self.MAXHITS: row[self.status_col]="False"
                        if row[self.status_col]=="False": continue # string comparison
                        status=self.process_row(row)
                        if status: naccepted+=1

        def process_row(self,row):
                # alignment quality
                qcov=float(row[self.qcov_col])
                scov=float(row[self.scov_col])
                pide=float(row[self.pide_col])
                lali=float(row[self.lali_col])
                # strict // permissive filtering
                if not self.FILTER_PERMISSIVE:
                       status =  qcov >= self.QCOVCUTOFF and \
                                scov >= self.SCOVCUTOFF and \
                                pide >= self.MINPIDECUTOFF and \
                                pide <= self.MAXPIDECUTOFF
                else:
                       status =  (qcov >= self.QCOVCUTOFF or \
                                scov >= self.SCOVCUTOFF or \
                                lali > self.MINLALI) and \
                                pide >= self.MINPIDECUTOFF \
                                and pide <= self.MAXPIDECUTOFF
                row[self.status_col]=str(status)
                return(status)

