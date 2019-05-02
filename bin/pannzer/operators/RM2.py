from myoperator import BlockOperator
import sys
from PannzerFunctions import logmod

class RM2(BlockOperator):
        """
        Calculate regression model 2 score.

        RM2 = 0.59 + 0.53 * logmod(avg(valencia_wordscore)) + 0.02 * avg(jaccard_wordscore) + 0.0004 * GSZ

        Reorders clusters in cluster_data by RM2.

        Creates cluster_data column 'RM2'
        Inputs: data columns 'clusid','valencia_wordscore','jaccard_wordscore'
                cluster_data columns 'clusid','cluster_GSZ'
        """
        def __init__(self, glob):
                sys.stderr.write("# Init RM2\n")
                self.glob = glob
                [self.data,self.cluster_data]=glob.use_sheets(["data","cluster"])
                [self.clusid_col1,self.val_col,self.jac_col]=self.data.use_columns(['clusid','valencia_wordscore','jaccard_wordscore'])
                [self.rm2_col,self.valavg_col,self.jacavg_col,self.clusid_col2,self.gsz_col]=self.cluster_data.use_columns(['RM2','val_avg','jac_avg','clusid','cluster_GSZ'])

        def process(self,block):
                if self.data.nrows==0: return
                # calculate average wordscores per cluster
                val_summa={}
                jac_summa={}
                n={}
                for row in block:
                        clusid=row[self.clusid_col1]
                        if not clusid in val_summa:
                                n[clusid]=0
                                val_summa[clusid]=0.0
                                jac_summa[clusid]=0.0
                        n[clusid]+=1
                        val_summa[clusid]+=float(row[self.val_col])
                        jac_summa[clusid]+=float(row[self.jac_col])
                val_avg={}
                jac_avg={}
                for clusid in n.keys():
                        val_avg[clusid]=val_summa[clusid]/n[clusid]
                        jac_avg[clusid]=jac_summa[clusid]/n[clusid]
                # save RM2 per cluster
                for row in self.cluster_data.block:
                        clusid=row[self.clusid_col2]
                        if not clusid in val_avg:
                                x=0.0
                                y=0.0
                        else:
                                x=logmod(val_avg[clusid])
                                y=jac_avg[clusid]
                        gsz=float(row[self.gsz_col])
                        rm2=0.59 + 0.53 * x + 0.02 * y + 0.0004 * gsz
                        #print "#rm2",x,y,gsz,0.59*x,0.02*y,0.0004*gsz,rm2
                        row[self.rm2_col]=str(rm2) # data is string
                        row[self.valavg_col]=str(x)
                        row[self.jacavg_col]=str(y)
                # reverse sort clusters by RM2
                self.cluster_data.sort_block(self.rm2_col,reverse=True)

