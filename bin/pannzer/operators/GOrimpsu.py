from myoperator import RowOperator
import sys

class GOrimpsu(RowOperator):
        """
        Add GO class and count vectors to data.

        Creates data columns 'GOclass','GOclass_count'
        Inputs: data columns 'spid','status'
        """
        def __init__(self, glob):
                sys.stderr.write("# Init GOrimpsu\n")
                self.glob=glob
                [self.go_col,self.gocount_col,self.spid_col,self.status_col]=glob.use_sheet("data").use_columns(['GOclass','GOclass_count','spid','status'])
                glob.use_online_dictionaries(["GOIDELIC","GODICT"])

        def process(self,row,verbose=False):
                # no GO terms imported for hits removed by Filter
                go=""
                cnt=""
                if row[self.status_col]!="False":
                  spid=row[self.spid_col]
                  try:
                        (db,accession,pid)=spid.split("|") # tr|E0SNF6|E0SNF6_DICD3
                        if accession in self.glob.GOdict:
                                go=self.glob.GOdict[accession].replace(","," ").split() # slow function, therefore here and not during initialization
                                for goid in go:
                                        cnt+=self.glob.GOcounts[goid]+" " # lookup count of goid
                  except:
                          if verbose: sys.stderr.write("# Warning: no GOcounts for spid = %s\n" %(spid))
                row[self.go_col]=" ".join(go).rstrip() # data is string
                row[self.gocount_col]=cnt.rstrip() # data is string

