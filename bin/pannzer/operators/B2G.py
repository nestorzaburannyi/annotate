from myoperator import RowOperator
import sys

class B2G(RowOperator):
        """
        Test adding evidence weights for Blast2GO function.

        Creates data columns 'GOclass','GOclass_count','evidence_code_weight'
        Inputs: data columns 'spid'
        """
        def __init__(self, glob):
                sys.stderr.write("# Init GOrimpsu with evidencecode\n")
                self.glob=glob
                [self.go_col,self.gocount_col,self.spid_col,self.status_col,self.ecw_col]=glob.use_sheet("data").use_columns(['GOclass','GOclass_count','spid','status','evidence_code_weight'])
                glob.use_online_dictionaries(["GOIDELIC","GODICT","ECWEIGHT"])

        def process(self,row,verbose=True):
                # no GO terms imported for hits removed by Filter
                go=""
                cnt=""
                ecw=""
                if row[self.status_col]!="False":
                  spid=row[self.spid_col]
                  try:
                        (db,accession,pid)=spid.split("|") # tr|E0SNF6|E0SNF6_DICD3
                        if accession in self.glob.GOdict:
                                go=self.glob.GOdict[accession].replace(","," ").split() # slow function, therefore here and not during initialization
                                for goid in go:
                                        cnt+=self.glob.GOcounts[goid]+" " # lookup count of goid
                                ecw=self.glob.GOdict_weights[accession].replace(","," ").split()
                  except:
                          if verbose: sys.stderr.write("# Warning: no GOcounts for spid = %s\n" %(spid))
                row[self.go_col]=" ".join(go).rstrip() # data is string
                row[self.gocount_col]=cnt.rstrip() # data is string
                row[self.ecw_col]=" ".join(ecw).rstrip()

