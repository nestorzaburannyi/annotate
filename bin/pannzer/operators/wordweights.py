from myoperator import RowOperator
from PannzerFunctions import FormFactor
import sys

class wordweights(RowOperator):
        """
        Adds FormFactor, total TFIDF, word and TFIDF vectors.

        Creates data column 'FF'
        Input: data column 'desc'
        """
        def __init__(self,glob):
                sys.stderr.write("# Init wordweights\n")
                [self.ff_col,self.cleandesc_col,self.qpid_col]=glob.use_sheet("data").use_columns(["FF","cleandesc","qpid"])
                self.f=FormFactor()
                [self.a]=glob.use_operators(['tfidfvector'])

        def process(self,row):
                self.a.process(row) # desc -> cleandesc, wordvector, tfidfvector
                cleandesc=row[self.cleandesc_col].strip()
                ff=str(self.f.formfactor(cleandesc))
                row[self.ff_col]=ff

