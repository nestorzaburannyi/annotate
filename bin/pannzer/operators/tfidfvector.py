import sys
from myoperator import RowOperator

# Clone a rowwise operator. Class name and file name (operators/classname.py) must be identical.
class tfidfvector(RowOperator) :
    """
Input: qpid, desc
Output: word, wordcount, termidf
"""

    def __init__(self,glob):
        self.glob=glob
        # use online dictionary. Object handles in glob are hardcoded
        [self.a,self.b]=self.glob.use_operators(['Cleandesc','TFIDF'])

    def process(self,row):
        self.a.process(row) # add cleandesc
        self.b.process(row) # add trdif

