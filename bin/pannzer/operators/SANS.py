import sys
from myoperator import RowOperator

# Clone a rowwise operator. Class name and file name (operators/classname.py) must be identical.
class SANS(RowOperator) :
    """This operator don't do nothing."""
    def __init__(self,glob):
        pass

    def process(self,row):
        pass
