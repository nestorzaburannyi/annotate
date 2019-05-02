import sys
import os
from abc import ABCMeta, abstractmethod

class Operator(object) :
    __metaclass__ = ABCMeta

    @abstractmethod
    def process(self) : pass
    def finalise(self): pass

class RowOperator(Operator) : pass
class BlockOperator(Operator) : pass
class TextOperator(Operator): pass
class OperatorError(Exception) : pass

def initialise_operators(path):
    p = os.path.abspath(sys.argv[0]) + os.sep + path
    if not p in sys.path :
        sys.path.insert(0, path)

def get_operator(operator_name) :
    try :
        m = __import__(operator_name, globals(), locals(), operator_name)
    except ImportError as ie :
        raise OperatorError("could not import module %s (%s)" % (operator_name, str(ie)))

    try :
        c = getattr(m, operator_name)
    except AttributeError as ae :
        raise OperatorError("module loaded, but could not find operator called %s (%s)" % (operator_name, str(ae)))

    return c
