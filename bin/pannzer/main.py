from __future__ import absolute_import
from __future__ import print_function
from . import myoperator
import globs

if __name__ == '__main__' :
    glob=globs.globs()
    operator_dir = 'operators'
    operators = ['taxon','lineage','taxonassignment']

    myoperator.initialise_operators(operator_dir)

    for operator_name in operators :
        try :
            Op = myoperator.get_operator(operator_name)
        except myoperator.OperatorError as oe :
            print("Error", oe)
            continue

        x = Op(glob)

        if isinstance(x, myoperator.RowOperator) :
            x.process("row")
        elif isinstance(x, myoperator.BlockOperator) :
            x.process("block")
	else:
	    x.process("unknown")
