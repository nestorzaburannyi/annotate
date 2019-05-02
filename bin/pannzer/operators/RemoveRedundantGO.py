from myoperator import BlockOperator
import sys

class RemoveRedundantGO(BlockOperator):
        """
        Set status of redundant GO classes to False. Redundant classes are a parent or child of higher ranking class.

        Sets row_status to False if redundant. Redundant GO classes will be hidden from output.

        Inputs: goclass_data column 'goid'. Block is assumed sorted.
        Outputs: row_status attribute of goclass_data spreadsheet
        """
        def __init__(self, glob):
                sys.stderr.write("# Init RedundantGO\n")
                self.glob=glob
                self.goclass_data=glob.use_sheet("goclass")
                glob.use_online_dictionaries(["GOIDELIC"])
                [self.goid_col]=self.goclass_data.use_columns(["goid"])

        def process(self,block):
                accepted={}
                rejected={}
                for i in range(0,self.goclass_data.nrows):
                        self.goclass_data.row_status[i]=False # overwrite accepted
                for i in range(0,self.goclass_data.nrows):
                        goid=self.goclass_data.block[i][self.goid_col]
                        # already rejected as parent of previously accepted
                        if goid in rejected:
                                continue
                        if not goid in self.glob.GOparents: continue # must be root
                        # reject if child of previously accepted
                        for p in self.glob.GOparents[goid]:
                                if p in accepted:
                                    rejected[goid]=True
                                    break
                        if goid in rejected:
                                continue
                        # accept goid, reject parents
                        accepted[goid]=True
                        self.goclass_data.row_status[i]=True
                        for p in self.glob.GOparents[goid]:
                                rejected[p]=True

