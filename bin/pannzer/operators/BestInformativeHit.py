from myoperator import BlockOperator

class BestInformativeHit(BlockOperator):
        """
        Uses best-hit method to predict DE and GO.

        Select first hit with noninformative words from filtered hitlist.
        Select first hit with GO annotations from filtered hitlist.

        Creates summary table columns 'qpid', 'bits', 'desc' (copied from selected row of data)
        Creates goclass table columns 'qpid','bits','ontology','goid','desc' (copied from selected row of data)

        Inputs: data columns 'status', 'bits', 'qpid', 'desc', 'spid'
        """
        def __init__(self, glob):
                self.glob=glob
                [self.data,self.summary_data,self.goclass_data]=glob.use_sheets(["data","summary","goclass"])
                [self.qpid_col1,self.bits_col1,self.desc_col1,self.status_col1,self.spid_col1,self.ff_col1,self.goclass_col]=self.data.use_columns(['qpid','bits','desc','status','spid','FF','GOclass'])
                self.summary_data.use_columns(['qpid','bits','desc'])
                self.goclass_data.use_columns(['qpid','bits','ontology','goid','desc'])
                glob.use_online_dictionaries(['GOIDELIC'])
                [self.b,self.c,self.d]=glob.use_operators(['Filter','FF','GOrimpsu'])
                # parameter nicknames
                self.FFCUTOFF=glob.param['PANZ_FFCUTOFF']

        def process(self,block):
                # get query
                try:
                        qpid=block[0][self.qpid_col1]
                except:
                        return
                # preprocess
                for row in block:
                        self.c.process(row) # calculate form factor
                        self.d.process(row) # associate GO terms
                self.b.process(block) # set status(qcov,scov,pide|MAXHITS)
                bestbits="0.0"
                bestdesc='No informative hits'
                # output only one row per block
                for row in block:
                        if row[self.status_col1]=="False": continue
                        if row[self.spid_col1]==qpid : continue
                        if float(row[self.ff_col1])<self.FFCUTOFF: continue
                        # first inforative hit
                        bestbits=row[self.bits_col1]
                        bestdesc=row[self.desc_col1]
                        break
                datarow=[qpid,bestbits,bestdesc]
                self.summary_data.append_row(datarow)
                # lift first GO annotation
                bestgo=""
                bestgobits="0.0"
                bestdesc=""
                for row in block:
                        if row[self.status_col1]=="False": continue
                        if row[self.spid_col1]==qpid : continue
                        if row[self.goclass_col]!="":
                                bestgo=row[self.goclass_col]
                                bestgobits=row[self.bits_col1]
                                break
                # one goid per row
                for goid in bestgo.split(' '):
                        if not goid in self.glob.godesc: continue
                        datarow=[qpid,bestgobits,self.glob.ontology[goid],goid,self.glob.godesc[goid]]
                        self.goclass_data.append_row(datarow)
