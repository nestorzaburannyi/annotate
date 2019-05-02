from myoperator import BlockOperator

class BestInformativeHit_DE(BlockOperator):
        """
        Select first hit without noninformative words from filtered hitlist.

        Creates cluster_data columns 'qpid', 'bits', 'desc' (copied from selected row of data)

        Inputs: data columns 'status', 'bits', 'qpid', 'desc', 'spid'
        """
        def __init__(self, glob):
                [self.data,self.summary_data]=glob.use_sheets(["data","summary"])
                [self.qpid_col1,self.bits_col1,self.desc_col1,self.status_col1,self.spid_col1,self.ff_col1]=self.data.use_columns(['qpid','bits','desc','status','spid','FF'])
                [self.qpid_col2,self.bits_col2,self.desc_col2]=self.summary_data.use_columns(['qpid','bits','desc'])
                [self.b,self.c]=glob.use_operators(['Filter','FF'])
                # parameter nicknames
                self.FFCUTOFF=glob.param['PANZ_FFCUTOFF']

        def process(self,block):
                # get query
                try:
                        qpid=block[0][self.qpid_col1]
                except:
                        return
                # preprocess
                for row in block: self.c.process(row) # calculate form factor
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
                self.summary_data.append_row(row)
                row[self.qpid_col2]=qpid
                row[self.bits_col2]=bestbits
                row[self.desc_col2]=bestdesc

