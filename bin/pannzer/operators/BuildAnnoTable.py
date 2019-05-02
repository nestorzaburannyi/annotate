from myoperator import BlockOperator
import sys
from PannzerFunctions import FormFactor,DE_PPV_euk,DE_PPV_bac,GO_PPV,GO_argot_PPV,GO_jac_PPV,GO_hyge_PPV,GO_slow_hyge_PPV

class BuildAnnoTable(BlockOperator):
        """
        Save one prediction per row

        Creates anno_data columns 'type','score', 'PPV', 'id', 'desc'
        type := original_DE | DE | BP | MF | CC | original_DE
        score is RM2, RM3 or bits
        description is DE or GO prediction

        Inputs: cluster_data and goclass_data
        """
        def __init__(self, glob):
                sys.stderr.write("# Init BuildAnnoTable\n")
                self.glob = glob
                self.f=FormFactor()
                [self.data,self.cluster_data,self.goclass_data,self.anno_data]=glob.use_sheets(["data","cluster","goclass","anno"])
                self.annocols=self.anno_data.use_columns(["qpid","type","score","PPV","id","desc"])
                [self.qpid_col]=self.data.use_columns(['qpid'])
                [self.rm2_col,self.desc_col1,self.genename_col]=self.cluster_data.use_columns(['RM2','desc','genename'])
                [self.ontology_col,self.goid_col,self.desc_col2]=self.goclass_data.use_columns(['ontology','goid','desc'])
                glob.use_online_dictionaries(['GOIDELIC'])
                # parameter nicknames
                self.predictors=[x for x in glob.param['PANZ_PREDICTOR'].upper().split(',') if x != "DE"]
                self.predictor_cols={}
                for predictor in self.predictors:
                        self.predictor_cols[predictor]=self.goclass_data.use_columns([predictor+'_score',predictor+'_PPV',predictor+'_rank'])
                self.BESTCLUSTER=glob.param['PANZ_BESTCLUSTER']
                [self.a]=glob.use_operators(["RemoveRedundantGO"])

        def process(self,block):
                #sys.stderr.write("# This is BuildAnnoTable.process data: %i cluster: %i goclass: %i\n" %(self.data.nrows,self.cluster_data.nrows,self.goclass_data.nrows))
                if self.data.nrows==0: return # no hits
                qpid=block[0][self.qpid_col]
                # original description
                desc=self.glob.QUERYDESC
                ff=self.f.formfactor(desc)
                self.anno_data.append_row([qpid,'original_DE',self.glob.QUERYKINGDOM,'n.d.',str(ff),desc])
                # query sequence
                self.anno_data.append_row([qpid,'qseq','n.d.','n.d.','n.d.',self.glob.QUERYSEQUENCE])
                # accepted DE predictions
                for i in range(0,self.cluster_data.nrows):
                        if not self.cluster_data.row_status[i]: continue
                        row=self.cluster_data.block[i]
                        score=float(row[self.rm2_col])
                        desc=row[self.desc_col1]
                        ff=self.f.formfactor(desc)
                        if self.glob.QUERYKINGDOM == 'euk':
                                ppv=DE_PPV_euk(score)
                        else:
                                ppv=DE_PPV_bac(score)
                        self.anno_data.append_row([qpid,"DE",str(score),str(ppv),str(ff),desc])
                        gn=row[self.genename_col]
                        if gn != "": self.anno_data.append_row([qpid,"GN","n.d.","n.d.","n.d.",gn])
                        if self.BESTCLUSTER: break # best only
                # nonredundant GO predictions per ontology
                for predictor in self.predictors: self.go_anno(qpid,predictor)

        def go_anno(self,qpid,predictor):
                score_col=self.predictor_cols[predictor][0]
                ppv_col=self.predictor_cols[predictor][1]
                rank_col=self.predictor_cols[predictor][2]
                label=predictor
                # sort block
                self.goclass_data.sort_block(rank_col)
                # remove redundant
                if self.glob.param['PANZ_FILTER_OUTPUT']: self.a.process(self.goclass_data.block)
                # output rows
                for i in range(0,self.goclass_data.nrows):
                        if not self.goclass_data.row_status[i]: continue
                        row=self.goclass_data.block[i]
                        score=float(row[score_col])
                        ppv=row[ppv_col]
                        if score<=0: continue
                        goid=row[self.goid_col]
                        desc=row[self.desc_col2]
                        ec=None
                        kegg=None
                        if goid in self.glob.EC: ec=self.glob.EC[goid]
                        if goid in self.glob.KEGG: kegg=self.glob.KEGG[goid]
                        if ec == '.': ec=None
                        if kegg == '.': kegg=None
                        ontology=row[self.ontology_col]
                        self.anno_data.append_row([qpid,'%s_%s' %(ontology,label),str(score),str(ppv),goid,desc])
                        if ec: self.anno_data.append_row([qpid,'EC_%s' %label,str(score),str(ppv),ec,'GO:%s' %goid])
                        if kegg: self.anno_data.append_row([qpid,'KEGG_%s' %label,str(score),str(ppv),kegg,'GO:%s' %goid])


