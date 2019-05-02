from myoperator import BlockOperator
import sys

class Pannzer(BlockOperator):
        """
        Generate description (DE) or Gene Ontology (GO) predictions using Pannzer scoring functions.
        Results are saved to cluster_data (DE predictions) or goclass_data (GO predictions) and
        anno_data (concise results).

        Parameter METHOD is "DE", "GO" or "DEGO". If parameter BESTCLUSTER is True, GO predictions
        uses only the best description-cluster. If parameter FILTER_OUTPUT is True, empty clusters
        are removed from DE output and redundant GO classes are removed from GO output.
        """
        def __init__(self, glob,verbose=True):
                sys.stderr.write("# Init Pannzer\n")
                self.glob=glob
                # define nicknames for parameters
                self.FILTER_OUTPUT=glob.param['PANZ_FILTER_OUTPUT']
                self.BESTCLUSTER=glob.param['PANZ_BESTCLUSTER']
                x=glob.param['PANZ_PREDICTOR'].upper()
                self.do_DE= "DE" in x
                self.do_B2GO = ("B2GO" in x)
                x = [a for a in x.split(',') if ( a != "DE" and a != "B2GO" ) ]
                self.do_GO= len(x)>0
                sys.stderr.write('# len(x) %d\n# do_B2GO ' %len(x) + str(self.do_B2GO) + "\n")
                self.do_clustering=(self.do_DE or self.BESTCLUSTER)
                sys.stderr.write("# bestcluster: "+str(self.BESTCLUSTER)+" do_DE: "+str(self.do_DE)+" do_GO: "+str(self.do_GO)+"\n")
                #
                [self.data,self.cluster_data]=glob.use_sheets(['data','cluster'])
                [self.a,self.d1,self.d,self.d2,self.b,self.c]=glob.use_operators(['Taxonomy','FF','Filter','DE_status','RM1','Cleandesc'])
                # DE-clustering operations: TFIDF+Clustering
                [self.e,self.f]=glob.use_operators(['TFIDF','Clustering'])
                if self.do_DE:
                        # DE-prediction operations: Wordscores+Cluster_GSZ+RM2+output_DE
                        [self.g,self.h,self.i,self.k]=glob.use_operators(['Wordscores','Cluster_GSZ','RM2','output_DE'])
                # GO-prediction operations: GOrimpsRM3+RedundantGO
                [self.l,self.m,self.n]=glob.use_operators(['GOrimpsu','RM3','BLAST2GO2'])
                [self.o]=glob.use_operators(['BuildAnnoTable'])
                # column indices
                [self.rm1sum_col,self.clustersize_col,self.rm2_col,self.clusid_col2]=self.cluster_data.use_columns(['cluster_RM1sum','cluster_size','RM2','clusid'])
                [self.clusid_col1,self.status_col,self.bits_col,self.desc_col,self.qseq_col,self.status_col]=self.data.use_columns(['clusid','DE_status','bits','desc','qseq','status'])

        def process(self,block):
                # sort SSRL
                if len(block)==0: return
                try:
                  block.sort(key = lambda x: float(x[self.bits_col]), reverse=True)
                except: # patch corrupted transmission
                  for row in block:
                        try: 
                          float(row[self.bits_col])
                        except:
                          row[self.bits_col]="0.0"
                  block.sort(key = lambda x: float(x[self.bits_col]), reverse=True)
                for row in block: self.d1.process(row) # add form factor
                self.d.process(block)  # Filter status
                for row in block: self.d2.process(row) # add DE_status = FF and status
                self.a.process(block) # add taxdist
                self.set_QueryData() # capture self.glob.QUERYSPECIES etc.
                for row in block: # common-to-all operations
                        self.b.process(row) # add RM1
                        self.c.process(row) # add cleandesc
                if self.do_clustering: # generate DE clusters
                        for row in block: self.e.process(row) # add word,wordcount,termidf
                        self.f.process(block) # add clusid
                if self.do_DE: # generate DE predictions
                        self.g.process(block) # add valencia_wordscore, jaccard_wordscore
                        self.h.process(block) # add cluster_GSZ
                        self.i.process(block) # add RM2
                        self.k.process(block) # select DE predictions
                        if self.FILTER_OUTPUT: self.removeBadClusters() # hide bad DEs from output
                if self.do_GO: # generate GO predictions
                        if self.BESTCLUSTER: self.select_bestcluster() # set status=False except best cluster
                        for row in block: self.l.process(row) # GOrimpsu
                        if self.do_GO: self.m.process(block) # generate GO predictions
                if self.do_GO or self.do_DE:
                        self.o.process(block) # fill summary table
                if self.do_B2GO: self.n.process(block) # generate BLAST2GO predictions

        def removeBadClusters(self):
                for i in range(0,self.cluster_data.nrows):
                        if self.cluster_data.block[i][self.rm1sum_col]=="0" or \
                           self.cluster_data.block[i][self.clustersize_col]=="0":
                                self.cluster_data.row_status[i]=False # data is string

        def select_bestcluster(self):
                bestcluster="0"
                bestrm2=0.0
                for row in self.cluster_data.block:
                        rm2=float(row[self.rm2_col])
                        if rm2>bestrm2:
                                bestrm2=rm2
                                bestcluster=row[self.clusid_col2]
                # set status=False except best cluster
                for row in self.data.block:
                        row[self.status_col]= str( row[self.clusid_col1] == bestcluster )

        def selectBestDE(self):
                if self.data.nrows<1: return
                self.cluster_data.row_status[0]=True
                for i in range(1,self.cluster_data.nrows):
                        self.cluster_data.row_status[i]=False

        def set_QueryData(self):
                """
                Set QUERYSPECIES, QUERYDESC, QUERYSEQUENCE,QUERYKINGDOM
                using taxonomy
                """
                # initialize
                self.glob.QUERYSEQUENCE=''
                self.glob.QUERYDESC=''
                self.glob.QUERYKINGDOM='bac'
                if "Eukaryota" in self.a.querylineage: self.glob.QUERYKINGDOM='euk'
                # column indices
                [desc_col,qseq_col,isquery_col,status_col]=self.data.get_col_index(['desc','qseq','isquery','status'])
                for row in self.data.block:
                    if row[isquery_col]=="1":
                        row[status_col]="False" # exclude query from prediction
                        self.glob.QUERYDESC=row[desc_col].strip()
                        self.glob.QUERYSEQUENCE=row[qseq_col].strip()
                        #sys.stderr.write("# Query species = %s, desc = %s, lseq = %i, kingdom = %s\n" %(glob.QUERYSPECIES, glob.QUERYDESC,len(glob.QUERYSEQUENCE),glob.QUERYKINGDOM))
                        break

