from myoperator import BlockOperator

class DEcluster(BlockOperator):
        """
        Cluster descriptions.

        Input: 'desc'
        Outputs: 'cleandesc','FF','clusid'
        """
        def __init__(self, glob):
                self.glob=glob
                [self.d1,self.c,self.e,self.f]=glob.use_operators(['FF','Cleandesc','TFIDF','Clustering'])

        def process(self,block):
                for row in block:
                        self.d1.process(row) # add form factor
                        self.c.process(row) # add cleandesc
                        self.e.process(row) # add word,wordcount,termidf
                self.f.process(block) # add clusid

