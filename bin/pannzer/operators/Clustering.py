from __future__ import print_function
from myoperator import BlockOperator
import sys
try:
  import fastcluster
  import numpy
except ImportError as ie:
  print(str(ie),file=sys.stderr)
  sys.exit(1)

class Clustering(BlockOperator):
        """
        Add unique description labels to data.
        Cluster unique descriptions using idf-weighted cosine distance
        Add cluster labels to data.

        Creates data columns 'clusid', 'descid'
        Inputs: data columns 'cleandesc','termidf','word'
        """
        def __init__(self, glob):
                sys.stderr.write("# Init Clustering\n")
                self.data=glob.use_sheet("data")
                [self.termidf_col,self.cleandesc_col,self.word_col,self.clusid_col,self.descid_col]=self.data.use_columns(['termidf','cleandesc','word','clusid','descid'])
                # parameter
                self.CLUSTERING_CUTOFF=glob.param['PANZ_CLUSTERING_CUTOFF']

        def process(self,block):
                nrows=len(block)
                if nrows==0: return
                nu=self.unique_descriptions(block)
                # special case: only one unique description in block
                if nu==1:
                        for row in block: row[self.clusid_col]="1"
                        return
                # * average linkage of unique descriptions is not equivalent
                #       with clustering all rows *

                # create X(n,d) word vectors for fastcluster
                D=self.wordspace(block)
                if False:
                        colsums=[]
                        for x in D[0]: colsums.append(0.0)
                        rowsums=[]
                        for row in D:
                                s=0.0
                                i=0
                                for x in row:
                                        s+=x
                                        colsums[i]+=x
                                        i+=1
                                rowsums.append(s)
                        print("colsums",colsums)
                        print("rowsums",rowsums)
                # create hierarchical clustering Y
                try:
                        Y=fastcluster.linkage(D,method='average',metric='cosine',preserve_input=False)
                        del D
                        # cut tree and traverse subtrees to assign cluster identifiers
                        roots=self.cuttree(Y)
                        clusid=self.getclusid(nrows,roots,Y)
                except:
                        sys.stderr.write("# Warning: fastcluster returned error. Setting single cluster\n")
                        # put all hits in one cluster
                        clusid=numpy.empty((nrows),dtype=numpy.int)
                        for i in range(0,nrows): clusid[i]=0
                # add clusid to block data
                for i in range(0,nrows):
                        block[i][self.clusid_col]=str(clusid[i])

        def unique_descriptions(self,block):
                """fills 'descid' column"""
                # determine unique cleaned descriptions
                x={}
                for row in block: x[row[self.cleandesc_col]]=1
                descid=0
                descmap={}
                for a in x.keys():
                        descmap[a]=str(descid)
                        descid+=1
                for row in block:
                        row[self.descid_col]=descmap[row[self.cleandesc_col]]
                return(descid)

        def wordspace(self,block):
                """Returns feature vector matrix D(data.nrows,nword).
                Features are all words that occur in block.
                """
                # map words
                nword=0
                wordcol={}
                for row in block:
                        for word in row[self.word_col].split():
                                if not word in wordcol:
                                        wordcol[word]=nword
                                        nword+=1
                # create feature vectors
                nrows=len(block)
                D=numpy.empty((nrows,nword),dtype=numpy.double)
                for i in range(0,nrows):
                        for j in range(0,nword): D[i,j]=0.0
                        idflist=block[i][self.termidf_col].split()
                        wlist=block[i][self.word_col].split()
                        for k in range(0,len(wlist)):
                                if idflist[k]=='n.d.': idflist[k]="0.0001" # HACK
                                D[i,wordcol[wlist[k]]]=idflist[k]
                return(D)

        def cuttree(self,Y):
                # cut tree at dsmcutoff: subtree root nodes are anything farther
                roots=[]
                i=len(Y)
                while i>0:
                        i=i-1
                        if Y[i,2]>self.CLUSTERING_CUTOFF:
                                roots.append(int(Y[i,0]))
                                roots.append(int(Y[i,1]))
                        else:
                                break
                return(roots)

        def getclusid(self,nrows,roots,Y):
                # traverse subtrees
                clusid=numpy.empty((nrows),dtype=numpy.int)
                for i in range(0,nrows): clusid[i]=0
                label=0
                stack=[]
                for x in roots:
                        label=label+1
                        stack.append(x)
                        while len(stack)>0:
                                x=stack.pop()
                                if x<nrows:
                                        clusid[x]=label
                                elif x>=nrows:
                                        stack.append(int(Y[x-nrows,0]))
                                        stack.append(int(Y[x-nrows,1]))
                return(clusid)

