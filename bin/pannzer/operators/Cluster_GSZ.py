from myoperator import BlockOperator
import sys
import numpy
import GSZ
from PannzerFunctions import sampleStats

class Cluster_GSZ(BlockOperator):
        """
        Calculate GSZ score: calculateGSZscore(sampleSize, score, sampleScoreMean, sampleScoreVariance, sizeInBackGroud, backGroundSize)

        sampleSize = size of filtered hitlist
        score = sum of RM1 in cluster
        sampleScoreMean = mean of RM1 in filtered hitlist
        sampleScoreVariance = variance of RM1 in filtered hitlist
        sizeInBackGround = sum of DE-counts of DEs occurring in cluster
        backGroundSize =nprot

        Creates cluster_data columns 'qpid','clusid''cluster_GSZ','cluster_RM1sum','cluster_desccount','cluster_size'
        Inputs: data columns 'RM1','DE_status','cleandesc','clusid','bits'
        """

        def __init__(self, glob):
                sys.stderr.write("# Init Cluster_GSZ\n")
                self.glob = glob
                [self.data,self.cluster_data]=glob.use_sheets(["data","cluster"])
                [self.rm1_col,self.status_col,self.cleandesc_col,self.clusid_col,self.bits_col,self.qpid_col]=self.data.use_columns(['RM1','DE_status','cleandesc','clusid','bits','qpid'])
                self.outcol=self.cluster_data.use_columns(['qpid','clusid','cluster_GSZ','cluster_RM1sum','cluster_size','cluster_desccount'])
                glob.use_online_dictionaries(['WORDCOUNT','DESCCOUNT','NPROT'])
                self.cluster_data.hide_from_output(['clusid'])

        def process(self,block):
                nrows=len(block)
                if nrows==0: return # no hits
                qpid=block[0][self.qpid_col]
                # sample is hitlist: sampleSize, sampleScoreMean, sampleScoreVariance
                (sampleSize,sampleScoreMean,sampleScoreVariance)=sampleStats(block,self.status_col,self.rm1_col)
                # per cluster: score and sizeInBackGround
                sizeInBackGround={}
                de_seen={}
                score={}
                bitsum={}
                size={}
                for row in block:
                        clusid=row[self.clusid_col]
                        if not clusid in score:
                                sizeInBackGround[clusid]=0
                                score[clusid]=0.0
                                bitsum[clusid]=0.0
                                size[clusid]=0
                        if row[self.status_col] == 'False': continue
                        rm1=float(row[self.rm1_col])
                        score[clusid]+=rm1
                        if row[self.bits_col]=="n.d.": continue
                        bits=float(row[self.bits_col])
                        bitsum[clusid]+=bits
                        size[clusid]+=1
                        # description counts
                        cleandesc=row[self.cleandesc_col]
                        if not cleandesc in de_seen:
                                if cleandesc in self.glob.desccounts:
                                        cnt=self.glob.desccounts[cleandesc]
                                        if not cnt: cnt="1"
                                        sizeInBackGround[clusid]+=int(cnt)
                                else:
                                        sizeInBackGround[clusid]+=1 # description missing from counts but found in hitlist
                                        sys.stderr.write("# WARNING: cleandesc not in desccounts: %s %i %s\n" %(clusid,sizeInBackGround[clusid],cleandesc))
                                de_seen[cleandesc]=True
                # calculate GSZ per cluster
                gsz={}
                for clusid in score:
                        # reject root class
                        if size[clusid] == sizeInBackGround[clusid]: continue
                        # sanity checks: cluster size <= bg <= nprot
                        if size[clusid] > sizeInBackGround[clusid]:
                                sys.stderr.write("# ERROR: Cluster_GSZ sample > background %s %d %d %i %i\n" %(clusid,score[clusid],bitsum[clusid],size[clusid],sizeInBackGround[clusid]))
                        # GSZ returns zero if sizeInBackGround<2
                        gsz=GSZ.calculateGSZscore(sampleSize,score[clusid],sampleScoreMean,sampleScoreVariance,sizeInBackGround[clusid],self.glob.nprot)
                        if numpy.isnan(gsz) or numpy.isinf(gsz): gsz=0.0
                        # save qpid,clusid,GSZ, cluster_RM1sum, cluster_desccount, cluster_size per cluster
                        self.cluster_data.append_row([row])
                        self.cluster_data.block[-1][self.outcol[0]]=qpid
                        self.cluster_data.block[-1][self.outcol[1]]=clusid
                        self.cluster_data.block[-1][self.outcol[2]]=str(gsz)
                        self.cluster_data.block[-1][self.outcol[3]]=str(score[clusid])
                        self.cluster_data.block[-1][self.outcol[4]]=str(size[clusid])
                        self.cluster_data.block[-1][self.outcol[5]]=str(sizeInBackGround[clusid])

