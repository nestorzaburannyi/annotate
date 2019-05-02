from myoperator import BlockOperator
import sys, numpy, math
from PannzerFunctions import sampleStats,logmod
from Hypergeometric import calculate_p_value_for_hypergeometric
import GSZ

class RM3(BlockOperator):
        """
        Computes GO prediction scores for all predictors in glob.param['PANZ_PREDICTOR'] list.

        Uses private variables and dictionaries for scoring functions.
        Scoring functions are named this.PREDICTGOR_score.
        Scoring functions are selected by setting glob.param['PANZ_PREDICTOR'].

        Creates goclass_data columns qpid,ontology,goid,desc.
        Creates goclass_data columns predictor_score, predictor_PPV, predictor_rank.

        Creates data column 'GO_status' to include at most MAXHITS hits with GO annotation.

        Inputs: data columns 'RM1', 'status', 'GOclass', 'GOclass_count','bits'
        """

        def __init__(self, glob):
                sys.stderr.write("# Init RM3\n")
                self.glob=glob
                self.hyper_p_cache={}
                [self.data,self.goclass_data]=glob.use_sheets(['data','goclass'])
                [self.status_col,self.rm1_col,self.go_col,self.gocount_col,self.bits_col,self.qpid_col]=self.data.use_columns(['GO_status','RM1','GOclass','GOclass_count','bits','qpid'])
                self.goclass_data.use_columns(["qpid","ontology","goid","desc"])
                [self.ontology_col,self.rank_col]=self.goclass_data.get_col_index(['ontology','rank'])
                glob.use_online_dictionaries(["GOIDELIC"])
                # nickname parameters
                self.MAXHITS=glob.param['PANZ_MAXHITS']
                self.MAXCACHEKEYS=glob.param['HYGE_MAXCACHEKEYS']
                # predictors append (score,ppv,rank) to goclass_data
                self.predictors=[x for x in glob.param['PANZ_PREDICTOR'].upper().split(',') if x != "DE"]
                self.predictor_cols={} # [0]=score_col, [1]=ppv_col, [2]=rank_col
                self.predictor_funcs=[]
                for predictor in self.predictors:
                        self.predictor_cols[predictor]=self.goclass_data.use_columns([predictor+'_score',predictor+'_PPV',predictor+'_rank'])
                        func=getattr(self,predictor+'_score',None)
                        self.predictor_funcs.append(func)

        def process(self,block):
                if self.data.nrows == 0: return
                qpid=block[0][self.qpid_col]
                # GO_status is True is sbjct has GO annotations. Hitlist filtered here.
                n=0
                for row in block:
                        if row[self.status_col]=="False": continue # already rejected by Filter
                        if n>=self.MAXHITS: continue
                        if len(row[self.go_col])>0:
                                row[self.status_col]="True"
                                n+=1
                        else:
                                row[self.status_col]="False"
                self.sampleSize=n
                # store useful attributes in private dictionaries (per block)
                self.score={} # sum of RM1-scores per GOclass
                self.bitsum={} # sum of bits per GOclass
                self.sizeInBackGround={} # GOclass_count per GOclass [could take directly from GOIDELIC]
                self.goclasscount={} # for jaccard term
                self.totalbitsum=1.0
                self.sampleSizeOntology={"BP":0,"CC":0,"MF":0,"":0}
                for row in block:
                        if row[self.status_col] == 'False': continue
                        goidlist=row[self.go_col].rstrip().split()
                        gocntlist=row[self.gocount_col].rstrip().split()
                        rm1=float(row[self.rm1_col])
                        if len(goidlist) != len(gocntlist):
                                continue
                        for goid in goidlist:
                                if not goid in self.sizeInBackGround:
                                        self.score[goid]=0.0
                                        self.bitsum[goid]=0.0
                                        self.sizeInBackGround[goid]=int(self.glob.GOcounts[goid]) #int(gocntlist[i])
                                        self.goclasscount[goid]=0
                                self.score[goid]+=rm1
                                if row[self.bits_col]=="n.d.": continue
                                self.bitsum[goid]+=float(row[self.bits_col])
                                self.goclasscount[goid]+=1
                                self.sampleSizeOntology[self.glob.ontology[goid]]+=1
                        if row[self.bits_col]=="n.d.": continue
                        self.totalbitsum+=float(row[self.bits_col])
                # calculacte GSZ, RM3 per goid and save in goclass_data
                for goid in self.score:
                        # reject root classs
                        if self.sizeInBackGround[goid] == self.glob.rootcount[self.glob.ontology[goid]]: continue
                        if self.sizeInBackGround[goid]<1: self.sizeInBackGround[goid]=1
                        # useful private variables
                        self.jsc=self.goclasscount[goid]/float(self.sizeInBackGround[goid]+self.sampleSize-self.goclasscount[goid])
                        self.ontology=self.glob.ontology[goid]
                        self.rootcount=self.glob.rootcount[self.ontology]
                        # assuming fixed column order
                        datarow=[qpid,self.ontology,goid,self.glob.godesc[goid]]
                        # for each scoring function: append (score,ppv,rank) to datarow
                        for func in self.predictor_funcs:
                                if func is None: continue
                                score,ppv=func(goid)
                                # assuming fixed column order
                                datarow.append(str(score))
                                datarow.append(str(ppv))
                                datarow.append("0") # rank not known yet
                        # save in summary table
                        self.goclass_data.append_row(datarow)
                # reverse sort goclass_data block on PPV
                for predictor in self.predictors:
                        score_col=self.predictor_cols[predictor][1]
                        rank_col=self.predictor_cols[predictor][2]
                        self.rank_score(score_col,rank_col, True)
                # flush big cache
                if len(self.hyper_p_cache) > self.MAXCACHEKEYS: self.hyper_p_cache={}

        def rank_score(self, rm3_col, rank_col, sort_reverse):
                # reverse sort goclass_data block on RM3
                self.goclass_data.sort_block(rm3_col,reverse=sort_reverse)
                # overwrite rank per ontology in sorted block
                ontology_rank={}
                for row in self.goclass_data.block:
                        ontology=row[self.ontology_col]
                        if not ontology in ontology_rank: ontology_rank[ontology]=0
                        ontology_rank[ontology]+=1
                        row[rank_col]=str(ontology_rank[ontology])

        def RM3_score(self,goid):
                """Calculates regression model 3 score.

        RM3 = 0.36 + 0.01 * log(JSC) + 0.02 * logmod(GSZ) - 3.67e-5 * sqrt(|GO|) - 0.007 * sqrt(|SSRL|)

        where
                JSC = number of GOid instances in hitlist / number of GOid instances in database

                GSZ uses RM1 scores summed over GOid instances. Using bitsum instead of RM1sum does not change ranks, in practice.

                This is the Pannzer regression model RM3 from CAFA2. Returns (raw score, PPV)."""
                # sample is hitlist: sampleSize, sampleScoreMean, sampleScoreVariance
                (sampleSize,sampleScoreMean,sampleScoreVariance)=sampleStats(self.data.block,self.status_col,self.rm1_col)
                gsz=GSZ.calculateGSZscore(sampleSize,self.score[goid],sampleScoreMean,sampleScoreVariance,self.sizeInBackGround[goid],self.rootcount)
                if numpy.isinf(gsz) or numpy.isnan(gsz): gsz=0.0
                rm3=0.36 + 0.01 * numpy.log(self.jsc) + 0.02 * logmod(gsz) - 3.67e-5 * numpy.sqrt(self.sizeInBackGround[goid]) - 0.007 * numpy.sqrt(sampleSize) # sampleSize == n
                if rm3<0: rm3=0.0
                ppv=GO_PPV(rm3)
                return(rm3,ppv)

        def ARGOT_score(self,goid):
                "Bitscore-weighted information content. Returns (raw score, PPV)."
                argot=self.bitsum[goid]/self.totalbitsum*numpy.log( float(self.rootcount) / self.sizeInBackGround[goid])/numpy.log(2.0)
                ppv=GO_argot_PPV(argot)
                return(argot,ppv)

        def JAC_score(self,goid):
                "Bitscore-weighted Jaccard coefficient of similarity. Returns (raw score, PPV)."
                jac=self.bitsum[goid]/self.totalbitsum*self.jsc
                ppv=GO_jac_PPV(jac)
                return(jac,ppv)

        def HYGE_score(self,goid):
                "Probability from hypergeometric distribution. Returns (raw score, PPV)."
                hyge=calculate_p_value_for_hypergeometric(self.goclasscount[goid],self.rootcount,self.sampleSizeOntology[self.ontology],self.sizeInBackGround[goid],hyper_p_cache=self.hyper_p_cache)
                if hyge <= 0: return(0.0,0.0)
                hyge=-math.log(hyge)
                ppv=GO_slow_hyge_PPV(hyge)
                return(hyge,ppv)

        def BG_score(self,goid):
                "Background frequency of GO classes in sequence neighborhood. Returns (raw score, PPV)."
                bg=float(self.sizeInBackGround[goid])/self.rootcount
                ppv=bg
                return(bg,ppv)

        def SANS_score(self,goid):
                "Frequency of GO classes in sequence neighborhood. Returns (raw score, PPV)."
                obs=float(self.goclasscount[goid])/self.sampleSize
                ppv=min(0.99,obs)
                return(obs,ppv)

def GO_PPV(x):
                "argument is float(RM3)"
                if x<0.466: return(0.1205+1.6223*x)
                return(min(1.0,-0.5582+5.6307*x-5.4758*x*x))

def GO_argot_PPV(x):
                "argument is float(RM3_argot)"
                if x < 0: return(0.0)
                z=math.sqrt(x)
                if z<3.709538: return(0.3047 + 0.1452*z)
                return(min(1.0,0.67302 + 0.04591*z))

def GO_jac_PPV(x):
                "argument is float(RM3_jac)"
                if x <= 0: return(0.0)
                z=math.log(x)
                if z >= -3.527961: return(min(1.0,0.89756+0.01835*z))
                return(min(1.0,max(0.0,0.97768+0.04106*z)))

def GO_hyge_PPV(x):
                "argument is float(RM3_hyge) where RM3_hyge is -log(hypergeometric_pvalue)"
                # this is fitted to "fast hyge" with nprot, sampleSize
                # "slow hyge" is more accurate, i.e. reaches higher PPV values
                if x<1: x=1
                return(min(1.0,0.32608+0.05392*math.log(x)))

def GO_slow_hyge_PPV(x):
                "argument is float(RM3_hyge) where RM3_hyge is -log(hypergeometric_pvalue)"
                if x<1: x=1
                z=math.log(x)
                if z>0.6155868: return(min(1.0,0.51099+0.06829*z))
                return(min(1.0,max(0.0,0.2557+0.483*z)))

