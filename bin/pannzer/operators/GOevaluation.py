# cut -f 2,5,7 goa_uniprot_noIEA.gaf

from myoperator import BlockOperator
from Read_and_Print import read_dict_data
from PannzerFunctions import Propagate
import sys,re,math

class GOevaluation(BlockOperator):
    """
Input: qpid, goid, PPV_RM3, PPV_argot, PPV_jac, PPV_hyge
Output: metric (Fmax, wFmax, Jacc, wJacc, Smin), PPV_RM3, PPV_argot, PPV_jac, PPV_hyge

Propagates GO predictions to parents, compares to truth (goa_truth_propagated) and computes GO-evaluation metrics
    """
    def __init__(self,glob):
        # set parameters
        self.glob=glob
        self.glob.use_online_dictionaries(['GOIDELIC'])
        truthfile=glob.param['eval_TRUTH']
        self.in_truth=read_dict_data(truthfile,'qpid','propagated',header=True,UseColNames=True)
        # write blockwise summmary
        [self.data,self.propagated_data,self.summary_data]=glob.use_sheets(["data","propagated","summary"])
        # column indices in GAF file
        self.scorefunctions=glob.param['eval_SCOREFUNCTIONS'].split()
        self.nscores=len(self.scorefunctions)
        tmp=self.data.use_columns(["qpid","goid"]+self.scorefunctions)
        self.qpid_col1=tmp[0]
        self.goid_col1=tmp[1]
        self.score_cols1=tmp[2:]
        self.propagated_data.use_columns(["qpid","acc","goid","in_truth","ontology"]+self.scorefunctions)
        self.precision=4
        self.scores={}
        self.summary_data.use_columns(["ontology","metric"]+self.scorefunctions)
        # cafa-evaluation
        self.T={}
        self.nP={}
        self.wT={}
        self.ntarget={}
        for ontology in ('BP','MF','CC'):
                self.T[ontology]=0
                self.nP[ontology]=0
                self.wT[ontology]=0.0
                self.ntarget[ontology]=0
        self.P={} # p[scorefunction]{ontology}[tau]
        self.wP={}
        self.TP={}
        self.wTP={}
        for x in self.scorefunctions:
                self.P[x]={} # ontology hash
                self.wP[x]={} # ontology hash
                self.TP[x]={} # ontology hash
                self.wTP[x]={} # ontology hash
                for ontology in ('BP','MF','CC'):
                        self.P[x][ontology]={} # tau hash
                        self.wP[x][ontology]={} # tau hash
                        self.TP[x][ontology]={} # tau hash
                        self.wTP[x][ontology]={} # tau hash
        self.target_not_in_truth={}

    def process(self,block):
        try:
                qpid=block[0][self.qpid_col1]
        except:
                return
        try:
                (db,acc,pid)=qpid.split("|") # tr|E0SNF6|E0SNF6_DICD3
        except:
                acc=qpid
        try:
                truth=self.in_truth[acc].split(',')
        except:
                self.target_not_in_truth[acc]=1
                return
        # accumulate truth counts [ only including proteins which have predictions ... ]
        x={}
        for goid in truth:
                if not goid in self.glob.ontology: continue
                ontology=self.glob.ontology[goid]
                self.T[ontology]+=1
                x[ontology]=1 # targets in ontology
                self.wT[ontology]+=float(self.glob.IC[goid])
        for ontology in x.keys(): self.ntarget[ontology]+=1 # block protein has truth in ontology
        # copy scores to hash
        self.scores={} # initialize block
        for row in block:
                goid=row[self.goid_col1]
                goid=re.sub(r'GO:','',goid)
                self.scores[goid]=[]
                # append truncated scores
                for col in self.score_cols1: self.scores[goid].append(row[col][0:self.precision])
        # propagate
        tmp=list(self.scores.keys())
        for goid in tmp:
                if not goid in self.glob.GOparents: continue # root
                for p in self.glob.GOparents[goid]:
                        if not p in self.scores:
                                self.scores[p]=[]
                                for col in range(0,self.nscores): self.scores[p].append(0.0)
        # propagate scores from direct assignments to parental lineage
        for goid in tmp:
                if not goid in self.glob.GOparents: continue # root
                # test parental lineages
                for p in self.glob.GOparents[goid]:
                        if p=='': continue
                        if p=='n.d.': continue
                        for col in range(0,self.nscores):
                                x=float(self.scores[goid][col])
                                if x > float(self.scores[p][col]): self.scores[p][col]=x
        # write propagated scores to propagated_data
        for goid in self.scores.keys():
                if goid in truth:
                        in_truth="1"
                else:
                        in_truth="0"
                if goid not in self.glob.ontology: continue
                ontology=self.glob.ontology[goid]
                if ontology=='': continue
                self.nP[ontology]+=1
                datarow=[qpid,acc,goid,in_truth,ontology]+self.scores[goid]
                self.propagated_data.append_row(datarow)
                ic=float(self.glob.IC[goid])
                for i in range(0,self.nscores):
                        x=self.scorefunctions[i]
                        tau=self.scores[goid][i]
                        if not tau in self.P[x][ontology]:
                                self.P[x][ontology][tau]=0
                                self.wP[x][ontology][tau]=0.0
                                self.TP[x][ontology][tau]=0
                                self.wTP[x][ontology][tau]=0.0
                        self.P[x][ontology][tau]+=1
                        self.wP[x][ontology][tau]+=ic
                        if in_truth=="1":
                                self.TP[x][ontology][tau]+=1
                                self.wTP[x][ontology][tau]+=ic

    def finalise(self):
        for ontology in ('BP','MF','CC'):
          for metric in ('Fmax','wFmax','Jacc','wJacc','Smin'):
            datarow=[ontology,metric]
            datarow1=[ontology,'avg-ru']
            datarow2=[ontology,'avg-mi']
            for scorefun in self.scorefunctions:
                tmp=list(self.P[scorefun][ontology].keys())
                tmp.sort(key=float)
                tmp.reverse()
                if metric == "Smin":
                        opt=float('inf')
                else:
                        opt=0.0
                opt_ru=0.0
                opt_mi=0.0
                P=0.0
                wP=0.0
                TP=0.0
                wTP=0.0
                for tau in tmp:
                        P+=self.P[scorefun][ontology][tau]
                        wP+=self.wP[scorefun][ontology][tau]
                        TP+=self.TP[scorefun][ontology][tau]
                        wTP+=self.wTP[scorefun][ontology][tau]
                        if metric == "Smin":
                                mi=wP-wTP
                                ru=self.wT[ontology]-wTP
                                S=math.sqrt(ru*ru+mi*mi)
                                if S<opt:
                                        opt=S
                                        opt_ru=ru
                                        opt_mi=mi
                                continue
                        if metric == "Fmax" and P>0.0:
                                S=2*TP/float(self.T[ontology]+P)
                        elif metric == "wFmax" and wP>0.0:
                                S=2*wTP/(self.wT[ontology]+wP)
                        elif metric == "Jacc" and P>0.0:
                                S=TP/float(self.T[ontology]+P-TP)
                        elif metric == "wJacc" and wP>0.0:
                                S=wTP/(self.wT[ontology]+wP-wTP)
                        if S>opt: opt=S
                if self.T[ontology]==0:
                        opt=0.0 # empty True set -> predict nothing!; correction for Smin
                elif metric == 'Smin' and self.ntarget[ontology]>0:
                        opt=opt/self.ntarget[ontology]
                        opt_ru=opt_ru/self.ntarget[ontology]
                        opt_mi=opt_mi/self.ntarget[ontology]
                # normalize Smin
                datarow.append(str(opt))
                if metric == 'Smin':
                        datarow1.append(str(opt_ru))
                        datarow2.append(str(opt_mi))
            if metric == 'Smin': # output ru, mi components
                self.summary_data.append_row(datarow1)
                self.summary_data.append_row(datarow2)
            self.summary_data.append_row(datarow)
          sys.stderr.write("# ontology %s: %i target proteins, %i nodes in True set, %i total nodes in Predicted set, information content of T set: %f\n" %(ontology,self.ntarget[ontology],self.T[ontology],self.nP[ontology],self.wT[ontology]))
        # output
        self.summary_data.output(result=True)
        sys.stderr.write("# %i targets not found in truth\n" %len(self.target_not_in_truth))

