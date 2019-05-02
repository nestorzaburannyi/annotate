from __future__ import print_function
from myoperator import RowOperator
from PannzerFunctions import FormFactor,Cleaner
import sys,math
from Read_and_Print import read_dict_data

class DE_evaluation(RowOperator):
        """
        DE evaluation. Reads reference into a dictionary. Reference must be prepared beforehand with wordweights operator. Predictions are DE.out by Pannzer (qpid, desc).

        Prepare reference of truth:

grep '^>' human_subset.fasta| perl -pe 's/ /\t/' | perl -pe 's/^>//' | python runsanspanz.py -f tab -m wordweights -c "qpid desc" -k 1000 > human_subset_de.tab

        Calculate DSM (description similarity measure) of predictions DE.out to reference:

python runsanspanz.py -m DSM -i DE.out -f tab --eval_TRUTH human_subset_de.tab -o ",DE.eval"

        Input: qpid, desc
        Output: qpid, desc, correct_desc, FF, informative/noninformative classes, DSM

        """
        def __init__(self,glob):
                sys.stderr.write("# Init DSM\n")
                self.glob=glob
                [self.data,self.summary_data]=glob.use_sheets(["data","summary"])
                [self.a]=self.glob.use_operators(["wordweights"])
                self.summary_cols=self.summary_data.use_columns(["qpid","FF_pred","FF_true","class_pred","class_true","TFIDF_pred","TFIDF_true","DSM","desc","correct_desc"])
                [self.qpid_col,self.desc_col,self.word_col,self.termidf_col,self.vectorlength_col,self.FF_col]=self.data.use_columns(["qpid","desc","word","termidf","vector_length","FF"])
                self.ff_cutoff=float(glob.param['PANZ_FFCUTOFF'])
                self.trueword=read_dict_data(glob.param['eval_TRUTH'],'qpid','word',header=True,UseColNames=True)
                self.truevectorlength=read_dict_data(glob.param['eval_TRUTH'],'qpid','vector_length',header=True,UseColNames=True)
                self.trueff=read_dict_data(glob.param['eval_TRUTH'],'qpid','FF',header=True,UseColNames=True)
                self.truedesc=read_dict_data(glob.param['eval_TRUTH'],'qpid','desc',header=True,UseColNames=True)

        def process(self,row):
                # link to reference
                qpid=row[self.qpid_col]
                if not qpid in self.trueword:
                        print("# Warning: %s not in reference\n" %qpid, file=sys.stderr)
                        return # no reference
                self.a.process(row) # add FF, TFIDF variables
                ff=row[self.FF_col]
                correct_ff=self.trueff[qpid]
                desc=row[self.desc_col]
                wordlist2=self.trueword[qpid]
                dsm_yy=float(self.truevectorlength[qpid])
                wordlist1=row[self.word_col]
                termidf1=row[self.termidf_col]
                dsm_xx=float(row[self.vectorlength_col])
                correct_desc=self.truedesc[qpid]
                # termidf-weighted cosine similarity
                dsm=self.DSM(wordlist1,termidf1,wordlist2,dsm_xx,dsm_yy)
                # write to summary table - fixed order
                datarow=[qpid,ff,correct_ff,self.informative(ff),self.informative(correct_ff),row[self.vectorlength_col],self.truevectorlength[qpid],str(dsm),desc,correct_desc]
                self.summary_data.append_row(datarow)

        def DSM(self,wordlist1,termidf1,wordlist2,l1,l2):
                " Calculate cosine similarity of termidf vectors between two strings."
                if l1==0 or l2==0: return(0.0)
                tmp={}
                x=wordlist1.split()
                z=termidf1.split()
                for i in range(0,len(x)): tmp[x[i]]=float(z[i])
                s=0.0
                for w in wordlist2.split():
                        if w in tmp: s+=tmp[w]*tmp[w]
                return(s/l1/l2)

        def informative(self,ff):
                if float(ff)<self.ff_cutoff: return('uncharacterized')
                return('informative')

