from myoperator import RowOperator
import math,sys

class TFIDF(RowOperator):
        """
        Generate cleaned word vectors and respective count vector
        and term idf vector.

        term idf = log(total number of terms in the dataset /
                       number of documents where terms appears)
        # we assume any term occurs just once in any document

        Creates data columns 'word','wordcount','termidf'
        Inputs: data column 'cleandesc'
        """
        def __init__(self, glob):
                sys.stderr.write("# Init TFIDF\n")
                self.glob = glob
                # define nicknames for column indices
                [self.cleandesc_col,self.word_col,self.wordcount_col,self.termidf_col,self.vectorlength_col]=self.glob.use_sheet("data").use_columns(['desc','word','wordcount','termidf','vector_length'])
                # use online dictionary. Object handles in glob are hardcoded
                self.glob.use_online_dictionaries(["WORDCOUNT"])

        def process(self,row, verbose=False):
                cleandesc=row[self.cleandesc_col]
                # create word, wordcount, termidf vectors
                tmp=cleandesc.upper().split(" ")
                words=[]
                counts=[]
                termidf=[]
                tmp.sort()
                ssq=0.0
                for word in tmp:
                        if not word in self.glob.wordcounts:
                                if verbose: sys.stderr.write("# Warning: unknown word %s\n%s\n" %(word,tmp))
                                continue
                        words.append(word)
                        cnt=self.glob.wordcounts[word]
                        if not cnt: cnt="1"
                        counts.append(str(cnt))
                        # PK's script uses nwordtotal instead of nprot
                        x=math.log(self.glob.nwordtotal/float(cnt))
                        ssq+=x*x
                        termidf.append(str(x))
                row[self.word_col]=" ".join(words)
                row[self.wordcount_col]=" ".join(counts)
                row[self.termidf_col]=" ".join(termidf)
                row[self.vectorlength_col]=str(math.sqrt(ssq))
