from myoperator import BlockOperator
import sys

class Wordscores(BlockOperator):
        """
        Calculate
                valencia_wordscore = bitsum-weighted word counts (average over words in description)
                jaccard_wordscore = jaccard-weighted valencia_wordscore (average over words in description)

        Creates data columns 'jaccard_wordscore','valencia_wordscore'
        Inputs: data columns 'bits','status','word','wordcount'
        """
        def __init__(self, glob):
                sys.stderr.write("# Init Wordscores\n")
                self.glob = glob
                [self.jac_col,self.val_col,self.bits_col,self.status_col,self.word_col,self.wordcount_col]=glob.use_sheet("data").use_columns(['jaccard_wordscore','valencia_wordscore','bits','DE_status','word','wordcount'])
                glob.use_online_dictionaries(["WORDCOUNT"])
                self.JACCARD_MINCOUNT=glob.param['PANZ_JACCARD_MINCOUNT']
                # bad words in description
                self.UNINFORMATIVE =['HYPOTHETICAL', 'UNCHARACTERIZED', 'PUTATIVE', \
                'CONTIG', 'PREDICTED', 'PROBABLE', 'FRAGMENT', 'GENOME', \
                'PROTEIN', 'CHROMOSOME', 'POSSIBLE', 'SIMILAR', 'PROTEINS', \
                'HOMOLOG', 'POSSIBLE', 'CONSERVED', 'HOMOLOGOUS', 'COMPLETE', \
                'SHOTGUN', 'CDNA', 'FAMILY']

        def process(self,block):
                # valencia_wordscore = bitsum-weighted word counts
                # jaccard_wordscore = jaccard-weighted valencia_wordscore
                word_bitsum={}
                word_count={}
                total_bits=1.0 # offset to avoid division by zero
                for row in block:
                        if row[self.status_col] == "False": continue # skip filtered rows
                        bits=float(row[self.bits_col]) # each word weighted by bitscore
                        for word in row[self.word_col].split():
                                if word in self.UNINFORMATIVE: continue # skip bad word
                                if not word in word_bitsum:
                                        word_bitsum[word]=0.0
                                        word_count[word]=0
                                word_bitsum[word]+=bits
                                word_count[word]+=1.0
                                total_bits+=bits
                # normalize
                for word in word_bitsum: word_bitsum[word]=word_bitsum[word]/total_bits
                for word in word_count:
                        cnt=1
                        if word in self.glob.wordcounts: cnt=self.glob.wordcounts[word]
                        if cnt<self.JACCARD_MINCOUNT: cnt=self.JACCARD_MINCOUNT
                        word_count[word]=word_count[word]/cnt # jaccard
                # calculate average word scores over description
                for row in block:
                        val=0.0
                        jac=0.0
                        de_count=0
                        if row[self.status_col] == 'True': # exclude filtered rows
                                for word in row[self.word_col].split():
                                        if word in self.UNINFORMATIVE: continue # skip bad word
                                        val+=word_bitsum[word] # normalized
                                        jac+=word_count[word]*word_bitsum[word] # jaccard*valencia
                                        de_count+=1
                        # word average
                        if de_count<1: de_count=1
                        row[self.val_col]=str(val/de_count) # average over words in description
                        row[self.jac_col]=str(jac/de_count)

