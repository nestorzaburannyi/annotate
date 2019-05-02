import re,sys
import numpy,math

def Propagate(goidlist,GOparentsDict):
    """
    Input list of GO identifiers
    Return list inclusing all GO parents
    """
    parents={}
    stack=[]
    for goid in goidlist: stack.append(goid)
    while len(stack)>0:
        goid=stack.pop()
        parents[goid]=1
        if not goid in GOparentsDict: continue
        tmp=GOparentsDict[goid]
        if tmp=='': continue # root
        for p in tmp.split(','):
                if not p in parents: stack.append(p)
    tmp=list(parents.keys())
    tmp.sort()
    return(tmp)

def Cleaner(description, remove_abbr = False):
    """
    Removes unwanted characters, spaces, locus tags and one-word descriptions
    withs letters and numbers. Returns upper case version of the desription.
    """
    description = description.strip()

    #LH remove uniprot tags OS=... GN=... PE=...
    description=re.sub(r' \w{2}=.*$',"",description)

    #remove unwanted characters
    #"""
    for a in ["(", ")", "[", "]", "{", "}"]:
        description = description.replace(a, "")
    #"""
    description = re.compile(r"[\"/=\\:;,\|%\.#!<>\?_'\^\+\-]").sub(" ", description)
    description = re.compile(r"[\(\)\[\]\{\}]").sub("", description)

    description = re.sub(r'-{2,10}', "-", description)   #if description has "----" converts to "-"
    #replace these characters with space
    #"""
    for a in ["\"", "/", "=", "\\", ":", ";", ",", "|", "%", ".", "#", "!", ">", "<", "?", "_", "'", "-"]:
        description = description.replace(a, " ")
    #"""
            ##### CLEANING LOCUS TAGS ##########
    clean_de = re.sub(r'\w+_\w+\s?','',description)
    if len(clean_de) == 0:
        clean_de = "Hypothetical protein" #>> uniprot uses "uncharacterized protein"
    ##### CLEAN ONE WORD DESCRIPTIONS WITH LETTERS AND NUMBERS ######
    split_array = clean_de.split()
    pattern_digit = re.compile(r'[0-9]')
    pattern_dash = re.compile(r"-")
    if len(split_array) == 1:
        if pattern_digit.search(clean_de) and not re.search("-", clean_de):
            clean_de = "Hypothetical protein"
            description = re.sub(r'\s{2,10}', " ", clean_de)  #if description has "    " converts to " "
    ######### REMOVE WORDs WITH NUMBERS AND LETTERS ## FOR TFIDF CALCULATIONS
    else:
      if remove_abbr:
        list_out = []
        for tmp_wrd in split_array:
          if not ( pattern_digit.search(tmp_wrd) and re.search(r'[a-zA-Z]', tmp_wrd)) :
             list_out.append( tmp_wrd )
          #
        description = " ".join(list_out)
    return description.upper().strip()


def createCorrectedTaxDistances(maxpathlength=70):
    """
    Creates and return a list of corrected taxonomic distances.
    Values for this list is calculated from a regression curve created by Petri T&ouml;r&ouml;nen
    Curve is shifted upwards to ensure nonnegativity and the tail is set to zero.
    """
    weights = numpy.array([-0.0704054300, 0.0006559828, 0.4615545000])
    tmp=[]
    for i in range(0,maxpathlength):
        x=sum(numpy.array([i, i**2, numpy.log(i+1)]) * weights)
        tmp.append(x)
    tmp_min = min(tmp)
    for i in range(0,maxpathlength):
        tmp[i]=tmp[i]-tmp_min
    min_index = tmp.index(min(tmp))
    for i in range(min_index, len(tmp)):
        tmp[i] = 0.0
    return tmp

def logmod(x):
        """
        Petri Toronen's squashing function

                        1+log(x), if x>1
        logmod(x)=      x, if -1<=x<=1
                        -1-log(-x), if x<-1
        """
        if x>1: return(1+numpy.log(x))
        if x<-1: return(-1-numpy.log(-x))
        return(x)

def sampleStats(block,status_col,rm1_col):
        """
        Calculates size, mean and average over a column of data matrix needed to calculate GSZ.

        Inputs:
                block = array of data arrays
                status_col = index of status column; used to filter data
                rm1_col = index of RM1 column; used as score in GSZ

        Outputs:
                sampleSize = size of filtered hitlist
                sampleScoreMean = mean of RM1 in filtered hitlist
                sampleScoreVariance = variance of RM1 in filtered hitlist
        """
        sampleSize=0
        summa=0.0
        sumsq=0.0
        for row in block:
                if row[status_col] == 'False': continue
                sampleSize+=1
                try:
                        rm1=float(row[rm1_col])
                except:
                        rm1=0.0
                summa+=rm1
                sumsq+=rm1*rm1
        # GSZ sample is the whole hitlist
        if sampleSize<1: return(0,0.0,0.0) # can't divide by zero
        sampleScoreMean=summa/sampleSize
        sampleScoreVariance=(sumsq-summa*summa/sampleSize)/sampleSize

        return(sampleSize,sampleScoreMean,sampleScoreVariance)

#############################
# Special function (class)  #
# FormFactor.formfactor     #
#############################

class FormFactor:
        """Informative description penalized for black-listed terms and non-alphabetic characters
        """

        def __init__(self):
                self.blacklist={'UNCHARACTERI': 60,'WGS PROJECT': 44, '_': 60, 'GENOMIC': 60,
                                'WHOLE GENOME SHOTGUN': 80,  'CONTIG': 60, 'SCAFFOLD': 60,
                                 'PARTIAL': 14, 'HYPOTHETICAL': 24, 'CDNA': 16, 'WGS': 16,
                                 'UNKNOWN': 14,'-LIKE': 10, 'ORF': 6, 'COMPLETE ': 60,
                                'PREDICTED PROTEIN': 60, 'LOCUS': 60, 'SEQUENCE': 60,
                                'HOMOLOG': 3, 'ISOFORM': 4,
                                'PUTATIVE': 12, 'SIMILAR TO': 10,
                                'FRAGMENT': 8, 'PROBABLE': 8,
                                'CHROMOSOME': 40, 'POSSIBLE': 8,
                                'CONSERVED': 9, 'FAMILY': 6, 'RELATED': 10}
                self.badpattern=re.compile(r'[\w\d]{3}\-*\d')
                self.bad1=re.compile(r'[A-Z]\d')
                self.bad2=re.compile(r'\d[A-Z]')
                self.bad3=re.compile(r'^PROTEIN')
                self.bad4=re.compile(r'PROTEIN$')
                self.goodpattern=re.compile(r'\w\-[\d\,\'\)\(]+\-')

        def formfactor(self,desc):
                if not desc: return(0.0)
                ntotal=len(desc)
                if ntotal<1: return(0.0)
                undesired=3.0 # filter gene name; pseudocount to disfavour single-word descriptions
                test=desc.upper()
                # undesirables = non-alphabetic characters
                for c in test:
                        if c == ' ' or c == '-' or c == "'":
                                pass # no penalty for spaces in multiword description
                        elif c<'A' or c>'Z':
                                undesired+=2.0
                # mask letter-number combinations
                for word in test.split():
                        if self.bad1.search(word): undesired+=4
                        if self.bad2.search(word): undesired+=4
                if self.bad3.search(test): undesired+=12
                if self.bad4.search(test): undesired+=12

                # undesirables = black-listed terms
                for term in self.blacklist:
                        if term in test: undesired+=self.blacklist[term]
                # calculate form factor
                ff=1.0-undesired/ntotal
                #print "' Form Factor of %s is %f" %(desc,ff)
                return ff

############################
# PPV estimation functions #
############################

def DE_PPV_euk(x):
        "argument is float(RM2)"
        if x < 0.705: return(max(0.0, -1.253 + 2.295*x ))
        if x < 3.05855: return(min(1.0,-0.05199 + 0.66827*x - 0.10945*x*x))
        return(0.9680758)

def DE_PPV_bac(x):
        "argument is float(RM2)"
        if x < 1.322: return(max(0.0, -0.5057 + 1.4713*x - 0.3757*x*x))
        return(min(1.0, 0.77258 - 0.01382*x + 0.01632*x*x))

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
        return(max(0.0,0.97768+0.04106*z))

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
        return(max(0.0,0.2557+0.483*z))

