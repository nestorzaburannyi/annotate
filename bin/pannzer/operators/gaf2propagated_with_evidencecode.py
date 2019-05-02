# cut -f 2,5,7 goa_uniprot_noIEA.gaf

from myoperator import BlockOperator
from Read_and_Print import read_dict_data
import sys,re,math

class gaf2propagated_with_evidencecode(BlockOperator):
    """
Input: qpid, goid (multiline block)
Output: qpid, propagated_goids (one line per qpid)

IDA     1           350,497     inferred from direct assay
IMP     1           223,662     inferred from mutant phenotype
IGI     1            41,964     inferred from genetic interaction
IPI     1           207,391     inferred from physical interaction
IEP     1            24,128     inferred from expression pattern
TAS     0.9         146,953     traceable author statement
NAS     0.8          18,123     non-traceable author statement
IC      0.9           7,919     inferred by curator
ISS     0.8         280,161     inferred from sequence or structural similarity
RCA     0.8           4,752     inferred from reviewd computational analysis
IEA     0.7     405,214,644     inferred from electonic annotation
ND      0.5         226,319     no biological data available
NR      0.5               0     not recorded

IBA     0.8       1,464,064     inferred from biological aspect of ancestor
IGC     0.7             444     inferred from genomic context
IRD     0.7              33     inferred from rapid divergence
IKR     0.8             273     inferred from key residues
ISM     0.8           5,370     inferred from sequence model
EXP     1             3,719     inferred from experiment
ISO     0.8         103,925     inferred from sequence orthology
ISA     0.8           9,602     inferred from sequence alignment

High-throughput equivalents of I?? codes IMP, IGI, IDA, IEP, IPI
 HMP
 HGI
 HDA
 HEP
 HPI
    """
    def __init__(self,glob):
        # set parameters
        self.glob=glob
        # propagate assigned GOids to parents
        obotabfile='obo.tab'
        self.mapped_to=read_dict_data(obotabfile,'goid','mapped_to',header=True,UseColNames=True) # goid, mapped_to
        self.directparents=read_dict_data(obotabfile,'goid','directParentSet',header=True,UseColNames=True) # goid,parentset
        # write blockwise summmary
        [self.data,self.godict_data]=glob.use_sheets(["data","godict"])
        # column indices in GAF file
        [self.acc_col1,self.goid_col1,self.ec_col1]=self.data.use_columns(["qpid","goid","evidence_code"])
        self.godict_data.use_columns(["qpid","propagated","evidence_weights"])
        self.evidence_weight={'IDA': 1.0, 'IMP': 1.0, 'IGI': 1.0, 'IPI': 1.0, 'IEP': 1.0, 'TAS': 0.9, 'NAS': 0.8, 'IC': 0.9, 'ISS': 0.8, 'RCA': 0.8, 'IEA': 0.7, 'ND': 0.5, 'NR': 0.5, 'IBA': 0.8, 'IGC': 0.7, 'IRD': 0.7, 'IKR': 0.8, 'ISM': 0.8, 'EXP': 1.0, 'ISO': 0.8, 'ISA': 0.8, 'HDA': 1.0, 'HMP': 1.0, 'HGI': 1.0, 'HPI': 1.0, 'HEP': 1.0 }

    def process(self,block):
        try:
                acc=block[0][self.acc_col1]
        except:
                return
        if acc[0]=='!': return # gaf comment line
        goidlist=[]
        x={}
        GOweight={} # initialize for each block
        for row in block:
                goid=row[self.goid_col1]
                goid=re.sub(r'GO:','',goid)
                try:
                        w=self.evidence_weight[row[self.ec_col1]]
                except: # don't trust unknown evidence codes
                        w=0.5
                        sys.stderr.write("# Unknown evidence code: %s\n" %row[self.ec_col1])
                GOweight[goid]=w
                # direct assignment + parents
                x[goid]=1
        # propagate to all parents, max evidence weight
        stack=x.keys()
        parents={}
        while len(stack)>0:
            goid=stack.pop()
            w=GOweight[goid]
            parents[goid]=1
            if not goid in self.directparents: continue
            tmp=self.directparents[goid]
            if tmp=='': continue # root
            for p in tmp.split(','):
                    if not p in parents: stack.append(p)
                    # propagate maximum evidence weight
                    if not p in GOweight:
                        GOweight[p]=w
                    elif w > GOweight[p]:
                        GOweight[p]=w
        propagated=parents.keys()
        propagated.sort()
        # created matching list of evidenceweights
        ecw=[]
        for goid in propagated: ecw.append(str(GOweight[goid])) #ecw+=str(GOweight[goid])+','
        #for goid in propagated: ecw+=str(GOweight[goid])+','
        # write to godict spreadsheet
        self.godict_data.append_row([acc,','.join(propagated),','.join(ecw)])

