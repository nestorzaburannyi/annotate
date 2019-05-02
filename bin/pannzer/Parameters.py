from __future__ import absolute_import
from __future__ import print_function
import sys, os
import argparse
from Read_and_Print import read_dict_data, read_dict_counts, read_dict_GOdict #, read_dict_PHR
import SpreadSheet
import myoperator
from config import BugsyConfig
import socket

############################################################################
# GLOBALS-CLASS : Data, dictionaries and parameters visible to all classes #
############################################################################
class WorkSpace:
        """
        This collects all the variables that need to be visible to
        other objects.

        If you need a new global variable, add it as self.param['new_global_variable']
        """
        def __init__(self, configFile = None, verbose = False):
          # operators
          pathname = os.path.dirname(sys.argv[0])
          fullpath = os.path.abspath(pathname)
          operator_dir =  fullpath + os.sep + 'operators'
          myoperator.initialise_operators(operator_dir)
          # set parameter defaults, parse command-line arguments
          self.get_parameters(configFile=configFile)
          if verbose: self.print_parameters()
          self.check_parameter_values()
          # block-specific data --> overwrite from block if defined there
          self.QUERYSPECIES=self.param['input_QUERYSPECIES']
          self.QUERYSEQUENCE=None
          self.QUERYDESC=None
          self.QUERYKINGDOM='bac'
          # manage spreadsheets
          self.sheet_names={} # key=name, value=object handle
          self.nsheet=0
          self.sheets=[]
          self.use_sheet("data") # sheet[0] in Runner
          # dictionaries
          self.dictlist=[]
          self.PHR=None         # key = accession number, value = description
          self.lineage=None     # key=species, value=lienage; used in Taxonomy
          self.taxid=None       # key=species, value=taxid
          self.wordcounts=None  # key=word, value=count; used in TFIDF
          self.desccounts=None  # key=DE, value=count; used in Cluster_GSZ
          self.GOcounts=None    # key=goid, value=count; used in GOrimpsu
          self.GOdict=None      # key=accession_number, value=list of goids including parents; used in GOrimpsu, GO_evaluation
          self.GOdict_weights=None # GO-weights for Blast2GO emulator
          self.GOdict_noIEA=None # key=accession number, value=list of goids includeing parents; used in GO_Evaluation
          self.GOparents=None   # key=goid, value=csv parent terms; used in  GO_evaluation
          self.godesc=None      # key=goid, value=description; used in RM3
          self.ontology=None    # key=goid, value=BP|MF|CC|null; used in RM3

        def use_operators(self,list,operator_dir='operators'):
          "Return array of operator object handles based on name"
          handles=[]
          myoperator.initialise_operators(operator_dir)
          for name in list:
                Op=myoperator.get_operator(name)
                x=Op(self) # instantiated as __init__(glob)
                handles.append(x)
          return(handles)

        def use_sheets(self,list):
          "Return an array of SpreadSheet object handles based on name"
          handles=[]
          for name in list:
                handles.append(self.use_sheet(name))
          return(handles)

        def use_sheet(self,name):
          "Return a SpreadSheet object handle based on name"
          if not name in self.sheet_names:
                        # create new SpreadSheet object
                        self.sheets.append(SpreadSheet.SpreadSheet())
                        self.sheet_names[name]=self.sheets[self.nsheet]
                        self.nsheet+=1
          return(self.sheet_names[name])

        def use_online_dictionaries(self,list):
          "add to dictlist for lazyRunner"
          for x in list:
                if not x in self.dictlist:
                        self.dictlist.append(x)
                        sys.stderr.write("# Added online dictionary %s\n" %x)

        def load_lineage(self,fn):
          if not self.lineage: # load dictionary
                sys.stderr.write("# Loading taxonomy %s\n" %fn)
                self.lineage=read_dict_data(fn,'Scientific name','Lineage',True,True,"\t",False,True)
                sys.stderr.write("# %i species in taxonomy\n" %len(self.lineage))
          return(len(self.lineage))

        def load_taxid(self,fn):
          if not self.taxid: # load dictionary
                sys.stderr.write("# Loading taxid %s\n" %fn)
                self.taxid=read_dict_data(fn,'Scientific name','Taxon',True,True,"\t",False,True)
                sys.stderr.write("# %i species in taxonomy\n" %len(self.lineage))
          return(len(self.taxid))

        def load_wordcounts(self,fn):
          if not self.wordcounts: # load dictionary
                sys.stderr.write("# Loading word count %s\n" %fn)
                self.wordcounts=read_dict_counts(fn,1,0)
                sys.stderr.write("# nwordtotal = %i\n" %len(self.wordcounts))
          return(len(self.wordcounts))

        def load_desccounts(self,fn):
          if not self.desccounts: # load dictionary
                sys.stderr.write("# Loading description count %s\n" %fn)
                self.desccounts=read_dict_counts(fn,1,0)
          return(len(self.desccounts))

#        def load_PHR(self,fn):
#          if not self.PHR: # load dictionary
#                sys.stderr.write("# Loading PHR dictionary %s\n" %fn)
#                self.PHR=read_dict_PHR(fn)
#                sys.stderr.write("# %i keys in PHR\n" %len(self.PHR))
#          return(len(self.PHR))

        def load_GOdict(self,fn):
          if not self.GOdict: # load dictionary
                sys.stderr.write("# Loading GO dictionary %s\n" %fn)
                #self.GOdict=read_dict_GOdict(fn)
                self.GOdict=read_dict_data(fn,0,1,False,False,"\t",True,False)
                self.GOdict_weights=read_dict_data(fn,0,2,False,False,"\t",True,False)
                sys.stderr.write("# %i keys in GOdict\n" %len(self.GOdict))
          return(len(self.GOdict))

        def load_GOdict_noIEA(self,fn):
          if not self.GOdict_noIEA: # load dictionary
                sys.stderr.write("# Loading GO dictionary %s\n" %fn)
                self.GOdict_noIEA=read_dict_data(fn,0,1,False,False,"\t",True,False)
                sys.stderr.write("# %i keys in GOdict_noIEA\n" %len(self.GOdict_noIEA))
          return(len(self.GOdict_noIEA))

        def load_goidelic(self,fn):
                "load all small dictionaries that use goid as key: GOparents, godesc, ontology, GOcounts, EC, KEGG, IC"
                sys.stderr.write("# Loading GO-id dictionaries from %s\n" %fn)
                fh=open(fn,'r')
                lines=fh.readlines()
                fh.close()
                self.load_goidelic_data(lines)

        def load_goidelic_data(self,lines):
                self.GOcounts={}
                self.ontology={}
                self.godesc={}
                self.GOparents={}
                self.rootcount={}
                self.EC={}
                self.KEGG={}
                self.IC={}
                for line in lines:
                    try:
                        (cnt,ont,goid,desc,parlist,ec,kegg,ic)=line.split("\t")
                        self.GOcounts[goid]=cnt
                        self.ontology[goid]=ont
                        self.godesc[goid]=desc
                        if parlist != '': self.GOparents[goid]=parlist.rstrip().split(",")
                        if ec != '': self.EC[goid]=ec
                        if kegg != '': self.KEGG[goid]=kegg
                        self.IC[goid]=ic.rstrip() # last field on line
                    except:
                        pass
                        #print >> sys.stderr, "# ERROR reading goidelic data: ",line
                # size of ontology-assigned database
                for goid in self.GOcounts.keys():
                        try:
                                cnt=int(self.GOcounts[goid])
                                ontology=self.ontology[goid]
                                if not ontology in self.rootcount: self.rootcount[ontology]=0
                                if cnt>self.rootcount[ontology]: self.rootcount[ontology]=cnt
                        except:
                                #print >> sys.stderr, "# ERROR in goidelic:",goid
                                continue # header gives ValueError
                for ontology in self.rootcount: sys.stderr.write("# ontology %s count %i\n" %(ontology,self.rootcount[ontology]))

        def load_nprot(self,fn):
                fh=open(fn,'r')
                x=fh.readline()
                fh.close()
                return(int(x))

        def load_nwordtotal(self,fn):
                sys.stdout.write("# this is load_nwordtotal %s\n" %(fn))
                fh=open(fn,'r')
                x=fh.readline()
                fh.close()
                sys.stdout.write("# got nwordtotal %s\n" %x)
                i=int(x)
                return(i)

        def get_parameters(self,configFile=None):
                """Set default values of global parameters.
                Convention: for compatibility with ConfigParser, parameter name is of the form section_param, where section is input|DATA|PANZ|CONN.
                Convention: all parameter names must be uppercase, because ConfigParser returns lowercased names!
                This sets default values to all parameters that can be read from the command line.
                Parameters can be read from config file or written to config file.
                """
                # get command-line arguments
                parser=argparse.ArgumentParser(description="Parse arguments from command line.")
                parser.add_argument('-i','--input_FILE',help='input file name',default='--',required=False)
                parser.add_argument('-o','--input_OUTFILES',help='comma-separated list of output file names',default="--",required=False)
                parser.add_argument('-f','--input_FORMAT',help='FASTA|tab. Input format can be FASTA-formatted sequences or tabular data',default="FASTA",required=False)
                parser.add_argument('-c','--input_COLNAMES',help='column names (string)',default=None,required=False)
                parser.add_argument('-b','--input_BLOCK_COLUMN_NAME',help='column name to identify query block',default='qpid',required=False)
                parser.add_argument('-s','--input_QUERYSPECIES',help='query species for Pannzer',default="auto",required=False)
                parser.add_argument('-a','--input_PARAM_FILE_NAME',help='Configuration file',default=configFile,required=False)
                parser.add_argument('-k','--input_CHUNK',type=int,help='integer. Chunk size of data stream (to make fewer client-server calls)',default=100,required=False)
                parser.add_argument('--input_GODICT_noIEA',help='file with table of non-IEA GOA assignments',default=None,required=False) # must read file, not online
                parser.add_argument('-m','--input_OPERATOR',help='name of operator: Pannzer|BestInformativeHit|SANS|taxonassignment|FF|gaf2propagated|GOevaluation|etc.',default="Pannzer",required=False)
                parser.add_argument('--input_LIVEDATA',help='Progress meter',default=None,required=False)
                # GO evaluation
                parser.add_argument('--eval_OBOTAB',help='GO hierachy for propagation',default='obo.tab',required=False)
                parser.add_argument('--eval_TRUTH',help='propagated reference of truth in GO evaluation',default="goa_truth_propagated",required=False)
                parser.add_argument('--eval_SCOREFUNCTIONS',help='score functions to be evaluated',default="RM3_PPV ARGOT_PPV JAC_PPV HYGE_PPV SANS_PPV",required=False)
                # data files: DATADIR/dictfile
                parser.add_argument('-d','--DATA_DIR',help='directory where Uniprot + GO data reside',default="/data/uniprot",required=False)
                # Pannzer parameters
                parser.add_argument('--PANZ_CLUSTERING_CUTOFF',type=float,help='float. Cosine distance cutoff of DE clustering in Pannzer',default=0.7,required=False)
                parser.add_argument('--PANZ_FILTER_PERMISSIVE',action='store_true',help='Flag for permissive coverage filtering in Pannzer (True: qcov or scov; False: qcov and scov)',default=False,required=False)
                parser.add_argument('--PANZ_MINLALI',type=int,help='integer. Minimum alignment length if permissive coverage filtering is True',default=100,required=False)
                parser.add_argument('--PANZ_QCOVCUTOFF',type=float,help='float. Query coverage cutoff in Pannzer',default=0.6,required=False)
                parser.add_argument('--PANZ_SCOVCUTOFF',type=float,help='float. Sbjct coverage cutoff in Pannzer',default=0.6,required=False)
                parser.add_argument('--PANZ_MINPIDECUTOFF',type=float,help='float. Minimum sequence identity in Pannzer',default=0.4,required=False)
                parser.add_argument('--PANZ_MAXPIDECUTOFF',type=float,help='float. Maximum sequence identity in Pannzer',default=1.0,required=False)
                parser.add_argument('--PANZ_FILTER_OUTPUT',action='store_true',help='remove redundant GO terms',default=False,required=False)
                parser.add_argument('--PANZ_BESTCLUSTER',action='store_true',help='output only one DE prediction per query',default=False,required=False)
                parser.add_argument('--PANZ_MAXHITS',type=int,help='integer. Maximum number of sequence hits in Pannzer',default=100,required=False)
                parser.add_argument('--PANZ_FFCUTOFF',type=float,help='float. Form Factor cutoff for informative/non-informative hits',default=0.2,required=False)
                parser.add_argument('--PANZ_REMOVE_ABBR',action='store_true',help='remove abbreviations in Cleandesc function',default=False,required=False)
                parser.add_argument('--PANZ_PREDICTOR',help='CSV list of Pannzer scoring functions (DE,RM3,ARGOT,JAC,HYGE)',default='DE,ARGOT',required=False)
                # client-server parameters
                lhost=socket.gethostname()
                parser.add_argument('-S','--CONN_SANSPORT',type=int,help='integer. Port number of uniprot SANS-server',default=54321,required=False)
                parser.add_argument('-T','--CONN_SANSHOST',help='SANS-server host',default='localhost',required=False)
                parser.add_argument('-H','--CONN_HOSTNAME',help='DictServer host',default=lhost,required=False)
                parser.add_argument('-P','--CONN_PORTNO',type=int,help='integer. Port number of DictServer',default=50002,required=False)
                parser.add_argument('-R','--CONN_REMOTE',action='store_true',help='Flag for HTTP requests to access remote online servers',default=False,required=False)
                # BLAST2GO parameter
                parser.add_argument('-B2G','--B2G_THRESH',type=float,help='Threshold for BLAST2GO for accepting predictions.',default=55)
                # SANS parameters
                parser.add_argument('-u','--SANS_SSEQ',action='store_true',help='output sbjct sequence from SANS',default=False,required=False)
                parser.add_argument('-v','--SANS_RANGES',action='store_true',help='output alignment ranges from SANS',default=False,required=False)
                parser.add_argument('--SANS_H',type=int,help='integer. Number of hits from SANS',default=100,required=False)
                parser.add_argument('--SANS_HX',type=int,help='integer. Width of suffix array neighborhood in SANS',default=200,required=False)
                parser.add_argument('--SANS_R',type=int,help='integer. Number of rejects in SANS',default=100,required=False)
                parser.add_argument('--SANS_VOTELIST_SIZE',type=int,help='integer. Size of vote list in SANS',default=200,required=False)
                parser.add_argument('--SANS_PROTOCOL',type=int,help='integer. SANS protocol',default=1,required=False)
                # taxonassignment
                parser.add_argument('--TAXON_DEPTH',help='Number of taxonomic levels output by taxonassigment',default=2,required=False)
                # args is a namespace; self.param is a dictionary
                args = parser.parse_args()
                self.param=vars(args)
                # non-command-line defaults:
                # - DictServer data files hardcoded in update script
                self.param['DATA_GODICT']='godict.txt' # 'godict.noIEA.txt'
                self.param['DATA_GOIDELIC']='mergeGO.out' # 'mergeGO.noIEA.out'
                self.param['DATA_NPROT']='nprot'
                self.param['DATA_NWORDTOTAL']='nwordtotal'
                self.param['DATA_WORDCOUNTS']='uniprot.word.uc.counts'
                self.param['DATA_DESCCOUNTS']='uniprot.desc.uc.counts'
                self.param['DATA_TAXONOMY']='taxonomy-all.tab'
                self.param['DATA_PHR']='uniprot.phr'
                # - constants
                self.param['PANZ_JACCARD_MINCOUNT']=10
                self.param['HYGE_MAXCACHEKEYS']=10000
                # override defaults if there is a config file
                cf=self.param['input_PARAM_FILE_NAME']
                self.readConfigFile(cf)

        def readConfigFile(self,filename):
                """Read config file and fill self.param dictionary"""
                if filename is None: return
                sys.stderr.write('# readConfigFile %s\n' %filename)
                bc=BugsyConfig(filename)
                for section in bc.config.sections():
                        tmp=bc.config.items(section) # returns array of tuples [('coffee',1),('tea',2)]
                        for (name,value) in tmp:
                                key="%s_%s" %(section,name.upper())
                                value=bc.get(section,name) # value is typed
                                self.param[key]=value
                                sys.stderr.write('# set %s = %s\n' %(key,str(value)))

        def writeConfigFile(self,filename):
                """Write current parameter configuration to a ConfigParser config file.
                Sections are mandatory in ConfigParser. Sections are split from parameter name as in
                PANZ_MAXHITS -> section = PANZ, parameter name = MAXHITS, value = self.param['PANZ_MAXHITS']
                """
                # open file
                fh=open(filename,'w')
                # sort self.param alphabetically
                tmp=self.param.keys()
                tmp.sort()
                old_section=''
                for key in tmp:
                # split section,parametername from self.param key
                        (section,name)=key.split('_',1)
                        # write new section header
                        if section != old_section:
                                fh.write("\n[%s]\n\n" %section)
                                old_section=section
                        # write value
                        fh.write("%s = %s\n" %(name,str(self.param[key])))
                fh.close()

        def print_parameters(self):
                tmp=self.param.keys()
                tmp.sort()
                for key in tmp:
                        print("# %s:" %key, self.param[key], file=sys.stderr)

        def check_parameter_values(self):
                PANZ_CLUSTERING_CUTOFF=self.param['PANZ_CLUSTERING_CUTOFF']
                PANZ_QCOVCUTOFF=self.param['PANZ_QCOVCUTOFF']
                PANZ_SCOVCUTOFF=self.param['PANZ_SCOVCUTOFF']
                PANZ_MINPIDECUTOFF=self.param['PANZ_MINPIDECUTOFF']
                PANZ_MAXPIDECUTOFF=self.param['PANZ_MAXPIDECUTOFF']
                # value must be between zero and one
                ok=True
                if PANZ_CLUSTERING_CUTOFF<0 or PANZ_CLUSTERING_CUTOFF>1:
                        print("# ERROR: CLUSTERING_CUTOFF out of range [0,1]",PANZ_CLUSTERING_CUTOFF)
                        ok=False
                if PANZ_QCOVCUTOFF<0 or PANZ_QCOVCUTOFF>1:
                        print("# ERROR: QCOVCUTOFF out of range [0,1]",PANZ_QCOVCUTOFF)
                        ok=False
                if PANZ_SCOVCUTOFF<0 or PANZ_SCOVCUTOFF>1:
                        print("# ERROR: SCOVCUTOFF out of rang [0,1]",PANZ_SCOVCUTOFF)
                        ok=False
                if PANZ_MINPIDECUTOFF<0 or PANZ_MINPIDECUTOFF>1:
                        print("# ERROR: MINPIDECUTOFF out of range [0,1]",PANZ_MINPIDECUTOFF)
                        ok=False
                if PANZ_MAXPIDECUTOFF<0 or PANZ_MAXPIDECUTOFF>1:
                        print("# ERROR: MAXPIDECUTOFF out of range [0,1]",PANZ_MAXPIDECUTOFF)
                        ok=False
                if PANZ_MINPIDECUTOFF>PANZ_MAXPIDECUTOFF:
                        print("# ERROR: MINPIDECUTOFF>MAXPIDECUTOFF",PANZ_MINPIDECUTOFF, PANZ_MAXPIDECUTOFF)
                        ok=False
                if ok: return
                sys.exit()

