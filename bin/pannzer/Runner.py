from __future__ import absolute_import
from __future__ import print_function
import sys
import requests
import myoperator
import Parameters, DictServer, XMLParser
from PannzerFunctions import Cleaner

class Runner:
        def __init__(self, glob, operator_name=None, CHUNK=100, liveData=None, MAXRES=64000, PACKETSIZE=100000):
                """
glob = object handle containing spreadsheets,dictionaries,parameters
CHUNK = number of entries (query sequences, blocks) in buffer. Buffer is used for lazy dictionary lookup
liveData = name of status file (number of processed queries). None implies no output
                """
                self.pythonversion=int(sys.version_info[0])
                print("# Python version: ",self.pythonversion, file=sys.stderr)
                self.glob=glob
                self.MAXRES=MAXRES # maximum size of SANSparallel query string
                self.liveData=liveData
                self.sentquery=0
                self.CHUNK=CHUNK
                self.PACKETSIZE=PACKETSIZE # guard against long FASTA sequences
                self.colnames=None
                self.block_column_index=0
                # initialize operator
                self.operator_name=operator_name
                if not self.operator_name: self.operator_name=self.glob.param['input_OPERATOR']
                [self.myoperator]=self.glob.use_operators([operator_name])
                self.do_lazy=(len(self.glob.dictlist)>1 or (len(self.glob.dictlist)==1 and self.glob.dictlist[0]!='GOIDELIC')) # glob.dictlist defined in operators' __init__
                self.linewise=False
                if isinstance(self.myoperator, myoperator.RowOperator):
                        self.blockwise=False
                elif isinstance(self.myoperator, myoperator.BlockOperator):
                        self.blockwise=True
                elif isinstance(self.myoperator, myoperator.TextOperator):
                        self.linewise=True
                else:
                        sys.stderr.write("# Invalid operator %s. Exiting!\n" %operator_name)
                        sys.exit()
                # parameter nicknames
                self.REMOTE=self.glob.param['CONN_REMOTE']
                self.PORTNO=self.glob.param['CONN_PORTNO']
                self.H=self.glob.param['SANS_H']
                self.HX=self.glob.param['SANS_HX']
                self.R=self.glob.param['SANS_R']
                self.VOTELIST_SIZE=self.glob.param['SANS_VOTELIST_SIZE']
                self.SANSPROTOCOL=self.glob.param['SANS_PROTOCOL']
                self.SSEQ=self.glob.param['SANS_SSEQ']
                self.RANGES=self.glob.param['SANS_RANGES']
                self.SANSHOST=self.glob.param['CONN_SANSHOST']
                self.SANSPORT=self.glob.param['CONN_SANSPORT']
                self.DICTHOST=self.glob.param['CONN_HOSTNAME']
                self.DICTPORT=self.glob.param['CONN_PORTNO']

        def lazyRunner(self, infile, output_files=['--'], input_format="FASTA", colnames=None, queryspecies="auto", block_column_index=0, block_column_name=None):
                """
Main function to process data streams with SANSPANZ operators.

infile = input file name; None implies STDIN
output_files = output file names for each spreadsheet [data, DE predictions, GO predictions, annotations]; None implies no output
input_format = format of input file: "FASTA" = sequences in FASTA format, "tab" = tab-separated values
colnames = column names of input data stream. None implies column header is first non-comment line of input; automatic if input_format is FASTA
queryspecies = "auto" implies parsing from Uniprot entry's OS= tag, None implies input data has isquery/species columns, must be supplied if input_format is "FASTA"
                """
                self.have_colnames=False
                if colnames:
                        self.colnames=colnames.split()
                        self.colmap=self.glob.sheets[0].use_columns(self.colnames)
                        self.have_colnames=True
                        print("# Received colnames:",self.colnames,self.colmap, file=sys.stderr)
                self.blocking_initialized=(block_column_name is None) # block_column_index is overwritten in first process_chunk if block_column_name is defined
                self.block_column_name=block_column_name
                # load small dictionaries: GOIDELIC, including rootcount
                if "GOIDELIC" in self.glob.dictlist:
                        if self.REMOTE:
                                tmp=DictServer.DICTquery("GOIDELIC",self.pythonversion,REMOTE=self.REMOTE,HOSTNAME=self.DICTHOST,PORTNO=self.DICTPORT).split("\n")
                                self.glob.load_goidelic_data(tmp)
                        else:
                                fn=self.glob.param['DATA_DIR']+'/'+self.glob.param['DATA_GOIDELIC']
                                self.glob.load_goidelic(fn)
                        print('# loaded GOIDELIC remote=',self.REMOTE,len(self.glob.GOcounts),len(self.glob.GOparents),len(self.glob.ontology),len(self.glob.godesc),len(self.glob.EC),len(self.glob.
KEGG), file=sys.stderr)
                # chunk header is known if FASTA or colnames argument is defined; otherwise it must be captured from input data stream
                if input_format=="FASTA" and not self.linewise:
                        # get colnames from XMLparser with empty message
                        self.xml=XMLParser.XMLParser(queryspecies)
                        self.colnames=self.xml.stream("",True,False,self.SSEQ,self.RANGES)[0].rstrip().split("\t")
                        self.colmap=self.glob.sheets[0].use_columns(self.colnames)
                        self.have_colnames=True
                # open IO channels
                data_in=self.open_IO_channels(infile,output_files)
                for sheet in self.glob.sheets: sheet.output(header=True)
                # fill buffer
                self.olduid='?'
                self.iquery=0
                packet=''
                while True:
                        line=data_in.readline()
                        if not line: break
                        if line[0]=='#': continue # ignore comment lines
                         # shortcut: process line by line
                        if self.linewise:
                                self.myoperator.process(line)
                                continue
                        # CHUNK counter
                        if self.test_newentry(input_format,line):
                                if (self.iquery % self.CHUNK == 0) or (input_format=="FASTA" and not self.linewise and len(packet)>self.PACKETSIZE):
                                #if self.iquery % self.CHUNK == 0:
                                        self.process_chunk(input_format,packet)
                                        packet=''
                                self.iquery+=1
                        # append line to packet
                        packet+=line
                # last buffer
                self.process_chunk(input_format,packet) # linewise has empty packet, returns immediately
                self.myoperator.finalise()
                # close IO channels
                self.close_IO_channels()

        def process_chunk(self,input_format,packet):
                if len(packet)<1: return
                # split to lines (removes \n)
                if input_format=="FASTA":
                        lines=self.SANSquery(packet)[0].split("\n")
                elif input_format=="tab":
                        lines=packet.split("\n")
                else:
                        sys.stderr.write("ERROR: unknown input_format %s\n" %input_format)
                        return
                # add input colnames if not defined in operator
                if not self.have_colnames:
                        self.colnames=lines[0].split("\t") # comments were removed on data input
                        self.colmap=self.glob.sheets[0].use_columns(self.colnames)
                        startrow=1 # data starts after header row
                        self.have_colnames=True
                else:
                        startrow=0 # no header row
                if (not self.blocking_initialized) and self.block_column_name: # overwrite self.block_column_index
                        print("find %s in " %self.block_column_name,self.colnames, file=sys.stderr)
                        self.block_column_index=self.colnames.index(self.block_column_name)
                        self.blocking_initialized=True
                        sys.stderr.write("# block_column_index = %i name = %s\n" %(self.block_column_index,self.block_column_name))
                # capture dictionary keys from input data
                if self.do_lazy: self.load_private_online_dictionaries(lines[startrow:])
                # send one block or chunk-of-rows to operator
                for sheet in self.glob.sheets: sheet.empty_block()
                olduid='?'
                for line in lines[startrow:]:
                        if not line: continue
                        row=line.split("\t")
                        if self.blockwise:
                                if row[self.block_column_index] != olduid:
                                        olduid=row[self.block_column_index]
                                        self.myoperator.process(self.glob.sheets[0].block)
                                        # output result
                                        for sheet in self.glob.sheets:
                                                sheet.output(result=True)
                                                sheet.empty_block()
                        datarow=[]
                        self.glob.sheets[0].append_row(datarow) # empty, n.d.
                        datarow=self.glob.sheets[0].block[-1]
                        try:
                                for i in range(0,len(self.colmap)):
                                        ix=self.colmap[i]
                                        datarow[ix]=row[i]
                        except:
                                sys.stderr.write("# Input error: %s\n" %line)
                                self.glob.sheets[0].block.pop() # remove last entry
                                continue
                # last block
                if self.blockwise:
                        self.myoperator.process(self.glob.sheets[0].block)
                else:
                        for row in self.glob.sheets[0].block: self.myoperator.process(row)
                # output result
                for sheet in self.glob.sheets: sheet.output(result=True)
                # update liveData
                if self.liveData:
                        fh=open(self.liveData,"w")
                        fh.write("%i\n" %self.iquery)
                        fh.close()

        def load_private_online_dictionaries(self,lines):
                  uacc=[]
                  uspecies=[]
                  udesc=[]
                  uword=[]
                  if "GODICT" in self.glob.dictlist:
                        try:
                                spid_col=self.colnames.index('spid')
                        except:
                                spid_col=self.colnames.index('qpid')
                        uspid=self.catch_unique(spid_col,lines)
                        # grab accession: spid is db|accession|pid
                        x={}
                        for spid in uspid:
                                try:
                                        acc=spid.split("|")[1]
                                        x[acc]=1
                                except:
                                        sys.stderr.write("# Warning: no accession number from %s\n" %spid)
                        uacc=x.keys()
                  if "LINEAGE" in self.glob.dictlist or "TAXID" in self.glob.dictlist:
                        species_col=self.colnames.index('species') # use spreadsheet's colnames
                        uspecies=self.catch_unique(species_col,lines)
                  if "DESCCOUNT" in self.glob.dictlist or  "WORDCOUNT" in self.glob.dictlist:
                        desc_col=self.colnames.index('desc')
                        x={}
                        tmp=self.catch_unique(desc_col,lines)
                        # call Cleaner function
                        for desc in tmp: x[Cleaner(desc.upper())]=1 ## operators.Cleaner !!!
                        udesc=x.keys()
                        if "WORDCOUNT" in self.glob.dictlist:
                                x={}
                                for desc in udesc:
                                        for word in desc.split(): x[word.upper()]=1
                                uword=x.keys()
                  # compose DictServer message
                  msg=""
                  for key in uspecies: msg+="LINEAGE"+"\t"+key.upper()+"\n"
                  for key in uspecies: msg+="TAXID"+"\t"+key.upper()+"\n"
                  for key in uacc: msg+="GODICT"+"\t"+key.upper()+"\n"
                  for key in uword: msg+="WORDCOUNT"+"\t"+key.upper()+"\n"
                  for key in udesc: msg+="DESCCOUNT"+"\t"+key.upper()+"\n"
                  if "WORDCOUNT" in self.glob.dictlist: msg+="NWORDTOTAL\n"
                  if "DESCCOUNT" in self.glob.dictlist: msg+="NPROT\n"
                  if "ECWEIGHT" in self.glob.dictlist: # evidence code weights for Blast2GO
                        for key in uacc: msg+="ECWEIGHT"+"\t"+key.upper()+"\n"
                  # send request to DictServer
                  tmp=DictServer.DICTquery(msg,self.pythonversion,REMOTE=self.REMOTE,HOSTNAME=self.DICTHOST,PORTNO=self.PORTNO).split("\n")
                  # initialize dictionaries
                  self.glob.GOdict={}
                  self.glob.GOdict_weights={}
                  self.glob.desccounts={}
                  self.glob.lineage={}
                  self.glob.taxid={}
                  self.glob.wordcounts={}
                  # load values to dictionaries (keys occur in chunk)
                  for row in tmp:
                        if not row: continue
                        (table,key,value)=row.split("\t")
                        if table == "GODICT":
                                self.glob.GOdict[key]=value
                        elif table == "DESCCOUNT":
                                self.glob.desccounts[key]=int(value)
                        elif table == "LINEAGE":
                                self.glob.lineage[key]=value
                        elif table == "WORDCOUNT":
                                self.glob.wordcounts[key]=int(value)
                        elif table == "TAXID":
                                self.glob.taxid[key]=value
                        elif table == "NPROT":
                                self.glob.nprot=int(value)
                        elif table == "NWORDTOTAL":
                                self.glob.nwordtotal=int(value)
                        elif table == "ECWEIGHT":
                                self.glob.GOdict_weights[key]=value
                        else:
                                sys.stderr.write("# unknown table: %s\n" %str(row))
                  sys.stderr.write("# Dictionary sizes: GOdict %i, lineage %i, taxid %i, desccounts %i, wordcounts %i\n" %(len(self.glob.GOdict), len(self.glob.lineage), len(self.glob.taxid), len(self.glob.desccounts), len(self.glob.wordcounts)))

        def catch_unique(self,target_col,lines):
                "Returns list of unique keys in target_col of lines. Lines has rows of tab-separated data. Keys are uppercase."
                tmp={}
                for line in lines:
                    if not line: continue # skip empty line
                    try:
                        key=line.split("\t")[target_col].upper()
                        tmp[key]=1
                    except:
                        sys.stderr.write("# Warning: column %i not found on line %s\n" %(target_col,line))
                return(tmp.keys())

        def SANSquery(self,message,SANSURL="http://ekhidna2.biocenter.helsinki.fi/cgi-bin/sans/sans.cgi"):
                # send chunk of query sequences to SANSparallel server
                sys.stderr.write("# Calling SANSparallel, message size is %i bytes\n" %len(message))
                if self.REMOTE:
                        values={'query': message, 'mode': 'raw', 'db': 'uniprot', 'H': '100', 'protocol': '0'}
                        # try 3 times
                        itried=0
                        while itried<2:
                          try :
                                r=requests.post(SANSURL,data=values)
                                break
                          except requests.exceptions.RequestException as rerr :
                                print("Error,", str(rerr), file=sys.stderr)
                                itried+=1
                                #sys.exit(1)

                        print("# Result size is %i bytes" % len(r.text), file=sys.stderr)
                        if r.text == "" :
                                print("Error: calling SANSparallel remotely gave empty result", file=sys.stderr)
                                sys.exit(1)
                        #tmp=r.text
                        if self.pythonversion == 2:
                                tmp=r.text.encode('ascii','ignore')
                        else:
                                tmp=r.text
                        #print("# SANS result type is ",type(tmp),file=sys.stderr)
                else: # local server
                        # extract hdr,seq from FASTA
                        qbuffer=''
                        hdr=''
                        seq=''
                        for line in message.split("\n"):
                            if not line: continue
                            line=line.replace("\"","")
                            if line[0]=='>':
                                        self.sentquery+=1
                                        qbuffer+="%i %i %i 20 2 11 1.0 %i %i %i \"%s\" \"%s\" </QUERY>\n" %(self.sentquery, self.H, self.HX, self.R, self.VOTELIST_SIZE, self.SANSPROTOCOL,seq[0:self.MAXRES],hdr[0:self.MAXRES])
                                        hdr=line.strip()
                                        seq=''
                            else:
                                        seq+=line.strip().upper().replace(" ","")
                        # last entry
                        qbuffer+="%i %i %i 20 2 11 1.0 %i %i %i \"%s\" \"%s\" </QUERY>\n" %(self.sentquery, self.H, self.HX, self.R, self.VOTELIST_SIZE, self.SANSPROTOCOL,seq[0:self.MAXRES],hdr[0:self.MAXRES])
                        #sys.stderr.write("# qbuffer %s" %qbuffer)
                        # send request to SANSparallel server in proper format
                        tmp=DictServer.Generic_client(qbuffer,self.pythonversion,HOSTNAME=self.SANSHOST,PORTNO=self.SANSPORT)
                        #print("SANS result type is ",type(tmp), file=sys.stderr)
                        #sys.stderr.write(tmp)
                sys.stderr.write("# SANSparallel returned %i bytes\n" %len(tmp))
                # format conversion to tabular bytestream
                data=self.xml.stream(tmp.split("\n"),header=False,bracket=False,output_sseq=self.SSEQ,output_ranges=self.RANGES)
                return(data[0],data[1]) # SANS tabular, metadata

        def test_newentry(self,input_format,line):
                if input_format=="FASTA":
                        try:
                                if line[0]==">":
                                        return(True)
                                else:
                                        return(False)
                        except:
                                return(False)
                if input_format=="tab":
                        row=line.split("\t")
                        try:
                                uid=row[self.block_column_index]
                                if uid != self.olduid:
                                        self.olduid=uid
                                        return(True)
                                else:
                                        return(False)
                        except:
                                return(False)
                return(False)

        def open_IO_channels(self,infile,OUT_ARRAY):
                """
In SANSPANZ, the following sheets are used:
FileOut = name of data spreadsheet's output file (default = no output)
OUT_DE = name of cluster_data spreadsheet's output file (default = no output)
OUT_GO = name of goclass_data spreadsheet's output file (default = no output)
OUT_ANNO = name of anno_data spreadsheet's output file (default = no output)
                """
                data_in=sys.stdin
                if not infile == "--": data_in=open(infile,"r")
                # Redirect output to files; default = sys.stdout
                for i in range(0,self.glob.nsheet):
                        self.glob.sheets[i].connection=None
                        try:
                                x=OUT_ARRAY[i]
                        except:
                                x=self.operator_name+'.out_'+str(i)
                                sys.stderr.write('# set output file %i to %s\n' %(i,x))
                        if x == "--":
                                self.glob.sheets[i].fh=sys.stdout
                        elif x:
                                self.glob.sheets[i].fh=open(x,"w")
                        else:
                                self.glob.sheets[i].fh=None
                return(data_in)

        def close_IO_channels(self):
                # close output files
                for  i in range(0,self.glob.nsheet):
                        if self.glob.sheets[i].fh: self.glob.sheets[i].fh.close()

if __name__=="__main__":
        # glob is object handle containing spreadsheets,dictionaries,parameters
        glob=Parameters.WorkSpace(configFile=None)
        # liveData is update every CHUNK queries
        z=Runner.Runner(glob, method=glob.param['input_OPERATOR'], CHUNK=glob.param['input_CHUNK'], liveData=glob.param['input_LIVEDATA'])
        # input_format='FASTA'|'tab'
        # if colnames as argument then data has no header rows
        # queryspecies "auto" implies queries are Uniprot entries, None implies queryspecies is defined in input data and input_format is "tab", otherwise given by --queryspecies option on command line
        z.lazyRunner(glob.param['input_FILE'], ['--','DE.out','GO.out','anno.out'], input_format="FASTA", colnames=None, queryspecies=glob.param['input_QUERYSPECIES'])
