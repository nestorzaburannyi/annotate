from __future__ import print_function
# SANS tags:
#       QUERY nid= vote_cutoff= LSEQ=
#       SBJCT VOTE= TUPS= PIDE= LALI= BITS= EVALUE= DIAG= LSEQ=
# computed:
#       qpid, spid
#       qcov,scov = LALI/LSEQ
#       desc, species, gene
#       rank
#       qseq, seq

import re,sys

class XMLParser:
    def __init__(self,queryspecies="auto"):
            """ If queryspecies is "auto" then assume query comes from Uniprot and species is taken from OS= tag in header.
            """
            self.autospecies=(queryspecies=="auto") # Uniprot query
            self.queryspecies=queryspecies
            self.iquery=0
            self.init_entry()
            self.init_sbjct()
            self.inquery=0
            self.insbjct=0

    def stream(self,indata,header=False, bracket=False, output_sseq=False, output_ranges=False):
        self.init_entry()
        #sys.stderr.write("# data size = %i; output_ranges = %s\n" %(len(indata),output_ranges))
        # header line
        result=""
        if bracket: result='#{\n'
        if header:
                result += "\t".join('nid isquery qpid spid qcov scov bits pide lali desc species qseq vote genename'.split())
        if output_ranges:
                result += "\t"
                result += "\t".join('qfrom qto qseqlen sfrom sto sseqlen'.split())
        result += "\n"
        if bracket: result += "#}\n"
        iline=0
        metadata=[]
        sseq="n.d." # default: no subjct sequence output
        infasta=0
        for line in indata:
                #print('#stream: ', line, file=sys.stderr)
                iline+=1
                if not line: continue
                if line[0]=='#': continue
                line=line.strip()
                if len(line)<=1: continue
                # FASTA part
                if line[0]=='>':
                        try:
                                pid,hdr=line[1:].split(" ",1)
                        except:
                                pid=line[1:]
                                hdr=line[1:]
                        # extract species from hdr OS=... GN=...
                        if self.inquery==1 and self.autospecies:
                                try:
                                        self.queryspecies=re.sub(r' \w{2}=.*$',"",hdr[hdr.index(r"OS=")+3:])
                                except:
                                        self.queryspecies="unknown"
                        if self.insbjct:
                                try:
                                        genename=re.sub(r' \w{2}=.*$',"",hdr[hdr.index(r"GN=")+3:])
                                except:
                                        genename=""
                                try:
                                        species=re.sub(r' \w{2}=.*$',"",hdr[hdr.index(r"OS=")+3:])
                                except:
                                        species="unknown"
                        # remove uniprot tags from desc OS=... GN=... PE=...
                        desc=re.sub(r' \w{2}=.*$',"",hdr)
                        infasta=1
                elif infasta:
                        # sequence is next line from header line
                        seq=line.rstrip()
                        if self.inquery==1:
                                # print query row
                                self.qpid=pid
                                self.qseqlen=self.lseq
                                genename=""
                                result += "\t".join([str(self.iquery),"1",self.qpid,self.qpid,"1.00","1.00","0.0","1.00",str(self.qseqlen),desc,self.queryspecies,seq,"0",genename])
                                if output_ranges:
                                        result += "\t"
                                        result += "\t".join(["1",str(self.qseqlen),str(self.qseqlen),"1",str(self.qseqlen),str(self.qseqlen)])
                        elif self.insbjct==1:
                                # compute qcov,scov
                                qcov=str(self.lali/float(1++self.qseqlen))
                                scov=str(self.lali/float(1+self.lseq))
                                # print sbjct row
                                if not output_sseq: seq="" # suppress sequence output
                                result += "\t".join([str(self.iquery),"0",self.qpid,pid,qcov,scov,str(self.bits),str(self.pide),str(self.lali),desc,species,seq,self.vote,genename])
                                if output_ranges:
                                        result += "\t"
                                        result += "\t".join([str(self.qfrom),str(self.qto),str(self.qseqlen),str(self.sfrom),str(self.sto),str(self.lseq)])
                        result += "\n"
                        infasta=0
                else: # parse XML tags
                        x=line.split()
                        i=0
                        self.ok=1
                        while i<len(x):
                          try:
                                if re.search('<QUERY',x[i]):
                                        self.init_entry()
                                        self.inquery=1
                                        self.insbjct=0
                                        self.iquery+=1
                                        #print("iquery = %i iline = %i" %(self.iquery,iline), file=sys.stderr)
                                        if bracket: result += "#{\n"
                                elif x[i]=='<SBJCT':
                                        self.init_sbjct()
                                        self.inquery=0
                                        self.insbjct=1
                                elif x[i]=='LALI=':
                                        self.lali=int(x[i+1])
                                elif x[i]=='VOTE=':
                                        self.vote=str(x[i+1])
                                elif x[i]=='PIDE=':
                                        self.pide=float(x[i+1])
                                elif x[i]=='BITS=':
                                        self.bits=float(x[i+1])
                                elif x[i]=='LSEQ=':
                                        self.lseq=int(x[i+1])
                                elif x[i]=='QFROM=':
                                        self.qfrom=int(x[i+1])
                                elif x[i]=='QTO=':
                                        self.qto=int(x[i+1])
                                elif x[i]=='SFROM=':
                                        self.sfrom=int(x[i+1])
                                elif x[i]=='STO=':
                                        self.sto=int(x[i+1])
                                elif x[i]=='</SBJCT>':
                                        self.insbjct=0
                                elif x[i]=='</QUERY>':
                                        self.inquery=0
                                        if bracket: result += "#}\n"
                                elif x[i]=='nid=':
                                        self.nid=x[i+1] # this is not output, self.iquery is
                                elif x[i]=='<DATABASE=':
                                        metadata=[ x[i+1],x[i+3],x[i+5] ] # database version; letters, sequences
                          except:
                                self.ok=0
                                break
                          i+=1
                        if self.ok==0: print("# XML parse error: ", self.iquery, line, file=sys.stderr)
        #result += "\n#}\n"
        #print("# metadata: ",metadata, file=sys.stderr)
        return(result,metadata)

    def init_entry(self):
        self.qseqlen=1
        self.qpid=''

    def init_sbjct(self):
        self.lseq=1
        self.pid=""
        self.bits=0.0
        self.pide=0.0
        self.lali=0
        self.vote="0"
        self.qfrom=1
        self.qto=0
        self.sfrom=1
        self.sto=0

if __name__=="__main__":
        #queryspecies=sys.argv[1]
        p=XMLParser()
        indata=sys.stdin.readlines()
        result=p.stream(indata,header=True,output_ranges=True,output_sseq=False)
        print(result[0])
        print(result[1],file=sys.stderr)
