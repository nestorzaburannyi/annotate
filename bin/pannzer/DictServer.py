from __future__ import absolute_import
from __future__ import print_function
import sys,re
try:
  import SocketServer
  import threading
  import signal
  import socket
  import time
  import requests
except ImportError as ie:
  print(str(ie),file=sys.stderr)
  sys.exit(1)
import Parameters


class DictServerHandler(SocketServer.BaseRequestHandler) :
    def handle(self) :
        #print >> stderr, "--> RequestHandler called"
        message = ''
        while True :
            indata = self.request.recv(self.server.BUFSIZ)
            if indata :
                if self.server.stream :
                    self.server.func(message, self.request)
                else :
                    message += indata

            else :
                if not self.server.stream :
                    result = self.server.func(message, self.request)

                self.request.shutdown(0)
                break

class ThreadedDictServer(SocketServer.ThreadingMixIn, SocketServer.TCPServer, object) :
    def __init__(self, server_address, request_handler, func, stream=False, BUFSIZ=1024) :
        super(ThreadedDictServer, self).__init__(server_address, request_handler)
        self.func = func
        self.stream = stream
        self.BUFSIZ = BUFSIZ

class DictServer(object):
    """
    This is a special class fo SANSPANZ prediction server. Input/output is through socket.
    It only outputs glob.anno_data sheet. Inputs should be generated using sansparser.pl.
    """
    def __init__(self, glob):
        "glob is a globs object containing object handles to spreadheets, dictionaries and parameters."
        # load dictionaries
        self.glob = glob
        base=glob.param['DATA_DIR']+'/'
        fn=base+glob.param['DATA_PHR']
        #self.glob.load_PHR(fn) # key = accession number, value = description
        fn=base+glob.param['DATA_WORDCOUNTS']
        self.glob.load_wordcounts(fn)
        fn=base+glob.param['DATA_NWORDTOTAL']
        self.glob.nwordtotal=self.glob.load_nwordtotal(fn)
        fn=base+glob.param['DATA_DESCCOUNTS']
        self.glob.load_desccounts(fn)
        fn=base+glob.param['DATA_NPROT']
        self.glob.nprot=self.glob.load_nprot(fn)
        fn=base+glob.param['DATA_GODICT']
        self.glob.load_GOdict(fn)
        fn=base+glob.param['DATA_TAXONOMY']
        self.glob.load_lineage(fn)
        self.glob.load_taxid(fn)
        fn=base+glob.param['DATA_GOIDELIC']
        self.glob.load_goidelic(fn)

        print("# Dictionaries loaded", file=sys.stderr)
        print("# nwordtotal=%i" % self.glob.nwordtotal, file=sys.stderr)
        print("# nprot=%i" % self.glob.nprot, file=sys.stderr)

    def lookup_key_values(self, message, connection):
        """
        Inputs: (table\tkey) tuples separated by tabs and newlines
        Outputs: (table\tkey\tvalue) tuples separated by tabs and newlines
        """
        # build result
        result=''
        # do dictionary lookups
        for row in message.split("\n"):
            data=re.split("\W", row, 1) #row.split("\t")
            #print >> stderr, data
            dicti=data[0].upper()
            if dicti == "GOIDELIC":
                self.download_goidelic(connection)
                return
            if dicti == "NPROT":
                key=dicti
                value=self.glob.nprot
            elif dicti == "NWORDTOTAL":
                key=dicti
                value=self.glob.nwordtotal
            else:
                try:
                    key=data[1].upper()
                    if dicti=='GODICT':
                        value=self.glob.GOdict[key].replace("\t",",") # join parent classes to direct classes list
                    elif dicti=='DESCCOUNT':
                        value=self.glob.desccounts[key]
                    elif dicti=='LINEAGE':
                        value=self.glob.lineage[key]
                    elif dicti=="TAXID":
                        value=self.glob.taxid[key]
                    elif dicti=='WORDCOUNT':
                        value=self.glob.wordcounts[key]
                    elif dicti=='GOCOUNT':
                        value=self.glob.GOcounts[key]
                    elif dicti=='GODESC':
                        value=self.glob.godesc[key]
                    elif dicti=='ONTOLOGY':
                        value=self.glob.ontology[key]
                    elif dicti=='GOPARENTS':
                        value=self.glob.GOparents[key]
                    elif dicti=='PHR':
                        value=self.glob.PHR[key]
                    elif dicti=='ECWEIGHT':
                        value=self.glob.GOdict_weights[key]
                    else:
                        value="unknown dictionary %s" % dicti
                except:
                    continue

            result+="\t".join([dicti,key,str(value)])+"\n"

        # send result to client
        connection.sendall(result)

    def download_goidelic(self,connection):
        "send copy of mergeGO.out to remote user via socket connection"
        result=""
        for key in self.glob.GOcounts.keys():
            ontology=''
            desc=''
            plist=''
            ec=''
            kegg=''
            if key in self.glob.ontology: ontology=self.glob.ontology[key]
            if key in self.glob.godesc: desc=self.glob.godesc[key]
            if key in self.glob.GOparents: plist=",".join(self.glob.GOparents[key])
            if key in self.glob.EC: ec=self.glob.EC[key]
            if key in self.glob.KEGG: kegg=self.glob.KEGG[key]
            result+="\t".join([self.glob.GOcounts[key],ontology,key,desc,plist,ec,kegg,self.glob.IC[key]])+"\n"
        connection.sendall(result)

def Generic_client(message, pythonversion, HOSTNAME='localhost', PORTNO=50002, BUFSIZE=1000, verbose=False):
    "Sends message to PORTNO, returns result"
    # Create a TCP/IP socketi
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    # Connect the socket to the port where the server is listening
    server_address = (HOSTNAME, PORTNO)

    if verbose:
        print("# CLIENT server_address =", server_address, file=sys.stderr)
        print("CLIENT> connecting to %s port %s" % server_address, file=sys.stderr)

    sock.settimeout(600)
    #keyboardinterrupt=False
    exitaftercleanup = False

    try :
        sock.connect(server_address)
    except :
        print("CLIENT> connection failed, terminating...", file=sys.stderr)
        sys.exit(1)

    result=''
    try:
        if verbose:
            print('sending "%s"' % message, file=sys.stderr)
        if pythonversion == 3:
            sock.sendall(message.encode('utf-8')) # send bytes
        else:
            sock.sendall(message)
        sock.shutdown(1)

        # Receive the data in small chunks and retransmit it
        amount_received=0
        while True:
            try:
                if pythonversion == 3:
                    data = sock.recv(BUFSIZE).decode('ascii','ignore')
                else:
                    data = sock.recv(BUFSIZE)
                amount_received += len(data)

                if verbose:
                    print('CLIENT> so far received %i' % amount_received, file=sys.stderr)

                if data:
                    result+=data
                else:
                    if verbose:
                        print('CLIENT> no more data from', server_address, file=sys.stderr)
                    break

            except socket.timeout :
                print("CLIENT> socket timed out", file=sys.stderr)
                sock.close()
                sys.exit(1)

            except IOError as ioe : # socket.error is child of IOError from 2.6
                print('CLIENT> socket error (%s)' % str(ioe), file=sys.stderr)
                sock.close()
                sys.exit(1)
    finally:
        sock.close()

    return result

def DICTquery(message,pythonversion,REMOTE=False,HOSTNAME='localhost',PORTNO=50002,DICTURL="http://ekhidna2.biocenter.helsinki.fi/cgi-bin/sans/DictServer.cgi"):
        if len(message)<1: return("")
        sys.stderr.write("# Calling DictServer, REMOTE=%s message size is %i bytes\n" %(REMOTE,len(message)))
        if REMOTE:
                values={'query': message}
                r=requests.post(DICTURL,data=values)
                sys.stderr.write("# Result size is %i bytes\n" %len(r.text))
                tmp=r.text
        else: # local server
                tmp=Generic_client(message,pythonversion,HOSTNAME=HOSTNAME,PORTNO=PORTNO)
        sys.stderr.write("# DictServer returned %i bytes\n" %len(tmp))
        return(tmp)

def Generic_server(function, hostname='localhost', portno=50002, stream=False, verbose=False, BUFSIZE=1000):
        """
        Starts the server part of client-server on host (hostname) at port (portno).
        The server runs until killed. Function processes inputs and sends the result to connection.
        If stream is True, function must chop the stream to smaller chunks for processing.
        If stream is False, all inputs are processed in one chunk by function.

        DictServer: Generic_server(DictServer.lookup_key_values, portno=50002, stream=False)
        sanspanz-server: Generic_server(chop_FASTA_stream, portno=500001, stream=True)
        tab-server: Generic_server(chop_bracket_stream, portno=500001, stream=True)
        """
        server = ThreadedDictServer((hostname, portno), DictServerHandler, function, stream, BUFSIZE)
        server.serve_forever()

def main() :

    glob = Parameters.WorkSpace()
    d = DictServer(glob)

    print("--> Starting DictServer...", file=sys.stderr)
    #Generic_server(d.lookup_key_values, portno=glob.param.PORTNO, stream=False)
    hostname = glob.param['CONN_HOSTNAME']
    port = glob.param['CONN_PORTNO']
    server = ThreadedDictServer((hostname, port), DictServerHandler, d.lookup_key_values, False, 1024)

    server_thread = threading.Thread(target=server.serve_forever)
    server_thread.daemon = True
    server_thread.start()

    signal.pause()

    print("--> Shutting down DictServer...", file=sys.stderr)
    server.shutdown()
    server.server_close()

    print("--> Done!", file=sys.stderr)

    return 0

if __name__ == '__main__' :
    try :
        main()

    except KeyboardInterrupt :
        print("\n\n\nKilled by User (please wait for exit)...\n", file=sys.stderr)
