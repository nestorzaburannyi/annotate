from __future__ import absolute_import
from __future__ import print_function
import Parameters, Runner
import sys
#sys.path.append('/home/luholm/AAI-profiler')

if __name__ == '__main__' :
    example="""
This is a customized script that performs AAI-profiling starting from previously generated SANS results.

Usage:
        python %s -i sanstabfile -d HIDFOLER [-R]

Outputs:
        HIDFOLDER/tax.tab = onesided AAI profiles per species
        HIDFOLDER/tax2.tab = bidirectional AAI profiles per species
        HIDFOLRER/kt.txt = onesided krona input excluding top species
        HIDFOLDER/kt_self.txt = onesided krona input including top species

Shell script for AAI-profiler:
        python runsanspanz.py -m SANS -i inputfastafile -o /data/tmp/HID/sans.tab
        python AAI-profiler-tab.py -i /data/tmp/HID/sans.tab -d HIDFOLDER
        perl tax2scatter.pl "onesided AAI-profiler results" onesided HIDFOLDER 0.5 20 0 1 0 1 < HIDFOLDER/tax.tab > ~/public_html/AAI.html
        perl tax2scatter.pl "bidirectional AAI-profiler results" bidirectional HIDFOLDER 0.5 20 0 1 0 1 < HIDFOLDER/tax2.tab > ~/public_html/AAI2.html
        /usr/local/bin/ktImportText -u http://krona.sourceforge.net -o ~/public_html/krona.html HIDFOLDER/kt.txt
        /usr/local/bin/ktImportText -u http://krona.sourceforge.net -o ~/public_html/krona_self.html HIDFOLDER/kt_self.txt

        """ %sys.argv[0]
    if len(sys.argv) < 5 :
        print(example, file=sys.stderr)
        sys.exit(1)
    try :
        glob=Parameters.WorkSpace()
        # use remote access to SANSparallel and DictServer
        USEREMOTE=glob.param['CONN_REMOTE']
        MAXHITS=2000 # include all hits
        # use remote access to SANSparallel and DictServer
        glob.param['CONN_REMOTE'] = USEREMOTE
        print('# Remote access is',USEREMOTE, file=sys.stderr)
        # output files
        sanstabfile=glob.param['input_FILE']
        prefix=glob.param['DATA_DIR']
        taxtabfile=prefix+'/tax.tab'
        taxtabfile2=prefix+'/tax2.tab'
        kttxtfile=prefix+'/kt.txt'
        ktselftxtfile=prefix+'/kt_self.txt'
        # step 0: run SANS
        #print >> sys.stderr, '# Running SANS over the Internet. Progress monitored in ',liveData
        #z=Runner.Runner(glob,operator_name='SANS',CHUNK=100,liveData=liveData)
        #z.lazyRunner(glob.param['input_FILE'],[sanstabfile],input_format="FASTA")
        # * onesided AAI *
        # step 1: collect AAI profiles per species
        glob.param['PANZ_MAXHITS']=MAXHITS
        glob.param['input_QUERYSPECIES']='species'
        glob.kt_format=False
        print('# Running AAI. MAXHITS is ',glob.param['PANZ_MAXHITS'], file=sys.stderr)
        z=Runner.Runner(glob,operator_name='AAI',liveData=None)
        z.lazyRunner(sanstabfile, [None,taxtabfile], input_format='tab', block_column_name='qpid' )
        # step 2: get best hit from sans.tab including top species
        # re-initialize workspace to get clean sheet headers
        glob=Parameters.WorkSpace()
        glob.param['CONN_REMOTE'] = USEREMOTE
        glob.param['PANZ_MAXHITS']=1
        glob.param['input_QUERYSPECIES']='species'
        glob.kt_format=True
        print('# Extracting best hit. All species included. MAXHITS is ',glob.param['PANZ_MAXHITS'], file=sys.stderr)
        z=Runner.Runner(glob,operator_name='AAI',CHUNK=1000,liveData=None)
        z.lazyRunner(sanstabfile, [None,ktselftxtfile], input_format='tab', block_column_name='qpid' )
        # step 3: get best hit from sans.tab excluding top species
        glob=Parameters.WorkSpace()
        glob.param['CONN_REMOTE'] = USEREMOTE
        glob.param['PANZ_MAXHITS']=1
        glob.kt_format=True
        fh=open(taxtabfile,'r')
        species_col=fh.readline().rstrip().split("\t").index('species')
        glob.param['input_QUERYSPECIES']=fh.readline().split("\t")[species_col].rstrip()
        fh.close()
        print('# Extracting best hit. %s excluded. MAXHITS is ' %glob.param['input_QUERYSPECIES'],glob.param['PANZ_MAXHITS'], file=sys.stderr)
        z=Runner.Runner(glob,operator_name='AAI',CHUNK=1000,liveData=None)
        z.lazyRunner(sanstabfile, [None,kttxtfile], input_format='tab', block_column_name='qpid' )
        # * bidirectional AAI *
        glob=Parameters.WorkSpace()
        glob.param['PANZ_MAXHITS']=MAXHITS
        glob.param['input_QUERYSPECIES']='species'
        glob.kt_format=False
        print('# Running AAI2. MAXHITS is ',glob.param['PANZ_MAXHITS'], file=sys.stderr)
        z=Runner.Runner(glob,operator_name='AAI2',liveData=None)
        z.lazyRunner(sanstabfile, [None,taxtabfile2], input_format='tab', block_column_name='qpid' )

    except KeyboardInterrupt :
        print("CLIENT> killed by user...\n", file=sys.stderr)
        sys.exit(1)

