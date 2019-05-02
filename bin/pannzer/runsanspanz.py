from __future__ import absolute_import
from __future__ import print_function

import Parameters, Runner
import sys
#sys.path.append('/home/luholm/PANNZER2')


if __name__ == '__main__' :
    text="""
 python runtest.py -h
 python runtest.py -m taxonassignment -i sans.tab -f tab -o ",--" 2> e
 python runtest.py -m lineage -i sans.tab -f tab  2> e
 python runtest.py -m BestInformativeHit -i sans.tab -f tab -o ',--'
 python runtest.py -m Pannzer -i sans.tab -f tab 2> e
 python runtest.py -m SANS -f FASTA  < /data/liisa/test.fasta  > sans.tab
 python runtest.py -m Pannzer -f FASTA -i /data/liisa/test.fasta  -o 'sans.tab,DE.out,GO.out,anno.out' --PANZ_FILTER_PERMISSIVE 
 cut -f 2,5,7 /data/uniprot/goa_uniprot_noIEA.gaf | python runtest.py -m gaf2tab -f tab -b "acc" -c "acc goid" -o ",goa_truth_propagated"

 python runtest.py -m obo -i /data/uniprot/go-basic.obo -f txt -o 'obo.tab,'
 cut -f 2,5 /data/uniprot/goa_uniprot_all.gaf | python runtest.py -m gaf2propagated -f tab -c "qpid goid" -o ",godict.tab" 2> e
 cut -f 2 godict.tab | sort | uniq -c | perl -pe 's/^\s+//' | perl -pe 's/ /\t/' > godict_counts
 python runtest.py -m BayesIC -i godict_counts -f tab -c 'qpid propagated' -o ',obo_with_ic.tab' 

 python runtest.py -m GOevaluation -f tab -o ',,--' -i GO.out

If anyone desperately wants to use BLAST:

 perl runblast.pl /data/liisa/test.fasta /data/liisa/uniprot.fasta 40 > blast.tab 
 python runtest.py -m BestInformativeHit -i blast.tab -f tab -o ',--'
 python runtest.py -m Pannzer -f tab -i blast.tab  -o ',DE.out,GO.out,anno.out'
    """
    example="""
PANNZER2 function prediction: 
        python %s -m Pannzer -s "Pusa hispida saimensis" < querysequences.fasta
	- query species given by -s argument
        - output goes to STDOUT and files Pannzer.out_1, Pannzer.out_2, Pannzer.out_3

	python %s -m Pannzer -i querysequences.fasta -o 'sans.tab,DE.out,GO.out,anno.out'
	- query species parsed from Uniprot header
	- output goes to files sans.tab,DE.out,GO.out,anno.out

More options:
	python %s -h

More info:
	http://ekhidna2.biocenter.helsinki.fi/sanspanz

""" %(sys.argv[0],sys.argv[0],sys.argv[0])
    
    if len(sys.argv) == 1: 
        print(example, file=sys.stderr)
        sys.exit(1)

    try :
        glob=Parameters.WorkSpace()
        z=Runner.Runner(glob,operator_name=glob.param['input_OPERATOR'],CHUNK=glob.param['input_CHUNK'],liveData=glob.param['input_LIVEDATA'])
        z.lazyRunner(glob.param['input_FILE'], glob.param['input_OUTFILES'].split(","), queryspecies=glob.param['input_QUERYSPECIES'],colnames=glob.param['input_COLNAMES'], block_column_name=glob.param['input_BLOCK_COLUMN_NAME'], input_format=glob.param['input_FORMAT'] )

    except KeyboardInterrupt :
        print("CLIENT> killed by user...\n", file=sys.stderr)
        sys.exit(1)

