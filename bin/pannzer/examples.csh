# function prediction example
python runsanspanz.py -R -s "Macaca mulatta" < testdata/querysequences.fasta
echo "results: Pannzer.out_1 Pannzer.out_2 Pannzer.out_3"
# convert Pannzer result to HTML
perl anno2html.pl "test ARGOT" argot < Pannzer.out_3 > predictions.html
echo "result: predictions.html"
# named output files
python runsanspanz.py -R -m Pannzer -s "Macaca mulatta" -i testdata/querysequences.fasta -o ",DE.out,GO.out,anno.out"
echo "results: DE.out GO.out anno.out"
# save homology search result
python runsanspanz.py -R -m SANS -s "Macaca mulatta" -i testdata/querysequences.fasta -o sans.tab
python runsanspanz.py -R -m Pannzer -i sans.tab -f tab -o ",DE.out1,GO.out1,anno.out1" 
echo "results: sans.tab DE.out1 GO.out1 anno.out1"
# select Pannzer predictors
python runsanspanz.py -R -m Pannzer -i sans.tab -f tab -o ",DE.out2,GO.out2,anno.out2" --PANZ_PREDICTOR "DE,RM3,ARGOT,JAC,HYGE,SANS"
echo "results: DE.out2 GO.out2 anno.out2"
# relax filtering criteria
python runsanspanz.py -R --PANZ_FILTER_PERMISSIVE -m Pannzer -i sans.tab -f tab -o ",DE.out3,GO.out3,anno.out3"
echo "results: DE.out3 GO.out3 anno.out3"
# using BLAST 
#perl runblast.pl testdata/querysequences.fasta testdata/uniprot.fasta 40 > blast.tab
python runsanspanz.py -R -m BestInformativeHit -i testdata/blast.tab -f tab -o ',--' 2> err
python runsanspanz.py -R -m Pannzer -i testdata/blast.tab -f tab -o ',DE.out4,GO.out4,anno.out4'
echo "results: DE.out4 GO.out4 anno.out4"
# GO evaluation
# wget http://www.berkeleybop.org/ontologies/go/go-basic.obo
# python runsanspanz.py -m obo -i /data/uniprot/go-basic.obo -o 'obo.tab,'
python runsanspanz.py -R -m gaf2propagated -f tab -i testdata/gotest_truth.gaf --eval_OBOTAB 'testdata/obo.tab' -o ",gotest_truth_propagated" -c "qpid goid"
python runsanspanz.py -R -m GOevaluation -f tab -o ',,eval1' -i testdata/GO.out --eval_SCOREFUNCTIONS "RM3_PPV ARGOT_PPV JAC_PPV HYGE_PPV" --eval_TRUTH testdata/gotest_truth_propagated
echo "result: eval1"
# naive predictions and evaluation 
python runsanspanz.py -R -m naive -i testdata/target_identifiers.list -o ',naive_predictions.out' -c 'qpid' -f tab
python runsanspanz.py -R -m GOevaluation -i testdata/naive_predictions.out -f tab -o ',,eval2' --eval_SCOREFUNCTIONS "frequency" --eval_TRUTH testdata/gotest_truth_propagated  2> err 
echo "result: eval2"
# compare to example results  in testresults/

