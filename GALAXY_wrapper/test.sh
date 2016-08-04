cp run_enrichment.py ../
cp run_ranked_enrichment.py ../

#OPTIONS="-m cProfile"
#OPTIONS="-m memory_profiler"
echo 'test GSEA: p=100, rank100.txt'
time python $OPTIONS ../run_ranked_enrichment.py  GSEA -o test1.html -i test-data/rank100.txt -g test-data/goslim_generic.obo -a test-data/goslim.gaf -f html -p 100 --seed 1234
diff test-data/outputs/test1.html test1.html | head

echo 'test Ranked parent-child: intersection, rank100.txt'
time python $OPTIONS ../run_ranked_enrichment.py parent-child -o test2.html -i test-data/rank100.txt -g test-data/goslim_generic.obo -a test-data/goslim.gaf -f html -m intersection 
diff test-data/outputs/test2.html test2.html | head

echo 'test Parent-child: union, Bonferroni and BH_FDR corrections, list1000.txt'
time python $OPTIONS ../run_enrichment.py parent-child -o test3.csv -i test-data/list1000.txt -g test-data/goslim_generic.obo -a test-data/goslim.gaf -f tabular -c bonferroni bh_fdr -m union 
diff test-data/outputs/test3.csv test3.csv | head

echo 'test term-for-term: Bonferroni correction, list1000.txt'
time python $OPTIONS ../run_enrichment.py term-for-term -o test4.csv -i test-data/list1000.txt -g test-data/goslim_generic.obo -a test-data/goslim.gaf -f tabular -c bonferroni
diff test-data/outputs/test4.csv test4.csv | head




