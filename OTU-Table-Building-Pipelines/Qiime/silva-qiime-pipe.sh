#input 1 should be the input file
#input 2 should be the the output dir

pick_otus.py -r /fslgroup/fslg_jena_magic/compute/SILVA_132_QIIME_release/rep_set/rep_set_all/97/silva132_97.fna -i $1 -m uclust_ref -o $2

make_otu_table.py -i $2/$2_otus.txt -t /fslgroup/fslg_jena_magic/compute/SILVA_132_QIIME_release/taxonomy/taxonomy_all/99/consensus_taxonomy_7_levels.txt -o $2/$2-97.biom

biom convert -i $2/$2-97.biom -o $2/$2-97.otu_table_with_taxonomy.txt --to-tsv --header-key taxonomy --output-metadata-id "ConsensusLineage"


