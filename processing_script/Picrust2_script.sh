#!/bin/bash
#activating picrust2 environment
source activate picrust2
#place your fasta sequence and OTU biome table in the folder named data
#1st step is Sequence placement
place_seqs.py -s study_seqs_test.fasta -o placed_seqs.tre -p 4 --intermediate placement_working
#Kindly change the Number of processes to run in parallel depends on the system
#Hidden state prediction: Hidden-state prediction for 16S copy, E.C. numbers, and KO abundances per-genome can be run with these commands
hsp.py -i 16S -t placed_seqs.tre -o marker_nsti_predicted.tsv.gz -p 4 -n
hsp.py -i EC -t placed_seqs.tre -o EC_predicted.tsv.gz -p 4 
hsp.py -i KO -t placed_seqs.tre -o KO_predicted.tsv.gz -p 4
metagenome_pipeline.py -i test_input_sequence_abun.biom -m marker_nsti_predicted.tsv.gz -f EC_predicted.tsv.gz -o EC_metagenome_out --metagenome_contrib --strat_out
metagenome_pipeline.py -i test_input_sequence_abun.biom -m marker_nsti_predicted.tsv.gz  -f KO_predicted.tsv.gz -o KO_metagenome_out --metagenome_contrib --strat_out
#Infer pathway abundances
pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -o pathways_out --intermediate minpath_working  -p 4 --regroup_map --per_sequence_contrib 
#Add descriptions
add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO  -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC  -o pathways_out/path_abun_unstrat_descrip.tsv.gz


