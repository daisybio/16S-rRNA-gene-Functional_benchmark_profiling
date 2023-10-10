#!/bin/bash

#Input files:
#otu_table.tsv 
#taxtab.tsv
#Example Format for Taxatab.tsv
#Feature ID	Taxon	Confidence
#c28d0a4d9f1c7b39400278482dcb0d59	k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Prevotellaceae; g__Prevotella; s__copri	0.999999562
#bccda5b8e71776dfdfbd836b63e0ef33	k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__	0.999994775


#activate conda environment
#conda activate metgem_env


#Convert ASVs table into KO profiles and EC profiles
metgem markp -i "${PWD%*/*/*}"/Input_preparation/OTU_table.txt -t "${PWD%*/*/*}"/Input_preparation/Tax.txt -m k_core -o output_ko.tsv
metgem markp -i "${PWD%*/*/*}"/Input_preparation/OTU_table.txt -t "${PWD%*/*/*}"/Input_preparation/Tax.txt -m e_core -o output_ec.tsv
