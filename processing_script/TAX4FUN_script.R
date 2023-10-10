#Install the package in R using the command below
setwd("../Collaboration")
if (!require('Tax4Fun2'))install.packages(pkgs = "Tax4Fun2_1.1.5.tar.gz", repos = NULL, source = TRUE); library(Tax4Fun2)

#path_to_otus = "" > Specifiy the path to the fasta file containing your otu seqeunces
#path_to_otu_table = "" > Specifiy the path to the otu table (without taxonomic annotation)
#path_to_refernce_data = "" > Specifiy the path to the folder with the reference data
#path_to_temp_folder = "" > Specifiy the path to the tempory folder. Results will be stored here.
#database_mode = "Ref99NR" > Run either with 'Ref99NR' or 'Ref100NR' as database mode. The number refers to the clustering threshold used in uclust (99% and 100%, respectively)
#num_threads = 6 > Number of parallel threads for diamond
#use_force = TRUE > > Overwrite folder if exists (Default is FALSE)
#normalize_by_copy_number = TRUE > normalize by the average number of 16S rRNA copies calculated for each sequence in the reference database
#normalize_pathways = FALSE will affiliate the rel. abundance of each KO to each pathway it belongs to. By setting it to true, the rel. abundance is equally distributed to all pathways it was assigned to.)


#set the number of core
nc=4
#Input files
#provide the path for Otu file and reference data
otu_seq<-"C:/Users/selva/Documents/Cohort/16s_cnv_correction_databases-main/Collaboration/Custom Database/Tax4FUN/OTUs-Seqs.fasta"
ref_path<-"C:/Users/selva/Documents/Cohort/16s_cnv_correction_databases-main/Collaboration/Custom Database/Tax4FUN2/Tax4Fun2_ReferenceData_v2/"
otu_table_path<-"C:/Users/selva/Documents/Cohort/16s_cnv_correction_databases-main/Collaboration/Custom Database/Tax4FUN/seqtab_norm_EC_IRTREE.tsv"

buildReferenceData(path_to_working_directory = "C:/Users/selva/Documents/Cohort/16s_cnv_correction_databases-main/Collaboration/Custom Database/Tax4FUN/", use_force = FALSE, install_suggested_packages = TRUE)
testReferenceData(path_to_reference_data = "Tax4FUN2/Tax4Fun2_ReferenceData_v2/")
buildDependencies(path_to_reference_data = ref_path, install_suggested_packages = TRUE,use_force=TRUE)
#Run the reference blast
runRefBlast(path_to_otus = otu_seq, path_to_reference_data = ref_path, path_to_temp_folder = "C:/Users/selva/Documents/Cohort/16s_cnv_correction_databases-main/Collaboration/Custom Database/Tax4FUN/Kelp_Ref99NR", database_mode = "Ref99NR", use_force = T, num_threads = nc)

# Predicting functional profiles
makeFunctionalPrediction(path_to_otu_table = otu_table_path, path_to_reference_data = ref_path, path_to_temp_folder = "C:/Users/selva/Documents/Cohort/16s_cnv_correction_databases-main/Collaboration/Custom Database/Tax4FUN/Kelp_Ref99NR", database_mode = "Ref99NR", normalize_by_copy_number = TRUE, min_identity_to_reference = 0.97, normalize_pathways = TRUE)


##################Additional Analysis################################################
#Calculating (multi-)functional redundancy indices (experimental)
runRefBlast(path_to_otus = otu_seq, path_to_reference_data = ref_path, path_to_temp_folder = "Temp_Ref99NR", database_mode = "Ref99NR", use_force = T, num_threads = nc)

#Calculating FRIs
calculateFunctionalRedundancy(path_to_otu_table = otu_table_path, path_to_reference_data = ref_path, path_to_temp_folder = "Temp_Ref99NR", database_mode = "Ref99NR", min_identity_to_reference = 0.97)

#New in the latest pre-release (v1.1.6): prevalence_cutoff (see comment on pre-release)
calculateFunctionalRedundancy(path_to_otu_table = otu_table_path, path_to_reference_data = ref_path, path_to_temp_folder = "Temp_Ref99NR", database_mode = "Ref99NR", min_identity_to_reference = 0.97, prevalence_cutoff = 1.0)



#Output Files Needed for further analysis

#Output from from functional prediction analysis: functional_prediction.txt
#Output from from Pathway prediction analysis: pathway_prediction.txt
