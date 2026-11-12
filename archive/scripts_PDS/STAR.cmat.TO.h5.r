### read StarSoloResults
library("scCustomize")

ddir1 <- list.dirs("/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/In_house/snRNA_Nphs2mut_HeChen/Data/bam/STAR/STARpremRNA", recursive = F)
ddir2 <- list.dirs("/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/In_house/snRNA_Wt1hetdel_HeChen/Data/bam/STAR/STARpremRNA", recursive = F)[1:6]
dir3 <- list.dirs("/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/In_house/snRNA_Pdss2_Coq2_CemOezel/Data/bam/STAR/premRNA", recursive = F)

ddir <- dir3
lapply( seq(ddir) , function(ii){
  Create_10X_H5.starSolo(
    raw_data_file_path= paste0(ddir[ii], "/GeneFull/filtered"),
    source_type = "10X",
    save_file_path="/data/public/tpadvits/PROJECTS/PodocytePJ/PDS_submission/GEOsubmission_PDSmanuscript_snRNAseq/StarSolo/",
    save_name = paste0( sub(".*_(.*?)_.*", "\\1", basename(ddir[ii])), "_filtered_feature_bc_matrix" )
  )
})
