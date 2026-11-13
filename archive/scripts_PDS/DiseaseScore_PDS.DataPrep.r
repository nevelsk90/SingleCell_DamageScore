# ###################################################### #
#  prepare SCSN podocyte data for PDS analysis, select TFs, calculate PDS #
# ###################################################### #
mallinfo::mallinfo()
mallinfo::malloc.trim()

options(connectionObserver = NULL)
library(Matrix)
library( org.Mm.eg.db)
library(ggplot2 )
library(reshape2)
library(plyr )
library(Seurat )
library(viridis )
library(ggthemes )
library(AUCell )
library(ggpubr )
library(biomaRt )
library(GSEABase )
mart_homo <- useMart( "ensembl",dataset="hsapiens_gene_ensembl" , host="www.ensembl.org")
mart_mouse <- useMart( "ensembl",dataset="mmusculus_gene_ensembl")
tx2gene <- getBM( attributes=c('ensembl_gene_id', 'external_gene_name',"entrezgene_id"),  mart = mart_mouse)
tx2prot <-  getBM( attributes=c("uniprotswissprot","uniprot_gn_id", 'external_gene_name',"mgi_symbol"),  mart = mart_mouse)
#### load functions
# generate damage signature
source("/cellfile/PROJECTS/damage_scores/podocyte_damage_score_TimPadvitski/Output/00_legacy_unassigned/from_F/PODOCYTE_nonWD/DiseaseScore/cell-damage-score/Universal_Damage_Signature_Tim_Paula_Version.R")
# calculate damage score
source("/cellfile/PROJECTS/damage_scores/podocyte_damage_score_TimPadvitski/Code/archive/scripts_PDS/AUCell_script.r")
source("/cellfile/datapublic/tpadvits/PROJECTS/utility.SCRIPTS/usefulRfunc.r")
source( "/cellfile/datapublic/tpadvits/PROJECTS/utility.SCRIPTS/TFA.estimate_functions.r")

setwd("/data/user/tpadvits/PROJECTS/PDS")
DS_all <- read.table( header = T, "DS_all.20.09.2023.tsv")

### load  KFO and GSE146912 podocyte data
  {
  listPodo <- c("/cellfile/PROJECTS/podocyte_snRNAseq_models/snrnaseq_processing_TimPadvitski/Output/11_Nphs2mut/from_WD/Seurat/decontX.allcells_Seur_podo.rda",
                "/cellfile/PROJECTS/podocyte_snRNAseq_models/snrnaseq_processing_TimPadvitski/Output/10_Wt1hetdel/from_WD/Seurat/decontX.allcells_Seur_podo.rda",
                "/cellfile/PROJECTS/podocyte_snRNAseq_models/snrnaseq_processing_TimPadvitski/Output/14_Cem/from_WD/Seurat/decontX.allcells_Seur_podo.pdss2.rda" ,
                "/cellfile/PROJECTS/podocyte_snRNAseq_models/snrnaseq_processing_TimPadvitski/Output/12_Lmx1b/from_WD/Seurat/decontX.allcells_Seur_podo.rda",
                "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/Public/GSE146912_PDSvalidation/Data/processed/Seurat/decontX.allcells_Seur_podo.btbr.rda",
                "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/Public/GSE146912_PDSvalidation/Data/processed/Seurat/decontX.allcells_Seur_podo.cd2ap.rda",
                "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/Public/GSE146912_PDSvalidation/Data/processed/Seurat/decontX.allcells_Seur_podo.doxo.rda",
                "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/Public/GSE146912_PDSvalidation/Data/processed/Seurat/decontX.allcells_Seur_podo.nephritD1.rda",
                "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/Public/GSE146912_PDSvalidation/Data/processed/Seurat/decontX.allcells_Seur_podo.nephritD5.rda")
    listSCSN <- lapply(seq(listPodo), function(ii){
    datt <- readRDS(listPodo[ii])
  })
  names(listSCSN )<- c("Nphs2","Wt1","Pdss2","Lmx1b","btbr","cd2ap","doxo","nephr.D1","nephr.D5")
  listSCSN$Nphs2$gtype[listSCSN$Nphs2$sample=="143489"] <- "wt"
  listSCSN$Nphs2$group <- paste(listSCSN$Nphs2$gtype, listSCSN$Nphs2$age, sep = "_")
  listSCSN$Lmx1b$group <- paste(listSCSN$Lmx1b$gtype, listSCSN$Lmx1b$age, sep = "_")
  
  # unificate control names
  listSCSN <- lapply( seq(listSCSN), function(ii){
    newSeu <- listSCSN[[ii]]
    newSeu$gtypeDE <- ifelse( newSeu$gtype %in% c( "wt", "CD2AP_WT","Normal" ,"BTBR_ob_plus" ) ,
                              "control", "experimental")
    return(newSeu)
  })
  names(listSCSN )<- c("Nphs2","Wt1","Pdss2","Lmx1b","btbr","cd2ap","doxo","nephr.D1","nephr.D5")
  saveRDS( listSCSN , file="/cellfile/PROJECTS/damage_scores/podocyte_damage_score_TimPadvitski/Output/00_legacy_unassigned/SCSN_RNAseq/listSCSN.22.12.23.rda")
  
  }

### get list of genes detected in all experiments
  {
  allPodoGenes.inter <- Reduce(intersect, lapply( seq(listSCSN), function(ii)
    {
    expr <-  if( !is.null(listSCSN[[ii]]@assays$RNA) ) {
      as.matrix( listSCSN[[ii]]@assays$RNA@data )
    }  else as.matrix( listSCSN[[ii]]@assays$data@data )
    
    expr <- expr[ rowSums( round(expr) > 0 ) > ncol(expr)*0.01 , ]

    rownames(expr)
  } ) )
  
  allPodoGenes.uni <- Reduce( union,  
                              lapply( seq(listSCSN) ,
                                      function(ii){
                                        # all genes expressed in 1%of cells in at least one group
                                        iind <- rowSums( 
                                          sapply( unique(listSCSN[[ii]]$group), 
                                                  function(XX){
                                                    rowSums(round(listSCSN[[ii]]@assays$RNA@counts[,listSCSN[[ii]]$group==XX])>0) > 
                                                      0.05*ncol(listSCSN[[ii]]@assays$RNA@counts[,listSCSN[[ii]]$group==XX])
                                                  }) )
                                        print( summary( iind>0) )
                                        rownames(listSCSN[[ii]])[iind  >0 ]
                                      }))
  
  allPodoGenes_mean <- Reduce(cbind, lapply( seq(listSCSN), function(ii)
  {
    expr <-  if( !is.null(listKFO[[ii]]@assays$RNA) ) {
      as.matrix( listKFO[[ii]]@assays$RNA@data )
    }  else as.matrix( listKFO[[ii]]@assays$data@data )
    
    expr <- rowMeans( expr[ allPodoGenes , ] )
    
    return(expr)
  } ) )
  colnames(allPodoGenes_mean) <- names(listSCSN)
  
  saveRDS(allPodoGenes.uni, file="/cellfile/PROJECTS/damage_scores/podocyte_damage_score_TimPadvitski/Output/00_legacy_unassigned/from_F/PODOCYTE_nonWD/DiseaseScore/PDS_pseudotime/SCSN_allPodoGenes.uni.rda")
  
  saveRDS(allPodoGenes, file="/cellfile/PROJECTS/damage_scores/podocyte_damage_score_TimPadvitski/Output/00_legacy_unassigned/from_F/PODOCYTE_nonWD/DiseaseScore/PDS_pseudotime/SCSN_allPodoGenes.rda")
  saveRDS(allPodoGenes_mean, file="/cellfile/PROJECTS/damage_scores/podocyte_damage_score_TimPadvitski/Output/00_legacy_unassigned/from_F/PODOCYTE_nonWD/DiseaseScore/PDS_pseudotime/SCSN_allPodoGenes_mean.rda")
  
  write.table(allPodoGenes_mean, sep = "\t", quote = F, col.names = NA, file="/cellfile/PROJECTS/damage_scores/podocyte_damage_score_TimPadvitski/Output/00_legacy_unassigned/from_F/PODOCYTE_nonWD/DiseaseScore/PDS_pseudotime/SCSN_allPodoGenes_mean.tsv")
}

### downsample to 1 K per study
  {
  
  listSCSN_1K <- lapply( seq( listSCSN ), function(ii)
  {
    print(ii)
    newSeu <- listSCSN[[ii]]
    
   
    # balance samples
    Idents(newSeu)<- newSeu$sample
    newSeu <- subset( newSeu , downsample=1000 )
    
    # balance groups
    Idents(newSeu)<- newSeu$group
    newSeu <- subset( newSeu , downsample=1000 )
    
    
    #balance gtypes
    Idents(newSeu)<- newSeu$gtypeDE
    newSeu <- subset( newSeu, downsample=1000 )
    
    
    return(newSeu)
  })
  
  names(listSCSN_1K ) <- c("Nphs2","Wt1","Pdss2","Lmx1b","btbr","cd2ap","doxo","nephr.D1","nephr.D5")
    # # add  control to Adriamicin
  # SCSN_PDSlist_1K[[6]] <- merge(SCSN_PDSlist_1K[[6]], subset( SCSN_PDSlist_1K[[7]], subset= gtypeDE=="control"))
  
  
  saveRDS( listSCSN_1K , file="/cellfile/PROJECTS/damage_scores/podocyte_damage_score_TimPadvitski/Output/00_legacy_unassigned/from_F/PODOCYTES/disease.score/listSCSN_1K.22.12.23.rda")
  
  #### calculate PDS
  {
    listSCSN_1K.PDS <- lapply( seq(listSCSN_1K), function(ii)
    {
      print(ii)
      newSeu <- listSCSN_1K[[ii]]
      
      
      ## exclude genes non.expressed in more than certain percentage of cells
      expr <- newSeu@assays$RNA@counts
      expr <- expr[ rowSums( round(expr) > 0 ) > 0 , ]
      
      # calculate damage signatures
      set.seed(42)
      newSeu@meta.data$PDS <- DS_calc.func( exprMatrices = expr , 
                                            DSignature = DS_all , 
                                            ntop = 42 ,wghtd = T,
                                            ceilThrsh =  0.05 )
      # adjust ceilThrsh based on a total number of genes in the matrix!
      # the top should include ~ 1K genes
      
      
      return(newSeu)
    })
    
    names(listSCSN_1K.PDS) <- names(listSCSN_1K)
    listSCSN_1K.PDS_merged <- merge( listSCSN_1K.PDS[[1]], listSCSN_1K.PDS[2:length(listSCSN_1K.PDS)],
                                     add.cell.ids=  names(listSCSN_1K) )
    listSCSN_1K.PDS_merged$experiment <- sub("_.*","",colnames(listSCSN_1K.PDS_merged))
    
    saveRDS( listSCSN_1K.PDS , file="/cellfile/PROJECTS/damage_scores/podocyte_damage_score_TimPadvitski/Output/00_legacy_unassigned/from_F/PODOCYTES/disease.score/listSCSN_1K.22.12.23.rda")
    
    # plot
    VlnPlot( listSCSN_1K.PDS_merged, features = c("PDS","PDS.015" ,"PDS.03"), group.by = "experiment" , 
             split.by = "gtypeDE")

  }
  
  
}

### calculate PDS for all cells
  {

  ### calculate PDS for podocytes
  # # load PDS signatures
  # cv.lasso_listALL <- readRDS( file="/cellfile/PROJECTS/damage_scores/podocyte_damage_score_TimPadvitski/Output/00_legacy_unassigned/from_F/PODOCYTE_nonWD/DiseaseScore/Development_Validation/ML/cv.lasso_listALL.rda")
  # geneSS<- cv.lasso_listALL[c(2,7)]
  
  # dont subsample, only balance genotypes
    
  listSCSN.PDS <- lapply( seq(listSCSN), function(ii){
    print(ii)
    newSeu <- listSCSN[[ii]]


    ## exclude genes non.expressed in more than certain percentage of cells
    expr <- newSeu@assays$RNA@counts
    expr <- expr[ rowSums( round(expr) > 0 ) > 0 , ]

    # calculate damage signatures
    newSeu@meta.data$PDS <- DS_calc.func( exprMatrices = expr , 
                           DSignature = DS_all , 
                           ntop = 42 ,wghtd = T,
                           ceilThrsh =  0.05 )
      
    # adjust ceilThrsh based on a total number of genes in the matrix!
    # the top should include ~ 1K genes

  return(newSeu)
})

  # names(listSCSN.PDS) <- names(listSCSN )
  # # combine nephritis and  to Adriamicin
  # listSCSN.PDS[[6]] <- merge( listSCSN.PDS[[6]], subset( listSCSN.PDS[[7]],
  #                                                        subset=gtypeDE=="control"))
  ## save
  saveRDS( listSCSN.PDS , file="/cellfile/PROJECTS/damage_scores/podocyte_damage_score_TimPadvitski/Output/00_legacy_unassigned/from_F/PODOCYTES/disease.score/listSCSN.PDSlist_29.09.2023.rda")
  
 

}

### list of all mouse TFs
TFmouse_MAT <- read.table( header = T, sep = "\t", fill = T , "/media/tim_nevelsk/WD_tim/ANNOTATIONS/CISBP/Mus_musculus_2021_01_19_5 55_pm/TF_Information.txt")
TFmouse_all <- unique( TFmouse_MAT$TF_Name)


### select TFs based on expression in podocytes
  {
  
  ### get summary stat for each TF for each study wt and ctrl
  SCSNdata_TF_stat <- lapply( seq(listSCSN), function(ii)
  {
    print(ii)
    datt <- listSCSN[[ii]]@assays$RNA@counts

    ctrl <- listSCSN[[ii]]@assays$RNA@counts[ , listSCSN[[ii]]$gtypeDE=="control"]
    xprmnt <- listSCSN[[ii]]@assays$RNA@counts[ , listSCSN[[ii]]$gtypeDE!="control"]
    
    # filternonexpressed genes
    expr.ctrl <- tryCatch( ctrl[ rowSums(ctrl>0)>ncol(ctrl)*0.01 &
                                   rownames( ctrl) %in% TFmouse_all, ], error = function(e) NA)
    expr.xprmnt <- tryCatch( xprmnt[ rowSums(xprmnt>0)>ncol(xprmnt)*0.01 & 
                                       rownames( xprmnt) %in% TFmouse_all, ], error = function(e) NA)
    
    # subset expression to 513 TFs
    if(nrow(expr.ctrl>0)) XX1 <- t( apply( expr.ctrl , 1, summary) ) else XX1 <-NA
    if(nrow(expr.xprmnt>0)) XX2 <- t( apply( expr.xprmnt , 1, summary) )  else XX2 <-NA
    
    # print(dim(XX))
    return(list(XX1,XX2))
  } )
  SCSNdata_TF_stat <- Reduce( c , SCSNdata_TF_stat )
  SCSNdata_TF_stat <- SCSNdata_TF_stat[!is.na(SCSNdata_TF_stat)]
  ### select TFs for each study separately using mean and median, then take uniom
  
  # median above 0 and/or mean greater than 1
  TF_MeanMed <- Reduce( union , lapply( seq(SCSNdata_TF_stat) , function(ii){
    union( rownames(SCSNdata_TF_stat[[ii]])[which( SCSNdata_TF_stat[[ii]][,"Median"]>0)] , 
           rownames(SCSNdata_TF_stat[[ii]])[which( SCSNdata_TF_stat[[ii]][,"Mean"]>1)])
    
  }))
  TF_MeanMed <- c(TF_MeanMed, "Zfp423")
  saveRDS(TF_MeanMed, file="/cellfile/PROJECTS/damage_scores/podocyte_damage_score_TimPadvitski/Output/00_legacy_unassigned/from_F/PODOCYTE_nonWD/DiseaseScore/PDS_pseudotime/TFreg/TFA/TFs_SNSC.podo.MeanMed.08.01.24.rda")

  ### make a meme file with motifs 
  # read all meme mouse motifs
  mouseTFmm <- universalmotif::read_meme("/media/tim_nevelsk/WD_tim/ANNOTATIONS/CISBP/CISBP2core_mouse28.03.20.meme",
                            skip =0)
  # give names
  nname <- sapply( mouseTFmm, function(x) x["name"])
  names(mouseTFmm) <-sub( ".*__", "", nname)
  # capture all motifs expressed in podocytes (some motifs shared between TFs)
  mouseTFmm_sep <- strsplit( names(mouseTFmm) , ":")
  names(mouseTFmm_sep) <- names(mouseTFmm)
  mouseTFmm_sep <- setNames( rep(names(mouseTFmm_sep), lengths(mouseTFmm_sep)),
                             unlist(mouseTFmm_sep, use.names=F))
  nnames <- unique( mouseTFmm_sep[ grep( paste(TF_MeanMed,collapse="$|^") ,  names(mouseTFmm_sep) , ignore.case = F)] )
  mouseTFmm.podo <- mouseTFmm[ nnames ] 
  universalmotif::write_meme( mouseTFmm.podo, overwrite = TRUE,
                              "/media/tim_nevelsk/WD_tim/ANNOTATIONS/CISBP/CISBP2mouse_podoTF_08.01.24.meme")
  

  
  }

# ### add SCT corrected counts
# library( Seurat)
# listSCSN.PDS.sct <- listSCSN.PDS
# listSCSN.PDS.sct <- lapply( listSCSN.PDS.sct , SCTransform , 
#                             method = "glmGamPoi" , 
#                             residual.features = allPodoGenes )
# 
# saveRDS(listSCSN.PDS.sct , file="/cellfile/PROJECTS/damage_scores/podocyte_damage_score_TimPadvitski/Output/00_legacy_unassigned/from_F/PODOCYTES/disease.score/listSCSN.PDSlist.SCT_23.12.2022.rda" )

# ### study correlation of mean counts with gene lvlvs
#   {
#   ### snRNAseq KFO data
#   { 
#     library(Seurat)
#     ### load data after decontX
#     X.decontX <- readRDS("/cellfile/PROJECTS/podocyte_snRNAseq_models/snrnaseq_processing_TimPadvitski/Output/10_Wt1hetdel/from_WD/Seurat/decontX.allcells_Seur_podo.rda")
#     X.decontX <- NormalizeData(X.decontX)
#     # extract sample names
#     ccnames1 <- gsub(".*_|-.*" , "", colnames(X.decontX) )
#     ccnames11 <- gsub("_.*|SID" , "", colnames(X.decontX) )
#     
#     ### load data without ambient RNA filtering
#     X <- readRDS("/cellfile/PROJECTS/podocyte_snRNAseq_models/snrnaseq_processing_TimPadvitski/Output/10_Wt1hetdel/from_WD/Seurat/snRNAseq_Wt1hd_premRNA.4w12w25w.Seurat.Norm.rda")
#     # celect podocytes of 12-24 weeks age
#     ccnames22 <- sub("__.*" , "", colnames(X) )
#     X <- subset( X , cells= colnames(X)[ ccnames22 %in% ccnames11] ) 
#     ccnames2 <- sub(".*_" , "", colnames(X) )
#     X <- subset( X , cells = colnames(X)[ ccnames2 %in% ccnames1] ) 
#     
#     X <- NormalizeData(X)
#     X$sample <- X.decontX$sample[ match( colnames(X), sub( ".*_","" ,colnames(X.decontX)) ) ]
#     
#     ## SCT
#     X.decontX.sct <- X.decontX
#     # X.decontX.sct@assays$RNA@counts <- round(X.decontX.sct@assays$RNA@counts)
#     X.decontX.sct <- SCTransform(X.decontX, method = "glmGamPoi" ,  return.only.var.genes=F)
#     
#     X.sct <- X
#     X.sct <- SCTransform(X , method = "glmGamPoi" , return.only.var.genes=F)
#     
#     
#     gEne <- "Wt1"
#     X.tab <- t(rbind(
#       col.Means=colMeans(X.sct@assays$RNA@counts[ allPodoGenes, ]),
#       col.Sums=colSums(X.sct@assays$RNA@counts[ allPodoGenes, ]),
#       
#       X.sct@assays$RNA@counts[gEne,],
#       X.sct@assays$RNA@data[gEne,] ,
#       X.sct@assays$SCT@counts[gEne,],
#       X.sct@assays$SCT@scale.data[gEne,]))
#     
#     X.decontX.tab <- t(rbind(
#       col.Means=colMeans(X.decontX.sct@assays$RNA@counts[ allPodoGenes, ]),
#       col.Sums=colSums(X.decontX.sct@assays$RNA@counts[ allPodoGenes, ]),
#       
#       X.decontX.sct@assays$RNA@counts[gEne,],
#       X.decontX.sct@assays$RNA@data[gEne,] ,
#       X.decontX.sct@assays$SCT@counts[gEne,],
#       X.decontX.sct@assays$SCT@scale.data[gEne,]))
#     
#     colnames(X.tab) <- colnames(X.decontX.tab) <- 
#       c( "colMeans_RNA.counts", "colSums_RNA.counts", 
#          "RNA.counts" , "RNA.data" ,"SCT.counts" , "SCT.scaleData" )
#     
#     id <- 2
#     xx <- rbind( cor( X.tab[ X.sct$gtype=="wtype",id], 
#                       X.tab[ X.sct$gtype=="wtype",], method="spearman"),
#                  cor( X.tab[ X.sct$gtype!="wtype",id], 
#                       X.tab[ X.sct$gtype!="wtype",], method="spearman"),
#                  cor( X.tab[ X.sct$orig.ident =="121173",id], 
#                       X.tab[ X.sct$orig.ident=="121173",], method="spearman"),
#                  cor( X.tab[ X.sct$orig.ident=="121174",id], 
#                       X.tab[ X.sct$orig.ident=="121174",], method="spearman"),
#                  cor( X.tab[ X.sct$orig.ident=="121175",id], 
#                       X.tab[ X.sct$orig.ident=="121175",], method="spearman"),
#                  cor( X.tab[ X.sct$orig.ident =="121176",id], 
#                       X.tab[ X.sct$orig.ident=="121176",], method="spearman"),
#                  cor( X.tab[ X.sct$orig.ident=="121177",id], 
#                       X.tab[ X.sct$orig.ident=="121177",], method="spearman"),
#                  cor( X.tab[ X.sct$orig.ident=="121178",id], 
#                       X.tab[ X.sct$orig.ident=="121178",], method="spearman"),
#                  cor( X.decontX.tab[ X.decontX.sct$gtype=="wt",id], 
#                       X.decontX.tab[ X.decontX.sct$gtype=="wt",], method="spearman"),
#                  cor( X.decontX.tab[ X.decontX.sct$gtype!="wt",id], 
#                       X.decontX.tab[ X.decontX.sct$gtype!="wt",], method="spearman"),
#                  cor( X.decontX.tab[ X.decontX.sct$sample=="SID121173",id], 
#                       X.decontX.tab[ X.decontX.sct$sample=="SID121173",], method="spearman"),
#                  cor( X.decontX.tab[ X.decontX.sct$sample=="SID121174",id], 
#                       X.decontX.tab[ X.decontX.sct$sample=="SID121174",], method="spearman"),
#                  cor( X.decontX.tab[ X.decontX.sct$sample=="SID121175",id], 
#                       X.decontX.tab[ X.decontX.sct$sample=="SID121175",], method="spearman"),
#                  cor( X.decontX.tab[ X.decontX.sct$sample=="SID121176",id], 
#                       X.decontX.tab[ X.decontX.sct$sample=="SID121176",], method="spearman"),
#                  cor( X.decontX.tab[ X.decontX.sct$sample=="SID121177",id], 
#                       X.decontX.tab[ X.decontX.sct$sample=="SID121177",], method="spearman"),
#                  cor( X.decontX.tab[ X.decontX.sct$sample=="SID121178",id], 
#                       X.decontX.tab[ X.decontX.sct$sample=="SID121178",], method="spearman")
#     )
#     rownames(xx) <- c( "noFilt_ctrl", "noFilt_exprmnt", 
#                        "noFilt_121173" , "noFilt_121174" ,"noFilt_121175" ,
#                        "noFilt_121176" , "noFilt_121177" ,"noFilt_121178" ,
#                        "decontX_ctrl", "decontX_exprmnt",
#                        "decontX_121173" , "decontX_121174" ,"decontX_121175" ,
#                        "decontX_121176" , "decontX_121177" ,"decontX_121178" )
#     library(gplots)
#     heatmap.2(xx[,-c(1:2)], col = gplots::redblue( 15 ), margins = c(12,12), trace = "none",
#               cellnote= round( xx[,-c(1:2)], digits = 2), notecex = 1.5, notecol = "black", 
#               Rowv = F, Colv = F, main = "Spearman cor. of Nphs1 and the mean expression,\ndiffernt transformations of Wt1het.del data")
#     
#     
#   }
#   
#   ### scRNAseq GSE146912 data
#   { 
#     library(Seurat)
#     X.decontX <- readRDS("/cellfile/PROJECTS/podocyte_snRNAseq_models/snrnaseq_processing_TimPadvitski/Output/10_Wt1hetdel/from_WD/Seurat/decontX.allcells_Seur_podo.rda")
#     X.decontX <- NormalizeData(X.decontX)
#     ccnames1 <- gsub(".*_|-.*" , "", colnames(X.decontX) )
#     ccnames11 <- gsub("_.*|SID" , "", colnames(X.decontX) )
#     
#     # 
#     X <- readRDS("/cellfile/PROJECTS/podocyte_snRNAseq_models/snrnaseq_processing_TimPadvitski/Output/10_Wt1hetdel/from_WD/Seurat/snRNAseq_Wt1hd_premRNA.4w12w25w.Seurat.Norm.rda")
#     ccnames22 <- sub("__.*" , "", colnames(X) )
#     X <- subset( X , cells= colnames(X)[ ccnames22 %in% ccnames11] ) 
#     
#     ccnames2 <- sub(".*_" , "", colnames(X) )
#     X <- subset( X , cells = colnames(X)[ ccnames2 %in% ccnames1] ) 
#     
#     X <- NormalizeData(X)
#     X$sample <- X.decontX$sample[ match( colnames(X), sub( ".*_","" ,colnames(X.decontX)) ) ]
#     
#     ## SCT
#     X.decontX.sct <- X.decontX
#     # X.decontX.sct@assays$RNA@counts <- round(X.decontX.sct@assays$RNA@counts)
#     X.decontX.sct <- SCTransform(X.decontX, method = "glmGamPoi" ,  return.only.var.genes=F)
#     
#     X.sct <- X
#     X.sct <- SCTransform(X , method = "glmGamPoi" , return.only.var.genes=F)
#     
#     
#     gEne <- "Wt1"
#     X.tab <- t(rbind(
#       col.Means=colMeans(X.sct@assays$RNA@counts[ allPodoGenes, ]),
#       col.Sums=colSums(X.sct@assays$RNA@counts[ allPodoGenes, ]),
#       
#       X.sct@assays$RNA@counts[gEne,],
#       X.sct@assays$RNA@data[gEne,] ,
#       X.sct@assays$SCT@counts[gEne,],
#       X.sct@assays$SCT@scale.data[gEne,]))
#     
#     X.decontX.tab <- t(rbind(
#       col.Means=colMeans(X.decontX.sct@assays$RNA@counts[ allPodoGenes, ]),
#       col.Sums=colSums(X.decontX.sct@assays$RNA@counts[ allPodoGenes, ]),
#       
#       X.decontX.sct@assays$RNA@counts[gEne,],
#       X.decontX.sct@assays$RNA@data[gEne,] ,
#       X.decontX.sct@assays$SCT@counts[gEne,],
#       X.decontX.sct@assays$SCT@scale.data[gEne,]))
#     
#     colnames(X.tab) <- colnames(X.decontX.tab) <- 
#       c( "colMeans_RNA.counts", "colSums_RNA.counts", 
#          "RNA.counts" , "RNA.data" ,"SCT.counts" , "SCT.scaleData" )
#     
#     id <- 2
#     xx <- rbind( cor( X.tab[ X.sct$gtype=="wtype",id], 
#                       X.tab[ X.sct$gtype=="wtype",], method="spearman"),
#                  cor( X.tab[ X.sct$gtype!="wtype",id], 
#                       X.tab[ X.sct$gtype!="wtype",], method="spearman"),
#                  cor( X.tab[ X.sct$orig.ident =="121173",id], 
#                       X.tab[ X.sct$orig.ident=="121173",], method="spearman"),
#                  cor( X.tab[ X.sct$orig.ident=="121174",id], 
#                       X.tab[ X.sct$orig.ident=="121174",], method="spearman"),
#                  cor( X.tab[ X.sct$orig.ident=="121175",id], 
#                       X.tab[ X.sct$orig.ident=="121175",], method="spearman"),
#                  cor( X.tab[ X.sct$orig.ident =="121176",id], 
#                       X.tab[ X.sct$orig.ident=="121176",], method="spearman"),
#                  cor( X.tab[ X.sct$orig.ident=="121177",id], 
#                       X.tab[ X.sct$orig.ident=="121177",], method="spearman"),
#                  cor( X.tab[ X.sct$orig.ident=="121178",id], 
#                       X.tab[ X.sct$orig.ident=="121178",], method="spearman"),
#                  cor( X.decontX.tab[ X.decontX.sct$gtype=="wt",id], 
#                       X.decontX.tab[ X.decontX.sct$gtype=="wt",], method="spearman"),
#                  cor( X.decontX.tab[ X.decontX.sct$gtype!="wt",id], 
#                       X.decontX.tab[ X.decontX.sct$gtype!="wt",], method="spearman"),
#                  cor( X.decontX.tab[ X.decontX.sct$sample=="SID121173",id], 
#                       X.decontX.tab[ X.decontX.sct$sample=="SID121173",], method="spearman"),
#                  cor( X.decontX.tab[ X.decontX.sct$sample=="SID121174",id], 
#                       X.decontX.tab[ X.decontX.sct$sample=="SID121174",], method="spearman"),
#                  cor( X.decontX.tab[ X.decontX.sct$sample=="SID121175",id], 
#                       X.decontX.tab[ X.decontX.sct$sample=="SID121175",], method="spearman"),
#                  cor( X.decontX.tab[ X.decontX.sct$sample=="SID121176",id], 
#                       X.decontX.tab[ X.decontX.sct$sample=="SID121176",], method="spearman"),
#                  cor( X.decontX.tab[ X.decontX.sct$sample=="SID121177",id], 
#                       X.decontX.tab[ X.decontX.sct$sample=="SID121177",], method="spearman"),
#                  cor( X.decontX.tab[ X.decontX.sct$sample=="SID121178",id], 
#                       X.decontX.tab[ X.decontX.sct$sample=="SID121178",], method="spearman")
#     )
#     rownames(xx) <- c( "noFilt_ctrl", "noFilt_exprmnt", 
#                        "noFilt_121173" , "noFilt_121174" ,"noFilt_121175" ,
#                        "noFilt_121176" , "noFilt_121177" ,"noFilt_121178" ,
#                        "decontX_ctrl", "decontX_exprmnt",
#                        "decontX_121173" , "decontX_121174" ,"decontX_121175" ,
#                        "decontX_121176" , "decontX_121177" ,"decontX_121178" )
#     library(gplots)
#     heatmap.2(xx[,-c(1:2)], col = gplots::redblue( 15 ), margins = c(12,12), trace = "none",
#               cellnote= round( xx[,-c(1:2)], digits = 2), notecex = 1.5, notecol = "black", 
#               Rowv = F, Colv = F, main = "Spearman cor. of Nphs1 and the mean expression,\ndiffernt transformations of Wt1het.del data")
#     
#     
#   }
#  
# }
