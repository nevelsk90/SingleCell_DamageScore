# ################################################ #
# ####### analysis of GSE146912 ####### #
# ################################################ #
library(ggplot2)
library(ggthemes)
library(plyr)
library(dplyr)
library(Matrix)
library(Seurat)
library(biomaRt)
library(SingleCellExperiment)
library(matrixStats)


source( "/media/tpadvits/WD_tim/PROJECTS/myCode/Read10X_STAR.r" )
source( "/media/tpadvits/WD_tim/PROJECTS/myCode/usefulRfunc.r" )

# select 14 podocyte samples
lmat <- list.files(full.names = T , pattern = ".txt.gz" , 
                   "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/raw_tabs/GSE146912_raw_counts", recursive = F)
PodoMarks.mouse <-  c( "Wt1", "Nphs1", "Nphs2", "Synpo","Robo2","Thsd7a","Nebl","Ptpro","Magi2","Podxl")

#### removal of ambient RNA , DecontX ####
### sic! make sure you start with count matrices, not libnorm data
  {
  library(celda)
  
  for( i in seq(lmat)) {
    print(i)
    
    # read raw
    sce.raw <- read.table( gzfile( lmat[i] ), header = T  )
    sce.raw[, "SYMBOL" ] <- make.unique(  sce.raw[, "SYMBOL" ] )
    rownames(sce.raw) <- sce.raw[, "SYMBOL" ]
    sce.raw <-  sce.raw[ , setdiff(colnames(sce.raw), c("SYMBOL","IGIS")) ]
    sce.raw <- as.sparse( sce.raw )
    saveRDS( sce.raw , paste( sub( ".txt.gz", "" , lmat[i] ),
                              "__raw.sparseM.rda", sep = "") )
    
    sce.decont <- decontX( sce.raw )
    saveRDS( sce.decont , paste( dirname(lmat[i]),"/decontX/" ,
                                 sub( ".txt.gz", "" ,basename(lmat[i])),
                                 "__decontX.rda", sep = "") )
    
    rm( sce.decont, sce.raw )
    gc()
  }
  
}


#### basic QCs for raw and decontX data ####
  {

  # read data create subsampled seurat obj
  seuList.sub <- lapply( seq(lmat), function(i){
    print(i)
    cmat <- readRDS( paste( sub( ".txt.gz", "" , lmat[i] ),
                               "__raw.sparseM.rda", sep = "") )
    
    cmat.decontX <- readRDS( paste( dirname(lmat[i]),"/decontX/" ,
                                    sub( ".txt.gz", "" ,basename(lmat[i])),
                                    "__decontX.rda", sep = ""))

    # create seurat objects
    ssc.raw <- CreateSeuratObject( counts = cmat[,
                                                 sample( colnames(cmat), 2000 )] )
    ssc.decontX <- CreateSeuratObject( counts = cmat.decontX$decontXcounts[,
                                                                           sample( colnames(cmat.decontX$decontXcounts), 2000 )] )

    # # pdf( height = 10, width = 12, file =paste( lldir[i] , "/outs/filtered_decontX.",
    #                                            basename(lldir[i]),".decontPlot.pdf", sep = "") )
    # g1 <-  plotDecontXContamination( cmat.decontX  , size = 0.2) + ggtitle( basename(lldir[i]))                          
    # print(g1)
    # dev.off()  
    
    return( list( ssc.raw, ssc.decontX ))
    
    rm( cmat.decontX, cmat.raw , ssc.raw , ssc.decontX )
    gc()
    
  })
  
  # merge 
  sce.raw.Seurat <- merge( seuList.sub[[1]][[1]], y= lapply(seuList.sub[2:length(seuList.sub)] ,`[[` ,1),
                           add.cell.ids= sub( ".txt.gz","", basename(lmat) ) )
  sce.decontX.Seurat <- merge( seuList.sub[[1]][[2]] , y= lapply(seuList.sub[2:length(seuList.sub)] ,`[[` ,2),
                               add.cell.ids= sub( ".txt.gz","", basename(lmat) )  )
  # add labels
  sce.raw.Seurat$sample <- sub( "_HIT.*" , "" , colnames(sce.raw.Seurat))
  sce.decontX.Seurat$sample <- sub( "_HIT.*" , "" , colnames(sce.decontX.Seurat))
  
  # normalise and estimate cell.cycle
  m.s.genes <- fun_homoTOmouse(cc.genes.updated.2019$s.genes)$MUS
  m.g2m.genes <- fun_homoTOmouse(cc.genes.updated.2019$g2m.genes)$MUS
  sce.raw.Seurat <- NormalizeData(sce.raw.Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  sce.decontX.Seurat <- NormalizeData(sce.decontX.Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  
  sce.raw.Seurat <- CellCycleScoring(object = sce.raw.Seurat, g2m.features = m.g2m.genes ,
                                     s.features = m.s.genes )
  sce.decontX.Seurat <- CellCycleScoring(object = sce.decontX.Seurat, g2m.features = m.g2m.genes,
                                         s.features = m.s.genes )
  # add logged counts
  sce.raw.Seurat$nCount.log_RNA <- log(sce.raw.Seurat$nCount_RNA)
  sce.decontX.Seurat$nCount.log_RNA <- log(sce.decontX.Seurat$nCount_RNA)
  
  sce.raw.Seurat <- PercentageFeatureSet(sce.raw.Seurat, "^mt-", col.name = "percent_mito")
  sce.decontX.Seurat <- PercentageFeatureSet(sce.decontX.Seurat, "^mt-", col.name = "percent_mito")
  sce.raw.Seurat <- PercentageFeatureSet(sce.raw.Seurat, "^Rp[sl]", col.name = "percent_ribo")
  sce.decontX.Seurat <- PercentageFeatureSet(sce.decontX.Seurat, "^Rp[sl]", col.name = "percent_ribo")
  
  # add metadata   
  annot_tab <- read.table(sep ="," ,header=T, "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/raw_tabs/SraRunTable.txt")
  # sce.raw.Seurat$age <- annot_tab$Age_weeks[match( sub("SID","", sce.raw.Seurat$sample), annot_tab$CCG_Sample_ID)]
  # sce.raw.Seurat$gtype <- annot_tab$Genotype[match( sub("SID","", sce.raw.Seurat$sample), annot_tab$CCG_Sample_ID)]
  # sce.raw.Seurat$group <- paste( sce.raw.Seurat$gtype, sce.raw.Seurat$age , sep = "_")
  # 
  # sce.decontX.Seurat$age <- annot_tab$Age_weeks[match( sub("SID","", sce.decontX.Seurat$sample), annot_tab$CCG_Sample_ID)]
  # sce.decontX.Seurat$gtype <- annot_tab$Genotype[match( sub("SID","", sce.decontX.Seurat$sample), annot_tab$CCG_Sample_ID)]
  # sce.decontX.Seurat$group <- paste(sce.decontX.Seurat$gtype, sce.decontX.Seurat$age , sep = "_")
  # 
  sce.raw.Seurat$group <- gsub('.{2}$', '', sce.raw.Seurat$sample)
  sce.raw.Seurat$group <- as.factor(  sce.raw.Seurat$group )
  levels( sce.raw.Seurat$group ) <-  c( "Adriamycin_d10__wk12" ,"BTBR_ob_ob__wk12", 
                                      "BTBR_ob_ob__wk21", "BTBR_ob_plus__wk12", "CD2AP_KO__wk3",
                                      "CD2AP_WT__wk3" , "Nephritis_d1__wk12" ,  
                                      "Nephritis_d5__wk12", "Normal__wk12")
  
  sce.raw.Seurat$age <- sub(".*__","", sce.raw.Seurat$group)
  sce.raw.Seurat$gtype <- sub("__.*","", sce.raw.Seurat$group)


  
  # plot
  feats <- c("nFeature_RNA", "nCount.log_RNA", "percent_mito", "percent_ribo" )
  maxx <- c( 10000 , 10, 5, 10 )
  minn <- c( 0, 5, 0, 0 )
  
  ggpl <- lapply( seq(feats), function(i){
    gg1 <- VlnPlot( sce.raw.Seurat , group.by = "sample", 
                    features = feats[i], pt.size = 0, ncol = 1 ) +
      ylim(c(minn[i],  maxx[i])) +NoLegend() 
    
    gg2 <- VlnPlot( sce.decontX.Seurat , group.by = "sample", 
                    features = feats[i], pt.size = 0, ncol = 1) +
      ylim(c(minn[i],  maxx[i])) +NoLegend() 
    
    gg3 <- cowplot::plot_grid( gg1, gg2  , ncol = 2)
    
    return(gg3)
    
  })
  cowplot::plot_grid( plotlist = ggpl , ncol = 1)
  
  
}

#### read raw count matrices, create a seurat obj, select podocytes ####  
{
  ll.rawC <- list.files( full.names = T, pattern = "raw.sparseM.rda" , path = "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/raw_tabs/GSE146912_raw_counts")
  annot.tab <- read.csv("/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/raw_tabs/SraRunTable.txt")
  
  
  # create Seurat objects
  seur.list <- lapply( seq(ll.rawC), function(ii){
    
    expMat <- readRDS(ll.rawC[ii])
    
    print(ii)
    # extract all cells with nonzero expression 
    datt <- expMat[ rowSums( expMat > 0 )> 0 , ]
    dim(datt)
    
    seu <- tryCatch( {
      seu<- CreateSeuratObject( datt ) 
      subset( seu, subset= nFeature_RNA>=250 & nFeature_RNA < 10000 )
    } , 
    error = function(e) NA )
    return(seu)
  })
  names(seur.list) <- c("GSM4409514","GSM4409515","GSM4409520","GSM4409521","GSM4409522",
                        "GSM4409523","GSM4409518","GSM4409519","GSM4409517","GSM4409516",
                        "GSM4409510","GSM4409511","GSM4409512","GSM4409513","GSM4409507")
  
  seur <- merge( seur.list[[1]],
                 y = seur.list[2:length(seur.list)],
                 add.cell.ids = names(seur.list)
  )
  
  
  seur$geoID <- sub("_.*","",colnames(seur))
  
  # filter by mito content 
  seur$percent.mt <- PercentageFeatureSet(object = seur, pattern = "^mt")
  VlnPlot(seur,features = "percent.mt", group.by="geoID")
  VlnPlot(seur,features = "nFeature_RNA", group.by="geoID")
  seur <- subset(seur, subset=percent.mt < 10 )
  saveRDS( seur, file="/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/Public/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/Seurat/GSE146912_allcells_rawSeur.rda") 
  
  
  # cluster and Dimred
  seur <- NormalizeData(seur)
  # ### find variable features
  seur  <- FindVariableFeatures( seur , selection.method = "vst" )
  ### scale the data
  seur <- ScaleData( seur )
  
  ### run PCA
  seur <- RunPCA(seur, features = VariableFeatures(object = seur))
  print( ElbowPlot(seur , ndims = 50 )) # choosing number of PCs to include 
  
  # cluster
  seur <- FindNeighbors( seur, dims = 1:30 )
  seur <- FindClusters( seur, resolution = 0.2 )
  
  # run UMAP
  seur <- RunUMAP( seur, dims = 1:30 )
  seur$PodoMarks <- colMeans(seur@assays$RNA@data[ PodoMarks.mouse,],
                             na.rm = T)
  seur$gtype <- annot.tab$disease_model[ 
    match(  seur$geoID , annot.tab$GEO_Accession)]
  # saveRDS(seur, file="/cellfile/DATASETS/Human/Kidney/Multi_omics/Public/KPMP/Data/from_WD_2026-10/snRNAseq/Seurat/KPMP_snRNAseq.AlbCr.CKD_dimred.27.09.23.rda")
  # subset for plotting
  seur.sub <- subset( seur , downsample= 250 )
  
  # make plots
  # seur <- readRDS("/cellfile/DATASETS/Human/Kidney/Multi_omics/Public/KPMP/Data/from_WD_2026-10/snRNAseq/Seurat/KPMP_snRNAseq.AlbCr_PodoMarks.rda")
  p1 <- DimPlot( seur.sub, label = T, shuffle = T)
  p2 <- DimPlot( seur.sub, shuffle = T,label = T, group.by = "geoID")
  
  p3 <- FeaturePlot( seur.sub, features = "PodoMarks")
  p4 <- DimPlot( seur.sub, shuffle = T, label = T, group.by = "gtype")
  
  cowplot::plot_grid(plotlist = list(p1,p2,p3,p4), nrow = 2)
  
  # find cluster markers
  seur.clustMarks <- FindAllMarkers( seur.sub )
  library(dplyr)
  top10 <- seur.clustMarks %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  # plot
  seur.sub <- ScaleData(seur.sub , features = rownames(seur.sub))
  DoHeatmap( seur.sub, features =top10$gene )
  
  # select podo by cluster
  seur.podo <- subset( seur, id=c(3,7))
  
  saveRDS( seur.podo, file="/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/Public/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/Seurat/GSE146912_podo_rawSeur.rda") 
  
  
  
}

#### read counts adjusted by decontX, create a seurat obj, select podocytes ####  
{
  ll.decontC <- list.files( full.names = T, pattern = "decontX.rda" , path = "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/raw_tabs/GSE146912_raw_counts/decontX/")
  annot.tab <- read.csv("/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/raw_tabs/SraRunTable.txt")
  
  
  # create Seurat ll.decontC
  seur.list <- lapply( seq(ll.decontC), function(ii){
    
    expMat <- readRDS(ll.decontC[ii])$decontXcounts
    
    print(ii)
    # extract all cells with nonzero expression 
    datt <- expMat[ rowSums( expMat > 0 )> 0 , ]
    dim(datt)
    
    seu <- tryCatch( {
      seu<- CreateSeuratObject( datt ) 
      subset( seu, subset= nFeature_RNA>=250 & nFeature_RNA < 10000 )
    } , 
    error = function(e) NA )
    return(seu)
  })
  names(seur.list) <- c("GSM4409514","GSM4409515","GSM4409520","GSM4409521","GSM4409522",
                        "GSM4409523","GSM4409518","GSM4409519","GSM4409517","GSM4409516",
                        "GSM4409510","GSM4409511","GSM4409512","GSM4409513","GSM4409507")
  
  seur <- merge( seur.list[[1]],
                 y = seur.list[2:length(seur.list)],
                 add.cell.ids = names(seur.list)
  )
  
  
  seur$geoID <- sub("_.*","",colnames(seur))
  
  # filter by mito content 
  seur$percent.mt <- PercentageFeatureSet(object = seur, pattern = "^mt")
  VlnPlot(seur,features = "percent.mt", group.by="geoID")
  VlnPlot(seur,features = "nFeature_RNA", group.by="geoID")
  seur <- subset(seur, subset=percent.mt < 10 )
  saveRDS( seur, file="/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/Public/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/Seurat/GSE146912_decontX.allcells_rawSeur.Rdata") 
  
  
  # cluster and Dimred
  seur <- NormalizeData(seur)
  # ### find variable features
  seur  <- FindVariableFeatures( seur , selection.method = "vst" )
  ### scale the data
  seur <- ScaleData( seur )
  
  ### run PCA
  seur <- RunPCA(seur, features = VariableFeatures(object = seur))
  print( ElbowPlot(seur , ndims = 50 )) # choosing number of PCs to include 
  
  # cluster
  seur <- FindNeighbors( seur, dims = 1:30 )
  seur <- FindClusters( seur, resolution = 0.2 )
  
  # run UMAP
  seur <- RunUMAP( seur, dims = 1:30 )
  seur$PodoMarks <- colMeans(seur@assays$RNA@data[ PodoMarks.mouse,],
                             na.rm = T)
  seur$gtype <- annot.tab$disease_model[ 
    match(  seur$geoID , annot.tab$GEO_Accession)]
  # saveRDS(seur, file="/cellfile/DATASETS/Human/Kidney/Multi_omics/Public/KPMP/Data/from_WD_2026-10/snRNAseq/Seurat/KPMP_snRNAseq.AlbCr.CKD_dimred.27.09.23.rda")
  # subset for plotting
  seur.sub <- subset( seur , downsample= 250 )
  
  # make plots
  # seur <- readRDS("/cellfile/DATASETS/Human/Kidney/Multi_omics/Public/KPMP/Data/from_WD_2026-10/snRNAseq/Seurat/KPMP_snRNAseq.AlbCr_PodoMarks.rda")
  p1 <- DimPlot( seur, label = T, shuffle = T)
  # p2 <- DimPlot( seur.sub, shuffle = T,label = T, group.by = "geoID")
  
  p3 <- FeaturePlot( seur, features = "PodoMarks")
  p4 <- DimPlot( seur, shuffle = T, label = T, group.by = "gtype")
  
  cowplot::plot_grid(plotlist = list(p1,p3,p4), nrow = 2)
  
  # find cluster markers
  seur.clustMarks <- FindAllMarkers( seur.sub )
  top10 <- seur.clustMarks %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  # plot
  seur.sub <- ScaleData(seur.sub , features = rownames(seur.sub))
  DoHeatmap( seur.sub, features =top10$gene )
  
  # select podo by cluster
  seur.podo2 <- subset( seur, id=c(2,4))
  
  saveRDS( seur.podo2, file="/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/Public/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/Seurat/GSE146912_decontX.podo_rawSeur.rda") 
  
  
  
}


#### find markers, subset podocytes ####
{
  sce.decontX.Seurat.filt <- readRDS(
           file = "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/Public/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/Seurat/GSE146912_decontX.allcells_Seur.scDblFilt.rda")
  
  ### find markers
  sce.decontX_subsMarks <- FindAllMarkers( sce_subs )
  
  saveRDS( sce.decontX_subsMarks , 
           file="/cellfile/PROJECTS/podocyte_snRNAseq_models/snrnaseq_processing_TimPadvitski/Output/10_Wt1hetdel/from_WD/Seurat/snRNAseq.Wt1het.del_decontX.allcells.scDblFilt_clustMarks.Rdata")
  
  top10.decontX <- sce.decontX_subsMarks %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  
  DoHeatmap( sce_subs, features = top10.decontX$gene) + NoLegend()
  
  ### select podocytes
  sce.decontX_Podo <- subset( sce.decontX.Seurat.filt , idents= c( 2, 4 ))
  
  # save
  # saveRDS( sce.decontX_Podo ,
  #           file = "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/Seurat/decontX.allcells_Seur_podo.rda" )
  sce.decontX_Podo <- readRDS( "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/Seurat/decontX.allcells_Seur_podo.rda")
  sce.decontX_Podo <-  subset( sce.decontX_Podo , percent_mito <= 10 )
  
  # by experiments
  saveRDS(  subset( sce.decontX_Podo , subset=gtype=="Adriamycin_d10" ) ,
           file = "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/Public/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/Seurat/decontX.allcells_Seur_podo.doxo.rda" )
  saveRDS(  subset( sce.decontX_Podo , subset=gtype== "Nephritis_d5" )  ,
           file = "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/Public/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/Seurat/decontX.allcells_Seur_podo.nephritD5.rda" )
  saveRDS(  subset( sce.decontX_Podo , subset=gtype %in% c( "Normal" ,"Nephritis_d1" ))  ,
            file = "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/Public/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/Seurat/decontX.allcells_Seur_podo.nephritD1.rda" )
  saveRDS(  subset( sce.decontX_Podo , cells=which( sce.decontX_Podo$batch==3 ) ) ,
           file = "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/Public/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/Seurat/decontX.allcells_Seur_podo.cd2ap.rda" )
  saveRDS( subset( sce.decontX_Podo , cells=which( sce.decontX_Podo$batch==2 ) ) ,
           file = "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/Public/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/Seurat/decontX.allcells_Seur_podo.btbr.rda" )
  
  
  ### save raw podocytes ( no decontX)
  sce_Podo <- readRDS("/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/Public/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/Seurat/GSE146912_podo_rawSeur.rda")
  sce_Podo <-  subset( sce_Podo , subset=percent.mt <= 10 )
  
  saveRDS(  subset( sce_Podo , subset=gtype%in%"doxorubicin" ) ,
            file = "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/Public/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/Seurat/allcells_Seur_podo.doxo.rda" )
  saveRDS(  subset( sce_Podo , subset=gtype== "nephritis d5" )  ,
            file = "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/Public/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/Seurat/allcells_Seur_podo.nephritD5.rda" )
  saveRDS(  subset( sce_Podo , subset=gtype %in% c( "control" ,"nephritis d1" ))  ,
            file = "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/Public/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/Seurat/allcells_Seur_podo.nephritD1.rda" )
  saveRDS(  subset( sce_Podo , subset=gtype %in% c( "BTBR ob/+" ,"BTBR ob/ob" ) ) ,
            file = "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/Public/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/Seurat/allcells_Seur_podo.btbr.rda" )
  saveRDS( subset( sce_Podo , subset=gtype %in% c( "CD2AP KO" ,"CD2AP WT" ) ) ,
           file = "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/Public/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/Seurat/allcells_Seur_podo.cd2ap.rda" )
  
  }

#### calculate and plot PDS ####
{
  source("/cellfile/PROJECTS/damage_scores/podocyte_damage_score_TimPadvitski/Code/archive/scripts_PDS/AUCell_script.r")
  
  ##  load Damage signature
  DS_all <- read.table( header = T, "/cellfile/PROJECTS/damage_scores/podocyte_damage_score_TimPadvitski/Output/00_legacy_unassigned/from_F/PODOCYTE_nonWD/DiseaseScore/DamageSignatures/DS_all.02.06.2022.tsv")
  
  # load podocyte data
  sce.decontX_Podo <-  readRDS( "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/Seurat/decontX.allcells_Seur_podo.rda"  )
  sce.decontX_Podo <-  subset( sce.decontX_Podo , subset=percent_mito <= 10 )
  Idents(sce.decontX_Podo) <- "group" 
  sce.decontX_Podo <- subset( sce.decontX_Podo, downsample=1000)
  
  # no ambient removal
  annot.tab <- read.csv("/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/raw_tabs/SraRunTable.txt")
  
  sce_Podo <- readRDS("/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/Public/scRNA_PDSvalidation_GEO/Data/sc/GSE146912_scRNAseq/Seurat/GSE146912_podo_rawSeur.rda")
  sce_Podo <-  subset( sce_Podo , subset=percent.mt <= 10 )
  sce_Podo$age <- annot.tab$Age[ match(sce_Podo$geoID,  annot.tab$Sample.Name)]
  sce_Podo$group <- paste( sce_Podo$gtype, sce_Podo$age,sep = "_")
  Idents(sce_Podo) <- "group" 
  sce_Podo <- subset( sce_Podo, downsample=1000)
  
  ### for decontX
  expr <- sce.decontX_Podo@assays$RNA@counts
  expr <- expr[ rowSums( round(expr) > 0 ) > ncol(expr)*0.05 , ]
  expr.raw <- sce_Podo@assays$RNA@counts
  expr.raw <- expr.raw[ rowSums( round(expr.raw) > 0 ) > ncol(expr)*0.05  , ]
  
  
  set.seed(42)
  # adjust ceilThrsh based on a total number of genes in the matrix!
  # the top should include ~ 1K genes
  sce.decontX_Podo$PDS.0025 <-   DS_calc.func( exprMatrices = expr , ceilThrsh = 0.025 , 
                                              DSignature= DS_all , ntop=42  )
  
  set.seed(42)
  # adjust ceilThrsh based on a total number of genes in the matrix!
  # the top should include ~ 1K genes
  sce.decontX_Podo$PDS.005 <-   DS_calc.func( exprMatrices = expr , ceilThrsh = 0.05 , 
                                     DSignature= DS_all , ntop=42  )
  set.seed(42)
    sce.decontX_Podo$PDS.015 <-   DS_calc.func( exprMatrices = expr , ceilThrsh = 0.15 , 
                                              DSignature= DS_all , ntop=42  )
  set.seed(42)
    sce.decontX_Podo$PDS.03 <-   DS_calc.func( exprMatrices = expr , ceilThrsh = 0.3 , 
                                              DSignature= DS_all , ntop=42  )
    set.seed(42)
    sce.decontX_Podo$PDS.005all <-   DS_calc.func( exprMatrices = expr , ceilThrsh = 0.05 , 
                                               DSignature= DS_all , ntop=nrow(DS_all)  )
    
    set.seed(42)
    sce_Podo$PDS.0025 <-   DS_calc.func( exprMatrices = expr.raw , ceilThrsh = 0.025 , 
                                        DSignature= DS_all , ntop=42  )
    set.seed(42)
    sce_Podo$PDS.005 <-   DS_calc.func( exprMatrices = expr.raw , ceilThrsh = 0.05 , 
                                                DSignature= DS_all , ntop=42  )
    set.seed(42)
    sce_Podo$PDS.015 <-   DS_calc.func( exprMatrices = expr.raw , ceilThrsh = 0.15 , 
                                                DSignature= DS_all , ntop=42  )
    set.seed(42)
    sce_Podo$PDS.03 <-   DS_calc.func( exprMatrices = expr.raw , ceilThrsh = 0.3 , 
                                               DSignature= DS_all , ntop=42  )
    
    set.seed(42)
    seur.podo2 $PDS.0025 <-   DS_calc.func( exprMatrices = expr , ceilThrsh = 0.025 , 
                                         DSignature= DS_all , ntop=42  )
    set.seed(42)
    seur.podo2 $PDS.005 <-   DS_calc.func( exprMatrices = expr , ceilThrsh = 0.05 , 
                                        DSignature= DS_all , ntop=42 ,wghtd = T )
    set.seed(42)
    seur.podo2 $PDS.015 <-   DS_calc.func( exprMatrices = expr , ceilThrsh = 0.15 , 
                                        DSignature= DS_all , ntop=42  )
    set.seed(42)
    seur.podo2 $PDS.03 <-   DS_calc.func( exprMatrices = expr , ceilThrsh = 0.3 , 
                                       DSignature= DS_all , ntop=42, wghtd = T  )
    
    
    # plot
    VlnPlot( sce.decontX_Podo, features = c("PDS.0025","PDS.005","PDS.015" ,"PDS.03"),
             split.by = "gtype")
    VlnPlot( sce_Podo, features = c("PDS.0025","PDS.005","PDS.015" ,"PDS.03"), 
             split.by = "gtype")
    
    ggl1 <- VlnPlot( sce.decontX_Podo, idents = c("BTBR_ob_ob__wk12","BTBR_ob_ob__wk21","BTBR_ob_plus__wk12") ,
             features = c("PDS.0025","PDS.005","PDS.015" ,"PDS.03"), 
             flip = T, combine = F,pt.size = 0)
    ggl2 <- VlnPlot( sce_Podo, idents =  c( "BTBR ob/ob_12 weeks",   "BTBR ob/ob_21 weeks",    "BTBR ob/+_12 weeks"), 
             features = c("PDS.0025","PDS.005","PDS.015" ,"PDS.03"), 
             group.by = "group" , flip = T, combine = F,pt.size = 0)
    cowplot::plot_grid(plotlist = c(ggl1,ggl2), nrow = 2)
    
    # gg1<- ggplot(data=sce.decontX_Podo@meta.data[
    #   sce.decontX_Podo$group%in% c("BTBR_ob_ob__wk12","BTBR_ob_ob__wk21","BTBR_ob_plus__wk12"), ], 
    #        aes(x=PDS.0025, color=group))+geom_density()
    gg2<- ggplot(data=sce.decontX_Podo@meta.data[
      sce.decontX_Podo$group%in% c("BTBR_ob_ob__wk12","BTBR_ob_ob__wk21","BTBR_ob_plus__wk12"), ], 
      aes(x=PDS.005, color=group))+geom_density()
    gg3<- ggplot(data=sce.decontX_Podo@meta.data[
      sce.decontX_Podo$group%in% c("BTBR_ob_ob__wk12","BTBR_ob_ob__wk21","BTBR_ob_plus__wk12"), ], 
      aes(x=PDS.005all, color=group))+geom_density()
    gg4<- ggplot(data=seur.podo2@meta.data[
      seur.podo2$group%in% c("BTBR_ob_ob__wk12","BTBR_ob_ob__wk21","BTBR_ob_plus__wk12"), ], 
      aes(x=PDS.005, color=group))+geom_density()
    # gg5<- ggplot(data=sce.decontX_Podo@meta.data[
    #   sce.decontX_Podo$group%in% c("BTBR_ob_ob__wk12","BTBR_ob_ob__wk21","BTBR_ob_plus__wk12"), ], 
    #   aes(x=PDS.03, color=group))+geom_density()
    cowplot::plot_grid(plotlist = list(gg2,gg3,gg4), nrow =3)
    
    
  ### correlation betwen PDS and proteinuria
  {
    library(ggpubr)
   
    
    datt <- sce.decontX_Podo.sub@meta.data
    datt <- datt[ which(datt$batch==1 & datt$gtype!= "Adriamycin_d10"), ]
    aggPDS <- aggregate( .~group , data= datt[ , c( "group", "PDS.42.005")], 
                         mean )
    
    aggPDS$AlbCrRatio <- c( 542.5532 ,2.340426,  0.01)
   
    
    ggplot2::ggplot( data = aggPDS, aes( x=PDS.42.005 , y=log(AlbCrRatio))) +
      geom_point(  size=6, color="salmon") +
      theme_bw() +  theme( text = element_text(size = 22)) + 
      geom_smooth(method='lm', se = FALSE) + ggtitle( "GSE146912 nephritis")+ 
      stat_cor( size=8 )  +  geom_text(aes(label = group  ), size=6, position = "dodge")
    
  }
  
  

  

}




