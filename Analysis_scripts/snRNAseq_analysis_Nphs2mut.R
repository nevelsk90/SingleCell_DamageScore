# ############################################### #
# ####### analysis of Nphs2 mut. snRNAseq ####### #
# ############################################### #
# release memory
mallinfo::mallinfo()
mallinfo::malloc.trim()

library(ggplot2)
library(ggthemes)
library(plyr)
library(dplyr)
library(Matrix)
library(Seurat)
library(biomaRt)
library(SingleCellExperiment)
library(matrixStats)

source("/home/tim_nevelsk/PROJECTS/myCode/Read10X_STAR.r")
source( "/home/tim_nevelsk/PROJECTS/myCode/usefulRfunc.r")

mart_mouse = useMart("ensembl",dataset="mmusculus_gene_ensembl", host="www.ensembl.org")
egid2gname <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),  mart = mart_mouse)

lldir <- list.dirs(  "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/STAR/STARpremRNA", recursive = F)[ -c( 13,14, 17)]

#### removal of ambient RNA with decontX ####
  {
  
  library(celda)
 
    
  # select samples with PDS from the same Nphs2mut. experiment ( = 14 samples)
  lldir <- list.dirs(  "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/STAR/STARpremRNA", recursive = F)[ -c( 13,14, 17)]

  # # DecontX
   
  for( i in seq(lldir)) {
    print(i)
    sce <- Read10X_STAR(data.dir = paste( lldir[i] , "/GeneFull/filtered/", 
                                          sep = ""),  gene.column =2)
    # read raw
    sce.raw <- Read10X_STAR(data.dir = paste( lldir[i] , "/GeneFull/raw/", 
                                              sep = ""), gene.column =2)
    
    sce.decont <- decontX( sce, background = sce.raw)
    saveRDS( sce.decont , paste( lldir[i] , "/GeneFull/filtered_decontX.", 
                                 basename(lldir[i]), ".rda", sep = "") )
     rm( sce, sce.decont, sce.raw )
    gc()
  }
     
  

}

#### basic QCs for raw and decontX data ####
  {
   
    seuList_sub <- lapply( seq(lldir), function(i){
      
      print(i)
      
      # read count data
      cmat.raw <- Read10X_STAR(data.dir = paste( lldir[i] , "/GeneFull/filtered/", 
                                                     sep = ""),  gene.column =2)
      # read count data
      cmat.decontX <- readRDS(paste( lldir[i] , "/GeneFull/filtered_decontX.", 
                                     basename(lldir[i]), ".rda", sep = ""))
      
      # create seurat objects
      seu.raw <- subset( CreateSeuratObject( counts = cmat.raw ), 
                         downsample= 2000) 
      seu.decontX <- subset( CreateSeuratObject( counts = cmat.decontX$decontXcounts ), 
                             downsample= 2000)
      
        
      return( list( seu.raw, seu.decontX ) )
      
      rm( cmat.decontX, cmat.raw, seu.raw ,seu.decontX)
      gc()
      
    })
    
    # merge 
    sce.raw.Seurat <- lapply( seuList_sub, `[[` ,1 )
    sce.raw.Seurat <- merge( sce.raw.Seurat[[1]] , y= sce.raw.Seurat[2:length(sce.raw.Seurat)] , 
                                 add.cell.ids= sub(".*_","",sub("_S.*","",basename(lldir)))    )
    
    sce.decontX.Seurat <-  lapply( seuList_sub,  `[[` ,2 )
    sce.decontX.Seurat <- merge( sce.decontX.Seurat[[1]] , 
                                 y= sce.decontX.Seurat[2:length(sce.decontX.Seurat)] , 
                                 add.cell.ids= sub(".*_","",sub("_S.*","",basename(lldir)))    )
    
    # filter nonexpressed genes (less than 25 cells)
    sce.raw.Seurat <- CreateSeuratObject( counts = sce.raw.Seurat@assays$RNA@counts[
      rowSums( sce.raw.Seurat@assays$RNA@counts>0)>25,] )
    sce.decontX.Seurat <- CreateSeuratObject( counts = sce.decontX.Seurat@assays$RNA@counts[
      rowSums( sce.decontX.Seurat@assays$RNA@counts>0)>25,] )
    
    
    ### estimate quality params
    sce.raw.Seurat <- PercentageFeatureSet(sce.raw.Seurat, "^mt-", col.name = "percent_mito")
    sce.raw.Seurat <- PercentageFeatureSet(sce.raw.Seurat, "^Rp[sl]", col.name = "percent_ribo")
    
    sce.decontX.Seurat <- PercentageFeatureSet(sce.decontX.Seurat, "^mt-", col.name = "percent_mito")
    sce.decontX.Seurat <- PercentageFeatureSet(sce.decontX.Seurat, "^Rp[sl]", col.name = "percent_ribo")
    
    
    # normalise and estimate cell.cycle
    m.s.genes <- fun_homoTOmouse(cc.genes.updated.2019$s.genes)$MUS
    m.g2m.genes <- fun_homoTOmouse(cc.genes.updated.2019$g2m.genes)$MUS
    
    sce.raw.Seurat <- NormalizeData(sce.raw.Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
    sce.raw.Seurat <- CellCycleScoring(object = sce.raw.Seurat, g2m.features = m.g2m.genes,
                                           s.features = m.s.genes )
    sce.decontX.Seurat <- NormalizeData(sce.decontX.Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
    sce.decontX.Seurat <- CellCycleScoring(object = sce.decontX.Seurat, g2m.features = m.g2m.genes,
                                           s.features = m.s.genes )
    
    # add logged counts
    sce.raw.Seurat$nCount.log_RNA <- log(sce.raw.Seurat$nCount_RNA)
    sce.decontX.Seurat$nCount.log_RNA <- log(sce.decontX.Seurat$nCount_RNA)
    
    # add metadata    
    annot_tab <- read.table(sep = "\t",header = T, "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/Sample_Names_KFO.csv")
   
    sce.raw.Seurat$sample <- sub( "_.*" , "" , colnames(sce.raw.Seurat))
    sce.decontX.Seurat$sample <- sub( "_.*" , "" , colnames(sce.decontX.Seurat))
    
    sce.raw.Seurat$age <- annot_tab$Age_weeks[match( sce.raw.Seurat$sample, annot_tab$CCG_Sample_ID)]
    sce.raw.Seurat$gtype <- annot_tab$Genotype[match( sce.raw.Seurat$sample, annot_tab$CCG_Sample_ID)]
    sce.raw.Seurat$group <- paste(sce.raw.Seurat$gtype, sce.raw.Seurat$age , sep = "_")
    
    sce.decontX.Seurat$age <- annot_tab$Age_weeks[match( sce.decontX.Seurat$sample, annot_tab$CCG_Sample_ID)]
    sce.decontX.Seurat$gtype <- annot_tab$Genotype[match( sce.decontX.Seurat$sample, annot_tab$CCG_Sample_ID)]
    sce.decontX.Seurat$group <- paste(sce.decontX.Seurat$gtype, sce.decontX.Seurat$age , sep = "_")
    

    ### plot QCs
    feats <- c( "nFeature_RNA", "nCount.log_RNA", "percent_mito", "percent_ribo" )
    maxx <- c(10000, 10, 5, 10)
    minn <- c( 0, 3, 0, 0 )
    ggpl <- lapply( seq(feats), function(i){
     
            gg1 <- VlnPlot( sce.raw.Seurat , group.by = "sample", 
                      features = feats[i], pt.size = 0, ncol = 1) + 
              ylim(c(minn[i],  maxx[i])) + NoLegend() 
            
            gg2 <- VlnPlot( sce.decontX.Seurat , group.by = "sample", 
                      features = feats[i], pt.size = 0, ncol = 1 ) +
              ylim(c(minn[i],  maxx[i])) + NoLegend() 
      
    gg3 <-  cowplot::plot_grid( gg1, gg2 , ncol = 2 )

      return(gg3)
      
    })
    cowplot::plot_grid( plotlist = ggpl , ncol = 1 )
    
     
  }
  
#### merge samples, create an object ####
{
  lldir <- list.dirs(  "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/STAR/STARpremRNA", recursive = F)[ -c( 13,14, 17)]
  
  # merge
  seuList <- lapply( seq(lldir), function(i){
    
    print(i)
    
    # read count data
    cmat.decontX <-  readRDS(paste( lldir[i] , "/GeneFull/filtered_decontX.", 
                                    basename(lldir[i]), ".rda", sep = ""))
    
    # create seurat objects
    seu.decontX <- subset( CreateSeuratObject( counts = cmat.decontX$decontXcounts ))
    
    
    return( seu.decontX ) 
    
    rm( cmat.decontX, seu.decontX)
    gc()
    
  })
  names(seuList) <- sub(".*_","",sub("_S.*","",basename(lldir)))
  
  # make a barplot of total cell number 
  datt <- data.frame(Ncells=sapply( seuList,ncol))
  datt$sample <- sub( "SID", "", rownames(datt))
  ggplot2::ggplot( data=datt, aes(x=sample, y=Ncells)) + 
    geom_bar(stat="identity") +theme_bw()+ 
    theme(text = element_text( size = 20 ), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  # subsample
  seuList <- lapply( seuList , subset , downsample=8000 )
  
  # merge
  sce.decontX.Seurat <- merge( seuList[[1]] , y= seuList[2:length(seuList)] ,
                               add.cell.ids= sub(".*_","",sub("_S.*","",basename(lldir)))   )
  
  # add labels
  sce.decontX.Seurat$sample <- sub( "_.*" , "" , colnames(sce.decontX.Seurat))
  
  # add metadata    
  annot_tab <- read.table(sep = "\t",header = T, "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/Sample_Names_KFO.csv")
  
  sce.decontX.Seurat$age <- annot_tab$Age_weeks[match( sce.decontX.Seurat$sample, annot_tab$CCG_Sample_ID)]
  sce.decontX.Seurat$gtype <- annot_tab$Genotype[match( sce.decontX.Seurat$sample, annot_tab$CCG_Sample_ID)]
  sce.decontX.Seurat$group <- paste(sce.decontX.Seurat$gtype, sce.decontX.Seurat$age , sep = "_")
  
  # save
  saveRDS( sce.decontX.Seurat ,
           file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/Nphs2mut_decontX.allcells_rawSeur.Rdata")
  
  
  
}


#### do the Dim reduction and clustering ####
  {
    sce.decontX.Seurat <- readRDS( file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/Nphs2mut_decontX.allcells_rawSeur.Rdata" )
    
    ### find variable features
    sce.decontX.Seurat <- NormalizeData( sce.decontX.Seurat )
    sce.decontX.Seurat <- FindVariableFeatures( sce.decontX.Seurat , 
                                                selection.method = "vst", 
                                                nfeatures = 1000)
    
    ### scale the data
    sce.decontX.Seurat <- ScaleData(sce.decontX.Seurat)
    
    ### run PCA
    sce.decontX.Seurat <- RunPCA( sce.decontX.Seurat , 
                                  features = VariableFeatures( object = sce.decontX.Seurat ) )
    ElbowPlot( sce.decontX.Seurat , ndims = 50 ) # choosing number of PCs to include 
    
    # DimPlot( sce, reduction = "pca", group.by = "gtype", order =  sample(colnames(sce)) )
    
    # cluster
    sce.decontX.Seurat  <- FindNeighbors(sce.decontX.Seurat , dims = 1:20)
    sce.decontX.Seurat  <- FindClusters(sce.decontX.Seurat , resolution = 0.1)
    
    # 
    # run UMAP
    sce.decontX.Seurat <- RunUMAP(sce.decontX.Seurat, dims = 1:20 )
    
    saveRDS( sce.decontX.Seurat,
          file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/Nphs2mut_decontX.allcells_Seur.Rdata")
    
    
    
   ### QC for clusters
      sce.decontX.Seurat$nCount.log_RNA <- log(sce.decontX.Seurat$nCount_RNA)
      
      sce.decontX.Seurat <- PercentageFeatureSet(sce.decontX.Seurat, "^mt-", col.name = "percent_mito")
      sce.decontX.Seurat <- PercentageFeatureSet(sce.decontX.Seurat, "^Rp[sl]", col.name = "percent_ribo")
      
      feats <- c( "nFeature_RNA", "nCount.log_RNA", "percent_mito", "percent_ribo" )
      maxx <- c(10000, 10, 5, 10)
      minn <- c( 0, 3, 0, 0 )
      ggpl <- lapply( seq(feats), function(i){
        
        
        gg2 <- VlnPlot( sce.decontX.Seurat , group.by = "seurat_clusters", 
                        features = feats[i], pt.size = 0, ncol = 1 ) +
          ylim(c(minn[i],  maxx[i])) + NoLegend() 
        return( gg2 )
      })
        
       cowplot::plot_grid( plotlist = ggpl, ncol = 2 )
   
   
  }

#### find doublets with scDblFinder ####
#### https://bioconductor.org/books/release/OSCA/doublet-detection.html
  {
    
  # datt <- subset( sce.decontX.Seurat , downsample=500 )
    library(scDblFinder)
    library(SingleCellExperiment)

    
      # Setting up the parameters for consistency with denoisePCA();
      # this can be changed depending on your feature selection scheme.
      sce_scExp <- scDblFinder(  as.SingleCellExperiment( sce.decontX.Seurat ) ,
                                clusters=sce.decontX.Seurat$seurat_clusters,
                                samples= "sample", verbose=T )
      

    # save
    saveRDS( sce_scExp@colData ,
             file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/Nphs2mut_decontX.allcells_scDblFinder.rda")

# 
    ### plot doublet Finder
    g1 <- scater::plotColData( sce_scExp, x="seurat_clusters", y="scDblFinder.score",
                               colour_by="seurat_clusters", point_size =0 )
    g2 <- scater::plotUMAP( sce_scExp, colour_by="scDblFinder.class", point_size =0  )
    g3 <- scater::plotUMAP( sce_scExp, colour_by="scDblFinder.score", point_size =0  )
    g4 <- scater::plotUMAP( sce_scExp, text_by="seurat_clusters" ,
                            colour_by="group", point_size =0  )

    cowplot::plot_grid( plotlist = list(g1,g2,g3,g4 ) , ncol = 2)

    
      
}

#### recluster after removing doublets ####
{
  
  sce.decontX.Seurat.filt <- subset( sce.decontX.Seurat , 
                                     cells = which(sce_scExp$scDblFinder.class =="singlet"))
  
  ### find variable features
  sce.decontX.Seurat.filt <- NormalizeData( sce.decontX.Seurat.filt )
  sce.decontX.Seurat.filt <- FindVariableFeatures( sce.decontX.Seurat.filt , 
                                              selection.method = "vst", 
                                              nfeatures = 1000)
  
  ### scale the data
  sce.decontX.Seurat.filt <- ScaleData( sce.decontX.Seurat.filt )
  
  ### run PCA
  sce.decontX.Seurat.filt <- RunPCA( sce.decontX.Seurat.filt , 
                                features = VariableFeatures( 
                                  object = sce.decontX.Seurat.filt ) )
  ElbowPlot( sce.decontX.Seurat.filt , ndims = 50 ) # choosing number of PCs to include 
  
  # DimPlot( sce, reduction = "pca", group.by = "gtype", order =  sample(colnames(sce)) )
  
  # cluster
  sce.decontX.Seurat.filt  <- FindNeighbors(sce.decontX.Seurat.filt , dims = 1:20)
  sce.decontX.Seurat.filt  <- FindClusters(sce.decontX.Seurat.filt , resolution = 0.1)
  
  # 
  # run UMAP
  sce.decontX.Seurat.filt <- RunUMAP(sce.decontX.Seurat.filt, dims = 1:20 )
  
  saveRDS( sce.decontX.Seurat.filt,
           file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/Nphs2mut_decontX.allcells_Seur.scDblFilt.rda")
  
  
  
  # plot of subsampled cells
  set.seed(42)
  sce_subs <- subset( sce.decontX.Seurat.filt ,  downsample = 500  ) 
  
  p1 <- DimPlot( sce_subs , reduction = "umap", shuffle = T , label = T)
  p2 <- DimPlot( sce_subs, reduction = "umap", group.by = "sample", shuffle = T  )
  p3 <- DimPlot( sce_subs, reduction = "umap", group.by = "age", shuffle = T  )
  p4 <- DimPlot( sce_subs, reduction = "umap", group.by = "gtype", shuffle = T  )
  p5 <- FeaturePlot( sce_subs, reduction = "umap", order = F,features = "nCount_RNA" , 
                     max.cutoff = 15000)
  p6 <-   FeaturePlot( sce_subs, reduction = "umap", order = F, features = "Wt1" , 
                       min.cutoff = 0)
  
  cowplot::plot_grid( p1,p2, p3, p4,p5,p6, ncol = 2)
  
  # QC for clusters
  feats <- c( "nFeature_RNA", "nCount.log_RNA", "percent_mito", "percent_ribo" )
  maxx <- c(10000, 10, 5, 10)
  minn <- c( 0, 3, 0, 0 )
  ggpl <- lapply( seq(feats), function(i){
    
    
    gg2 <- VlnPlot( sce.decontX.Seurat.filt , group.by = "seurat_clusters", 
                    features = feats[i], pt.size = 0, ncol = 1 ) +
      ylim(c(minn[i],  maxx[i])) + NoLegend() 
    return( gg2 )
  })
  
  cowplot::plot_grid( plotlist = ggpl, ncol = 2 )
  
  
  ### Violin plot for podocyte marks
  VlnPlot( sce.decontX.Seurat.filt, 
           features = c("Wt1", "Nphs1", "Nphs2", "Podxl", "Magi2",
                        "Mafb","Foxd1","Plekhh2"), 
           pt.size = 0, ncol = 2 ) 
  
  
  # try NGCSC clustering 
  {
    
    sce_subs <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN.PDSlist_23.12.2022.rda")[[1]]
    sce_subs$PDS.42.005 <- listSCSN.PDS[[1]]$PDS.42.005[ match( colnames(sce_subs) , colnames(listSCSN.PDS[[1]] ))]
    
    library(dplyr)
    library(Seurat)
    library( ggplot2 )
    sce_subs@meta.data %>% colnames()
    sce_subs@meta.data <- sce_subs@meta.data %>% select(-RNA_snn_res.0.1)
    
    sce_subs <- FindNeighbors(sce_subs)
    for (n in seq(.01, 1, .05)){
      print(n)
      sce_subs <- FindClusters(sce_subs, resolution = n)
    }
    
    # try NGCSC clustering 
    library(NGCS)
    sce_subs@meta.data %>% colnames()
    NGCS::NGCS(sce_subs@meta.data, prefix = "RNA_snn")
    NGCS_out <- NGCS::NGCS(sce_subs@meta.data, prefix = "RNA_snn",
                           metadata_column_name = group)
    # NGCS_out$plot_node_tree+
    #   facet_wrap("group")
    
    sce_subs$Wt1 <- sce_subs@assays$RNA@data["Wt1",]
    NGCS_out <- NGCS::NGCS(sce_subs@meta.data, prefix = "RNA_snn",
                           metadata_column_name = "Wt1", suggest_cut =T )
    
    sce_subs[["NGCS_clustering"]] <- NGCS_out$nodes_selected$cells_NonGlobalClustering$id 
    saveRDS(sce_subs , file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/nphs2.Podo.NGCS.rda")
    NGCS::NGCS(sce_subs@meta.data, prefix = "RNA_snn",
               metadata_column_name = "PDS.42.005",suggest_cut = T, no_labels=T )
    
    
    ## find markers
    Idents(sce_subs) <- sce_subs$NGCS_clustering
    sce_subs_Marks <- FindAllMarkers( sce_subs )
    
    p1 <- DimPlot( sce_subs , reduction = "umap", shuffle = T , label = T)
    p2 <- DimPlot( sce_subs, reduction = "umap", group.by = "sample", shuffle = T  )
    p3 <- DimPlot( sce_subs, reduction = "umap", group.by = "age", shuffle = T  )
    p4 <- DimPlot( sce_subs, reduction = "umap", group.by = "gtype", shuffle = T  )
    
    
    cowplot::plot_grid( p1, p2, p3, p4, ncol = 2)
  }
  
  
  
}




#### find cluster markers, annotate  cell-types ####
  {

    # prepare markers markers
    marks.hmphr.sn <- read.table(header = T, sep = ";","/media/tim_nevelsk/WD_tim/ANNOTATIONS/kidney/Humphreys2020_mouse_snRNA.csv")
    # marks.hmphr.sc <- read.table(header = T, sep = ";","/media/tim_nevelsk/WD_tim/ANNOTATIONS/kidney/Humphreys2022_mouse_scRNA.csv")
    marks.sztk.sc <- readxl::read_xlsx(sheet = 1, skip = 1,  "/media/tim_nevelsk/WD_tim/ANNOTATIONS/kidney/Susztak2021_mouse.scRNA.xlsx")
    marks.sztk.sc <- lapply( as.list(marks.sztk.sc) , function(X) X <- X[!is.na(X)] )
    names(marks.sztk.sc) <- sub( " ","_" , names(marks.sztk.sc))
    # convert long df to a list
    marks.hmphr.sn <- split(marks.hmphr.sn$gene, marks.hmphr.sn$celltype)
    # marks.hmphr.sc <- split(marks.hmphr.sc$gene, marks.hmphr.sc$celltype)
    podoMarks_db <- c(marks.hmphr.sn, marks.hmphr.sc, marks.sztk.sc)
    
    # podoMarks_db <- lapply( seq(podoMarks_db), function(ii) podoMarks_db[[ii]][1:100])
    names(podoMarks_db) <- c( paste(names(marks.hmphr.sn),"hmphr.sn",sep = "_"),
                              paste(names(marks.hmphr.sc),"hmphr.sc",sep = "_"),
                              paste(names(marks.sztk.sc),"sztk.sc",sep = "_")
    )
    
    # prepare gene sets
    geneSets.hmphr.sn <- GSEABase::GeneSetCollection( lapply( seq(marks.hmphr.sn) ,
                                                              function(ii) GSEABase::GeneSet( marks.hmphr.sn[[ii]] , 
                                                                                              setName= names(marks.hmphr.sn)[ii] ) ))
    geneSets.sztk.sc <- GSEABase::GeneSetCollection( lapply( seq(marks.sztk.sc) ,
                                                             function(ii) GSEABase::GeneSet( marks.sztk.sc[[ii]] , 
                                                                                             setName= names(marks.sztk.sc)[ii] ) ))
    
   
    
    
    ### find markers
    # subsample to balance cell number per sample
    Idents(nphs2_seu) <- "sample"
    set.seed(42)
    sce_subs <- subset( nphs2_seu , downsample= 5000)
    # 
    Idents(sce_subs) <- "seurat_clusters"
    # sce_subs.subs <- subset( sce_subs.subs , downsample= 500)
    sce.decontX_subsMarks <- FindAllMarkers( sce_subs  )
    
    saveRDS( sce.decontX_subsMarks , 
             file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/snRNAseq.Nphs2m_decontX.allcells.scDblFilt_clustMarks.rda")
    
    sce.decontX_subsMarks <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/snRNAseq.Nphs2mut_decontX.allcells.scDblFilt_clustMarks.rda")
  
  # ### select podocytes
  # sce.decontX_Podo <- subset( top3lfcpval.decontX , idents= c(2,3) )
  # 
  # 
  # saveRDS( sce.decontX_Podo ,
  #          file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/decontX.allcells_Seur_podo.rda" )
  
  ### heatmap of markers with dotplot
  # select top 3 per cluster: first select more and then remove non-unique
  top3.decontX <- sce.decontX_subsMarks %>% group_by(cluster) %>%
    slice_min(n = 100, order_by = p_val )
  top3.decontX <- top3.decontX[ !duplicated(top3.decontX$gene),]
  top3.decontX <-  slice_max( top3.decontX , n = 5, order_by = avg_log2FC)
  clustMarks <- top3.decontX$gene
  # # endothelial subtypes
  # clustMarks <- c("Nrp1", "Kdr", "Flrt2","Sncaip","Lef1",
  #                "Slc14a1", "Vwf", "Eln", "Ccl21a") 
  # # stromal subclust
  # clustMarks <- c("Vcan", "Spon1", "Bmpr1b", "Dapk2", "Dcn", "Notch3",
  #                 "Myh11", "Acta2", "Ren1", "Piezo2", "Gata3")
  expr <- cbind.data.frame( t( sce_subs@assays$RNA@data[ clustMarks, ] ) ,
                              cluster= sce_subs$seurat_clusters )
  prcnt <- cbind.data.frame( t( sce_subs@assays$RNA@counts[ clustMarks, ]>0 ) ,
                            cluster= sce_subs$seurat_clusters )
  dat1 <- reshape2::melt( 
    aggregate( .~cluster , data=expr , FUN=mean, na.rm=T))
  dat2 <- aggregate( .~cluster , data=prcnt , FUN=sum, na.rm=T)/table(sce_subs$seurat_clusters) 
  dat2$cluster<- as.factor( 0:13)
  dat2 <- reshape2::melt( dat2)
  toPlot <- cbind.data.frame( dat1 , prcntg = dat2$value )
  
  # save 
  pdf( width =  10 , height = 6 , file="KFO.snRNAseq.Nphs2mut_clstMrk_htmpDtplt.pdf")
  a<-dev.cur()
  png( width =  800 , height = 500 , file="KFO.snRNAseq.Nphs2mut_clstMrk_htmpDtplt.png")
  dev.control("enable")

  ggplot( data=toPlot, aes(x=variable , y = cluster , 
                           color = value   ,
              size = prcntg  )) + 
    geom_point() + scale_y_discrete(limits=rev) +
    scale_color_viridis_c(name = 'mean log(count+ 1)', limits=c(0, 3.5), 
                          na.value = "#FDE725FF") + 
    cowplot::theme_cowplot()+ 
    theme( axis.text.x = element_text(angle = 90, vjust = 0.5,
                                      hjust=1))+ scale_size(range = c(0,10))

  dev.copy(which=a)
  dev.off()
  dev.off()
  
  
  ### relate cluster markers with cell-tpe marks
  sce.decontX_subsMarks <- readRDS(file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/snRNAseq.Nphs2mut_decontX.allcells.scDblFilt_clustMarks.rda" )

  
  nphs2_clustCtype.gsea <- lapply( levels(sce.decontX_subsMarks$cluster),
                                 function(clID)
                                 {
                                   print(clID)
                                   clXmark_LFC <- sce.decontX_subsMarks[
                                     sce.decontX_subsMarks$cluster==clID ,]
                                   clXmark_LFC <- setNames(clXmark_LFC$avg_log2FC , nm = clXmark_LFC$gene)
                                   clXmark_LFC <- clXmark_LFC[ order(-clXmark_LFC)]
                                   X <- fgseaMultilevel( stats= clXmark_LFC , scoreType="pos",
                                                         pathways = podoMarks_db ,nPermSimple = 10000,
                                                         nproc=4 , eps = 1e-100)
                                   return(X)
                                 })
  
  
  nphs2_clustCtype.gsea_top3 <- Reduce( union,  lapply( seq(nphs2_clustCtype.gsea), function(ii){
    X <- nphs2_clustCtype.gsea[[ii]]
    XX <- slice_min(X, n = 3, order_by = pval)
    XX <- XX$pathway[XX$padj<0.1]
    return(XX)
  }))
  
  
  
  nphs2_clustCtype.gsea_p <- Reduce( cbind, lapply( seq(nphs2_clustCtype.gsea), function(ii)
    {
    X <- nphs2_clustCtype.gsea[[ii]]
    X$logP <- scale( -log10(X$pval), center = F)
    
    X <- X[  match( nphs2_clustCtype.gsea_top2,  X$pathway),]
    X$logP
    # X$clust <- levels(sce.decontX_subsMarks$cluster)[ii]
    # return(X)
  }))
  colnames(nphs2_clustCtype.gsea_p) <- levels(sce.decontX_subsMarks$cluster)
  rownames(nphs2_clustCtype.gsea_p)<- nphs2_clustCtype.gsea_top2
  nphs2_clustCtype.gsea_p[is.na(nphs2_clustCtype.gsea_p)]<-0
  
  
  gg0 <-  pheatmap::pheatmap(t (nphs2_clustCtype.gsea_p), 
                             color=viridis::viridis(6),
                             cluster_rows = F,
                             fontsize = 16)
  
  pdf( width =  8 , height = 6 , file="/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/snRNAseq_Nphs2/ctypeAnalysis/KFO.snRNAseq.nphs2m_gsea_ctypeVSclstr.heatmap.pdf")
  a<-dev.cur()
  png( width =  700 , height = 500 , file="/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/snRNAseq_Nphs2/ctypeAnalysis//KFO.snRNAseq.nphs2m_gsea_ctypeVSclstr.heatmap.png")
  dev.control("enable")
  print(gg0)
  
  dev.copy(which=a)
  dev.off()
  dev.off()
  
  
  ### assign cell-types
  ctypes <- setNames( c( "Endo", "Mes", "Pod", "Pod", "PT", "Mes", "PC",
                         "TAL", "Immune", "PEC", "JGA", "IC", "Pod?Mes", "Prolif" ), 
                      nm = levels( nphs2_seu$seurat_clusters ) )
  nphs2_seu@meta.data$ctype <- ctypes[ match( nphs2_seu@meta.data$seurat_clusters ,
                                              names(ctypes))]
  
  saveRDS( nphs2_seu , "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/Nphs2mut_decontX.allcells_Seur.scDblFilt.rda")
  
  
  
  
  }

#### visualise data ####

gg1 <- DimPlot( nphs2_seu, shuffle = T, group.by = "sample" )
gg2 <- DimPlot( nphs2_seu, shuffle = T, group.by = "group" )
gg3 <- DimPlot( nphs2_seu, shuffle = T, group.by = "ctype", label = T )
gg4 <- FeaturePlot( nphs2_seu,  features = "Wt1" )

ggll <-  cowplot::plot_grid(plotlist = list(gg1, gg2, gg3,gg4),nrow = 2)

pdf( width =  12 , height = 10 , file="/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/snRNAseq_Nphs2/KFO.snRNAseq.Nphs2m_UMAPs.pdf")
a<-dev.cur()
png( width =  1000 , height = 800 , file="/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/snRNAseq_Nphs2/KFO.snRNAseq.Nphs2m_UMAPs.png")
dev.control("enable")
print(ggll)

dev.copy(which=a)
dev.off()
dev.off()



#### DE and trajectory for podoytes ####
setwd("/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/snRNAseq_Nphs2/Pseudotime")
# load podocytes
nphs2_seu <-  readRDS( file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/Nphs2mut_decontX.allcells_Seur.scDblFilt.rda")
sce_clust_SCE <- subset( nphs2_seu , subset=gtype=="Nphs2" & ctype="Pod")
set.seed(42)
Idents(sce_clust_SCE)<- sce_clust_SCE$sample
sce_clust_SCE <- subset( sce_clust_SCE , downsample=1000 )
Idents(sce_clust_SCE)<- sce_clust_SCE$group
sce_clust_SCE <- subset( sce_clust_SCE , downsample=1000 )
sce_clust_SCE <- as.SingleCellExperiment(sce_clust_SCE)

# quick pseudotime using clusters
# https://bioconductor.org/books/3.14/OSCA.advanced/trajectory-analysis.html 
pseudo.all <- TSCAN::quickPseudotime( sce_clust_SCE, start = 4,
                                      clusters = sce_clust_SCE$age, 
                                      use.dimred="UMAP",  with.mnn=F)

# order
mnn.pseudo <- rowMeans(pseudo.all$ordering, na.rm=TRUE)

# sce_clust_SCE <- sce_clust_SCE[sample(sce_clust_SCE, size = 1930, replace = F )

# save 
setwd("/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/snRNAseq_Nphs2/Pseudotime/")
# suffle cell number 
pdf( width =  10 , height = 5 , file="KFO.snRNAseq.Nphs2mut_clstMrk_htmpDtplt.pdf")
a<-dev.cur()
png( width =  800 , height = 400 , file="KFO.snRNAseq.Nphs2mut_clstMrk_htmpDtplt.png")
dev.control("enable")

pp1 <- scater::plotUMAP(sce_clust_SCE, colour_by="group", point_size = 1.5   )+
  coord_cartesian(xlim = c(-14,-9),  ylim = c(1,6))
pp2 <- scater::plotUMAP(sce_clust_SCE, colour_by=I(mnn.pseudo),
                        text_by="age", text_colour="red", point_size = 1.5 ) +
  geom_line(data=pseudo.all$connected$UMAP, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge))+
  coord_cartesian(xlim = c(-14,-9),  ylim = c(1,6))
cowplot::plot_grid(pp1 , pp2)

dev.copy(which=a)
dev.off()
dev.off()

### test Changes along a trajectory
XX <- subset(sce_clust_SCE,  rowSums(round(sce_clust_SCE@assays@data$counts)>0)>(1930*0.01))
pseudo <- TSCAN::testPseudotime( XX ,
                                 pseudotime=mnn.pseudo  )
pseudo <- pseudo[order(pseudo$p.value),]	
saveRDS(pseudo , file="Nphs2mut_snRNAseq.DEtest.MSTpsdtime.rda")

### volcanoplot
library(ggrepel)
toPlot <- as.data.frame(pseudo[!is.na(pseudo$p.value) & pseudo$logFC!=0,])
toPlot$gLabel <- rownames(toPlot)
toPlot$gLabel <- ifelse(toPlot$FDR<thrsh & abs(toPlot$logFC)>=0.1, 
                        toPlot$gLabel , "")
toPlot$gColor <-  ifelse(toPlot$FDR<thrsh & abs(toPlot$logFC)>=0.1, 
                         "red" , "black")
toPlot$logFC[toPlot$logFC< -0.5] <- -0.5
toPlot$logFC[toPlot$logFC>0.5] <- 0.5
toPlot$FDR.m.log10 <- -log10(toPlot$FDR)
toPlot$FDR.m.log10[(toPlot$FDR.m.log10)>20] <- 20

pdf( width =  8 , height = 8 , file="KFO.snRNAseq.Nphs2mut.pseudoDE_volcano.pdf")
a<-dev.cur()
png( width =  800 , height = 800 , file="KFO.snRNAseq.Nphs2mut.pseudoDE_volcano.png")
dev.control("enable")
ggplot( toPlot , aes(x= logFC ,y= (FDR.m.log10) ) ) +
  geom_point( alpha=0.3 ,  color=toPlot$gColor) +
  geom_hline(yintercept=-log10(thrsh), linetype="dashed", color = "black")+
  geom_vline(xintercept=-0.1, linetype="dashed", color = "black")+
  geom_vline(xintercept=0.1, linetype="dashed", color = "black")+
  geom_text_repel( aes(label=gLabel), max.overlaps=20, cex=4 ,max.time	=50) +
  # xlab( "mean LFC of GO") + ylab( "-log10( Fisher.pval )") + 
  # ggtitle(paste( "2-way plot of GO terms\n",
  #                dataName, sep = "")) + 
  # coord_cartesian( xlim = c(-0.2, max(GOrobert_LFC$Wt1het.del_4w)) ,
  #                  ylim = c(-0.2, max(GOrobert_LFC$Wt1het.del_12w)) ) +
  theme_bw()+theme( text = element_text(size = 24))

dev.copy(which=a)
dev.off()
dev.off()


### DE analysis of the progression with DEseq2
  {
    ress <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/snRNAseq_Nphs2/Pseudotime/Nphs2mut_snRNAseq_DEseq2.LRT.prog.rda")
    
    ### volcanoplot
    toPlot <- as.data.frame(ress[!is.na(ress$pvalue) & ress$log2FoldChange!=0,])
    toPlot$gLabel <- rownames(toPlot)
    toPlot$gLabel <- ifelse(toPlot$padj<thrsh & 
                              abs(toPlot$log2FoldChange)>=0.5, 
                            toPlot$gLabel , "")
    toPlot$gColor <-  ifelse(toPlot$padj<thrsh & 
                               abs(toPlot$log2FoldChange)>=0.5, 
                             "red" , "black")
    toPlot$log2FoldChange[toPlot$log2FoldChange< -3] <- -3
    toPlot$log2FoldChange[toPlot$log2FoldChange> 3] <- 3
    toPlot$FDR.m.log10 <- -log10(toPlot$padj)
    toPlot$FDR.m.log10[(toPlot$FDR.m.log10)>20] <- 20
    
    pdf( width =  8 , height = 8 , file="KFO.snRNAseq.Nphs2mut.dprogDEseq2_volcano.pdf")
    a<-dev.cur()
    png( width =  800 , height = 800 , file="KFO.snRNAseq.Nphs2mut.dprogDEseq2_volcano.png")
    dev.control("enable")
    ggplot( toPlot , aes(x= log2FoldChange ,y= (FDR.m.log10) ) ) +
      geom_point( alpha=0.3 ,  color=toPlot$gColor) +
      geom_hline(yintercept=-log10(thrsh), linetype="dashed", color = "black")+
      geom_vline(xintercept=-0.5, linetype="dashed", color = "black")+
      geom_vline(xintercept=0.5, linetype="dashed", color = "black")+
      geom_text_repel( aes(label=gLabel), max.overlaps=20, cex=4 ,max.time	=50) +
      # xlab( "mean LFC of GO") + ylab( "-log10( Fisher.pval )") + 
      # ggtitle(paste( "2-way plot of GO terms\n",
      #                dataName, sep = "")) + 
      # coord_cartesian( xlim = c(-0.2, max(GOrobert_LFC$Wt1het.del_4w)) ,
      #                  ylim = c(-0.2, max(GOrobert_LFC$Wt1het.del_12w)) ) +
      theme_bw()+theme( text = element_text(size = 24))
    
    dev.copy(which=a)
    dev.off()
    dev.off()
    
}

#### perform GO and pathway annotation 
  {
  ### source functions and prepare annotations 
  # source("/home/tim_nevelsk/PROJECTS/myCode/func_analysis.R")  
  pseudo <- readRDS(file="Nphs2mut_snRNAseq.DEtest.MSTpsdtime.rda")
  ### use Robert's function for GO analysis
  {
    library(DBI)
    
    # run function on list of DE results

    lapply(seq(), function(ii){
      datt <- X[[ii]]
      datt <- datt[!is.na(datt$p.value),]
      universe= unique( entr2gName$entrezgene_id[match(  rownames(datt), 
                                                         entr2gName$external_gene_name)] ) 
      universe <- as.character( universe[!is.na(universe)] )
      
      
      # prepare gene set, convert ensembleIDs to entrezIDs
      geneset=rownames(datt)[ datt$FDR < thrsh & !is.na(datt$FDR) ]
      geneset <- unique(entr2gName$entrezgene_id[match( geneset,  entr2gName$external_gene_name)])
      geneset <- geneset[!is.na(geneset)]
      geneset <- as.character(geneset)
      print(length(intersect(geneset,colnames(gomatrix))))
      
      # apply Robert's function that given a sparse matrix of GO terms (columns = genes, rows = GO terms)
      # a geneset of interest and a background set of genes (universe)
      # carry out clustering with members diverging by at most cut_max genes, and do enrichment testing.
      # Note, multiplicity adjustment is performed for the representative terms only.
      # RobertGO.10 <- sf.clusterGoByGeneset( gomatrix, geneset, universe, 
      #                                       min.genes=5, cut_max = 10 )
      RobertGO.50 <- sf.clusterGoByGeneset( gomatrix, geneset, universe, 
                                            min.genes=5, cut_max = 50  )
      
      # RobertGO.10$results$Term <- goterms[ match( 
      #   RobertGO.10$results$GO.ID, names(goterms))]
      RobertGO.50$results$Term <- goterms[ match( 
        RobertGO.50$results$GO.ID, names(goterms))]
    })
     
      
    saveRDS(RobertGO.10, file = "KFO.Nphs2mut.pseudoDE_GOrobert.cut10.rda")
    saveRDS(RobertGO.50, file = "KFO.Nphs2mut.pseudoDE_GOrobert.cut50.rda")
    
    ### create a vector with average LFC for all  GO terms
    LFCvec <-Reduce( c, sapply( seq( go2gName ) , function(ii){
        # print(ii)
        ggenes <- go2gName[[ii]]
        
        datt <- pseudo[!is.na(pseudo$p.value),]
        
        if( sum(rownames(datt)%in%ggenes)==0 ) {
          # filter out GOs with less than 3 DE genes detected in bulkRNAseq
          0 } else if( nrow(datt[ rownames(datt)%in%ggenes, ]) < 3 ){
          0 }  else mean(datt$logFC[ rownames(datt)%in%ggenes ], na.rm = T)
      }))
    
    names(LFCvec) <- names(go2gName) 
    
    # plot
    pdf( width =  18 , height = 12 , file="KFO.snRNAseq.Nphs2mut.pseudoDE_2wayPlot.pdf")
    a<-dev.cur()
    png( width =  1200 , height = 800 , file="KFO.snRNAseq.Nphs2mut.pseudoDE_2wayPlot.png")
    dev.control("enable")
    GO_LFC.PVALplot( GOrbrt=RobertGO.50 , LFCvec = LFCvec,
                     dataName = "Nphs2mut. pseudotime DE")
    
    
    dev.copy(which=a)
    dev.off()
    dev.off()
  }
  
  ### GSEA
  {
    library(fgsea)
    
    # run GSEA with shrunk lfc
      datt <- pseudo
      datt <- datt[!is.na(datt$p.value),]
      datt <- setNames( datt[["logFC"]] , rownames(datt))
      datt <- datt[order(-datt)]                              
      fgseaREACT <- fgseaMultilevel(pathways = reactPath_gName, nproc=4,
                                    stats    = datt,eps = 1e-20,
                                    minSize  = 5  ,maxSize = 500 )
      fgseaKEGG <- fgseaMultilevel(pathways = kegg_pathsList, nproc = 4, 
                                   stats    = datt,eps = 1e-20,
                                   minSize  = 5 ,maxSize = 500 )


    
    # save results
    saveRDS( fgseaKEGG , file = "KFO.Nphs2mut.pseudoDE_fgsea.KEGG.rda")
    saveRDS( fgseaREACT , file = "KFO.Nphs2mut.pseudoDE_fgsea.REACT.rda")

  }
  
  ### SPIA 
  {
    pseudoSPIA <- pseudo
    colnames(pseudoSPIA) <- c("log2FoldChange","pvalue", "padj")
    SPIA_bulkRNAseq_strictBkg <-SPIAonList( inputSPIA = pseudoSPIA , 
                                            strict = T , ID.key= "SYMBOL",
                            shrunk = F , geneThrsh= thrsh ) 
    
    
    SPIA_strictBkg_KEGG <- SPIA_bulkRNAseq_strictBkg[[1]]
    SPIA_strictBkg_REACT <- SPIA_bulkRNAseq_strictBkg[[2]]
    
    # save 
    saveRDS( SPIA_strictBkg_KEGG , file="KFO.Nphs2prog_spia.KEGG.rda")
    saveRDS( SPIA_strictBkg_REACT , file="KFO.Nphs2prog_spia.REACT.rda")
  
    # PLOT
    ggplots1 <- spia2way_ggplot(  pathtype="REACTOME", 
                                  SPIA_list = list(SPIA_strictBkg_REACT), 
                                  combinemethod = "fisher", threshold = thrsh )
    ggplots2 <- spia2way_ggplot(  pathtype="KEGG", 
                                  SPIA_list = list(SPIA_strictBkg_KEGG), 
                                  combinemethod = "fisher", threshold = thrsh )
    
    gridExtra::grid.arrange( grobs = list(ggplots1[[1]], ggplots2[[1]]), ncol=2 )
  }
}

### plot annotations
  {
  datTOplot <- fgseaKEGG[fgseaKEGG$pval<thrsh ,]
  
  datTOplot <- rbind( slice_min(datTOplot,n = 10,order_by =NES ) ,
              slice_max(datTOplot,n = 10,order_by =NES ))
  # datTOplot <- rbind( fgseaKEGG ,fgseaREACT )
  # datTOplot$dbase <- c(rep("KEGG", nrow(fgseaKEGG)), 
  #                      rep("REACTOME", nrow(fgseaREACT)))

  # shortn the labels
  datTOplot$pathway <- sapply( datTOplot$pathway , str_trunc, width = 40)
  
  
  ## save
  pdf( width =  12 , height = 8 , file="KFO.snRNAseq.Nphs2mut.pseudoDE_fgseaKEGG.REACT_hist.pdf")
  a<-dev.cur()
  png( width =  1000 , height = 600 , file="KFO.snRNAseq.Nphs2mut.pseudoDE_fgseaKEGG.REACT_hist.png")
  dev.control("enable")
  # plot
  ggplot(datTOplot, aes(x=-log10(pval), 
                              y=reorder( pathway, -log10(pval) ,
                                         FUN = mean, na.rm=T ) , 
                              fill=NES )) + 
    xlab("-log10(Fisher)") + ylab("pathways") +
    # scale_size(range = c(0, 15))+
    scale_fill_distiller( palette ="Spectral" ,
                          limits = c(-(max(abs(datTOplot$NES))),
                                     max(abs(datTOplot$NES))),
                          type = "seq", direction = -1) + 
    geom_bar(stat = "identity") + labs(fill='NES') +
    # facet_grid( rows = vars(dbase), scales = "free_y", space = "free")+
    theme_bw() + theme( text = element_text(size = 22 ),
                        axis.text.x = element_text(size = 14)) 
  dev.copy(which=a)
  dev.off()
  dev.off()
}




#### calculate and plot PDS ####
  {
    source("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/cell-damage-score/AUCell_script.r")
    
    ##  load Damage signature
    DS_all <- read.table( header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/platforms/DamageSignatures/DS_all.02.06.2022.tsv")
    
    # all Nphs2 data
    sce.decontX.Seurat.filt <- readRDS( file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/Nphs2mut_decontX.allcells_Seur.scDblFilt.rda")
    # load podocyte data
    sce.decontX_Podo <-  readRDS( file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/decontX.allcells.scDblFind_Seur_podo.rda")
    
    
    ### select podocytes
    sce.decontX_Podo.sub <- sce.decontX_Podo
    Idents(sce.decontX_Podo.sub) <- sce.decontX_Podo.sub$sample
    sce.decontX_Podo.sub <- subset( sce.decontX_Podo.sub , downsample=500 )
    
    ### calculate for podocytes
      expr <- sce.decontX_Podo.sub@assays$RNA@counts
      expr <- expr[ rowSums(expr)>10 , ]
      
      # adjust ceilThrsh based on a total number of genes in the matrix!
      # the top should include ~ 1K genes
      sce.raw.PDS.005 <-   DS_calc.func( exprMatrices = expr , ceilThrsh = 0.05 , 
                                     DSignature= DS_all , ntop=42  )
   

      sce.decontX_Podo.sub$PDS.42.005 <- sce.raw.PDS.005
      # saveRDS( sce.decontX_Podo , file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/decontX.allcells.scDblFind_Seur_podo.rda")
      
      
    # ### plot distributions 
    #   datt <- sce.decontX_Podo.sub@meta.data
    #   datt$group <- factor(  datt$group, levels = c("wt_4", "wt_8", "wt_12",
    #                                                "Nphs2_4", "Nphs2_6", "Nphs2_8", "Nphs2_12" ))
    #   
    #  ggplot2::ggplot( datt , aes( x=PDS.42.005, color= group ) ) + 
    #     scale_color_colorblind() +  geom_density( size=1.5 ) + theme_bw() + 
    #     labs(x = "podocyte damage score 42mark_0.05thrsh" , title = "Nphs2mut podocytes") +
    #     theme( text = element_text( size = 20 )) + xlim(c(-0.8, 0.2))+
    #     geom_vline( data=ddply( datt , "group", summarise, 
    #                             grp.mean=mean( PDS.42.005 )), aes(xintercept=grp.mean, 
    #                                                        color=group), linetype="dashed" )
      
   
    

    ###plot using vaseplots
    {
      # test differnce between means
      my_comparisons <- compare_means(formula = PDS.50.Podo ~ gtype ,  
                                      data = XX  ) 
      
      # make a vase plot https://github.com/heike/ggboxplots/blob/master/R/stat-hdr.r
      library( ggboxplots)
      library( ggpubr )
      
      na.rm=T
      ggplot2::ggplot( XX , aes( y=PDS.50.Podo, 
                                 x = reorder( orig.ident , PDS.50.Podo, FUN = median) , color=gtype) ) +
        geom_vase(fill = "white", width = 1,lwd=1.5) + 
        scale_color_colorblind() + theme_bw() +
        # coord_cartesian(ylim = quantile( sceHomoPodo@meta.data$PDS.50, c(0.02, 0.98) ))+
        theme( text =  element_text(size=20) ,
               axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        labs( title="Nphs2 mut. snRNAseq data", x ="samples ID, ordered by median" ,
              subtitle= paste( "Wilcoxon test: ",my_comparisons$group1[1],
                               "VS",my_comparisons$group2[1],"p =",my_comparisons$p.format[1])) 
      
    }
    
    ### correlation betwen PDS and proteinuria
     {
       library(ggpubr)
       annot_tab <- read.table(sep = "\t",header = T, "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/Sample_Names_KFO.csv")
       annot_tab$group <- paste(annot_tab$Genotype,annot_tab$Age_weeks,sep = "_")
       
       
       datt <- listSCSN.PDS$Nphs2@meta.data
       
       aggPDS <- aggregate( .~sample , data= datt[ , c( "sample", "PDS.42")], 
                            mean )
       
       aggPDS$AlbCrRatio <- annot_tab$AlbCrRatio[ match( aggPDS$sample , 
                                                         annot_tab$CCG_Sample_ID)]
       aggPDS$group <- annot_tab$group[ match( aggPDS$sample , annot_tab$CCG_Sample_ID)]
       
       gg1 <- ggplot2::ggplot( data = aggPDS, aes( x=PDS.42 , y=log(AlbCrRatio))) +
         geom_point(  size=6, color="salmon") +
         theme_bw() +  theme( text = element_text(size = 22)) + 
         geom_smooth(method='lm', se = FALSE) + ggtitle( "Nphs2mut podocytes")+ 
         stat_cor( size=8 )  +  geom_text(aes(label = sample  ), size=6, position = "dodge")
       
     }
   

     ### add corelation between AlbCr, PDS and percentage of captured podocytes
     {
       XX <- sce.decontX.Seurat.filt
       # balance by sample
       Idents(XX) <- XX$sample
       XX <-  subset( XX , downsample=5000)
       # subset podocytes only
       XX.podo <-  subset( XX , subset=ctype=="Pod")
       
       # calculate 
       # add result to the aggregated sample tab
       aggPDS$podoPer<- table(XX.podo$sample)/5000
       
       gg2<- ggplot2::ggplot( data = aggPDS, aes( x=podoPer , y=log(AlbCrRatio))) +
         geom_point(size=6, color="salmon") +
         theme_bw() +  theme( text = element_text(size = 22)) + 
         geom_smooth(method='lm', se = FALSE)  + ggtitle( "Nphs2mut podocytes") +
         stat_cor( size=8 )  +  geom_text(aes(label = sample ), size=6, position = "dodge")
       
       cowplot::plot_grid( gg1, gg2 , ncol = 2)
       
        }
    
      
      ### plot podocytes
      sce_scExp <- as.SingleCellExperiment( sce.decontX_Podo.sub )
      gg1 <- scater::plotUMAP(sce_scExp, colour_by="PDS.42", point_size =1  )+ 
        xlim(c(-15 , -7.5))+ylim(c(-2.5,7))
      gg2 <- scater::plotUMAP(sce_scExp, colour_by="group", point_size =1  )+ 
        xlim(c(-15 , -7.5))+ylim(c(-2.5,7))
      
      cowplot::plot_grid( gg1, gg2 , ncol = 2)
      
  }

#### DE and trajectory analysis for all clusters ####
  {
  ### load subsampled data, no doublets removal
    sce.decontX_Podo <-  readRDS( file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/decontX.allcells.scDblFind_Seur_podo.rda")
    
  kidney_Nphs2mut_Seurat$Cell_Type <- factor( kidney_Nphs2mut_Seurat$seurat_clusters, levels = levels(Idents(kidney_Nphs2mut_Seurat)), 
                                              labels = c("endothelium","mesangium","podocytes","proximal tubules",
                                                         "(distal) connecting tubules, CDPC","podocytes","doublets: endothelium + mesangium",
                                                         "loop of Henle (ascending)", "immune cells","doublets: mesangium + podocytes",
                                                         "loop of Henle (descending)","doublets: endothelium + podocytes", 
                                                         "juxtaglomerular apparatus", "S3 proximal tubule","collecting duct intercalated cells",
                                                         "cycling cells", "doublets") )
  Idents(kidney_Nphs2mut_Seurat) <- kidney_Nphs2mut_Seurat$Cell_Type
  
  ###
  kidney_Nphs2mut_ctrVSmut.DE <- lapply( levels(kidney_Nphs2mut_Seurat$Cell_Type),
                                         function(ii) {
                                           
                                           w4_DE <- tryCatch( FindMarkers( kidney_Nphs2mut_Seurat, subset.ident = ii , group.by = 'groups',
                                                                           ident.1 = "4-5_Podocin.wt.wt", ident.2 = "4-5_Podocin.R231Q.A286V" ,
                                                                           min.cells.group	= 50 ), error = function(e){ return(NULL) } )
                                           if( !is.null(w4_DE) ) { 
                                             w4_DE$Gene.name <- rownames(w4_DE)
                                             w4_DE$age <- rep("4", nrow(w4_DE))
                                           } 
                                           
                                           w6_DE <- tryCatch( FindMarkers( kidney_Nphs2mut_Seurat, subset.ident = ii , group.by = 'groups',
                                                                           ident.1 ="6_Podocin.wt.wt", ident.2 = "6_Podocin.R231Q.A286V" ,
                                                                           min.cells.group	= 50 ), error = function(e){ return(NULL) } )
                                           if( !is.null(w6_DE) ) { 
                                             w6_DE$Gene.name <- rownames(w6_DE)
                                             w6_DE$age <- rep("6", nrow(w6_DE))
                                           }
                                           
                                           w8_DE <- tryCatch( FindMarkers( kidney_Nphs2mut_Seurat, subset.ident = ii , group.by = 'groups', 
                                                                           ident.1 ="8_Podocin.wt.wt", ident.2 = "8_Podocin.R231Q.A286V",
                                                                           min.cells.group	= 50 ), error = function(e){ return(NULL) } )
                                           if( !is.null(w8_DE) ) { 
                                             w8_DE$Gene.name <- rownames(w8_DE)
                                             w8_DE$age <- rep("8", nrow(w8_DE))}
                                           
                                           w12_DE <- tryCatch( FindMarkers( kidney_Nphs2mut_Seurat, subset.ident = ii , group.by = 'groups',
                                                                            ident.1 = "12_Podocin.wt.wt", ident.2 = "12_Podocin.R231Q.A286V" ,
                                                                            min.cells.group	= 50 ), error = function(e){ return(NULL) } )
                                           if( !is.null(w12_DE) ) { 
                                             w12_DE$Gene.name <- rownames(w12_DE)
                                             w12_DE$age <- rep("12", nrow(w12_DE))
                                           }
                                           
                                           
                                           DEtab <- rbind(w4_DE , w6_DE, w8_DE, w12_DE)
                                           return(DEtab)
                                           
                                         } )
  
  library(scater)
  kidney_Nphs2mut_mutProg.DE <- lapply( levels(kidney_Nphs2mut_Seurat$Cell_Type), 
                                        function(ii){
                                          sce_clust <- subset( kidney_Nphs2mut_Seurat , ident=ii, subset= gtype=="Podocin.R231Q.A286V")
                                          sce_clust_SCE <- as.SingleCellExperiment(sce_clust)
                                          
                                          pseudo.all <- TSCAN::quickPseudotime( sce_clust_SCE,
                                                                                clusters = sce_clust_SCE$age, 
                                                                                use.dimred="UMAP",  with.mnn=F)
                                          
                                          mnn.pseudo <- rowMeans(pseudo.all$ordering, na.rm=TRUE)
                                          
                                          ffile <- paste("/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/snRNAseq_podo/Seurat/premRNA/DE/",
                                                         "Nphs2_FSGStraj_clust_",ii, ".pdf",sep = "")
                                          cairo_pdf( ffile, height=5, width=10, fallback_resolution=1200 )
                                          
                                          pp1 <- scater::plotUMAP(sce_clust_SCE, colour_by="groups", point_size = 1   )
                                          pp2 <- scater::plotUMAP(sce_clust_SCE, colour_by=I(mnn.pseudo), text_by="age", text_colour="red") +
                                            geom_line(data=pseudo.all$connected$UMAP, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge))
                                          
                                          
                                          
                                          cowplot::plot_grid(pp1 , pp2)
                                          dev.off()
                                          
                                          # Changes along a trajectory
                                          pseudo <- TSCAN::testPseudotime( sce_clust_SCE , pseudotime=mnn.pseudo )
                                          pseudo <- pseudo[order(pseudo$p.value),]	
                                          
                                          return(pseudo)
                                          
                                        })
  
}


