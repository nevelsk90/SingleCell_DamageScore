#### spatial transcriptomics ####  
library( Matrix)
library(Seurat)
library( ggplot2)
library(ggthemes)
library(ggpubr)
library(RColorBrewer)
library(ggeasy)
# F to calculate damage score
source("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/cell-damage-score/AUCell_script.r")
source("/home/tim_nevelsk/PROJECTS/myCode/usefulRfunc.r")

# release memory
mallinfo::mallinfo()
mallinfo::malloc.trim()
gc()

# podocyte cell-type markers
PodoMarks.mouse <-  c( "Wt1", "Nphs1", "Nphs2", "Synpo","Robo2","Thsd7a","Nebl","Ptpro","Magi2","Podxl")

# podocyte damage signature
DS_all <- read.table( header = T,  file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")
# convert gene names to human orthologs
DS.HOMO <-  fun_homoTO.FROMmouse( gns = DS_all$gene_symbol, TO = F)
DS_all.HOMO <-  DS_all
DS_all.HOMO$HOMO <- unlist(DS.HOMO$HOMO)

## podocyte markers
kidneyMarks_sc.suszt <- read.csv(sep = "\t","/media/tim_nevelsk/WD_tim/ANNOTATIONS/kidney/Susztak2021_mouse_scRNA.csv")
kidneyMarks_sc.Humphr <- read.csv(sep = ";","/media/tim_nevelsk/WD_tim/ANNOTATIONS/kidney/Humphreys2022_mouse_scRNA.csv")

#### GSE182939 mouse data from Humphrey https://pubmed.ncbi.nlm.nih.gov/34853151/  ####
  {
  ll<- list.dirs("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE182939_RAW", recursive = F)
lldir <- paste0( ll,"/outs")

mm.kidney_sptl <- lapply(lldir, function(X){
  Load10X_Spatial( data.dir = X, filename="filtered_feature_bc_matrix.h5" )
})

# run normalisation and Dim reduction
mm.kidney_sptl.dimred <- lapply(seq(mm.kidney_sptl), 
                                function(ii){
                                  
                                  sce.decontX.Seurat <- mm.kidney_sptl[[ii]]
                                  
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
                                  print( ElbowPlot( sce.decontX.Seurat , ndims = 50 ) ) # choosing number of PCs to include 
                                  
                                  # DimPlot( sce, reduction = "pca", group.by = "gtype", order =  sample(colnames(sce)) )
                                  
                                  # cluster
                                  sce.decontX.Seurat  <- FindNeighbors(sce.decontX.Seurat , dims = 1:20)
                                  sce.decontX.Seurat  <- FindClusters(sce.decontX.Seurat , resolution = 1.5)
                                  
                                  # 
                                  # run UMAP
                                  sce.decontX.Seurat <- RunUMAP(sce.decontX.Seurat, dims = 1:20 )
                                  return(sce.decontX.Seurat)
                                })
names( mm.kidney_sptl ) <- names( mm.kidney_sptl.dimred )<- basename( ll)
saveRDS(mm.kidney_sptl.dimred, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE182939_RAW/mm.kidney_sptr.Seur_dimred.rda")


## annotate  podocytes
{ 
  # prepare gene sets
  marks.hmphr.sn <- read.table(header = T, sep = ";","/media/tim_nevelsk/WD_tim/ANNOTATIONS/kidney/Humphreys2020_mouse_snRNA.csv")
  marks.hmphr.sn <- split(marks.hmphr.sn$gene, marks.hmphr.sn$celltype)
  geneSets.hmphr.sn <- GSEABase::GeneSetCollection( lapply( seq(marks.hmphr.sn) ,
                                                            function(ii) GSEABase::GeneSet( marks.hmphr.sn[[ii]] , 
                                                                                            setName= names(marks.hmphr.sn)[ii] ) ))
  
  
  # plot podomarks
  mm.kidney_sptl.dimred <- lapply(seq( mm.kidney_sptl.dimred), function(ii){
    datt<-  mm.kidney_sptl.dimred[[ii]]
    datt$podoMarks <-  colMeans(datt@assays$Spatial@data[
      c("Wt1","Nphs1","Nphs2","Synpo","Robo2","Thsd7a","Nebl","Ptpro"),], na.rm = T)
    datt$condition <- sub( "_.*", "",names(mm.kidney_sptl.dimred)[ii])
    return(datt)
    })
  
  gglist <- lapply(seq( mm.kidney_sptl.dimred), function(ii){
    datt<-  mm.kidney_sptl.dimred[[ii]]
   
    gg1 <- DimPlot(datt, label = T,shuffle = T)+ggtitle(unique(datt$condition))
    gg2<- FeaturePlot( datt, features = "podoMarks" )
    return(list(gg1, gg2))
  })
  cowplot::plot_grid(plotlist  = Reduce(c,(gglist)) , nrow = 5)
  
  }


### calculate PDS
{
  DS_all <- read.table( header = T,  file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")
  mm.kidney_sptl.dimred <- readRDS( file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE182939_RAW/mm.kidney_sptr.Seur_dimred.rda")

  
  ### for all cells
  podoClID <- c(7,8,7,8,7)
  mm.kidney_sptl.dimred.podo <- lapply( seq(mm.kidney_sptl.dimred),function(ii){
    subset(mm.kidney_sptl.dimred[[ii]], id=podoClID[ii])
  })
  
  mm.kidney_sptl.dimred.podo <- merge(mm.kidney_sptl.dimred.podo[[1]],
                                      y = mm.kidney_sptl.dimred.podo[2:5])
  
  
  expMat.podo <-mm.kidney_sptl.dimred.podo@assays$Spatial@counts
  # expMat.podo <- seur.BTBR.podo@assays$Spatial@counts
  expMat.podo <- expMat.podo[ rowSums( round(expMat.podo)>0)> 10 ,]
  set.seed(42)
  mm.kidney_sptl.dimred.podo$PDS <- DS_calc.func( exprMatrices = expMat.podo,
                                                   ceilThrsh = 0.3 ,
                                                   DSignature = DS_all , 
                                                   ntop = 42, wghtd = T, progStat = T)
  ## plot PDS
  mm.kidney_sptl.dimred.podo@meta.data$condition <- 
    factor(mm.kidney_sptl.dimred.podo@meta.data$condition,
           levels = c("fsham","f4hr","f12hr","f2dps","f6wks"))
  # boxplot
  toPlot <- mm.kidney_sptl.dimred.podo@meta.data
  ggplot( data = toPlot, aes(y=PDS, x=condition))+ 
    geom_boxplot()+ 
    theme( text = element_text(size = 20)) +theme_bw() +
    theme( text=element_text(size=24))
  
}


# spatial plots
# VlnPlot( fsham_137 , features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
library("scales")
sfactor <- c(2.4,2.0,1.8,1.8,2.2)
gglist <- lapply( seq(mm.kidney_sptl.dimred.PDS42), function(ii){
  
  SpatialFeaturePlot( mm.kidney_sptl.dimred.PDS42[[ii]], 
                      features = "PDS.42podo" ,pt.size.factor = sfactor[ii]) +
    scale_fill_gradientn(colours = colorspace::diverge_hcl(8, "Blue-Red 2"), 
                         values = rescale(c(-0.085, -0.0245,0.025)),
                         guide = "colorbar", limits=c(-0.085, 0.025)) +
    theme(legend.position = "right")+
    ggtitle(names(mm.kidney_sptl.dimred.PDS42)[ii])
})
cowplot::plot_grid( plotlist = gglist)
# fraction of podocytes
Podo.frac <- sapply( seq(mm.kidney_sptl.dimred.PDS42),  function(ii){
  as.numeric(summary(mm.kidney_sptl.dimred.PDS42[[ii]]$Podo)[3]) /
    as.numeric(summary( mm.kidney_sptl.dimred.PDS42[[ii]]$Podo)[2] )
  
})
names(Podo.frac) <- sub( "_.*", "",names(mm.kidney_sptl.dimred.PDS42))
barplot(Podo.frac, main="fraction of annotated podocytes")
}

#### GSE171406 mouse data from KPMP human study https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10356613/#MOESM6 ####
  
    ll<- list.dirs(
    "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE171406_RAW", 
    recursive = F)[1:4]
  
  kidney.AKI_sptl <- lapply(ll, function(X){
    Load10X_Spatial( data.dir = X, filename="filtered_feature_bc_matrix.h5" )
  })
  
  # run normalisation and Dim reduction
  kidney.AKI_sptl.dimred <- lapply(seq(kidney.AKI_sptl), 
                                  function(ii){
                                    
                                    sce.decontX.Seurat <- kidney.AKI_sptl[[ii]]
                                    
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
                                    print( ElbowPlot( sce.decontX.Seurat , ndims = 50 ) ) # choosing number of PCs to include 
                                    
                                    # DimPlot( sce, reduction = "pca", group.by = "gtype", order =  sample(colnames(sce)) )
                                    
                                    # cluster
                                    sce.decontX.Seurat  <- FindNeighbors(sce.decontX.Seurat , dims = 1:20)
                                    sce.decontX.Seurat  <- FindClusters(sce.decontX.Seurat , resolution = 1)
                                    
                                    # 
                                    # run UMAP
                                    sce.decontX.Seurat <- RunUMAP(sce.decontX.Seurat, dims = 1:20 )
                                    
                                    # add mito and ribo genes  percentage
                                    sce.decontX.Seurat <- PercentageFeatureSet(sce.decontX.Seurat, "^mt-", col.name = "percent_mito")
                                    sce.decontX.Seurat <- PercentageFeatureSet(sce.decontX.Seurat, "^Rp[sl]", col.name = "percent_ribo")
                                    
                                    return(sce.decontX.Seurat)
                                  })
  names( kidney.AKI_sptl.dimred )<- basename( ll)
  saveRDS(kidney.AKI_sptl.dimred, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE171406_RAW/kidney.AKI_sptl.dimred.rda")
  
  
  ## annotate  podocytes
  { 
    # prepare gene sets
    marks.hmphr.sn <- read.table(header = T, sep = ";","/media/tim_nevelsk/WD_tim/ANNOTATIONS/kidney/Humphreys2020_mouse_snRNA.csv")
    marks.hmphr.sn <- split(marks.hmphr.sn$gene, marks.hmphr.sn$celltype)
    
    kidney.AKI_sptl.dimred.mm <- lapply( c(1,3,4), function(ii){
      datt<- kidney.AKI_sptl.dimred[[ii]]
      datt$podoMarks <-  colMeans(datt@assays$Spatial@data[
          c("Wt1","Nphs1","Nphs2","Synpo","Robo2","Thsd7a","Nebl","Ptpro"),], na.rm = T)
      datt$condition <-   names(kidney.AKI_sptl.dimred)[ii]
      datt  <- FindNeighbors(datt , dims = 1:10)
      datt  <- FindClusters(datt , resolution = 1)
      
      # 
      # run UMAP
      datt <- RunUMAP(datt, dims = 1:10 )
      
      return(datt)
  
    })

    ##3 plot 
    gglist <- lapply(1:3, function(ii){
      datt<- kidney.AKI_sptl.dimred.mm[[ii]]
      gg1 <- DimPlot(datt, label = T,shuffle = T)
      gg2<- FeaturePlot(datt, features = "podoMarks")
      return(list(gg1, gg2))
    })
        cowplot::plot_grid(plotlist  = Reduce(c,(gglist)) , nrow = 3)
        
 
    
    # plot UMAP for podocyte annotations
    
    pdf( width =  10 , height = 12 , file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Spatial.trnscrptm/GSE171406/GSE171406_AUCell.podo_plots.pdf")
    # a<-dev.cur()
    # png( width =  600 , height = 3200 , file="AUCell_sztk.sc21_hist.png")
    # dev.control("enable")
    par(mfrow=c(4,3))
 
    dev.off()
    
    
    
  }
  
  
  ### calculate PDS
  {
    DS_all <- read.table( header = T,  file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")
    kidney.AKI_sptl.dimred.mm <- readRDS( file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE182939_RAW/mm.kidney_sptr.Seur_dimred.rda")
    
    
    ### for all cells
    # resolution=3, ndim=10
    podoClID.mm <- c(17,13,16)
    # # resolution=1, ndim=10
    # podoClID.mm <- c(10,11,9)
    
    kidney.AKI_sptl.dimred.Podo <- lapply(1:3, function(ii){
     subset( kidney.AKI_sptl.dimred.mm[[ii]], id=podoClID.mm[[ii]])
    })
    kidney.AKI_sptl.dimred.Podo <- merge( kidney.AKI_sptl.dimred.Podo[[1]],
                                          kidney.AKI_sptl.dimred.Podo[2:3])
  
    
    expMat.podo <-kidney.AKI_sptl.dimred.Podo@assays$Spatial@counts
    # expMat.podo <- seur.BTBR.podo@assays$Spatial@counts
    expMat.podo <- expMat.podo[ rowSums( round(expMat.podo)>0)> 10 ,]
    kidney.AKI_sptl.dimred.Podo$PDS <- DS_calc.func( exprMatrices = expMat.podo,
                                        ceilThrsh = 0.3 ,
                                        DSignature = DS_all , 
                                        ntop = 42, wghtd = T, progStat = T)
    ## plot PDS
      # boxplot
    toPlot <- kidney.AKI_sptl.dimred.Podo@meta.data
    ggplot( data = toPlot, aes(y=PDS, x=condition))+ 
      geom_boxplot()+ 
      theme( text = element_text(size = 20)) +theme_bw() +
      theme( text=element_text(size=24))
    
  }
  
  
  # spatial plots
  # VlnPlot( fsham_137 , features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  library("scales")
  # sfactor <- c(2.4,2.0,1.8,1.8,2.2)
  gglist <- lapply( seq(hs.kidney_sptl.dimred.PDS42), function(ii){
    
    SpatialFeaturePlot( hs.kidney_sptl.dimred.PDS42[[ii]], 
                        # pt.size.factor = sfactor[ii] ,
                        features = "PDS.42podo" ) +
      # scale_fill_gradientn(colours = colorspace::diverge_hcl(8, "Blue-Red 2"), 
      #                      values = rescale(c(-0.085, -0.0245,0.025)),
      #                      guide = "colorbar", limits=c(-0.085, 0.025)) +
      theme(legend.position = "right")+
      ggtitle(names(hs.kidney_sptl.dimred.PDS42)[ii])
  })
  cowplot::plot_grid( plotlist = gglist)
  # fraction of podocytes
  Podo.frac <- sapply( seq(mm.kidney_sptl.dimred.PDS42),  function(ii){
    as.numeric(summary(mm.kidney_sptl.dimred.PDS42[[ii]]$Podo)[3]) /
      as.numeric(summary( mm.kidney_sptl.dimred.PDS42[[ii]]$Podo)[2] )
    
  })
  names(Podo.frac) <- sub( "_.*", "",names(mm.kidney_sptl.dimred.PDS42))
  barplot(Podo.frac, main="fraction of annotated podocytes")
  
#### GSE190094 mouse data from KPMP human study https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10356613/#MOESM6 ####

  ## annot.data
  annot.GSE190094 <- GEOquery::getGEO("GSE190094")
  annot.GSE190094 <- annot.GSE190094$GSE190094_series_matrix.txt.gz@phenoData@data
  
  ### format coordinate files to use with Python scripts
    {
    # read spatial info
    dir <- "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/beadLoc_qc/"
    ffile <- list.files( path = dir, pattern = "_qc.csv",full.names = T )
    
    # BTBR mouse samples are 13w.old
    lapply( which(annot.GSE190094$`age:ch1`=="13 week") , function(ii){
      datt <- read.csv( ffile[[ii]])
      colnames(datt) <- c("barcode",   "x",   "y", "cell_type")
      write.table(datt, sep = ",", col.names = NA,
                  file=sub("beadLoc_qc/","BTBR/beadLocPython_qc",
                           sub("Puck.*Bead","Bead",ffile[[ii]])))
      })
  }
  
  ### generate Seurat object 
    {
    
    ## read spatial info
    dir <- "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/beadLoc_qc/"
    ffile <- list.files( path = dir, pattern = "_qc.csv" , full.names = T )
    
    slideFile <- lapply(seq(ffile) , function(ii){
      ReadSlideSeq( coord.file = ffile[[ii]], assay = "Spatial")
      })
    names(slideFile) <- sub("_.*","",basename(ffile))
    # saveRDS(slideFile, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/GSE190094_beadslocation.qc.list.rda")
    slideFile <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/GSE190094_beadslocation.qc.list.rda")
    
    ###  create Seurat objects
    # read in count data from sparse matrices
    expr.qc <- readRDS(file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/GSE190094_DGE.qc.list.rda")
    # dims sapply(expr.qc, ncol)[c(1:5,45:49)]
    # 30344 24438 30182 29011 15592 12541 12533 22430 11005 28577 20424

    seur.list.BTBR <- lapply( which(annot.GSE190094$`age:ch1`=="13 week") , 
                              function(ii){
                                print(ii)
                                seur = CreateSeuratObject(counts = expr.qc[[ii]], assay="Spatial")
                                
                                seur@images$image  <- slideFile[[ii]]
                                
                                seur$cell_type <-  seur@images$image@coordinates$cell_type
                                seur$geoID <- annot.GSE190094$geo_accession[ii]
                                seur$sampleID <- annot.GSE190094$title[ii]
                                seur$gtype <- annot.GSE190094$`genotype:ch1`[ii]
                                return(seur)
                              })
    names(seur.list.BTBR) <- annot.GSE190094$geo_accession[annot.GSE190094$`age:ch1`=="13 week"]
    saveRDS( seur.list.BTBR , file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/GSE190094_BTBR.seur.rda")

    ### QC plots
    
    # barplot od number of the median fetures per sample
    nFeature_list <- Reduce( rbind, sapply( seq(seur.list.BTBR), function(ii){
      summary(seur.list.BTBR[[ii]]$nFeature_Spatial)
    }))
    rownames(nFeature_list) <- paste0( annot.GSE190094$geo_accession[1:49],"_",
                                      annot.GSE190094$`genotype:ch1`[1:49])
    par(mar=c(3,12,2,2))
    barplot(nFeature_list[,3], las=2 , horiz = T)
    abline(v=100, col="red")
    abline(v=150, col="red")
    abline(v=200, col="red")
    # samples with more than median of 200 features 
    seur.list.BTBR.test <- seur.list.BTBR[ which( nFeature_list[,3] > 200) ]
    
    # plot number of cells
    barplot( sapply(seur.list.BTBR , ncol) )
    
    }
  
  seur.list.BTBR <- readRDS( file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/GSE190094_BTBR.seur.rda")
  ### annotate podocytes manually using podo markers
  # filter clusteres with the average coverage < 100 features per bead
  # save filtered data and plot UMAPs here /media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/filtered.by.nFeatr  
  ### plot podocytes cluster before and after filterings
  lapply( 1:49, function(ii)
    {
    print(ii)
    sce.decontX.Seurat <-  seur.list.BTBR[[ii]]
    
    # add mito and ribo genes  percentage
    sce.decontX.Seurat <- PercentageFeatureSet(sce.decontX.Seurat, "^mt-", col.name = "percent_mito")
    sce.decontX.Seurat <- PercentageFeatureSet(sce.decontX.Seurat, "^Rp[sl]", col.name = "percent_ribo")
    # hist((sce.decontX.Seurat$nFeature_Spatial), xlim = c(0,2000), breaks = 50)
    ### find variable features
    sce.decontX.Seurat <- NormalizeData( sce.decontX.Seurat )
    
    sce.decontX.Seurat$PodoMarks <- colMeans(sce.decontX.Seurat@assays$Spatial@data[
      c("Wt1","Nphs1","Nphs2","Synpo","Robo2","Thsd7a","Nebl","Ptpro"),], na.rm = T)
    
    sce.decontX.Seurat <- FindVariableFeatures( sce.decontX.Seurat ,
                                                selection.method = "vst",
                                                nfeatures = 1000)
    
    ### scale the data
    sce.decontX.Seurat <- ScaleData(sce.decontX.Seurat)
    
    ### run PCA
    sce.decontX.Seurat <- RunPCA( sce.decontX.Seurat ,
                                  features = VariableFeatures( object = sce.decontX.Seurat ) )
    pp001 <-  ElbowPlot( sce.decontX.Seurat , ndims = 50 )  # choosing number of PCs to include
    
    # DimPlot( sce, reduction = "pca", group.by = "gtype", order =  sample(colnames(sce)) )
    
    # cluster
    sce.decontX.Seurat  <- FindNeighbors(sce.decontX.Seurat , dims = 1:30)
    sce.decontX.Seurat  <- FindClusters(sce.decontX.Seurat , resolution = 1)
    
    #
    # run UMAP
    sce.decontX.Seurat <- RunUMAP(sce.decontX.Seurat, dims = 1:30 )
    
    # add mito and ribo genes  percentage
    pp0 <- DimPlot( sce.decontX.Seurat, shuffle = T, label = T)
    pp1 <- FeaturePlot( sce.decontX.Seurat, features = "PodoMarks",
                        max.cutoff = 1 )
    pp2 <- VlnPlot(sce.decontX.Seurat,"nFeature_Spatial" )+
      theme(legend.position = 'none')
    
    ## filter 
    highCovCl <- aggregate( nFeature_Spatial~seurat_clusters , 
                            data=sce.decontX.Seurat@meta.data, FUN=median)
    
    if( sum( highCovCl$nFeature_Spatial>100)>0 ){
      sce.decontX.Seurat.filt <- subset( sce.decontX.Seurat, 
                                         subset=seurat_clusters %in% 
                                           highCovCl$seurat_clusters[
                                             highCovCl$nFeature_Spatial>100
                                           ])
      sce.decontX.Seurat.filt <- NormalizeData( sce.decontX.Seurat.filt )
      sce.decontX.Seurat.filt <- FindVariableFeatures( sce.decontX.Seurat.filt ,
                                                       selection.method = "vst",
                                                       nfeatures = 1000)
      
      ### scale the data
      sce.decontX.Seurat.filt <- ScaleData( sce.decontX.Seurat.filt )
      
      ### run PCA
      sce.decontX.Seurat.filt <- RunPCA( sce.decontX.Seurat.filt ,
                                         features = VariableFeatures( object = sce.decontX.Seurat.filt ) )
      pp002 <-  ElbowPlot( sce.decontX.Seurat.filt , ndims = 50 )  # choosing number of PCs to include
      
      # cluster
      sce.decontX.Seurat.filt  <- FindNeighbors(sce.decontX.Seurat.filt , dims = 1:20)
      sce.decontX.Seurat.filt  <- FindClusters(sce.decontX.Seurat.filt , resolution = 1)
      
      # run UMAP
      sce.decontX.Seurat.filt <- RunUMAP( sce.decontX.Seurat.filt , dims = 1:20 )
      
      # add mito and ribo genes  percentage
      pp3 <- DimPlot( sce.decontX.Seurat.filt ,shuffle = T , label = T)
      pp4 <- FeaturePlot( sce.decontX.Seurat.filt, features = "PodoMarks",
                          max.cutoff = 1 )
      pp5 <- VlnPlot(sce.decontX.Seurat.filt,"nFeature_Spatial")+
        theme(legend.position = 'none')
      
      
      
      ggl <- cowplot::plot_grid( plotlist = 
                                   list(pp0,pp1,pp2,pp3,pp4,pp5), nrow = 2)
      
      png( width = 2000 , height = 1000, file = paste0( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/filtered.by.nFeatr/",
                                                        unique(sce.decontX.Seurat.filt$geoID),"_",ii,".png"))
      print(ggl)
      dev.off()
      
      saveRDS(sce.decontX.Seurat.filt, file= paste0( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/filtered.by.nFeatr/",
                                                     unique(sce.decontX.Seurat.filt$geoID),"_",ii,".rda"))
      
    }
    
    
    
    # release memory
    mallinfo::mallinfo()
    mallinfo::malloc.trim()
    gc()
    
  })
  
  ### select 8 test samples with higher coverage
    {
    testIDs <- sub( "_.*", "", list.files(pattern = ".*BeadLocationsForR_qc",
                                          path = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/test8") )
    
    names(seur.list.BTBR) <-  annot.GSE190094$geo_accession[ 
      which(annot.GSE190094$`age:ch1`=="13 week")]
    
    seur.list.BTBR.test <-  lapply( seq(testIDs), function(ii){
      
      datt.glomAnnot <- glom_annot.list[[ii]]
      print(testIDs[ii])
      datt <-  seur.list.BTBR[[testIDs[ii]]]
      datt$cell_type.Glom <- datt.glomAnnot$cell_type[ match( colnames(datt),
                                                              datt.glomAnnot$barcode ) ]
      datt$cluster <- datt.glomAnnot$cluster[ match( colnames(datt),
                                                     datt.glomAnnot$barcode ) ]
      datt$podo.Glom <- ifelse(datt$cell_type.Glom=="Podocyte", "Podocyte", NA)
      datt$podo <- ifelse( datt$cell_type=="Podocyte", "Podocyte", NA )
      
      datt$PodoMarks <- colMeans(datt@assays$Spatial@data[
        c("Wt1","Nphs1","Nphs2","Synpo","Robo2","Thsd7a","Nebl","Ptpro"),],
        na.rm = T)
      
      return(datt)
    })
    
    ### calculate PDS
    {
      
      ### extract podocytes using author definition
      {
        ## provided by authors ( KNN filatrated)
        seur.list.BTBR.podo <- lapply( seq(seur.list.BTBR), function(ii){
          datt <- seur.list.BTBR[[ii]]
          subset(  datt, subset= cell_type=="Podocyte")
        })
        seur.BTBR.podo <- merge( seur.list.BTBR.podo[[2]], 
                                 y=seur.list.BTBR.podo[3:length(seur.list.BTBR.podo)])
        
        
        ### find variable features
        seur.BTBR.podo <- NormalizeData( seur.BTBR.podo )
        saveRDS( seur.BTBR.podo , file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/GSE190094_BTBR.Podocyte.seur.rda")
        
        # bio replicates
        seur.BTBR.podo$sampleBio <- stringr::str_sub(seur.BTBR.podo$sampleID, end = -2)
        
        # plots
        VlnPlot(seur.BTBR.podo, features = "Nphs2",group.by = "sampleBio")
        VlnPlot(seur.BTBR.podo, features = "nFeature_Spatial",group.by = "sampleBio")
        VlnPlot(seur.BTBR.podo, features = "nFeature_Spatial",group.by = "geoID")
        
        
      }
      
      
      DS_all <- read.table( header = T,  file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")
      
      # selected by authors
      seur.BTBR.podo <- readRDS( file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/GSE190094_BTBR.Podocyte.seur.rda")
      seur.BTBR.podo <- subset( seur.BTBR.podo, subset= geoID %in% c("GSM5713342","GSM5713369","GSM5713371", "GSM5713375","GSM5713376","GSM5713378"))
      
      # select based on clusters
      seur.BTBR.podo <- seur.list.BTBR.test.podo
      seur.BTBR.podo$sampleBio <- stringr::str_sub(seur.BTBR.podo$sampleID, end = -2)
      # 
      ### for all cells
      set.seed(42)
      # calculate PDS seperately for podocytes and nonpodocytes
      expMat.podo <-seur.BTBR.podo@assays$Spatial@counts
      # expMat.podo <- seur.BTBR.podo@assays$Spatial@counts
      expMat.podo <- expMat.podo[ rowSums( round(expMat.podo)>0)> 10 ,]
      seur.BTBR.podo$PDS <- DS_calc.func( exprMatrices = expMat.podo,
                                          ceilThrsh = 0.3 ,
                                          DSignature = DS_all , 
                                          ntop = 42, wghtd = T, progStat = T)
      
      # boxplot
      gg1 <- ggplot( data = seur.BTBR.podo@meta.data , aes(y=gtype, x=PDS))+ 
        geom_boxplot()+ 
        # facet_grid(cols = vars(Podo)) + 
        theme( text = element_text(size = 20)) +theme_bw() +
        theme( text=element_text(size=24))
      gg2 <- ggplot( data = seur.BTBR.podo@meta.data , aes(y=sampleBio, x=PDS))+ 
        geom_boxplot()+ 
        # facet_grid(cols = vars(Podo)) + 
        theme( text = element_text(size = 20)) +theme_bw() +
        theme( text=element_text(size=24))
      gg3 <- ggplot( data = seur.BTBR.podo@meta.data , aes(y=geoID, x=PDS))+ 
        geom_boxplot()+ 
        # facet_grid(cols = vars(Podo)) + 
        theme( text = element_text(size = 20)) +theme_bw() +
        theme( text=element_text(size=24))
      
      
      
      cowplot::plot_grid(plotlist = list(gg1,gg2,gg3), nrow = 1)
      
      
    }
  }
  
  ### analyse all samples with detectable podocyte cluster
    {
      testIDs <- c(2,11,16:19,21,22,25,26,30:47,49)
      test_geoID <- annot.GSE190094$geo_accession[ testIDs ]
      
      # list of all samples with detectable podocytes
      ll<- list.files(pattern = "rda",full.names = T,"/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/filtered.by.nFeatr")
      names(ll) <- sub( "_.*", "", basename(ll))
    
      ### filtr out clusters with average less than 100 features per bead and plot annot.BTBR.pdclust
      lapply( seq(test_geoID),
              function(ii){
                
                ID <- test_geoID[[ ii ]]
                # 
                print( paste( ii, "   " ,ID ))
                
                # load data
                datt <- readRDS(ll[ ID ])
                datt$PodoMarks <-  colMeans( datt@assays$Spatial@data[
                  rownames(datt@assays$Spatial@data) %in% PodoMarks.mouse, ] , na.rm = T )
                
                X <- aggregate( nFeature_Spatial ~ seurat_clusters,
                                data=datt@meta.data, FUN=mean, na.rm=T)
                seur <- subset( datt , subset= seurat_clusters%in% 
                                  X$seurat_clusters[ X$nFeature_Spatial>=100])                                    
                
                # ### find variable features
                seur  <- FindVariableFeatures( seur , selection.method = "vst" )
                ### scale the data
                seur <- ScaleData( seur )
                
                ### run PCA
                seur <- RunPCA(seur, features = VariableFeatures(object = seur))
                print( ElbowPlot(seur , ndims = 50 )) # choosing number of PCs to include 
                
                # cluster
                seur <- FindNeighbors( seur, dims = 1:20 )
                seur <- FindClusters( seur, resolution = 1 )
                
                # run UMAP
                seur <- RunUMAP( seur, dims = 1:20 )
                
                
                pp1 <- DimPlot( seur, shuffle = T, label = T)
                pp2 <- FeaturePlot( seur , features = "nFeature_Spatial", label = T)
                pp3 <- FeaturePlot( seur , features = "PodoMarks",max.cutoff = 2, label = T)
                
                ggl <- cowplot::plot_grid( plotlist = list(pp1,pp2, pp3), nrow = 1)
                
                
                png(height = 600, width = 1800, filename = paste0("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Spatial.trnscrptm/GSE190094/BTBR/podoClust/",
                                                                  test_geoID[ii],"_filt.Dimplot.png"))
                print(ggl)
                dev.off()
                
                saveRDS( seur , file= paste0("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/podoClust/",
                                             test_geoID[ii], "_filtSeur.rda"))
                gc()
              })
      
      ### plot podocyteMarker expr and cluster marker heatmaps
      ll.filt <- list.files(pattern = "rda",full.names = T,"/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/podoClust/filtSeur")
      names(ll.filt) <- sub( "_.*", "", basename(ll.filt))
      
      seur.BTBR_podoClustID <-  lapply( seq(ll.filt),
              function(ii){
                require(dplyr)
                
                ID <- names(ll.filt)[ ii ]
                # 
                print( paste( ii, "   " ,ID ))
                
                # load data
                seur <- readRDS(ll.filt[ ID ])
              
                # find podocyte cluster
                X <- aggregate( PodoMarks~ seurat_clusters, data=seur@meta.data, FUN=mean, na.rm=T )
                podoClID <- X$seurat_clusters[ which.max(X$PodoMarks) ]
                # seur.podo <- subset( seur , subset= seurat_clusters== podoClID )
                return( podoClID )
                # # plot podoMarks
                # pp0 <- VlnPlot( seur , features = c("PodoMarks","nFeature_Spatial"))
                # 
                # png(height = 400, width = 1000 , filename = paste0("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Spatial.trnscrptm/GSE190094/BTBR/podoClust/",
                #                                                 test_geoID[ii],"_filt.podoMarks_vln.png"))
                # print(pp0)
                # dev.off()
                # 
                # seur.sub <- subset( seur , downsample= 500 )
                # seur.sub <- ScaleData(seur.sub,
                #                       features = rownames(seur.sub)[
                #                         rowSums(seur.sub@assays$Spatial@counts>0)>10])
                # # find cluster markers
                # seur.clustMarks <- FindAllMarkers( seur.sub )
                # top10 <- seur.clustMarks %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
                # pp <- DoHeatmap( seur.sub, features =top10$gene )
                # 
                # pdf(height = 20, width = 20 , filename = paste0("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Spatial.trnscrptm/GSE190094/BTBR/podoClust/",
                #                                                   test_geoID[ii],"_filt.clustMarks_heatmap.pdf"))
                # print(pp)
                # dev.off()

                # return(seur.clustMarks)
              })
      saveRDS(seur.BTBR_podoClustID)
      
      ### extract podocytes
      seur.BTBR_podo <-  lapply( seq(ll.filt),
                                  function(ii , podoClID =seur.BTBR_podoClustID ){
                                          require(dplyr)
                                          
                                          ID <- names(ll.filt)[ ii ]
                                          # 
                                          print( paste( ii, "   " ,ID ))
                                          
                                          # load data
                                          seur <- readRDS(ll.filt[ ID ])
                                          
                                          # find podocyte cluster
                                          seur.podo <- subset( seur , subset= seurat_clusters == podoClID[[ii]])
                                          # seur.podo <- subset( seur , subset= seurat_clusters== podoClID )
                                          return( seur.podo )
                                          # # plot podoMarks
                                          # pp0 <- VlnPlot( seur , features = c("PodoMarks","nFeature_Spatial"))
                                          # 
                                          # png(height = 400, width = 1000 , filename = paste0("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Spatial.trnscrptm/GSE190094/BTBR/podoClust/",
                                          #                                                 test_geoID[ii],"_filt.podoMarks_vln.png"))
                                          # print(pp0)
                                          # dev.off()
                                          # 
                                          # seur.sub <- subset( seur , downsample= 500 )
                                          # seur.sub <- ScaleData(seur.sub,
                                          #                       features = rownames(seur.sub)[
                                          #                         rowSums(seur.sub@assays$Spatial@counts>0)>10])
                                          # # find cluster markers
                                          # seur.clustMarks <- FindAllMarkers( seur.sub )
                                          # top10 <- seur.clustMarks %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
                                          # pp <- DoHeatmap( seur.sub, features =top10$gene )
                                          # 
                                          # pdf(height = 20, width = 20 , filename = paste0("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Spatial.trnscrptm/GSE190094/BTBR/podoClust/",
                                          #                                                   test_geoID[ii],"_filt.clustMarks_heatmap.pdf"))
                                          # print(pp)
                                          # dev.off()
                                          
                                          # return(seur.clustMarks)
                                        })
      
      seur.BTBR_podo <- merge( seur.BTBR_podo[[1]], y = seur.BTBR_podo[2:length(seur.BTBR_podo)] ,
                               add.cell.ids = names(ll.filt) )
      
      seur.BTBR_podo$sampleID <- sub("_.*","",colnames(seur.BTBR_podo))
      seur.BTBR_podo <- subset( seur.BTBR_podo , subset= PodoMarks >0 )
      barplot(table( seur.BTBR_podo$sampleID[ seur.BTBR_podo$PodoMarks>0]))
      abline(h =100, col="red", lwd=1)
      
      saveRDS( seur.BTBR_podo , file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/podoClust/seur.BTBR_podo.rda")
      
      
      ### calculate PDS
      {
    
          set.seed(42)
          # calculate PDS seperately for podocytes and nonpodocytes
          expMat.podo  <-  seur.BTBR_podo@assays$Spatial@counts
          # expMat.podo <- seur.BTBR.podo@assays$Spatial@counts
          expMat.podo  <- expMat.podo[ rowSums( round(expMat.podo)>0)> 0 ,]
          seur.BTBR_podo$PDS <- DS_calc.func( exprMatrices = expMat.podo,
                                     ceilThrsh = 0.3 ,
                                     DSignature = DS_all , 
                                     ntop = 42, wghtd = T, progStat = T)
          

        
        
      ggplot( seur.BTBR_podo@meta.data, aes(x=geoID, color=gtype,
                                        y=PDS)) + 
          geom_boxplot()+geom_point( alpha=0.4)+
          theme(legend.position = "right")+ theme_bw()+ 
          scale_color_colorblind()+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
     
      
       }
      
    }

  ### plot KNN curated OG glom annotations
    {
      
     
      ### Spatial plot podocytes to decide on number of clusters per sample for KNN clustering
        {
         gglist <- lapply( seq(testIDs) , function(ii){
           datt <-  seur.list.BTBR[[testIDs[ii]]]
           datt$podo <- ifelse( datt$cell_type=="Podocyte", "Podocyte", NA )
           
           SpatialDimPlot( datt , alpha = 0.7,stroke = 0.1 , 
                           group.by = "podo") + 
             theme(legend.position = "right")+ 
             scale_fill_colorblind()+ggtitle( testIDs[ii])
         })
         ggn <- cowplot::plot_grid( plotlist = gglist , nrow  = 2)
         png(width = 2000, height = 800, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Spatial.trnscrptm/GSE190094/BTBR/test8/SpatialDimPlot_podoFilt.png",
         )
         print(ggn)
         dev.off()
       }
      
      ### Spatial plots  for glomerular annotations after KNN clustering
        {
        ll <-  list.files(pattern = ".*_glom_annot",full.names = T,
                          path = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/test8")
        glom_annot.list <- lapply(ll, read.csv, row.names=1)
        names(glom_annot.list) <- sub("_.*","",basename(ll))
        
        gglist <-  lapply( seq(glom_annot.list), function(ii){
          
          datt.glomAnnot <- glom_annot.list[[ii]]
          print(testIDs[ii])
          datt <-  seur.list.BTBR[[testIDs[ii]]]
          datt$cell_type.Glom <- datt.glomAnnot$cell_type[ match( colnames(datt),
                                                                  datt.glomAnnot$barcode ) ]
          datt$cluster <- datt.glomAnnot$cluster[ match( colnames(datt),
                                                         datt.glomAnnot$barcode ) ]
          datt$podo.Glom <- ifelse(datt$cell_type.Glom=="Podocyte", "Podocyte", NA)
          datt$podo <- ifelse( datt$cell_type=="Podocyte", "Podocyte", NA )
          
          gg2<- SpatialDimPlot( datt , alpha = 0.7,stroke = 0.1 , 
                                group.by = "cell_type.Glom") + 
            theme(legend.position = "right")+ 
            scale_fill_colorblind()+ggtitle( names(glom_annot.list)[ii])
          
          return(gg2)
          
        })
        ggn <- cowplot::plot_grid( plotlist = gglist , nrow  = 2)
        
        png(width = 2000, height = 800, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Spatial.trnscrptm/GSE190094/BTBR/test8/SpatialDimPlot_glomAnnot.png",
        )
        print(ggn)
        dev.off()
      }
    
      ### plot coverage in podocyte before and after aggregation   
        {
           # plot
           gglistDat <- lapply( seq(glom_annot.list),  function(ii){
             print(ii)
             
             datt.glomAnnot <- glom_annot.list[[ii]]
             print(testIDs[ii])
             datt <-  seur.list.BTBR.test[[ii]]
             
             datTOplot <- datt@meta.data[datt$cell_type.Glom=="Podocyte" &
                                           !is.na(datt$cell_type.Glom),]
             
             
             ## aggregate 
             datt0 <- subset(datt , subset= cell_type.Glom %in% 
                               c("Podocyte"))
             
             XX <-  datt0@assays$Spatial@counts
             XX <- t( XX[ rowSums(XX>0, na.rm = T)>0 , ] )
             datTOplotAgg <- aggregate( XX , by = list(datt0$cluster) , FUN=sum) 
             rownames(datTOplotAgg) <- paste0("clust", datTOplotAgg$Group.1)
             datTOplotAgg <- t( datTOplotAgg[,colnames(datTOplotAgg)!="Group.1"])
             datTOplotAgg <- data.frame( nFeature_aggPodoPerClust=colSums(datTOplotAgg>0))
             datTOplotAgg$geoID <- unique(datt0$geoID)
             
             return( list(datTOplot , datTOplotAgg))
             
           })
           datTOplot <- Reduce( rbind , lapply( gglistDat, "[[" , 1))
           datTOplot$gtype <- annot.GSE190094$`genotype:ch1`[
             match( datTOplot$geoID ,annot.GSE190094$geo_accession ) ]
           datTOplotAgg <- Reduce( rbind , lapply( gglistDat, "[[" , 2) )
           datTOplotAgg$gtype <- annot.GSE190094$`genotype:ch1`[
             match( datTOplotAgg$geoID ,annot.GSE190094$geo_accession ) ]
           
           
           gg1<- ggplot( datTOplot, aes(x=geoID, color=gtype,
                                        y=nFeature_Spatial)) + 
             geom_boxplot()+geom_point( alpha=0.5)+
             theme(legend.position = "right")+ theme_bw()+ 
             scale_color_colorblind()+
             theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
           gg2<- ggplot(  datTOplotAgg , aes(x=geoID, color=gtype,
                                             y=nFeature_aggPodoPerClust)) +
             geom_boxplot()+  geom_point( alpha=0.5)+
             theme(legend.position = "right")+theme_bw()+ 
             scale_color_colorblind()+ theme(axis.title.x=element_blank(),
                                             axis.text.x=element_blank(),
                                             axis.ticks.x=element_blank())
           
           ggn <- cowplot::plot_grid(plotlist = list(gg1,gg2), nrow  = 2)
           
           png(width = 600, height = 800, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Spatial.trnscrptm/GSE190094/BTBR/test8/nFeatures_glomAnnotVSaggPodoPerGlom.png",
           )
           print(ggn)
           dev.off()
           
         }

      ### aggregate podocytes over gloms to increase the coverage
        {
           
           
           PodoAggPerGlom <- lapply( seq(glom_annot.list),  function(ii){
             print(ii)
             
             datt.glomAnnot <- glom_annot.list[[ii]]
             print(testIDs[ii])
             datt <-  seur.list.BTBR.test[[ii]]
             
             ## aggregate 
             datt0 <- subset(datt , subset= cell_type.Glom %in% 
                               c("Podocyte"))
             
             XX <-  datt0@assays$Spatial@counts
             XX <- t( XX[ rowSums(XX>0, na.rm = T)>0 , ] )
             datTOplotAgg <- aggregate( XX , by = list(datt0$cluster) , FUN=sum) 
             rownames(datTOplotAgg) <- paste0("clust", datTOplotAgg$Group.1)
             datTOplotAgg <- t( datTOplotAgg[,colnames(datTOplotAgg)!="Group.1"])
             
             
             return( datTOplotAgg )
             
           })
           
           
           ### calculate PDS for  podocyte  aggregated per glom
           
           PodoAggPerGlom_PDS <- Reduce( rbind, lapply( seq(PodoAggPerGlom), 
                                                        function(ii)
                                                        {
                                                          set.seed(42)
                                                          # calculate PDS seperately for podocytes and nonpodocytes
                                                          expMat.podo  <-  PodoAggPerGlom[[ii]]
                                                          # expMat.podo <- seur.BTBR.podo@assays$Spatial@counts
                                                          expMat.podo  <- expMat.podo[ rowSums( round(expMat.podo)>0)> 1 ,]
                                                          PDS <- DS_calc.func( exprMatrices = expMat.podo,
                                                                               ceilThrsh = 0.1 ,
                                                                               DSignature = DS_all , 
                                                                               ntop = 42, wghtd = T, progStat = T)
                                                          PDS_dat <- data.frame( glomID = colnames(expMat.podo), 
                                                                                 PDS=PDS , 
                                                                                 geoID = testIDs[ii] ,
                                                                                 gtype= annot.GSE190094$`genotype:ch1`[
                                                                                   annot.GSE190094$geo_accession ==testIDs[ii] ])
                                                          return(PDS_dat)
                                                          
                                                        }))
           
           
           # boxplot
           gg3 <- ggplot( PodoAggPerGlom_PDS, aes(x=geoID, color=gtype,
                                                  y=PDS)) + 
             geom_boxplot()+geom_point( alpha=0.4)+
             theme(legend.position = "right")+ theme_bw()+ 
             scale_color_colorblind()+
             theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
           
           
           ggn <- cowplot::plot_grid(plotlist = list(gg1,gg2,gg3), nrow  = 3)
           
           png(width = 600, height = 1200, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Spatial.trnscrptm/GSE190094/BTBR/test8/nFeatures_vs_PDSaggGloms.png",
           )
           print(ggn)
           dev.off()
           
         }
         
      
    }
  
  ### annotate podocytes by Podo-markers, run KNN filtration and glom IDtion
    {
    ### save my podo annotation - an input for KNN filtration
      {
     
      ll.seur <- list.files( pattern = "rda" , full.names = T,
                        path = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/podoClust/filtSeur/")
      names(ll.seur) <- sub( "_.*", "", basename(ll.seur))
      # # select 8 test samples
      # ll <- ll[ sub( "_.*","",basename(ll)) %in% testIDs ] 
      # seur.list.BTBR.test <- lapply( ll, readRDS)
      
      # read in OG cell annotations 
      ll.anot <-list.files( pattern = "BeadLocationsForR" , full.names = T,
                       path = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/beadLocPython_qc" )
      names( ll.anot) <-  sub( "_.*", "", basename(ll.anot))
      ll.anot <- ll.anot[names(ll.seur) ]
      OG_ctypes <- lapply( ll.anot,  read.csv , row.names=1)
      
      ### save my podo annotation - an input for KNN filtration
      # # 8 test samples
      # podoClIDs <- c(12,17,14,15,13,17,11,11)
      # 29  samples with detectable podo clusters
      podoClIDs <- unlist( seur.BTBR_podoClustID) # 14 13 17 14 17 16  9 18 14 14 12 14 11 15 18 12 15 16 17 19 15 17 15 17 15 10 14 14 16
      # podoClIDs <- c(14, 13,17, 14, 17, 16 , 9, 18, 14 ,14 ,12, 14, 11 ,15 ,18 ,12 ,15, 16, 17 ,19, 15, 17, 15 ,17, 15 ,10, 14 ,14 ,16) -1
      
      lapply( seq( podoClIDs) , function(ii){
        
        # OG annotation file 
        og_annot <- OG_ctypes[[ii]]
        
        # data
        seur <- readRDS( ll[ii])
        
        # select cells from podocyte cluster 
        tim_annot <- colnames(  seur )[ seur$seurat_clusters==as.numeric(as.character(podoClIDs[ii])) ] 
                                       
      
        og_annot$timPodo_annot <- ifelse( og_annot$barcode%in% tim_annot,
                                          "PodoTim", NA)

          write.csv( og_annot[ !is.na( og_annot$timPodo_annot) , 1:3],
                   file=  paste0("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/podoClust/BeadLocationsPodoTim_2/",
                                 names(ll.anot)[ii],"_BeadLocationsPodoTim.csv"))
          # write.csv( og_annot[ , 1:3],
          #            file=  paste0("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/podoClust/ogBeadLocations/",
          #                          names(ll.anot)[ii],"_ogBeadLocations.csv"))
      })
      
      
    }
 
    ### read results of KNN filtration, adjust global annotation accordingly
      {
      ll <-  list.files( pattern = "Tim_eroded.csv" , full.names = T,
                              path = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/podoClust/BeadLocationsPodoTim/" )
      PodoTim_filtered <- lapply( ll , read.csv , row.names=1 )
      
      lapply( seq( PodoTim_filtered) , function(ii)
        {
        print(ii)
        
        og_annot <- OG_ctypes[[ii]]
        timFilt_annot <- PodoTim_filtered[[ii]]
        # remove og podocyte identities
        og_annot$cell_type <- ifelse( og_annot$cell_type =="Podocyte" , 
                                      "oldPodo", og_annot$cell_type)
        
        # provide Tims filtrated podo idents 
        og_annot$cell_type <- ifelse( og_annot$barcode %in% timFilt_annot$barcode , 
                                      "Podocyte", og_annot$cell_type)
        # og_annot <- og_annot[ !is.na(og_annot$cell_type),] 
        og_annot$barcode <- make.unique(og_annot$barcode)
        
        write.csv( og_annot ,
                   file=  sub( "BeadLocationsPodoTim_eroded","BeadLocationsAll_eroded", ll[[ii]]))
        
      })
      
      
    }
    
    ### read results of glom annotation
     ll <-list.files( pattern = "glom_annot" , full.names = T,
                     path = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/podoClust/BeadLocationsPodoTim/" )
    TimGlom_ctypes <- lapply( ll,  read.csv , row.names=1)
    names(TimGlom_ctypes) <- sub( "_.*","",basename(ll))
    
    ll.seur <- list.files( pattern = "rda" , full.names = T,
                      path = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/podoClust/filtSeur/")
    names(ll.seur) <- sub( "_.*", "", basename(ll.seur))
    ll.seur <- ll.seur[names(TimGlom_ctypes)]
    
    ### aggregate podocytes over gloms to increase the coverage, calculate PDS
      {
      
      
      TimPodoAggPerGlom <- lapply( seq(TimGlom_ctypes),  function(ii)
        {
        print(ii)
        
        datt.glomAnnot <- TimGlom_ctypes[[ii]]
        print(names(TimGlom_ctypes)[ii])
        datt <-  readRDS( ll.seur[ ii ] )
        
        ## aggregate 
        datt$cell_type.Glom <- datt.glomAnnot$cell_type[ 
          match( colnames(datt), datt.glomAnnot$barcode)]
        datt0 <- subset(datt , subset= cell_type.Glom %in% 
                          c("Podocyte"))
        datt0$glomClust <-  datt.glomAnnot$cluster[ 
          match( colnames(datt0), datt.glomAnnot$barcode) ]
        
        XX <-  datt0@assays$Spatial@counts
        XX <- t( XX[ rowSums(XX>0, na.rm = T)>0 , ] )
        datTOplotAgg <- aggregate( XX , by = list(datt0$glomClust) , FUN=sum ) 
        rownames(datTOplotAgg) <- paste0("clust", datTOplotAgg$Group.1)
        datTOplotAgg <- t( datTOplotAgg[,colnames(datTOplotAgg)!="Group.1"])
        
        
        return( datTOplotAgg )
        
      })
      
      
      ### calculate PDS for  podocyte  aggregated per glom
      TimPodoAggPerGlom_PDS <-lapply( 
        seq(TimPodoAggPerGlom),function(ii)
          {
          set.seed(42)
          # calculate PDS seperately for podocytes and nonpodocytes
          expMat.podo  <-  TimPodoAggPerGlom[[ii]]
          # expMat.podo <- seur.BTBR.podo@assays$Spatial@counts
          expMat.podo  <- expMat.podo[ rowSums( expMat.podo >0.00000000000)> 0 ,]
          PDS <- DS_calc.func( exprMatrices = expMat.podo,
                               ceilThrsh = 0.05,
                               DSignature = DS_all , 
                               ntop = 42, wghtd = T, progStat = T)
          PDS_dat <- data.frame( glomID = colnames(expMat.podo), 
                                 PDS=PDS , 
                                 geoID = names(TimGlom_ctypes)[ii] ,
                                 gtype= annot.GSE190094$`genotype:ch1`[
                                   annot.GSE190094$geo_accession ==
                                     names(TimGlom_ctypes)[ii] ])
          return(PDS_dat)
          
  })
      
      
      ### boxplot PDS
      toPlot <- Reduce( rbind, TimPodoAggPerGlom_PDS)
      gg3 <- ggplot( toPlot, aes(x=geoID, color=gtype,
                                             y=PDS)) + 
        geom_boxplot()+geom_point( alpha=0.4)+
        theme(legend.position = "right")+ theme_bw()+ 
        scale_color_colorblind()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      
      
      ggn <- cowplot::plot_grid(plotlist = list(gg1,gg2,gg3), nrow  = 3)
      
      png(width = 600, height = 1200, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Spatial.trnscrptm/GSE190094/BTBR/test8/nFeatures_vs_PDSaggGloms_TimAnnot.png",
      )
      print(ggn)
      dev.off()
      
    }
    

    ### relate size of glom, number of podocytes and PDS
      {
      # read glom stat 
      ll <-list.files( pattern = "_glom_coord" , full.names = T,
                       path = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/podoClust/BeadLocationsPodoTim/" )
        
        TimGlom_coord_list <- lapply( ll,  read.csv , row.names=1)
      names(TimGlom_coord_list) <- sub( "_.*","" , basename(ll))
      TimGlom_coord_list <- TimGlom_coord_list[ names(TimGlom_ctypes)]
      
      # combine glom stat with PDS
      TimGlom_coord <- Reduce( rbind , lapply( seq(TimGlom_coord_list),
                                               function(ii){
        print(ii)
        datt <- cbind( TimGlom_coord_list [[ii]], TimPodoAggPerGlom_PDS[[ii]])
        
        # add count of ctypes per glom
        datt <- cbind( datt, t(as.data.frame.matrix(table(TimGlom_ctypes[[ii]][
          ,c("cell_type", "cluster")]))))
        
      }))
      
    
      
      # calculate fraction of podocytes in all glom cells
      TimGlom_coord$PodoFrct <- rowSums(TimGlom_coord[,c("Podocyte","oldPodo")])/
        rowSums(TimGlom_coord[ , c("EC" , "MC" , "Podocyte","oldPodo")
      ])
      # binarise genotype, BTBR-wt/wt == 2
      TimGlom_coord$gtypeBin <- as.numeric(factor( as.factor( TimGlom_coord$gtype), 
                                        levels = c("BTBR-wt/wt","BTBR-ob/ob")))
      TimGlom_coord$normPodo <-  rowSums(TimGlom_coord[,c("Podocyte","oldPodo")])/( 3.14 * TimGlom_coord$radii^2)
      # calculate average N of features in podocytes per sample
      featuresPERglom <-  Reduce( rbind , lapply( seq(TimGlom_ctypes),  function(ii)
        {
        print(ii)
        
        datt.glomAnnot <- TimGlom_ctypes[[ii]]
        print(names(TimGlom_ctypes)[ii])
        datt <-  readRDS( ll.seur[ ii ] )
        
        ## aggregate 
        datt$cell_type.Glom <- datt.glomAnnot$cell_type[ 
          match( colnames(datt), datt.glomAnnot$barcode)]
        datt0 <- subset(datt , subset= cell_type.Glom %in% 
                          c("Podocyte"))
        datt0$glomClust <-  datt.glomAnnot$cluster[ 
          match( colnames(datt0), datt.glomAnnot$barcode) ]
        
        XX <-  datt0@meta.data
        datTOplotAgg <- aggregate( nFeature_Spatial~ glomClust, data = XX , FUN=mean ) 
        datTOplotAgg$geoID <- names(TimGlom_ctypes)[ii]
        return( datTOplotAgg )
        
      }))
      TimGlom_coord$nFeatures <- featuresPERglom$nFeature_Spatial
      TimGlom_coord$gtype <- factor( as.factor(TimGlom_coord$gtype), 
                                     levels = c("BTBR-wt/wt","BTBR-ob/ob"))
      # saveRDS(TimGlom_coord , file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/podoClust/TimGlom_coord.rda")
        
      ### plot coorHeatmap
      {
        toPlot <- TimGlom_coord[ ,c(
          "radii","PDS","Podocyte","gtypeBin","nFeatures","PodoFrct","normPodo")]
        corMatVal <- psych::corr.test( toPlot , method = "spearman", adjust = "fdr")
        corrplot::corrplot(corMatVal$r, p.mat  = corMatVal$p , 
                           col=rev(COL2('RdBu', 200)),method = 'color',
                           # addCoef.col = "black",
                           tl.cex = 2 ,
                           sig.level = c(0.01, 0.05),insig = 'label_sig')
        
        toPlot1 <- TimGlom_coord[TimGlom_coord$gtype!="BTBR-ob/ob",c(
          "radii","PDS","Podocyte","PodoFrct","nFeatures","normPodo")]
        corMatVal1 <- psych::corr.test( toPlot1 , method = "spearman", adjust = "fdr")
        corrplot::corrplot(corMatVal1$r, p.mat = corMatVal1$p , type = "upper",
                           col=rev(COL2('RdBu', 200)),method = 'color',
                           # addCoef.col = "black",  
                           tl.cex = 2, 
                           tl.pos = "full",title = "control",
                           sig.level = c(0.01, 0.05),insig = 'label_sig')
        
        toPlot2 <- TimGlom_coord[TimGlom_coord$gtype=="BTBR-ob/ob",c(
          "radii","PDS","Podocyte","PodoFrct","nFeatures","normPodo")]
        corMatVal2 <- psych::corr.test( toPlot2 , method = "spearman", adjust = "fdr")
        corrplot::corrplot( corMatVal2$r, p.mat = corMatVal2$p , 
                            # addCoef.col = "black", 
                            tl.cex = 2 , 
                            col=rev(COL2('RdBu', 200)),
                            method = 'color',
                            add = T , type = "lower", tl.pos = "n" ,
                            sig.level = c(0.01, 0.05),insig = 'label_sig',
                            title = "BTBR-ob/ob")
      }
     
      ### plot box- and scatter- plots for PDS and glom attributes
      {
        library(rstatix)
        toPlot <- TimGlom_coord[, !colnames(TimGlom_coord)%in%
                                  c("x","y","gtypeBin","EC","MC","PDS","oldPodo","PodoFrct",
                                    "nFeatures")]
        toPlot <- reshape2::melt( toPlot ,
                                  id.vars=c("glomID","geoID","gtype") )
        stat.test <- toPlot %>%
          group_by(variable) %>%
          wilcox_test(value ~ gtype) %>%
          adjust_pvalue( method = "fdr") 
        
        
        ## boxplot
        gg1 <- ggplot( data = toPlot , aes(x=gtype, y=value , color=gtype))+
          geom_boxplot(lwd=1.5)+ scale_color_colorblind()+
          # geom_dotplot(alpha=0.5)+ 
          theme_bw()+ 
          facet_grid( rows = vars(variable), scales="free")+
          theme( axis.text.x = element_blank(),
                 axis.title.x = element_blank(),
                 # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                 text=element_text( size=20),legend.position="bottom")+
          stat_pvalue_manual(stat.test, label = "p.adj", y.position = c(350,40,0.003),size = 6)
        # facet_grid(gtype~gtype, scales = "free_y")
        
        ### scatterplots
        toPlot <- reshape2::melt( TimGlom_coord[,c("glomID","PDS","gtype","normPodo")], 
                                  id.vars=c("glomID","gtype","PDS"))
        gg2 <- ggplot2::ggplot( data=toPlot, aes(y=PDS, x=value, color=gtype))+
          geom_point( size=5, alpha=0.5)+ scale_color_colorblind()+
          geom_smooth( method=lm ,se = FALSE, lwd=2) + theme_bw()+
          theme( text = element_text( size=20),
                 legend.position = "none")+
          stat_cor(method="spearman",size=8 )+
          facet_grid( cols = vars(variable), scales = "free" )
        
        gg <- cowplot::plot_grid( plotlist = list(gg2,gg1), nrow = 1, rel_widths = c(2,1))
        
        pdf( width = 8, height = 6 , file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Spatial.trnscrptm/GSE190094/BTBR/podoClust/Glom.metrics_scatter.bxplts.test.pdf")
        print(gg)
        dev.off()
        
      
      }
     
      ### spatial plot: gloms as circles colored by PDS
      {
        library(ggforce)
        toPlotID <- "GSM5713367"
        toPlot <- TimGlom_coord[TimGlom_coord$geoID== "GSM5713367",] 
        
        cc <- ggplot(toPlot, aes(x0=x, y0=y, r=radii, color=PDS, fill=PDS)) + 
          geom_circle(alpha=0.4) + coord_equal() + theme_classic()+
          easy_remove_axes()+
          scale_color_gradientn( colours = brewer.pal( n = 9 , name = "YlOrBr"))+
          scale_fill_gradientn( colours = brewer.pal( n = 9 , name = "YlOrBr"))
        
        seur <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/podoClust/filtSeur/GSM5713367_filtSeur.rda")
        seur$glom.KNN <-  TimGlom_ctypes$GSM5713367$cell_type[ match( colnames(seur),
                                                                      TimGlom_ctypes$GSM5713367$barcode)]
        seur$glom.KNN[!is.na(seur$glom.KNN)] <- "glom" 
        
        SpatialDimPlot(seur, group.by = "glom.KNN",stroke = 0.1)+scale_fill_colorblind()+
          cc
        
        }
      
    
    }
    }
  
  
  
#### KPMP patient 29.10282 ####
  
    kpmp.sptl_29.10282 <- Read10X("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/KPMP/29-10282/deb3c13e-b49f-4812-9932-3760a3f20dd8_expression_matrix")
    kpmp.sptl_29.10282.seur <- CreateSeuratObject(kpmp.sptl_29.10282, assay="Spatial")
    
    # run normalisation and Dim reduction
    kpmp.sptl_29.10282.seur <- NormalizeData(kpmp.sptl_29.10282.seur)
    kpmp.sptl_29.10282.seur$WT1NPHS12 <- colMeans(kpmp.sptl_29.10282.seur@assays$Spatial@layers$data[
      rownames(kpmp.sptl_29.10282.seur@assays$Spatial) %in% c("WT1", "NPHS1", "NPHS2"),])
    
    ### download barcode with coordinates from loupe browser
    coord.df <- read.csv(header = T, row.names=1, sep = ",","/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/KPMP/29-10282/lowres.KPMP/Spatial-Projection.csv")
    coord.df <- coord.df[colnames(kpmp.sptl_29.10282.seur),]
    colnames(coord.df) <- c("imagerow" ,"imagecol")
    kpmp.sptl_29.10282.seur@images$image =  new(
      Class = 'SlideSeq',
      assay = "Spatial",
      key = "image_",
      coordinates = coord.df
    )
    # add recent glom annotation by Martin
    glom.annot <- read.table(header = T, row.names=1, sep = ",","/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/KPMP/29-10282/Gloms_V2_15-11-23.csv")
    kpmp.sptl_29.10282.seur$glom.annot <- glom.annot$Gloms[ match( 
      colnames(kpmp.sptl_29.10282.seur), rownames(glom.annot) )]
    kpmp.sptl_29.10282.seur$glom.annot <- ifelse( kpmp.sptl_29.10282.seur$glom.annot=="",NA,kpmp.sptl_29.10282.seur$glom.annot)

    ### calcualte PDS
    expMat <- kpmp.sptl_29.10282.seur@assays$Spatial$counts
    expMat <- expMat[rowSums(expMat>0)>0,]
    
    ### barplots for various AUCell ranking thresholds
    # gglist <- lapply(seq(from=0.05,to=0.4,by=0.05), function(ii)
    #   {
    #   print(ii)
    #   set.seed(42)
    #   kpmp.sptl_29.10282.seur$PDS.42 <- DS_calc.func( exprMatrices = expMat ,
    #                                                   ceilThrsh = ii,
    #                                                   DSignature = DS_all.HOMO , wghtd = T,
    #                                                   geneIDname = "HOMO",ntop = 42 )
    #   # boxplots
    #   gg<- ggplot(kpmp.sptl_29.10282.seur@meta.data, 
    #          aes(x=glom.annot, y=PDS.42, fill=glom.annot))+ geom_boxplot()+theme_bw()+
    #     geom_point(size=3,alpha = 0.3)+
    #     theme( text=element_text(size = 20),
    #            axis.text.x = element_blank())+
    #     annotate(geom="text",label=paste0( "AUCell.thrsh=",ii), x=2, y=0)+
    #     theme(legend.position = "none")
    #   return(gg)
    #    })
    # 
    # cowplot::plot_grid(plotlist =  gglist, ncol = 4)
    
   # use 0.15 thrshld for 42 DS for all cells
    # one glom is missed (probably healthy)?
    set.seed(42)
    kpmp.sptl_29.10282.seur$PDS.42 <- DS_calc.func( exprMatrices = expMat ,
                                                    ceilThrsh = 0.1,
                                                    DSignature = DS_all.HOMO , 
                                                    wghtd = T,
                                                    geneIDname = "HOMO",
                                                    ntop = 42 )
    
    kpmp.sptl_29.10282.seur@meta.data$glomType <-  stringr::str_sub(  kpmp.sptl_29.10282.seur@meta.data$glom.annot, end = -2)
    
    saveRDS(kpmp.sptl_29.10282.seur, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/KPMP/29-10282/kpmp.sptl_29.10282.seur.rda")
    
     # boxplots
     ggplot(kpmp.sptl_29.10282.seur@meta.data, 
                aes(x=glomType, y=PDS.42, fill=glomType))+ 
       geom_boxplot()+theme_bw()+
      geom_point(size=3,alpha = 0.3)+
      theme( text=element_text(size = 20),
             axis.text.x = element_blank())+
      annotate(geom="text",label= "AUCell.thrsh=0.3", x=2, y=0)
    
    # dimensional plot
    {
      kpmp.sptl_29.10282.seur$glom.PDS <- ifelse(kpmp.sptl_29.10282.seur$glom.annot=="",
                                                 NA, kpmp.sptl_29.10282.seur$PDS.42)
      
      gg1<- SpatialFeaturePlot( kpmp.sptl_29.10282.seur, features = "PDS.42",
                                pt.size.factor = 2 ) +
        theme(legend.position = "bottom", text = element_text(size=20))+
        scale_fill_gradientn( colours = brewer.pal( n = 9 , name = "YlOrBr"))+
        
        scale_x_reverse()
      gg2<- SpatialFeaturePlot( kpmp.sptl_29.10282.seur, features = "glom.PDS",
                                pt.size.factor = 2 ) + 
        scale_fill_gradientn( colours = brewer.pal( n = 9 , name = "YlOrBr"))+
        theme(legend.position = "bottom", text = element_text(size=20)) +
        scale_x_reverse()
      
      gg3<- SpatialDimPlot( kpmp.sptl_29.10282.seur, group.by = "glom.annot",
                            pt.size.factor = 2) + 
        theme(legend.position = "bottom", text = element_text(size=20)) +
        guides(fill = guide_legend(override.aes = list(size=8)))+
        scale_x_reverse()
      gg4<- SpatialFeaturePlot( kpmp.sptl_29.10282.seur, features =  "WT1NPHS12",
                                pt.size.factor = 2) + 
        scale_fill_gradientn( colours = brewer.pal( n = 9 , name = "YlOrBr"))+
        theme(legend.position = "bottom", text = element_text(size=20))+
        scale_x_reverse()
        # scale_fill_viridis_b(option = "A",direction = -1,n.breaks=7,)
      
      cowplot::plot_grid(plotlist = list( gg3 , gg2 , gg4), ncol = 1)
    
      pd
    }
    



    
  
  