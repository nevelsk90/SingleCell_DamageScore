#=====### PDS analysis in human data #=========###
options( connectionObserver = NULL )

.libPaths(c("/home/tim_nevelsk/R/x86_64-pc-linux-gnu-library/4.0", "/media/tim_nevelsk/WD_tim/SOFT/R"))

library( biomaRt )
library( Matrix )
library( GSEABase )
library( ggplot2 )
library( Seurat )
library( org.Hs.eg.db )
library(AnnotationDbi )
library( ggpubr )
library( ggthemes )
# release memory
mallinfo::malloc.trim()
gc()

#### load DS ####
# setwd("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/HumanValidation")
# calculate damage score
source("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/cell-damage-score/AUCell_script.r")
source("/home/tim_nevelsk/PROJECTS/myCode/usefulRfunc.r")

 ### use all inclusive signature
  # # human gene annotation
  ensembl_human = useMart("ensembl",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
  tx2gene_human <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name',"entrezgene_id"),  mart = ensembl_human)

  # load damage signatures
  DS_all <- read.table( header = T,  file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")
  # convert gene names to human orthologs
  DS.HOMO <-  fun_homoTO.FROMmouse( gns = DS_all$gene_symbol, TO = F)
  DS_all.HOMO <-  DS_all
  DS_all.HOMO$HOMO <- unlist(DS.HOMO$HOMO)



#### NEPTUNE Marray: GSE104066 ####
  
    # do analysis with updated annotation
    XX <- GEOquery::getGEO( "GSE104066" )
    Marray_expr <- XX$GSE104066_series_matrix.txt.gz@assayData$exprs
    # convert rownames to gene symbols
    rownames(Marray_expr)  <- XX$GSE104066_series_matrix.txt.gz@featureData@data$ENTREZ_GENE_ID[ match( 
      rownames(Marray_expr)  , XX$GSE104066_series_matrix.txt.gz@featureData@data$ID)]
    rownames(Marray_expr) <- mapIds(org.Hs.eg.db, rownames(Marray_expr) , "SYMBOL", "ENTREZID")
    Marray_expr <- Marray_expr[ !is.na( rownames(Marray_expr)),]
    Marray_expr <- Marray_expr[log(rowMeans( Marray_expr) )>1, ]
    # calculate the DAMAGE score
    Marray_PDS  <- DS_calc.func( exprMatrices = Marray_expr,
                                 # DSignature = DS_all.HOMO ,
                                 DSignature = DS_all.HOMO ,
                                 geneIDname = "HOMO", 
                                 ceilThrsh = 0.05 , 
                                 ntop = 42, 
                                 wghtd = T )
    
    Marray_PDS <- cbind.data.frame( as.data.frame(Marray_PDS) , phenotype = as.factor( XX$GSE104066_series_matrix.txt.gz@phenoData@data$`diagnosis:ch1`))
    Marray_PDS$phenotype <- factor(Marray_PDS$phenotype ,   
                                   labels = c("FSGS", "IgAN", "Healthy"))
    Marray_PDS$phenotype <- relevel(  Marray_PDS$phenotype , ref = "Healthy" )
    Marray_PDS$gtypeDE <- ifelse(Marray_PDS$phenotype == "Healthy", "control","experimental")
    ### plot jitter
    library( ggpubr )
  
    # toPlot1 <- reshape2::melt( Marray_PDS)
    gg1 <- ggplot2::ggplot( Marray_PDS , aes( y=Marray_PDS, 
                                              x=phenotype, 
                                              color=phenotype,
                                              groups=phenotype)) +
      scale_color_colorblind() +  
      geom_boxplot()+
      geom_jitter(width = 0.1, size=6) + 
      ggtitle("GSE104066 Micro-array")+
      theme_bw() + theme(text = element_text(size = 24) , legend.position = "none") +
      stat_compare_means( size = 8, test="t.test", comparisons = list(c("Healthy","FSGS"),c("Healthy","IgAN")))
    gg1
  
  
#### NEPTUNE : GSE140989  ####
### Single cell transcriptomics identifies focal segmental glomerulosclerosis 
### remission endothelial biomarker. JCI Insight 2020 Mar 26;5(6). PMID: 32107344

    ### primary data processing
    {
      # XX <- read.table( file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/NEPTUNE/GSE140989_KPMP_Premiere_reference_SCT_Normalizeddata.txt", row.names = 1, header = T, sep = "\t")
      # # sort files in folders
      # ll <- list.files(full.names = T,"/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/NEPTUNE/scRNAseq")
      # dirnames <- unique( sub("-genes.*|-features.*|-matrix.*|-barcodes.*","",ll) )
      # lapply(dirnames, dir.create )
      # # modify file names
      # ll <- list.files(full.names = T, recursive = T, "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/NEPTUNE/scRNAseq")
      # lapply( ll , function(x) file.rename( x, paste( dirname(x),"/",sub(".*-","",x), sep = "")) )
      # read the data
      dirnames <- list.dirs(recursive = F, path = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/GSE140989_NEPTUNE/scRNAseq")
      names(dirnames) <- sub( "_.*", "", basename(dirnames) )
      datt <-  Read10X( data.dir = dirnames , strip.suffix = T)
      sce <- CreateSeuratObject( datt, min.cells = 20, min.features = 200 )
      
      # add mitochondrial percentage
      sce$percent.mt <- PercentageFeatureSet(object = sce, pattern = "^MT")
      sce <- subset(sce, subset = percent.mt < 25  )
      # add metadata
      sce$sample <- sub( "_.*", "", colnames(sce) )
      sce$gtype <- ifelse( sce$sample %in%c("GSM4191952","GSM4191953","GSM4191954"),"preperfusion" , 
                           ifelse( sce$sample %in%c("GSM4191960","GSM4191961","GSM4191962","GSM4191963","GSM4191964"),
                                   "allograft","tumor_nephrectomy") )
      saveRDS(sce, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/GSE140989_NEPTUNE/GSE140989_Human.NEPTUNE.scRNAseq_Seurat.rds")
      
    }
    
    ### basic data processing and Viz
    {
      sce <- readRDS( file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/GSE140989_NEPTUNE/scRNAseq/Seurat/GSE140989_Human.NEPTUNE.scRNAseq_Seurat.rds")
      
      ### normialise the data
      sce <- NormalizeData( sce, normalization.method = "LogNormalize", scale.factor = 10000)
      # ### find variable features
      sce  <- FindVariableFeatures( sce , selection.method = "vst", nfeatures = 2000)
      ### scale the data
      sce <- ScaleData(sce )
      
      ### run PCA
      sce <- RunPCA(sce, features = VariableFeatures(object = sce))
      ElbowPlot(sce , ndims = 50 ) # choosing number of PCs to include 
      # DimPlot( sce, reduction = "pca", group.by = "gtype", order =  sample(colnames(sce)) )
      
      # cluster
      sce <- FindNeighbors( sce, dims = 1:20 )
      sce <- FindClusters( sce, resolution = 0.1 )
      
      # run UMAP
      sce <- RunUMAP( sce, dims = 1:20 )
      # saveRDS(sce, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/NEPTUNE/GSE140989_Human.NEPTUNE.scRNAseq_SeuratClust.rds")
      
      # subset
      set.seed(42)
      sce_subs <- subset( sce ,  downsample = 1000  ) 
      # Idents(sce_subs) <- sce_subs@meta.data$clust
      p1 <- DimPlot( sce_subs, reduction = "umap", order =  sample(colnames(sce_subs)) , label = T)
      # p2 <- DimPlot( sce_subs, reduction = "umap", group.by = "sample", shuffle = T )
      # p3 <- DimPlot( sce_subs, reduction = "umap", group.by = "batch", order =  sample(colnames(sce_subs)) )
      p4 <- DimPlot( sce_subs, reduction = "umap", group.by = "gtype", shuffle = T )
      p5 <-   FeaturePlot( sce_subs, reduction = "umap", order = F,features = "nCount_RNA" , 
                           max.cutoff = 15000)
      p6 <-   FeaturePlot( sce_subs, reduction = "umap", order = F, features = "WT1" , 
                           min.cutoff = 0 , max.cutoff = 20)
      p7 <-   FeaturePlot( sce_subs, reduction = "umap", order = F, features = "NPHS2" , 
                           min.cutoff = 0 , max.cutoff = 20 )
      p8 <-   FeaturePlot( sce_subs, reduction = "umap", order = F, features = "MEIS2" , 
                           min.cutoff = 0 , max.cutoff = 20 )
      
      cowplot::plot_grid(p1,p4, p5,p6, p7, p8, ncol = 3)
      
      # subset podocytes
      sceHomoPodo <- subset(sce, ident=7)
      sceGloms <- subset(sce, ident=c(4, 6, 7), downsample = 500)
      
      DimPlot( sceHomoPodo, reduction = "umap", group.by = "gtype", shuffle = T ) +
        xlim(c(1,4))+ylim(c(-13,-10))
      saveRDS(sceHomoPodo , file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/NEPTUNE/Seurat/GSE140989_Human.NEPTUNE.scRNAseq_Seurat.Podocytes.rds")
      
      
    }
    
    ### (pseudo) bulk
    {
      ### read bulk RNAseq data
        ll <- list.files(full.names = T, path = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/GSE140989_NEPTUNE/bulk")
      GSE140989_bulk <- Reduce( cbind, lapply( ll , function(xx){
        read.table(gzfile(xx), row.names = 1)
      }))
      colnames(GSE140989_bulk) <- sub(  ".txt.gz","", basename(ll) )
      
      
      ### generate pseudo-bulk from snRNAseq
      sceHomoPodo <- readRDS( file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/GSE140989_NEPTUNE/scRNAseq/Seurat/GSE140989_Human.NEPTUNE.scRNAseq_Seurat.Podocytes.rds")
      sceHomoPodo_pbulk <-  cbind.data.frame( sample=(sceHomoPodo$sample), t(as.matrix(sceHomoPodo@assays$RNA@counts)))
      sceHomoPodo_pbulkM <- aggregate( .~sample , data = sceHomoPodo_pbulk,  FUN = sum )
      XX <- data.frame( row.names = colnames(sceHomoPodo_pbulkM)[-1], 
                        t(sceHomoPodo_pbulkM[,-1]) )
      colnames( XX ) <- sceHomoPodo_pbulkM$sample
      
      ### combine bulk and pbulk
      GSE140989_bulk$gName <- tx2gene_human$external_gene_name[ match( 
        rownames(GSE140989_bulk) , tx2gene_human$ensembl_gene_id )]
      GSE140989_both <- cbind.data.frame(  XX , GSE140989_bulk[ match(sub( "\\..*","", rownames(XX)), (GSE140989_bulk$gName)),] )
      GSE140989_both <- GSE140989_both[ colnames(GSE140989_both) !="gName"]
      GSE140989_both[is.na(GSE140989_both)] <- 0
      
      ### calculate the score add PDS as a metadata
      GSE140989_PDS  <- DS_calc.func( exprMatrices = GSE140989_both,
                                      DSignature =  DS_all.HOMO ,
                                      geneIDname = "HOMO",
                                      ntop = 42,
                                      wghtd = T)
      GSE140989_PDS <- cbind.data.frame( as.data.frame(GSE140989_PDS) ,
                                         phenotype = as.factor( 
                                           c( rep("control",ncol(XX)),
                                              rep("experiment",ncol(GSE140989_bulk)-1))))
      
      ### plot distributions
      library( ggpubr )

      gg2 <- ggplot2::ggplot( GSE140989_PDS , 
                              aes( y = GSE140989_PDS, 
                                   x = phenotype, 
                                   color= phenotype)) +
        scale_color_colorblind() +  
        geom_boxplot()+
        geom_jitter(width = 0.1, size=6) + 
        ggtitle("GSE140989 bulk/pseudo-bulk")+
        theme_bw() + theme(text = element_text(size = 24) , legend.position = "none") +
        stat_compare_means( size = 8)
      gg2
      
      
    }
    
  
  
#### GSE176465  ####
  ### ! there are very little podocytes for some patients
    {
    ll <- list.dirs( recursive = T, 
                     "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/GSE176465_scRNAseq/GSE176465_RAW")
    ll <- ll[ grepl( "GRCh38", ll) ]
    # # rename genes.tsv.gz to  features.tsv.gz to work with Seurat
    # old.names <- list.files( pattern = "genes", recursive = T, full.names = T,
    #                          path = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/GSE176465_scRNAseq/GSE176465_RAW")
    # new.names <- gsub("genes", "features", old.names)
    # file.rename(old.names, new.names)
    
    datt <- Read10X( data.dir = ll ,
                     gene.column = 2 ,
                     cell.column = 1 ,
                     unique.features = TRUE , 
                     strip.suffix = T )
    
    
    # analyse with Seurat
    seu <- CreateSeuratObject( datt , project ="GSE176465", )
    saveRDS(seu , file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/GSE176465_scRNAseq/GSE176465_raw.seurat.rda"   )   
    
    # process
    seu <- readRDS(file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/GSE176465_scRNAseq/GSE176465_raw.seurat.rda"  )
    # decide cutoff
    UMI.counts <- colSums( seu@assays$RNA@layers$counts)
    UMI.counts <- (UMI.counts[order(- UMI.counts )])
    plot( log10(1:length(UMI.counts)), log10(UMI.counts),   type="l")
    nFeature_RNA <- seu$nFeature_RNA
    nFeature_RNA <- nFeature_RNA[ order(-nFeature_RNA)]
    plot( log10(1:length( nFeature_RNA)), log10(nFeature_RNA),   type="l")
    
    # filter
    seu$UMI.counts <-  colSums( seu@assays$RNA@layers$counts)
    seu.filt <- subset(seu , subset= nFeature_RNA>200 )
    
    # analyse
    seu.filt <- NormalizeData( seu.filt )
    seu.filt <- FindVariableFeatures(seu.filt, verbose = T )
    seu.filt <- ScaleData( seu.filt ) 
    
    seu.filt <- RunPCA( seu.filt, verbose = T )
    ElbowPlot( seu.filt , ndims = 50)
    seu.filt <- FindNeighbors( seu.filt, dims = 1:20)
    seu.filt <- FindClusters( seu.filt, resolution = 0.5 )
    seu.filt <- RunUMAP( seu.filt, dims = 1:20 )

    ### check expression of podocytes
    Seurat::FeaturePlot( seu.filt , 
                         order=T , min.cutoff = 0,
                         features = "NPHS2" , 
                         label = T, cols=c("lightgray","red"))
   
    # add annotation 
    GSM <- sub( "_.*", "",
                basename(list.dirs( recursive = F ,  
    path = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/GSE176465_scRNAseq/GSE176465_RAW")) ) 
    seu.filt$GSM <- factor( sub("_.*","",colnames(seu.filt)) ,
                       levels= unique(sub("_.*","",colnames(seu.filt))) , 
                       labels = GSM )
    annot <- read.table( header = T, sep = ",",
                         "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/GSE176465_scRNAseq/SraRunTable.txt")
    PatientAnnot <- read.table( header = T, sep = "\t",
                                "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/GSE176465_scRNAseq/GSE176465_patientTOsample.csv")
    seu.filt$patient <- factor(  sub("_.*","",colnames(seu.filt)) ,
                           levels = unique(sub("_.*","",colnames(seu.filt)) ) ,
                           labels = PatientAnnot$Subject )
    seu.filt$biopsy_diagnosis <- factor(  sub("_.*","",colnames(seu.filt)) ,
                            levels = unique(sub("_.*","",colnames(seu.filt)) ) ,
                            labels = PatientAnnot$Subject )
    
    # add biopsy diagnose
    XX <- unique(annot[, c("biopsy_diagnosis", "Sample.Name")])
    seu.filt$biopsy_diagnosis <- XX$biopsy_diagnosis[
      match( seu.filt$GSM, XX$Sample.Name ) ]
    saveRDS( seu.filt  , file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/GSE176465_scRNAseq/GSE176465_seurat.rda")
    # 
    # 
    
  }
    
  ## plot PDS
    {
      seu.filt <- readRDS( file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/GSE176465_scRNAseq/GSE176465_seurat.rda" )
  
      Seurat::FeaturePlot( seu.filt ,
                       order=T , min.cutoff = 0,
                       features =  "WT1" ,
                       label = T, cols=c("lightgray","red"))

      seu.podo <- subset( seu.filt , subset= seurat_clusters == 9 )
      expMat <- seu.podo@assays$RNA@layers$counts
      rownames(expMat)<- rownames(seu.podo@assays$RNA)
      colnames(expMat)<- colnames(seu.podo@assays$RNA)
      expMat <- expMat[rowSums(expMat>0)> 0 ,] # >0 for scRNAseq, >10 for snRNAseq

      ### calcualte PDS, ceilThrsh must be >0.05 for HUMAN scRNAseq and snRNAseq
      set.seed(42)
      seu.podo$PDS  <- DS_calc.func( exprMatrices = expMat ,
                                                      ceilThrsh = 0.1 ,
                                                      DSignature = DS_all.HOMO,
                                                      wghtd = T,
                                                      progStat = T,
                                                      geneIDname = "HOMO",
                                                      ntop = 42)
      
      # aggregate
      seu.podo_pbulkM <- aggregate( .~patient+ biopsy_diagnosis+GSM,
                                    data = seu.podo@meta.data[
                                      , c("PDS","patient","biopsy_diagnosis","GSM")],
                                    FUN = mean )

      gg4 <- ggplot2::ggplot(  seu.podo_pbulkM  , aes(
        y= PDS, x=biopsy_diagnosis, color=patient)) +
        scale_color_colorblind() +
        geom_jitter(width = 0.1, size=6) +
        ggtitle("GSE176465 pseudo-bulk")+
        theme_bw() + theme(text = element_text(size = 24) ,
                           legend.position = "none")
         # stat_summary(fun.y= median, fun.ymin=median, fun.ymax=median, geom = "crossbar", width = .2, color = "red")
        # stat_compare_means( size = 8)
      gg4

      # correlate with UPCR
      seu.podo_pbulk.sbj <- aggregate( .~patient + biopsy_diagnosis,
                                          data = seu.podo@meta.data[
                                            , c("PDS","patient","biopsy_diagnosis")],
                                          FUN = mean )
      upcr.tab <- read.table(header = T, sep = "\t","/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/GSE176465_scRNAseq/patient_data.csv")
      seu.podo_pbulk.sbj$UPCR <- as.numeric( upcr.tab$UPCR..at.urine.collection.[
        match( seu.podo_pbulk.sbj$patient , upcr.tab$Study.participants) ])
      
      seu.podo_pbulk.sbj$UPCR.highest <- as.numeric( upcr.tab$UPCR..highest.recorded.[
        match( seu.podo_pbulk.sbj$patient , upcr.tab$Study.participants ) ])
      seu.podo_pbulk.sbj$
      
      cor.test( as.numeric(seu.podo_pbulk.sbj$UPCR),
                as.numeric(seu.podo_pbulk.sbj$PDS), method="spearman")
        
      ggplot( seu.podo_pbulk.sbj , aes( x = log10( UPCR ),
                                        y = PDS ))+
        geom_point(size=6, aes( color=biopsy_diagnosis))+theme_bw()+ 
        theme( text = element_text(size = 24), legend.position = "bottom")+ggtitle("GSE176465 pseudo-bulk")+
        geom_smooth( method = "lm", color="black",se = F)+
        stat_cor( method = "spearman", size=8)
        
      }

#### KPMP snRNAsea https://atlas.kpmp.org/explorer/ ####
  
    ### annotation 
      {
      KPMP.annotation <- read.csv( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/metadata/e20f1612-849c-488e-a06c-3049098a4bb1_20230914_OpenAccessClinicalData.csv")  
      { 
        KPMP.annotation$AlbCr  <- factor( as.factor(KPMP.annotation$Albuminuria..mg...Binned. ),
                                              labels = c(NA,30,1500,150,400,750))
        KPMP.annotation$AlbCr <- as.numeric(as.character(KPMP.annotation$AlbCr  ))
        
        KPMP.annotation$Proteinuria  <- factor( as.factor(KPMP.annotation$Proteinuria..mg...Binned. ),
                                                labels = c(NA,100,1500,350,750))
        KPMP.annotation$Proteinuria <- as.numeric(as.character( KPMP.annotation$Proteinuria  ))
        
        KPMP.annotation$Age  <- factor( as.factor(KPMP.annotation$Age..Years...Binned. ),
                                        labels = c(NA,15,25,35,45,55,65,75,85))
        KPMP.annotation$Age <- as.numeric(as.character( KPMP.annotation$Age  ))
        
        KPMP.annotation$A1c  <- factor( as.factor(KPMP.annotation$A1c......Binned. ),
                                        labels = c(NA,6,9,7,8))
        KPMP.annotation$A1c <- as.numeric(as.character( KPMP.annotation$A1c  ))
        
        KPMP.annotation$eGFR  <- as.numeric( sub( "-.*", "", KPMP.annotation$Baseline.eGFR..ml.min.1.73m2...Binned. )) + 5 
       
        KPMP.annotation$Diabetes.Duration..Years.  <- as.numeric( as.character( factor( 
          as.factor(KPMP.annotation$Diabetes.Duration..Years. ) , labels = c(0,2,12,17,22,27,32,37,42,7,52,57,62))))
        KPMP.annotation$Hypertension.Duration..Years.  <- as.numeric( as.character( factor( 
          as.factor( KPMP.annotation$Hypertension.Duration..Years.) , labels = c(0,2,12,17,22,27,32,37,42,47,7))))
        
        
         }
      
      ### original albinuria measurments
       eGFR_continuos_BC_AU <-  readxl::read_xlsx("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/metadata/OP-BCAU5812-PlasmaBiomarkerDataFile-2022.xlsx")
       AlbCre_continuos_BC_AU <-  readxl::read_xlsx( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/metadata/OP-BC_AU5812-UrineBiomarkerData-2022.xlsx")
       AlbCre_continuos_MSDQ <-  readxl::read_xlsx( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/metadata/OP-MSDQ120-UrineBiomarkerData-2022.xlsx")
       
       # calculate ACR for MSDQ samples
       AlbCre_continuos_MSDQ$ACR <- AlbCre_continuos_MSDQ$`Urine Albumin Concentration (mg/dL)`/
         (AlbCre_continuos_MSDQ$`Urine Creatinine Concentration (mg/dL)`*0.001)
       # folter BC_AU samples
       AlbCre_continuos_BC_AU <- AlbCre_continuos_BC_AU[
         AlbCre_continuos_BC_AU$`Sample Type`=="Spot Urine Supernatant" &
           AlbCre_continuos_BC_AU$`Participant Visit`=="Enrollment",]
       
       # combine in one dataframe
       AlbCre_continuos<- cbind.data.frame(AlbCre_continuos_BC_AU[,c(1:6,13)], 
                                           ACR_MSDQ= AlbCre_continuos_MSDQ$ACR[ 
                                             match( AlbCre_continuos_BC_AU$Participant.ID, 
                                                              AlbCre_continuos_MSDQ$Participant.ID)],
                                           eGFR.CrCys= eGFR_continuos_BC_AU$eGFRCrCys[
                                             match( AlbCre_continuos_BC_AU$Participant.ID,
                                                    eGFR_continuos_BC_AU$Participant.ID)] )
       AlbCre_continuos$BiopsUrin.lag <- (AlbCre_continuos$`Biopsy Time Interval Days Post Consent`+1)-
         (AlbCre_continuos$`Sample Collection Interval Days Post Consent`+1)
       AlbCre_continuos$ACR.mean <- rowMeans(cbind( AlbCre_continuos$`ACR (mg/g)`,
                                                AlbCre_continuos$ACR_MSDQ ), na.rm=T)
     
       ### update original annot table with continuous eGFR and AlbCr
        KPMP.annotation.contn <- cbind( KPMP.annotation , 
                                  AlbCre_continuos[ 
                match( KPMP.annotation$Participant.ID , 
                       AlbCre_continuos$Participant.ID ), 7:11 ])

      
      # ## which RNAseq datasets have albuminuria or proteinuria measurments
      # {
      #   IDintrst <- unique( KPMP.annotation$Participant.ID[(KPMP.annotation$Albuminuria..mg...Binned.)!=""] )
      #   IDintrst.albcr.cont  <- KPMP.annotation$Participant.ID[ !is.na(KPMP.annotation$AlbCr.cont)]
      #   
      #   KPMP.sn.avail <- read.table(sep = ",",header = T, fill=T, 
      #                               "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/snRNAseq/atlas_repository_filelist-20230926.csv")
      #   IDavail.sn <- unique(KPMP.sn.avail$Participant.ID[KPMP.sn.avail$Workflow.Type=="Expression Matrix"])
      #   
      #   KPMP.sc.avail <- read.table(sep = ",",header = T, fill=T, 
      #                               "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/snRNAseq/atlas_repository_filelist-20230926_scRNAseq.csv")
      #   IDavail.sc <- unique(KPMP.sc.avail$Participant.ID[KPMP.sc.avail$Workflow.Type=="Expression Matrix"])
      #   
      #   IDintrst.all <- union( IDintrst, IDintrst.albcr.cont)
      #   IDretrive <- data.frame( Participant.ID = IDintrst2, 
      #                            # older.sn = IDintrst %in% unique(sceHomoPodo$partID) ,
      #                            IDretrive.sn_AlbCr=IDintrst2 %in% intersect(IDintrst , IDavail.sn ),
      #                            IDretrive.sc_AlbCr=IDintrst2 %in%  intersect(IDintrst , IDavail.sc) , 
      #                            IDretrive.sn_AlbCr.Prt=IDintrst2 %in%  IDavail.sn ,
      #                            IDretrive.sc_AlbCr.Prt=IDintrst2 %in%  IDavail.sc )
      #  
      #   write.table(IDretrive, sep = "\t", row.names = F, quote = F,
      #               file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/AlbCre_patientID_snsc_AlbCr.Prt.tsv")
      #   
      # }
     
    }
 
    ### primary processing
      {
      
      podoMarks <- c("WT1","NPHS1","NPHS2","SYNPO","ROBO2","THSD7A","NEBL","PTPRO","MAGI2","PODXL")
      
      ### snRNAseq
        {
        ### read expression matrix and make a seurat obj
        lldir <- list.dirs( recursive = F, 
                            path = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/snRNAseq/expression.matrix")
        annot.mtx <- read.table( header = T,sep = ",","/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/snRNAseq/atlas_repository_filelist-20231204.csv")
        
        ### create seur obj
        file.issue <- list.files( path = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/snRNAseq/expression.matrix", recursive = T)
       file.issue[! basename( file.issue) %in% c("matrix.mtx.gz",
                                                  "features.tsv.gz",
                                                  "barcodes.tsv.gz")]
        
        seur.list <- lapply( seq(lldir), function(ii){
          print(ii)
          datt <- tryCatch(  Read10X(lldir[[ii]]), error = function(e) NA )
          if( is.list(datt)) datt <- datt$`Gene Expression`
          dim(datt)
          # extract all cells with nonzero expression of podoMarks
          datt <- tryCatch( {
            datt <- datt[ , colSums(datt[ podoMarks,]>0)>3 ]
            datt[ rowSums( datt > 0 )> 0 , ]
          } , error = function(e) NA )
          print( dim(datt))
          
          seu <- tryCatch( {
            seu<- CreateSeuratObject( datt ) 
            subset( seu, subset= nFeature_RNA>=500 & nFeature_RNA<10000 )
          } , 
            error = function(e) NA )
          return(seu)
        })
        # kpmp.ckdID <- KPMP.annotation$Participant.ID[KPMP.annotation$Tissue.Type=="CKD" ]
        names( seur.list ) <- annot.mtx$Participant.ID[ match( basename(lldir) ,
                                                             sub( ".zip", "", annot.mtx$File.Name ) ) ]
        names( seur.list ) <- make.unique( names(seur.list) )
         seur.list <- seur.list[!is.na(seur.list)]
        
        saveRDS(seur.list, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/snRNAseq/Seurat/KPMP_snRNAseq.AlbCr_PodoMark.list.08.12.23.rda") 

        # merge all selected
        seur.list <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/snRNAseq/Seurat/KPMP_snRNAseq.AlbCr_PodoMark.list.08.12.23.rda" )

         seur <- merge( seur.list[[1]],
                       y = seur.list[2:length(seur.list)],
                       add.cell.ids = names(seur.list)
                      )
         
        seur$partID <- sub("_.*","",colnames(seur))

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
        seur <- FindNeighbors( seur, dims = 1:10 )
        seur <- FindClusters( seur, resolution = 0.1 )
        
        # run UMAP
        seur <- RunUMAP( seur, dims = 1:10 )
        seur$PodoMarks <- colMeans(seur@assays$RNA@data[ podoMarks,],
          na.rm = T)
        seur$gtype <- KPMP.annotation$Tissue.Type[ 
          match(  seur$partID,KPMP.annotation$Participant.ID)]
        seur$Protocol<- KPMP.annotation$Protocol[ 
          match(  seur$partID,KPMP.annotation$Participant.ID)]
        # saveRDS(seur, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/snRNAseq/Seurat/KPMP_snRNAseq.AlbCr.CKD_dimred.27.09.23.rda")
        # subset for plotting
        seur.sub <- subset( seur , downsample= 500 )
        # find cluster markers
        seur.clustMarks <- FindAllMarkers( seur.sub )
        top10 <- seur.clustMarks %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
        DoHeatmap( seur.sub, features =top10$gene )
        
        # make plots
        # seur <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/snRNAseq/Seurat/KPMP_snRNAseq.AlbCr_PodoMarks.rda")
        p1 <- DimPlot(seur.sub, label = T, shuffle = T)
        p2 <- DimPlot(seur.sub, shuffle = T,label = T, group.by = "partID")
        
        p3 <- FeaturePlot(seur.sub, features = "PodoMarks")
        p4 <- DimPlot(seur.sub, shuffle = T,label = T, group.by = "gtype")

        cowplot::plot_grid(plotlist = list(p1,p2,p3,p4), nrow = 2)
        
        # select podo by cluster
        seur.podo <- subset( seur, id=5)
        saveRDS(seur.podo, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/snRNAseq/Seurat/KPMP_snRNAseq.AlbCr_Podo_04.12.23.rda") 
        
      }
    
      ### scRNAseq
        {
        
        # sc
        lldir.sc <- list.dirs( recursive = F, 
                               path = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/scRNAseq/expression.matrix")
          annot.mtx.sc <- read.table( header = T,sep = ",","/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/scRNAseq/atlas_repository_filelist-20230927_scRNAseq.csv")
        
        # first check if you can read all files, modify folder names if that's the problem
          file.issue <- list.files( path = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/scRNAseq/expression.matrix", recursive = T)
          file.issue <- file.issue[ ! basename( file.issue)%in%c("matrix.mtx.gz",
                                                    "features.tsv.gz",
                                                    "barcodes.tsv.gz")]
          
        ### create seur obj
        seur.sc.list <- lapply( seq(lldir.sc), function(ii)
          {
          print(ii)
          datt <- tryCatch(  Read10X(lldir.sc[[ii]]), error = function(e) NA )
          if( is.list(datt)) datt <- datt$`Gene Expression`
          # extract all cells with nonzero expression of podoMarks
          datt <- tryCatch( {
            datt <- datt[ , colSums(datt[ podoMarks,]>0)>3 ]
            datt[ rowSums( datt > 0 )> 0 , ]
          } , error = function(e) NA )
          print( dim(datt))
          
          seu <- tryCatch(  {   
            seu<- CreateSeuratObject( datt ) 
            subset( seu, subset= nFeature_RNA>=500 & nFeature_RNA<10000 )
          } , 
            error = function(e) NA )
          return(seu)
        })
               # kpmp.ckdID <- KPMP.annotation$Participant.ID[KPMP.annotation$Tissue.Type=="CKD" ]
        names( seur.sc.list ) <- annot.mtx.sc$Participant.ID[ match( basename(lldir.sc) ,
                                                             sub( ".zip", "", annot.mtx.sc$File.Name ) ) ]
        seur.sc.list <- seur.sc.list[!is.na(seur.sc.list)]
        
        saveRDS( seur.sc.list , file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/scRNAseq/Seurat/KPMP_scRNAseq.AlbCr_PodoMark.08.12.23.rda") 
        
        # merge all selected
        seur <- merge( seur.sc.list[[1]],
                       y = seur.sc.list[2:length(seur.sc.list)],
                       add.cell.ids = names(seur.sc.list)
        )
        seur$partID <- sub("_.*","",colnames(seur))
        seur <- subset( seur , subset= nFeature_RNA>500 & nFeature_RNA<10000)
        # saveRDS(seur , file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/scRNAseq/Seurat/KPMP_scRNAseq.AlbCr_PodoMarks.04.12.23.rda") 
        
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
        seur$PodoMarks <- colMeans(seur@assays$RNA@data[ podoMarks,], na.rm = T)
        saveRDS(seur, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/scRNAseq/Seurat/KPMP_scRNAseq.AlbCr_PodoMark.08.12.23.rda")
        
        # make plots
        # seur <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/scRNAseq/Seurat/KPMP_scRNAseq.AlbCr_PodoMarks.04.12.23.rda")
        p1 <- DimPlot(seur, label = T, shuffle = T)
        p2 <- DimPlot(seur, shuffle = T,label = T, group.by = "partID")
        p3 <- FeaturePlot(seur, features = "PodoMarks")

        cowplot::plot_grid(plotlist = list(p1,p2,p3), nrow = 1)
        
        # select podo by cluster
        seur.podo <- subset( seur, id=6 )
        saveRDS(seur.podo, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/scRNAseq/Seurat/KPMP_scRNAseq.AlbCr_Podo_08.12.23.rda") 
        
        }
      
      # ### extract.podocytes based on podo marks
      # seur.list <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/snRNAseq/Seurat/KPMP_snRNAseq.AlbCr_27.09.23.rda")
      # 
      # KPMP.snRNAseq_podo.Seur  <- lapply(seq(seur.list), function(ii){
      #   datt<- seur.list[[ii]]
      #   datt <- NormalizeData(datt)
      #   datt$WT1NPHS12 <- colMeans(datt@assays$RNA@data[c("WT1","NPHS1","NPHS2"),],na.rm = T)
      #   podo <- subset(datt, subset=WT1NPHS12>0.1)
      #   return(podo)
      # })
      # names(KPMP.snRNAseq_podo.Seur) <- annot.mtx$Participant.ID[match(
      #   basename(lldir), sub( ".zip", "", annot.mtx$File.Name)
      # )]
      # # combine in one obj
      # KPMP.snRNAseq_podo.Seur <- merge( KPMP.snRNAseq_podo.Seur[[1]],
      #                                   y=KPMP.snRNAseq_podo.Seur[2:length(KPMP.snRNAseq_podo.Seur)],
      #                                   add.cell.ids=names(KPMP.snRNAseq_podo.Seur)
      # )
      # KPMP.snRNAseq_podo.Seur$partID <- sub("_.*","",colnames(KPMP.snRNAseq_podo.Seur))
      # saveRDS(KPMP.snRNAseq_podo.Seur, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/snRNAseq/Seurat/KPMP_snRNAseq.AlbCr_podo_27.09.23.rda") 
      
    }
   
    ### calculate PDS for sc and snRNAseq
      {
    
     KPMP.snRNAseq_podo.Seur <- readRDS( file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/snRNAseq/Seurat/KPMP_snRNAseq.AlbCr_Podo_08.12.23.rda") 
     KPMP.scRNAseq_podo.Seur <- readRDS( file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/scRNAseq/Seurat/KPMP_scRNAseq.AlbCr_Podo_08.12.23.rda") 
     seqType.sn <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/HUMAN/KPMP/snRNAseq/Seurat/seur.sn.types.rda")
     seqType.sn$partID <- make.unique(seqType.sn$partID)
     # add seq type
     KPMP.snRNAseq_podo.Seur$seqType <- seqType.sn$seqType[ 
       match( KPMP.snRNAseq_podo.Seur$partID , seqType.sn$partID) ] 
     KPMP.scRNAseq_podo.Seur$seqType <- "scRNAseq"
     # fix nonexisting IDs
     KPMP.snRNAseq_podo.Seur$partID[
       KPMP.snRNAseq_podo.Seur$partID=="31-10210.1"
     ] <- "31-10210"

      # combine sc and sn
      KPMP.scsnRNAseq_podo.Seur <- merge( KPMP.snRNAseq_podo.Seur , KPMP.scRNAseq_podo.Seur)
      
     ## add  combination of patient ID and metadata
     KPMP.scsnRNAseq_podo.Seur$partIDsnsc <- paste(KPMP.scsnRNAseq_podo.Seur$partID, 
                                               KPMP.scsnRNAseq_podo.Seur$seqType, sep = "_")

     
     # calculate PDS
     expMat <- KPMP.scsnRNAseq_podo.Seur@assays$RNA@counts
      expMat <- expMat[rowSums(expMat>0)> 0 ,] # >0 for scRNAseq, >10 for snRNAseq
       
       ### calcualte PDS, ceilThrsh must be >0.05 for HUMAN scRNAseq and snRNAseq
       set.seed(42)
       KPMP.scsnRNAseq_podo.Seur$PDS  <- DS_calc.func( exprMatrices = expMat ,
                                                       ceilThrsh = 0.1 ,
                                                       DSignature = DS_all.HOMO, 
                                                       wghtd = T,
                                                       progStat = T,
                                                       geneIDname = "HOMO", 
                                                       ntop = 42)
       

       KPMP.scsnRNAseq_podo.Seur$ACR <- KPMP.annotation.contn$ACR.mean[
          match( KPMP.scsnRNAseq_podo.Seur$partID ,KPMP.annotation.contn$Participant.ID)]
    
       # aggregate
       XX <- KPMP.scsnRNAseq_podo.Seur@meta.data[ , c("PDS","partIDsnsc") ]
       # XX <- XX[ (XX$partID) %in% names(which( table(XX$partID)>10)), ]
       KPMP.PDS.agg <- aggregate(. ~ partIDsnsc, data=XX, FUN=mean, 
                                 na.action = na.pass)
       # plot 
       toPlot <- KPMP.PDS.agg
       toPlot$seqType <- sub( ".*_","",toPlot$partIDsnsc )
       toPlot$partID <- sub( "_.*", "", toPlot$partIDsnsc)
       
       toPlot <- cbind(toPlot, KPMP.annotation.contn[
         match( toPlot$partID ,KPMP.annotation.contn$Participant.ID),
         -c(7, 10:13)])
       saveRDS( toPlot , "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/HumanValidation/KPMP.all_scsnRNAseq_podo.PDS.aggregated.rda")
       
       ### remove samples with fewer than 3 cells
       X <- table( KPMP.scsnRNAseq_podo.Seur$partIDsnsc)
       toPlot.filt <- toPlot[ toPlot$partIDsnsc %in% names(X)[X>=3] &
                                # toPlot$Tissue.Type != "DM-R"  &
                                toPlot$Tissue.Type %in% c("CKD","AKI")  &
                                # toPlot$seqType=="snRNAseq" ,
                                # ,,"Healthy Reference"
                                # toPlot$BiopsUrin.lag >0 &
                                toPlot$Age < 75 ,  # filter samples younger than 75 for snRNAseq
       ]
       toPlot.filt$ACR.combo <- ifelse( is.na(toPlot.filt$ACR.mean) ,  toPlot.filt$AlbCr, toPlot.filt$ACR.mean)
       
    }
       
    ### dot plot for pbulk PDS
      {   
        
      # plot 
    toPlot <- aggregate(. ~ partIDsnsc, 
                        data= KPMP.scsnRNAseq_podo.Seur@meta.data[ 
                          , c("PDS","partIDsnsc") ] , 
                        FUN=mean, na.action = na.pass)
    toPlot$seqType <- sub( ".*_","",toPlot$partIDsnsc )
    toPlot$partID <- sub( "_.*", "", toPlot$partIDsnsc)
    
    toPlot <- cbind(toPlot, KPMP.annotation.contn[
      match( toPlot$partID ,KPMP.annotation.contn$Participant.ID),
      -c(7, 10:13)])
    # filter samples with fewer than 3 cells
    X <- table( KPMP.scsnRNAseq_podo.Seur$partIDsnsc)
    toPlot.filt <- toPlot
    toPlot.filt <- toPlot[ toPlot$partIDsnsc %in% names(X)[X>=3]  &
                             toPlot$Tissue.Type %in% c("Healthy Reference","CKD","AKI")  
                             # toPlot$seqType=="snRNAseq" ,
                             # ,,"Healthy Reference"
                             # toPlot$BiopsUrin.lag >0 &
                             # toPlot$Age < 75 
                           ,    ]
    
    toPlot.filt$Tissue.Type <- relevel( factor(toPlot.filt$Tissue.Type), ref = "Healthy Reference")
    
    gg3 <- ggplot2::ggplot(  toPlot.filt  , aes( 
      y= PDS, x=Tissue.Type, color=Tissue.Type)) +
      scale_color_colorblind() +        
      geom_boxplot()+
      geom_jitter(width = 0.1, size=6) + 
      ggtitle("KPMP pseudo-bulk")+
      theme_bw() + theme(text = element_text(size = 24) , 
                         legend.position = "none") +
      stat_compare_means( size = 8, 
                           method = "t.test",
                          comparisons = list(c("Healthy Reference", "CKD"),c("Healthy Reference","AKI")))
    gg3
  }
  
 
    ### correlation with clinical features
      {
    library(corrplot)
        toPlot.filt <-readRDS( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/HumanValidation/KPMP.CKD.AKI_scsnRNAseq_podo.PDS.aggregated.rda")
        
    ### all traits
    {         
      corMat <- Reduce( cbind.data.frame, 
                        lapply( 1:ncol(toPlot.filt),  function(ii){
                          X <- toPlot.filt[,ii]
                          if( !is.numeric(X)) { 
                            X <-  as.numeric( as.factor(X))
                            return(X)} else return(X)
                        }))
      rownames(corMat) <- toPlot.filt$partIDsnsc
      colnames(corMat) <- colnames(toPlot.filt)
      corMat <- as.data.frame( corMat[,apply(corMat, 2, sd,na.rm=T)>0] )
      
      XX <-  corMat[ corMat$seqType==1 , 
                     !colnames(corMat)%in% c("partIDsnsc","partID","Participant.ID","seqType")]
      XX <- as.data.frame( XX[,apply(XX, 2, sd,na.rm=T)>0] )
      corMatVal.sc.test <-  psych::corr.test( XX,
                                              method = "spearman", 
                                              use = "pairwise.complete.obs" )
      corMatVal.sc <- corMatVal.sc.test$r
      corMatVal.sc[is.na(corMatVal.sc)]<- 0
      
      
      XX2 <-  corMat[ corMat$seqType==2 , 
                      !colnames(corMat)%in% c("partIDsnsc","partID","Participant.ID","seqType")]
      XX2 <- as.data.frame( XX2[,apply(XX2, 2, sd,na.rm=T)>0] )
      corMatVal.sn.test <-  psych::corr.test( XX2,
                                              # "seqType","Tissue.Type")], 
                                              method = "spearman", 
                                              use = "pairwise.complete.obs" )
      corMatVal.sn <- corMatVal.sn.test$r
      corMatVal.sn[is.na(corMatVal.sn)]<- 0
      
      corMat.test <-  psych::corr.test( corMat[ ,!colnames(corMat)%in% c("partIDsnsc","partID",
                                                                         "Participant.ID")], 
                                        method = "spearman", adjust = "none",
                                        use = "pairwise.complete.obs" )
      corMat.test$r[is.na(corMat.test$r)]<- 0
      
      # make a plot
      corrplot::corrplot(corMatVal.sc, type = "upper", 
                         tl.pos = "full", 
                         order="original",
                         method = 'color',  col=rev(COL2('RdBu', 200)),
                         sig.level = 0.05,insig = 'label_sig',
                         p.mat =  corMatVal.sc.test$p, 
                         # title = "KPMP scRNAseq",
                         tl.cex = 1.5, cl.cex = 1.2, pch.cex = 1.5 )
      corrplot::corrplot(corMatVal.sn, type = "lower", tl.pos = "n", 
                         order="original",  col=rev(COL2('RdBu', 200)),
                         method = 'color',insig = 'label_sig',
                         sig.level = 0.05,p.mat =  corMatVal.sn.test$p, 
                         add = T,   
                         # title = "KPMP snRNAseq",
                         tl.cex = 1.5, cl.cex = 1.2, pch.cex = 1.5 )
      
      corrplot::corrplot(corMat.test$r,mar = c(0,0,0,0),
                         order="hclust",
                         col=rev(COL2('RdBu', 200)),
                         method = 'color', tl.pos = "l",
                         p.mat =corMat.test$p,
                         sig.level = 0.05,
                         # title = "KPMP sn.sc.RNAseq",
                         insig = 'label_sig',
                         tl.cex = 1.5, cl.cex = 1.2, pch.cex = 1.5 )
      
      
    }
    
    ### selected trait
    {
   
      toPlot.filt.red <- toPlot.filt[,c("PDS" ,
                                        "Age", "Sex", "Race", "Diabetes.Duration..Years.",
                                        "Hypertension.Duration..Years.", "eGFR.CrCys",
                                        "ACR.mean")]
      cor.test(toPlot.filt.red$PDS, toPlot.filt.red$ACR.mean, method = "spearman")
      
      corMat.red <- Reduce( cbind.data.frame, 
                            lapply( 1:ncol(toPlot.filt.red),  
                                    function(ii){
                                      X <- toPlot.filt.red[,ii]
                                      if( !is.numeric(X)) { 
                                        X <-  as.numeric( as.factor(X))
                                        return(X)} else return(X)
                                    }))
      rownames(corMat.red) <- toPlot.filt.red$partIDsnsc
      colnames(corMat.red) <- c("PDS" ,"Age", "Sex", "Race", "Diabetes Duration",
                                "Hypertension Duration", "eGFR Cr.Cys",
                                "ACR mg/g")
      corMat.red <- as.data.frame( corMat.red[,apply(corMat.red, 2, sd,na.rm=T)>0] )
      
      corMat.test <-  psych::corr.test( corMat.red,  method = "spearman", 
                                        use = "pairwise.complete.obs" , adjust = "none")
      # make a plot
      corrplot::corrplot(corMat.test$r,order="hclust", 
                         # type = "lower",
                         col= rev(COL2('RdBu', 200)),
                         tl.col = "black",
                         tl.pos = "l", 
                         method = 'color',
                         p.mat = corMat.test$p ,
                         sig.level = 0.05,
                         insig = 'label_sig',
                         tl.cex = 1.0, cl.cex = 1, pch.cex = 2)
      
    }
    
    ### boxplots
    {
      toPlot.boxplot <- toPlot[ toPlot$partIDsnsc %in% names(X)[X>=3] &
                                  
                                  #  toPlot$seqType=="snRNAseq" &
                                  # ,,"Healthy Reference"
                                  # toPlot$BiopsUrin.lag >0 &
                                  toPlot$Age != 75 ,  # filter samples younger than 75 for snRNAseq
      ]
      
      ggplot( data =toPlot.boxplot, aes(y=PDS, x=seqType, color=Tissue.Type))+
        geom_boxplot()
    }
        
     ### regression plot Alb.Cre.contn
      {
          ggplot(toPlot.filt ,
                 aes(x=PDS,
                     y=log(AlbCr),
                     color=seqType
                 ) )+ 
            geom_point( size=5, alpha=0.5)+ 
            geom_smooth( method=lm ,se = FALSE) + theme_bw()+
            theme( text = element_text( size=20))
          ggplot(toPlot.filt ,
                 aes(x=PDS,
                     y=log(ACR.mean),
                     color=seqType
                 ) )+ 
            geom_point( size=5, alpha=0.5)+ 
            geom_smooth( method=lm ,se = FALSE) + theme_bw()+
            theme( text = element_text( size=20))
        }
        
    
  }
       
  
  
  
  

#### save dotplots plots for multiple studies   ####
  gg <- cowplot::plot_grid( gg1, gg2 , gg3 , align = "h" , nrow = 1)
  # save the plot
  pdf(height = 8, width = 16, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/HumanValidation/PDS.42_Human.podo.jitter.v3.pdf")
  print(gg )
  dev.off()


    
    
    
    
    
    
    
    