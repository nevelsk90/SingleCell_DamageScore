###=====# PDS publication figure 3 code #=======###
# an input for each figure should be a data, 
# which requires no significant computations for plotting
# release memory
gc()
mallinfo::malloc.trim()
gc()


#### load DS and code ####

.libPaths(c("/home/tim_nevelsk/R/x86_64-pc-linux-gnu-library/4.0", 
            "/media/tim_nevelsk/WD_tim/SOFT/R"))
options( connectionObserver = NULL )

# library( org.Mm.eg.db )
library( ggplot2 )
library( cowplot)
library( reshape2 )
library( plyr )
library( Seurat )
library( viridis )
library( ggthemes )
library( AUCell )
library( ggpubr )
library( biomaRt )
library( RColorBrewer )
library( ggrepel )
library( scales )
library( ComplexHeatmap )

# setwd
setwd("/home/tim_nevelsk/PROJECTS/PODOCYTE/WRITING/PDS_manuscript/Figures/Figure4")

# source necessary code
source("https://raw.githubusercontent.com/PauUng/HepatocyteDamageScore/master/SharedFunctions.R")
source("https://raw.githubusercontent.com/nevelsk90/R_scripts/master/func_analysis.R")  

# load damage signatures
DS_all <- read.table( header = T,  file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")
# load 
allPodoGenes <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/SCSN_allPodoGenes.rda")
genecorrPDS.Sprmn_freq <- readRDS(  file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/genecorrPDS.Sprmn_freq.16.05.24.rda")
# expression data 

#### 1. visualising and analysing TRN ####
{

  
  ### load the prior

  ## ATACseq prior
  # ATACseq_tgenesM.TFtc <- readRDS( file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/ATACseq_TFtc.tgenesM_77.TF.rda")
  # ATACseq_tgenesM <- readRDS( file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/ATACseq_tgenesM.rda")
  ATACseq_tgenes <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/ATACseq/Podocytes/TOBIAS/FIMO/atac.podo_tobias.fimo110TF.cut1.p5e4_TFtc.rda")
  ATACseq_tgenes.S <- Reduce( cbind , lapply( ATACseq_tgenes , "[[", "score" ))
  ATACseq_tgenes.Q <- Reduce( cbind , lapply( ATACseq_tgenes , "[[", "qvalue" ))
  colnames(ATACseq_tgenes.S) <- colnames(ATACseq_tgenes.Q) <- names(ATACseq_tgenes)
  ATACseq_tgenesM <- ATACseq_tgenes.S * ( ATACseq_tgenes.Q <0.1)
  rownames(ATACseq_tgenesM) <- rownames(ATACseq_tgenes.S) <- rownames(ATACseq_tgenes.Q) <- rownames(ATACseq_tgenes[[1]])
  # filter for podocyte expressed
   ATACseq_tgenesM <- ATACseq_tgenesM[ rowSums(ATACseq_tgenesM)>0 , colSums(ATACseq_tgenesM)>0  ]

  ATACseq_tgenesM.podo <- ATACseq_tgenesM[ rownames( ATACseq_tgenesM ) %in% allPodoGenes ,  ]
  
  # # include Chip-seq results
  # # ## Wt1 and Tead1 Chipseq targets
  # WT1wt2014_TFtc <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/Chipseq/Chipseq.WT1wt2014_TFtc.rda")
  # WT1wt2014_TFtc.sig <- WT1wt2014_TFtc$score * ( WT1wt2014_TFtc$qvalue <0.3)
  # Tead1wt2018_TFtc <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/Chipseq/Chipseq.Tead1wt2018_TFtc.rda")
  # Tead1wt2018_TFtc.sig <-  Tead1wt2018_TFtc$score * ( Tead1wt2018_TFtc$qvalue <0.3)
  # identical(rownames(WT1wt2014_TFtc), rownames(ATACseq_tgenesM))
  # ATAC.Chip <- ATACseq_tgenesM
  # ATAC.Chip[ , "Cluster_1__Wt1"] <- (ATAC.Chip[ , "Cluster_1__Wt1"] + WT1wt2014_TFtc.sig)/2
  # ATAC.Chip[ , "Cluster_5__Tead1"] <- (ATAC.Chip[ , "Cluster_5__Tead1"] + Tead1wt2018_TFtc.sig)/2
  # ATACseq_tgenesM <- ATAC.Chip[ rowSums(ATAC.Chip)>0 , colSums(ATAC.Chip)>0  ]
  # ATACseq_tgenesM.podo <- ATACseq_tgenesM[ rownames( ATACseq_tgenesM ) %in% allPodoGenes ,  ]

  
  ### black-white heatmap of the whole TRN
  {
    toPlot <- ATACseq_tgenesM.podo
    # toPlot <- toPlot[ rowSums(toPlot)>0 , colSums(toPlot)>0  ]
    # toPlot <- toPlot[ rownames( toPlot ) %in% allPodoGenes ,  ]
    toPlot <- apply(toPlot , 2, rank,ties.method= "min")
    
    
    gg1 <- pheatmap::pheatmap(toPlot, treeheight_row = 100, show_rownames = F,
                              clustering_method = "ward.D2",
                              col= brewer.pal(9, "Greys"), treeheight_col = 0)
    
    
    
    # 
    
    pdf( width =  5 , height = 6 , file="TRN.14.TFtcCut0.1_blwh.heatmap.pdf")
    a<-dev.cur()
    png( width =  400 , height = 400 , file="TRN.14.TFtcCut0.1_blwh.heatmap.png")
    dev.control("enable")
    print( gg1 )
    dev.copy(which=a)
    dev.off()
    dev.off()
    
    
    
    
  }
  
  ### annotate TRN target gene clusters 
  {
    ### cluster 
    annotdf <-  as.dendrogram(gg1$tree_row)
    TRN_geneClust <- sort(cutree( dend, h=0.97e5))
    
    cols <- rainbow(16)
    names(cols) <- unique(TRN_geneClust)
    annR <- as.data.frame(TRN_geneClust)
    annR$TRN_geneClust <- as.factor( annR $TRN_geneClust)
    gg2 <- pheatmap::pheatmap(toPlot, treeheight_row = 100, show_rownames = F,
                              clustering_method = "ward.D2",
                              col= brewer.pal(9, "Greys"), treeheight_col = 0,
                              annotation_row = annR ,
                              cutree_rows = 16,
                              annotation_colors= list(TRN_geneClust = cols) )
    
    pdf( width =  8 , height = 8 , file="Supl.Fig4/TRN.14.TFtcCut0.1_blwh.heatmap.clust.pdf")
    a<-dev.cur()
    png( width =  600 , height = 600 , file="Supl.Fig4/TRN.14.TFtcCut0.1_blwh.heatmap.clust.png")
    dev.control("enable")
    print( gg2 )
    dev.copy(which=a)
    dev.off()
    dev.off()
    
    ###  annotate clusts
    ### use Robert's function for GO analysis
    {
      library(DBI)
      
      # run function on list of DE results
      GOrobert_res.cut10 <- lapply( unique(TRN_geneClust), function(ii)
      {
        print(ii)
        
        # datt <- names(TRN_geneClust)[TRN_geneClust==ii]
        # define background gene set
        universe= unique( entr2gName$entrezgene_id) 
        universe <- as.character( universe[!is.na(universe)] )
        
        
        # prepare gene set, convert ensembleIDs to entrezIDs
        geneset= names(TRN_geneClust)[TRN_geneClust==ii]
        geneset <- unique(entr2gName$entrezgene_id[match( geneset,  entr2gName$external_gene_name)])
        geneset <- geneset[!is.na(geneset)]
        geneset <- as.character(geneset)
        print(length(intersect(geneset,colnames(gomatrix))))
        
        # apply Robert's function that given a sparse matrix of GO terms (columns = genes, rows = GO terms)
        # a geneset of interest and a background set of genes (universe)
        # carry out clustering with members diverging by at most cut_max genes, and do enrichment testing.
        # Note, multiplicity adjustment is performed for the representative terms only.
        RobertGO <- sf.clusterGoByGeneset( gomatrix, geneset, universe, 
                                           min.genes=5, cut_max = 10 )
        
        return(RobertGO)
      })
      names(GOrobert_res.cut10) <- paste0( "clust_", unique(TRN_geneClust))
      saveRDS(GOrobert_res.cut10, file = "Supl.Fig4/TRNpodo.16clust_GOrobert.cut10.rda")
      # cut 50
      GOrobert_res.cut50 <- lapply( unique(TRN_geneClust), function(ii)
      {
        print(ii)
        
        # datt <- names(TRN_geneClust)[TRN_geneClust==ii]
        # define background gene set
        universe= unique( entr2gName$entrezgene_id) 
        universe <- as.character( universe[!is.na(universe)] )
        
        
        # prepare gene set, convert ensembleIDs to entrezIDs
        geneset= names(TRN_geneClust)[TRN_geneClust==ii]
        geneset <- unique(entr2gName$entrezgene_id[match( geneset,  entr2gName$external_gene_name)])
        geneset <- geneset[!is.na(geneset)]
        geneset <- as.character(geneset)
        print(length(intersect(geneset,colnames(gomatrix))))
        
        # apply Robert's function that given a sparse matrix of GO terms (columns = genes, rows = GO terms)
        # a geneset of interest and a background set of genes (universe)
        # carry out clustering with members diverging by at most cut_max genes, and do enrichment testing.
        # Note, multiplicity adjustment is performed for the representative terms only.
        RobertGO <- sf.clusterGoByGeneset( gomatrix, geneset, universe, 
                                           min.genes=5, cut_max = 50 )
        
        return(RobertGO)
      })
      names(GOrobert_res.cut50) <- paste0( "clust_", unique(TRN_geneClust))
      saveRDS(GOrobert_res.cut50, file = "Supl.Fig4/TRNpodo.16clust_GOrobert.cut50.rda")
      
      
      ### top functions per model
      # get a union of top N func in each dataset, select from primary terms only
      GOrobert_res <- GOrobert_res.cut50
      iids <- seq( GOrobert_res)
      thrsh <- 0.1
      
      topN_GOrb  <- Reduce( union , lapply( iids , function( ii,N=5 )
      {
        print(ii)
        datt <- GOrobert_res[[ii]]
        # select from primary terms only
        datt <- subset(datt$results,  Is.primary==TRUE & Primary.Fisher.adj < thrsh )
        datt <- datt[ order(datt$Fisher) ,]
        print(datt$GO.ID[ 1:min(N, nrow(datt)) ])
        # return N top primar terms
        return(datt$GO.ID[ 1:min(N, nrow(datt)) ])
      }))
      topN_GOrb <- topN_GOrb[!is.na(topN_GOrb)]
      
      gg<- GOrobert_barplot_clustered( iid = iids ,
                                       GOtoPlot = topN_GOrb,
                                       datGO = GOrobert_res,
                                       heatmap = T, 
                                       datOrder = names(GOrobert_res)[
                                         as.numeric(sub(".*_","",pp$tree_col$labels[pp$tree_col$order]))  
                                       ])
      
      ### cluster functions
      ss <- reshape(gg$data[,c( "Term","dataset","log10.Fisher")],
                    idvar = "Term", timevar = "dataset", direction = "wide")
      rownames(ss) <- ss$Term
      ss <- ss[,-1]
      pp<-  pheatmap::pheatmap(ss, cluster_rows = F, scale = "column")
      
      
      ## plot
      pdf( width =  16 , height = 14 , file="Supl.Fig4/TRNpodo.16clust_GOrobert.cut50.top5.pdf")
      a<-dev.cur()
      png( width =  1000 , height = 1000 , file="Supl.Fig4/TRNpodo.16clust_GOrobert.cut50.top5.png")
      dev.control("enable")
      
      print(gg )
      
      dev.copy(which=a)
      dev.off()
      dev.off()
      
      
      
      
    }
    
    
    ### check enrichment of relevant paths
    {
      pthsPlt <- read.table( sep = "\t","/home/tim_nevelsk/PROJECTS/PODOCYTE/WRITING/PDS_manuscript/Figures/Figure3/KFO_select.paths.txt")
      
      pthsPlt <- pthsPlt[!(pthsPlt$V1 %in% c( "Muscarinic acetylcholine receptors *REACT",
                                              "Metabolism *REACT")),]
      pthsPlt <-sub(" \\*" ,"__",pthsPlt)
      
      
      ### fisher test
      {
        PopSize <- 20000
        
        TRN_geneClust_pathAnnot.hyper<- Reduce( cbind,  lapply( 
          unique(TRN_geneClust), function(ii)
          {
            
            print(ii)
            list2<-length( names(TRN_geneClust)[TRN_geneClust==ii] )
            
            setNames( sapply( seq(pthsPlt), function(jj){
              
              pathGns <- pathDB[[ pthsPlt[[jj]] ]]
              
              list1 <-  length(  pathGns  )
              
              overlap <- length( intersect(  names(TRN_geneClust)[TRN_geneClust==ii],
                                             pathGns ) )
              
              phyper( overlap -1 , list2 , PopSize , list1  , lower.tail= FALSE)
            }) , nm = pthsPlt)
            
            
          } ) )
        
        
        colnames(TRN_geneClust_pathAnnot.hyper) <- 1:16
        
        
        
        # TRN_geneClust_pathAnnot.hyper.adj <- apply(TRN_geneClust_pathAnnot.hyper, 2, p.adjust , "fdr")
        # toPlot <- TRN_geneClust_pathAnnot.hyper.adj[ rowSums(TRN_geneClust_pathAnnot.hyper.adj<0.1)>0,]
        toPlot <- apply( TRN_geneClust_pathAnnot.hyper, 2, rank)
        
        sigTab <- TRN_geneClust_pathAnnot.hyper.adj
        sigTab <- ifelse( sigTab < 0.1, "*" ,"")
        gg3 <-   pheatmap::pheatmap(  log(toPlot) , cluster_cols = F, 
                                      display_numbers = sigTab,  
                                      fontsize_number=20, fontsize = 15,
                                      number_color = "turquoise",
                                      scale = "none" , color = rev(brewer.pal(name="Reds", n=9)))
        
        pdf( width =  12 , height = 10 , file="Supl.Fig4/TRNpodo_16TGclust_selP.Fisher.RankHeatmap.pdf")
        a<-dev.cur()
        png( width =  800 , height = 600 , file="Supl.Fig4/TRNpodo_16TGclust_selP.Fisher.RankHeatmap.png")
        dev.control("enable")
        print( gg3 )
        dev.copy(which=a)
        dev.off()
        dev.off()
      }
    }
    
    
    
    
    
  }
  
  
  ### check enrichment of target gene sets in selected.paths
  {
    pthsPlt <- read.table( sep = "\t","/home/tim_nevelsk/PROJECTS/PODOCYTE/WRITING/PDS_manuscript/Figures/Figure3/KFO_select.paths.txt")
    
    pthsPlt <- pthsPlt[!(pthsPlt$V1 %in% c( "Muscarinic acetylcholine receptors *REACT",
                                            "Metabolism *REACT")),]
    pthsPlt <-sub(" \\*" ,"__",pthsPlt)
    
    
    ### fisher test
    {
      PopSize <- 50000
      
      TRN_pathAnnot.hyper<- Reduce( cbind,  lapply( 
        1:ncol(ATACseq_tgenesM.podo), function(ii)
        {
          
          print(ii)
          
          tgenes <- rownames(ATACseq_tgenesM.podo)[ ATACseq_tgenesM.podo[,ii]>0 ]
          
          print(length(tgenes))
          list2<-length( tgenes )
          print(paste0(colnames(ATACseq_tgenesM.podo)[ii],"_",list2))
          pp <- setNames( sapply( seq(pthsPlt), function(jj){
            
            pathGns <- pathDB[[ pthsPlt[[jj]] ]]
            
            list1 <-  length(  pathGns  )
            
            overlap <- length( intersect(  tgenes,
                                           pathGns ) )
            
            phyper( overlap -1 , list2 , PopSize , list1  , lower.tail= FALSE)
          }) , nm = pthsPlt)
          
          # print(pp)
          return(pp)
          
          
        } ) )
      
      
      
      colnames(TRN_pathAnnot.hyper) <- colnames(ATACseq_tgenesM.podo)
      
      
      
      TRN_pathAnnot.hyper.adj <- apply(TRN_pathAnnot.hyper, 2, p.adjust , "fdr")
      # toPlot <- TRN_geneClust_pathAnnot.hyper.adj[ rowSums(TRN_geneClust_pathAnnot.hyper.adj<0.1)>0,]
      toPlot <- apply( TRN_pathAnnot.hyper, 2, rank)
      
      sigTab <- TRN_pathAnnot.hyper.adj
      sigTab <- ifelse( sigTab < 0.05, "*" ,"")
      gg3 <-   pheatmap::pheatmap(  log(toPlot) , cluster_cols = F, 
                                    display_numbers = sigTab,  
                                    fontsize_number=20, fontsize = 15,
                                    number_color = "turquoise",
                                    scale = "none" , color = rev(brewer.pal(name="Reds", n=9)))
      
      pdf( width =  12 , height = 10 , file="Supl.Fig4/TRNpodo_16TGclust_selP.Fisher.RankHeatmap.pdf")
      a<-dev.cur()
      png( width =  800 , height = 600 , file="Supl.Fig4/TRNpodo_16TGclust_selP.Fisher.RankHeatmap.png")
      dev.control("enable")
      print( gg3 )
      dev.copy(which=a)
      dev.off()
      dev.off()
    }
    
    #### Fgsea 
    {
      library(fgsea)
      # run GSEA with shrunk lfc
      TRN_fgsea_res <- lapply( 1:ncol(ATACseq_tgenes.S) , function(ii){
        print(ii)
        
        TRNnet <- podoGenesFACS
        datt <- setNames( ATACseq_tgenes.S[ , ii] , rownames(ATACseq_tgenes.S))
        datt <- datt[order(-datt)]   
        fgsea_res <- fgseaMultilevel(pathways = pathDB[ pthsPlt], nproc=4,
                                     stats    = datt , eps = 1e-20,
                                     minSize  = 3  , maxSize = 1000 )
        fgsea_res$testName <- colnames(ATACseq_tgenes.S)[ii]
        
        return(fgsea_res)
      })
      
      # extrct results for different databases
      names(TRN_fgsea_res) <-  sub( ".*__","",colnames(ATACseq_tgenes.S))
      
      # ### plot 
      # iids <- seq(TRN_fgsea_res)
      # fgsea_res <- TRN_fgsea_res
      # topN_fgsea  <- Reduce( union , lapply( iids , function( ii,N=3 )
      #   {
      #   print(ii)
      #   datt <- fgsea_res[[ii]]
      #   # select from primary terms only
      #   datt <- subset(datt,  pval < thrsh )
      #   datt <- datt[ order(datt$pval) ,]
      #   # return N top primar terms
      #   return(datt$pathway[ 1:min(N, nrow(datt)) ])
      # }))
      # topN_fgsea <- topN_fgsea[!topN_fgsea %in% c("salmon", "chartreuse3","Slit.Diaphr","Actin.Ctsklt.fig1f")]
      # 
      # gg1 <- fsgsea_barplot_clustered( iid=iids,
      #                                  pathToPlot= topN_fgsea,
      #                                  dat.fgsea = fgsea_res ,
      #                                  datOrder=names(TRN_fgsea_res) ,
      #                                  labelTRUNK=50 )
      # 
      # gg1
      
      TRN_fgsea_pval <- Reduce( cbind, lapply( seq(TRN_fgsea_res), function(ii){
        datt <- setNames( TRN_fgsea_res[[ii]]$pval , nm =  TRN_fgsea_res[[ii]]$pathway)
        datt <- datt[ match (pthsPlt,  names(datt))]
        names(datt) <- pthsPlt
        datt[is.na(datt)]<- 1
        return(datt)
      }))
      colnames(TRN_fgsea_pval) <- names(TRN_fgsea_res)
      
      
      sigTab <- TRN_fgsea_pval
      sigTab <- ifelse( sigTab < 0.05, "*" ,"")
      
      # to plot
      toPlot <- apply( TRN_fgsea_pval, 2, rank)
      gg3 <-   pheatmap::pheatmap(  log(toPlot) , cluster_cols = F, 
                                    display_numbers = sigTab,  
                                    fontsize_number=20, fontsize = 15,
                                    number_color = "turquoise",
                                    scale = "none" , color = rev(brewer.pal(name="Reds", n=9)))
      
      ### plot fraction
      toPlot <-     Reduce( cbind,  lapply( 
        1:ncol(ATACseq_tgenesM.podo), function(ii)
        {
          
          print(ii)
          
          tgenes <- rownames(ATACseq_tgenesM.podo)[ ATACseq_tgenesM.podo[,ii]>0 ]
          print( length(tgenes))
          
          pp <- setNames( sapply( seq(pthsPlt), function(jj){
            
            pathGns <- pathDB[[ pthsPlt[[jj]] ]]
            
            
            overlapFr <- length( intersect(  tgenes ,
                                             pathGns ) )/length(pathGns)
            
            
          }) , nm = pthsPlt)
          
          # print(pp)
          return(pp)
          
          
        } ) )
      toPlot[is.na(toPlot)] <- 0
      
      
      # toPlot <- apply( toPlot, 2, rank)
      colnames(toPlot) <-  sub( ".*__","",colnames(ATACseq_tgenes.S))
      toPlot <- toPlot[rowSums(toPlot)>0,]
      
      sigTab <- TRN_fgsea_pval
      sigTab <- ifelse( sigTab < 0.05, "*" ,"")
      gg4 <-   pheatmap::pheatmap(  (toPlot) , cluster_cols = F, 
                                    display_numbers = sigTab[ rowSums(TRN_fgsea_pval)!=14,],  
                                    fontsize_number=20, fontsize = 15,
                                    number_color = "turquoise",
                                    scale = "none" , color = (brewer.pal(name="Reds", n=9)))
      
      pdf( width =  12 , height = 8 , file="Supl.Fig4/TRNpodo_TFtgenes_OvlpFr.fgseaStar.Heatmap.pdf")
      a<-dev.cur()
      png( width =  800 , height = 500 , file="Supl.Fig4/TRNpodo_TFtgenes_OvlpFr.fgseaStar.Heatmap.png")
      dev.control("enable")
      print( gg4 )
      dev.copy(which=a)
      dev.off()
      dev.off()
      
    }
    
    ###  annotate  unique targets
    ### use Robert's function for GO analysis
    {
      library(DBI)
      source("/home/tim_nevelsk/PROJECTS/myCode/func_analysis.R")  
      
      
      # run function on list of DE results
      GOrobert_res.cut10 <- lapply( 1:ncol(ATACseq_tgenesM.podo), function(ii)
      {
        print(ii)
        
        # datt <- names(TRN_geneClust)[TRN_geneClust==ii]
        # define background gene set
        universe= unique( entr2gName$entrezgene_id) 
        universe <- as.character( universe[!is.na(universe)] )
        
        
        # prepare gene set, convert ensembleIDs to entrezIDs
        geneset=  rownames(ATACseq_tgenesM.podo)[ ATACseq_tgenesM.podo[,ii]>0 & 
                                                    rowSums(ATACseq_tgenesM.podo>0)==1  ]
        
        geneset <- unique(entr2gName$entrezgene_id[match( geneset,  entr2gName$external_gene_name)])
        geneset <- geneset[!is.na(geneset)]
        geneset <- as.character(geneset)
        print(length(intersect(geneset,colnames(gomatrix))))
        
        # apply Robert's function that given a sparse matrix of GO terms (columns = genes, rows = GO terms)
        # a geneset of interest and a background set of genes (universe)
        # carry out clustering with members diverging by at most cut_max genes, and do enrichment testing.
        # Note, multiplicity adjustment is performed for the representative terms only.
        RobertGO <- sf.clusterGoByGeneset( gomatrix, geneset, universe, 
                                           min.genes=5, cut_max = 10 )
        
        return(RobertGO)
      })
      names(GOrobert_res.cut10) <- paste0( "clust_", unique(TRN_geneClust))
      saveRDS(GOrobert_res.cut10, file = "Supl.Fig4/TRNpodo.uniqueTG_GOrobert.cut10.rda")
      # cut 50
      GOrobert_res.cut50 <- lapply( unique(TRN_geneClust), function(ii)
      {
        print(ii)
        
        # datt <- names(TRN_geneClust)[TRN_geneClust==ii]
        # define background gene set
        universe= unique( entr2gName$entrezgene_id) 
        universe <- as.character( universe[!is.na(universe)] )
        
        
        # prepare gene set, convert ensembleIDs to entrezIDs
        geneset= names(TRN_geneClust)[TRN_geneClust==ii]
        geneset <- unique(entr2gName$entrezgene_id[match( geneset,  entr2gName$external_gene_name)])
        geneset <- geneset[!is.na(geneset)]
        geneset <- as.character(geneset)
        print(length(intersect(geneset,colnames(gomatrix))))
        
        # apply Robert's function that given a sparse matrix of GO terms (columns = genes, rows = GO terms)
        # a geneset of interest and a background set of genes (universe)
        # carry out clustering with members diverging by at most cut_max genes, and do enrichment testing.
        # Note, multiplicity adjustment is performed for the representative terms only.
        RobertGO <- sf.clusterGoByGeneset( gomatrix, geneset, universe, 
                                           min.genes=5, cut_max = 50 )
        
        return(RobertGO)
      })
      names(GOrobert_res.cut50) <- paste0( "clust_", unique(TRN_geneClust))
      saveRDS(GOrobert_res.cut50, file = "Supl.Fig4/TRNpodo.16clust_GOrobert.cut50.rda")
      
      
      ### top functions per model
      # get a union of top N func in each dataset, select from primary terms only
      GOrobert_res <- GOrobert_res.cut50
      iids <- seq( GOrobert_res)
      thrsh <- 0.1
      
      topN_GOrb  <- Reduce( union , lapply( iids , function( ii,N=5 )
      {
        print(ii)
        datt <- GOrobert_res[[ii]]
        # select from primary terms only
        datt <- subset(datt$results,  Is.primary==TRUE & Primary.Fisher.adj < thrsh )
        datt <- datt[ order(datt$Fisher) ,]
        print(datt$GO.ID[ 1:min(N, nrow(datt)) ])
        # return N top primar terms
        return(datt$GO.ID[ 1:min(N, nrow(datt)) ])
      }))
      topN_GOrb <- topN_GOrb[!is.na(topN_GOrb)]
      
      gg<- GOrobert_barplot_clustered( iid = iids ,
                                       GOtoPlot = topN_GOrb,
                                       datGO = GOrobert_res,
                                       heatmap = T, 
                                       datOrder = names(GOrobert_res)[
                                         as.numeric(sub(".*_","",pp$tree_col$labels[pp$tree_col$order]))  
                                       ])
      
      ### cluster functions
      ss <- reshape(gg$data[,c( "Term","dataset","log10.Fisher")],
                    idvar = "Term", timevar = "dataset", direction = "wide")
      rownames(ss) <- ss$Term
      ss <- ss[,-1]
      pp<-  pheatmap::pheatmap(ss, cluster_rows = F, scale = "column")
      
      
      ## plot
      pdf( width =  16 , height = 14 , file="Supl.Fig4/TRNpodo.16clust_GOrobert.cut50.top5.pdf")
      a<-dev.cur()
      png( width =  1000 , height = 1000 , file="Supl.Fig4/TRNpodo.16clust_GOrobert.cut50.top5.png")
      dev.control("enable")
      
      print(gg )
      
      dev.copy(which=a)
      dev.off()
      dev.off()
      
      
      
      
    }
    
  }
  
  
  ### plot a network view
  {
    library(igraph)
    TF_MeanMed <- readRDS(  file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/TFA/TFs_SNSC.podo.MeanMed.08.01.24.rda")
    
    PodoTFnet <-ATACseq_tgenesM.podo[ (rownames(ATACseq_tgenesM.podo) %in% TF_MeanMed) & 
                                        rownames(ATACseq_tgenesM.podo) %in%  rownames(genecorrPDS.Sprmn_freq)[
                                          genecorrPDS.Sprmn_freq$total>=7], ]
    
    PodoTFnet <- PodoTFnet[ , colSums(PodoTFnet!=0)!=0 ]
    PodoTFnet <- PodoTFnet[ rowSums(PodoTFnet!=0)!=0,]
    
    edge_list <- reshape2::melt(t(PodoTFnet))
    edge_list <- edge_list[edge_list$value!=0,]
    toPlot <- igraph::graph_from_edgelist( as.matrix( edge_list[,1:2] ), directed = T )
    
    # 
    
    plot( toPlot, vertex.size=12, directed=T,
          vertex.color= ifelse( names(V(toPlot)) %in% colnames(PodoTFnet), "orange","skyblue"),
          # layout=layout.circle,
          edge.arrow.size=.2,
          vertex.label.cex=1.2,
          vertex.label.dist=0.1
    ) 
    
    pdf( width =  8 , height = 8 , file="TRN.14_TFs.PDScorCut4.pdf")
    a<-dev.cur()
    png( width =  600 , height = 600 , file="TRN.14_TFs.PDScorCut4.png")
    dev.control("enable")
    print( plot( toPlot, vertex.size=6, directed=T,
                 vertex.color= ifelse( names(V(toPlot)) %in% colnames(PodoTFnet), "orange","skyblue"),
                 # layout=layout.spring,
                 edge.arrow.size=.2,
                 vertex.label.cex=1.0,
                 vertex.label.dist=0.1
    ) )
    dev.copy(which=a)
    dev.off()
    dev.off()
  }
  
  
}


### pathway Fisher test
 {
  PopSize <- length( union( allPodoGenes , Reduce( union, pathDB)))
  

  genecorrPDS_pathAnnot.hyper<- Reduce( cbind,  lapply( seq(genecorrPDS.list) , function(ii)
    {
    
    print(names(gsets[ii]))
    
    # prepare gene set, convert ensembleIDs to entrezIDs
    # geneset <- names(gsets)[cclust==nname]
    geneset <- genecorrPDS.list[[ ii ]]
    list21 <- length(geneset )
    # geneset <- unique(datt$symbol[ datt$fromOverlappingOrNearest == "Overlapping"])
    # test all genes, cis and trans regulated genes
    
    XX <- sapply( seq(pathDB), function(jj){
      list1 <-  length( pathDB[[jj]] )
      overlap <- length( intersect( geneset ,
                                    pathDB[[jj]] ) )

      if(overlap>1) phyper( overlap -1 , list21 , PopSize , list1  , lower.tail= FALSE) else NA
    }) 
    names(XX) <- names(pathDB)

    return(XX)
  } ) )
  
  
  # adjust p-values
  genecorrPDS_pathAnnot.hyper.adj <- apply( genecorrPDS_pathAnnot.hyper, 2, p.adjust , "fdr")
  colnames(genecorrPDS_pathAnnot.hyper.adj) <- colnames(genecorrPDS_pathAnnot.hyper) <- 
    c(  "Damage_signature", "corrPDS", "NOcorrPDS")
  
  
  topN_fisher  <- Reduce( union , sapply( seq(genecorrPDS.list) , function( ii,N=12 )
    {
    print(ii)
    datt <- genecorrPDS_pathAnnot.hyper.adj[,ii]
    # select from primary terms only
    datt <- datt[ datt<0.05 &  !is.na(datt)]
    datt <- genecorrPDS_pathAnnot.hyper.adj[names(datt),ii]
    datt <- datt[order(datt)]
    # return N top primar terms
    return( names( datt[ 1:min(N, length(datt)) ]) )
  }))
  topN_fisher <- topN_fisher[ ! topN_fisher %in% c(
    "mediumpurple1","turquoise3","salmon", "chartreuse3","SARS-CoV-2 Infection__REACT",
    "SARS-CoV Infections__REACT","Axon guidance__KEGG")]
 
  # cluster 
  toPlot <- genecorrPDS_pathAnnot.hyper[ topN_fisher , ]
  # significance labels
  siglbl <- genecorrPDS_pathAnnot.hyper.adj[ rownames(toPlot),]
  siglbl <- ifelse(siglbl <0.05, "*","")
  
  
  siglbl[is.na(siglbl)] <- ""
  # toPlot[ toPlot >=0.05 ] <- NA
  toPlot <-  apply( apply(toPlot, 2, rank,na.last='keep'), 2, rank_normalize)
  rownames(toPlot) <- sapply( sub( "__.*","",rownames(toPlot)) , str_trunc, 50)
  toPlot[is.na(toPlot) ]<- 1.1
  better_col_palette <- viridis::viridis(9,direction = -1)
  
 
  ggl <- pheatmap::pheatmap( toPlot , cluster_cols = F , 
                             display_numbers = siglbl, number_color = "red",
                             fontsize_number = 18 ,
                             color = better_col_palette , fontsize = 16)
  
  
  
  pdf( width =  8 , height = 8 , file="Supl.Fig4/PDScor.6more_pathEnrich.Fisher_barplot.v3.pdf")
  a<-dev.cur()
  png( width =  600 , height = 600 , file="Supl.Fig4/PDScor.6more_pathEnrich.Fisher_barplot.v3.png")
  dev.control("enable")
  print(ggl)
  dev.copy(which=a)
  dev.off()
  dev.off()
  
}

#### 2. combining TRN with PDS to find interesting TFs #### 

### heatmap of intersect of gSets and TFtargets
  {

### load podocyte gene sets
    {
       # PMID:36307401 tab 1 and 2
    slitDiaphr <-  read.table(sep = "\t", header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/cytoscapeVIZ/Podo_areas/PPIsource/PMID:36307401_tab1tab2.slitDiaphr.csv")$GeneName
  
    # Schell et al. PMID:33514561 fig1f
    actCtsklt_fig1f <-read.table(sep = "\t", header = T, fill = T,  "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/cytoscapeVIZ/Podo_areas/C.shell/ActnCtskl_node.csv")
    actCtsklt_fig1f <- unique( actCtsklt_fig1f$gene_symbol)
    # # PMID:28536193 , figure 1E in Schell et al. 2017
    fcladhsn <-  read.table(sep = "\t", header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/cytoscapeVIZ/Podo_areas/PPIsource/PMID:28536193_SupplTab4_focalAdhesion.csv")$Gene.names...primary.
    
    # matrisome PMID:33761352, SupplTab7
    Mtrxsome <- read.table(sep = "\t", fill = T , header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/cytoscapeVIZ/Podo_areas/PPIsource/PMID:33761352_SupplTab7.fltrdMtxsome.csv")$Gene_names
    Mtrxsome <- unlist( strsplit(Mtrxsome, split = ";") )
    Mtrxsome <- unique(unlist( fun_homoTO.FROMmouse(Mtrxsome)$MUS) )
    
    # circadian genes
    circ.genes.bulk0.01 <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/circ.genes.bulk0.01.rda")
    
    # minimal gene set
   GSet <- list( slitDiaphr, actCtsklt_fig1f, 
                                          fcladhsn, Mtrxsome )
    names(GSet) <- c("Slit.Diaphr", "Actin.Ctsklt",
                             "Focal.Adhesion", "Mtrxsome" )
    

    }
    
    
### select genes that correlate with PDS 
    {
      ## decide in how many settings (gtype*study) a gene should correlate with PDS
      ## to be used for the heatmap
      Nn <- 7
      
      genecorrPDS.Sprmn_select <- genecorrPDS.Sprmn_freq[
        genecorrPDS.Sprmn_freq$total>=Nn & 
          rownames(genecorrPDS.Sprmn_freq) %in% allPodoGenes, ]
      
      genecorrPDS <- rownames(genecorrPDS.Sprmn_select)[
        !rownames(genecorrPDS.Sprmn_select)%in% DS_all$gene_symbol[1:42] ]
      
      ### save table 
      genecorrPDS_tab <- data.frame( gene_name=rownames(genecorrPDS.Sprmn_select), 
                                     damage_signature= ifelse( 
                                       rownames(genecorrPDS.Sprmn_select)%in% DS_all$gene_symbol[1:42] , "yes","no"))
      write.table(genecorrPDS_tab, file = "Genes.corrPDS_table.tsv", row.names = F )
      # genecorrPDS.ctrl <- rownames(genecorrPDS.Sprmn_select)[ 
      #   genecorrPDS.Sprmn_select$ctr.only >= genecorrPDS.Sprmn_select$exp.only &
      #     !rownames(genecorrPDS.Sprmn_select)%in% DS_all$gene_symbol[1:42] ]
      # 
      # genecorrPDS.xprmnt <- rownames(genecorrPDS.Sprmn_select)[ 
      #   genecorrPDS.Sprmn_select$ctr.only < genecorrPDS.Sprmn_select$exp.only &
      #     !rownames(genecorrPDS.Sprmn_select)%in% DS_all$gene_symbol[1:42]]
      # 
      # genecorrPDS.both <- rownames(genecorrPDS.Sprmn_select)[ 
      #   genecorrPDS.Sprmn_select$ctr.only == genecorrPDS.Sprmn_select$exp.only &
      #     !rownames(genecorrPDS.Sprmn_select)%in% DS_all$gene_symbol[1:42]]
      
      
      
      # genes that do not correlate with PDS
      genecorrPDS.Sprmn_r.cntrd <-  readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/SCSN_genecorrPDS_Sprmn.r.cntrd.rda")
      
      genecorrPDS.not <- rownames( genecorrPDS.Sprmn_r.cntrd[
        !rownames(genecorrPDS.Sprmn_r.cntrd)%in% DS_all$gene_symbol &
          rownames(genecorrPDS.Sprmn_r.cntrd) %in% allPodoGenes & 
          rownames(genecorrPDS.Sprmn_r.cntrd) %in% 
          rownames(genecorrPDS.Sprmn_freq)[
            genecorrPDS.Sprmn_freq$total==0
          ], ] )
      genecorrPDS.not <- sample(genecorrPDS.not , length(genecorrPDS ))
      
      # 
      genecorrPDS.list <- list ( DS_all$gene_symbol[1:42], genecorrPDS, genecorrPDS.not)
      names(genecorrPDS.list) <- c(  "Damage_signature", "genes.corPDS", "genes.NOcorPDS")
      
      ### combine with other podocyte gsets
      PodoPathGSet_new <- c( GSet ,   
                             circ.genes.bulk0.01=list(circ.genes.bulk0.01 ) ,
                             genecorrPDS.list  )
      
      PodoPathGSet_new <- lapply(PodoPathGSet_new , unique)
      # filter for podocyte expressed genes
      PodoPathGSet_new.podo <- lapply(seq(PodoPathGSet_new), function(ii){
        PodoPathGSet_new[[ii]][ PodoPathGSet_new[[ii]] %in% allPodoGenes]
      })
      names(PodoPathGSet_new.podo) <- names(PodoPathGSet_new)
      saveRDS(PodoPathGSet_new.podo, file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/PodoPathGSet_new.podo.rda")
    }
    
   
### test enrichment of TF targets in gSets 
    {
      PodoPathGSet_new.podo <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/PodoPathGSet_new.podo.rda")
      
      
      PodoPathGSet_PDS.test <- sapply( seq(PodoPathGSet_new.podo) , function( ii) 
      {
        print( names(PodoPathGSet_new)[ii])
        geneSet <- PodoPathGSet_new[[ii]]
        GRN <- ATACseq_tgenesM.podo
        print( paste0( length(geneSet), "_",length(geneSet[ geneSet %in% allPodoGenes])))
        geneSet <- geneSet[ geneSet %in% allPodoGenes]
        
        p.adjust( apply( GRN , 2, function(TF){
          TFtg <- names(TF)[TF>0] ## 
          # print( length(TFtg))
          overlap <- length( intersect( geneSet,  TFtg ) )
          list1 <- length( geneSet )
          list2 <- length( TFtg )
          popSize <- length( allPodoGenes )
          ph <- phyper(overlap-1, list1, popSize - list1, list2, lower.tail=FALSE )
          
          return(ph )
        }), method="fdr" )
        
      } )
      colnames(PodoPathGSet_PDS.test) <-names( PodoPathGSet_new ) 
      
      
      # PodoPathGSet_PDS.test[ PodoPathGSet_PDS.test >  0.05] <- NA
      PodoPathGSet_PDS.test <- PodoPathGSet_PDS.test[ rowSums(!is.na(
        PodoPathGSet_PDS.test))>1,
        colSums(!is.na(
          PodoPathGSet_PDS.test))>1]
      PodoPathGSet_PDS.test <- t(PodoPathGSet_PDS.test)
      siglbl <- round( PodoPathGSet_PDS.test , 4)
      siglbl <- ifelse( siglbl < 0.001, "***",
                        ifelse( siglbl< 0.01,"**",
                                ifelse( siglbl< 0.1,"*","")))
      
    }
   
    
### calculate and plot intersect between gSets and TFtargets
    {
      datt <- ATACseq_tgenesM.podo
      toPlot <- Reduce( rbind, lapply( PodoPathGSet_new , function(gSet){ 
        gSet <- gSet[ gSet%in% allPodoGenes ]
        apply( datt, 2, function(TF ){
          pathFr <- length( TF[ TF>0 & names(TF) %in% gSet]) /length(gSet)
          print( pathFr)
        })
      }))
      rownames(toPlot) <-  names( PodoPathGSet_new ) 
      # toPlot <- toPlot[rowSums(toPlot)>10,]
      # # scale
      # toPlot <-t(scale(t(toPlot),  center = F ))
      # toPlot <-  toPlot[,colSums(toPlot)>5]
      
      # # or rank
      # toPlot <- t( apply( toPlot ,1, rank, ties.method = "min" ) )
      # # col1 <-RColorBrewer::brewer.pal(10,"Paired")
      
      
      
      ### add labels form a hypergeometric test
      
      siglbl2 <- siglbl[ match( rownames(toPlot), rownames(siglbl)),
                         match( colnames(toPlot), colnames(siglbl))]
      
      # annotation for rows and columns
      annotation_col = data.frame(
        row.names = colnames(toPlot) ,
        tgene_N = ( colSums( ATACseq_tgenesM.podo[ ,colnames(toPlot)]>0 )) 
      )
      annotation_row = data.frame(
        row.names = rownames(toPlot) ,
        gSet_size =  ( Reduce( c, lapply( PodoPathGSet_new[ rownames(toPlot) ], function(datt){
          length( datt[ datt %in% allPodoGenes])
        }
        )) )
      )
      # colors for annotations
      ann_colors = list(
        
        gSet_size =  rev(grey.colors(2)),
        tgene_N =  rev( grey.colors(2) ) 
      )
      
      # specify a color pallete
      better_col_palette <- viridis::magma(30)
      colnames(toPlot) <- sub( ".*__", "", colnames(toPlot))
      
      ### make a plot
      gg1 <- ComplexHeatmap::pheatmap( toPlot , 
                                       color =  better_col_palette ,
                                       cluster_cols = T, 
                                       cluster_rows = F,  
                                       number_color = "lightgrey", 
                                       annotation_row = annotation_row ,
                                       annotation_col = annotation_col,
                                       annotation_colors = ann_colors ,
                                       annotation_legend = T,
                                       fontsize_number = 20 ,
                                       fontsize = 14, 
                                       display_numbers = siglbl2 
                                       
      )
      gg1
      
      pdf(height = 5, width = 9, file="podoPaths.TFprior_overlap.test_heatmap.rank.v2.pdf")
      ht <- draw(gg1) 
      dev.off()
      
      png(height = 400, width = 800, file="podoPaths.TFprior_overlap.test_heatmap.rank.v2.png")
      ht <- draw(gg1) 
      dev.off()
      
      ### find overlaps between gene sets and TF regulons 
      GG3 <-  newGOM(c( TFgsets,PodoPathGSet_new.podo),
                     genome.size= length( allPodoGenes) )
      
      drawHeatmap( GG3 , log.scale=T,
                   what="odds.ratio",
                   adj.p = T,
                   cutoff = 0.05)
      
      pdf(height = 14, width = 24, file="GRN.podoGsets_overlap.test_heatmap.odds.ratio.pdf")
      drawHeatmap( GG3 , log.scale=T,
                   what="odds.ratio",
                   adj.p = T,
                   cutoff = 0.05)
      dev.off()
    }
 
   
  
  }

### Suppl. fig motif enrichment of geneSet promoters
  { 
  
  # load gene sets
  PodoPathGSet_new.podo <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/PodoPathGSet_new.podo.rda")
  # load motifs
  motifFile <-  "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/ATACseq/Podocytes/TOBIAS/motif_clusters/motif_comparison_seqcor_Ward.cut.1_consensus_motifs.meme" 
  mmotifs <- universalmotif::read_meme( motifFile )
  
  motifFile_ind <- "/media/tim_nevelsk/WD_tim/ANNOTATIONS/CISBP/CISBP2mouse_podoTF_08.01.24.meme"
  mmotifs_ind <- universalmotif::read_meme( motifFile_ind )
  
  
   # mmotifs <- universalmotif::read_meme( "/media/tim_nevelsk/WD_tim/ANNOTATIONS/CISBP/CISBP2mouse_podoTF_08.01.24.meme" )
    names(mmotifs) <- sapply(seq(mmotifs) , function(ii) mmotifs[[ii]]@name )
    names(mmotifs_ind) <- sapply(seq(mmotifs_ind) , function(ii) mmotifs_ind[[ii]]@name )
    
  # mmotifs <- mmotifs[colnames(ATACseq_tgenesM.podo)]
  
  ### test enrichment
  library(memes)
  
  PodoPathGSet_new.podo_ame.15TF.runksum <- lapply( seq(PodoPathGSet_new.podo), 
                                            function(ii){
    
    ggenes <- PodoPathGSet_new.podo[[ii]]
    
    TF_motifEnrichTest( gset =   ggenes ,
                        bckgrGenes =  "all" ,
                        scoreMethod = "ranksum", 
                        motifDir = motifFile ,
                        outDir = paste0(
                          "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/AME_motifAnalysis/PodoPathGSet/15meta.podoMotifs.allBckgrnd.ranksum_PodoGsets/",
                          names(PodoPathGSet_new.podo)[ii]))
    
    
  } )
  
  PodoPathGSet_new.podo_ame.15TF <- lapply( seq(PodoPathGSet_new.podo), 
                                            function(ii){
                                             
                                             ggenes <- PodoPathGSet_new.podo[[ii]]
                                             
                                             TF_motifEnrichTest( gset =   ggenes ,
                                                                 bckgrGenes =  "all" ,
                                                                 scoreMethod = "fisher", 
                                                                 motifDir = motifFile ,
                                                                 outDir = paste0(
                                                                   "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/AME_motifAnalysis/PodoPathGSet/15meta.podoMotifs.allBckgrnd_PodoGsets/",
                                                                   names(PodoPathGSet_new.podo)[ii]))
                                             
                                             
                                           } )
  
  # PodoPathGSet_ame <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/AME_motifAnalysis/PodoPathGSet/15meta.podoMotifs.allBckgrnd.ranksum/PodoPathGSet_AMEmotifEnrich.rda")
  # PodoPathGSet_ame <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/AME_motifAnalysis/PodoPathGSet/15meta.podoMotifs.podoBckgrnd.ranksum/PodoPathGSet_AMEmotifEnrich.rda")
  PodoPathGSet_new.podo_ame <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/AME_motifAnalysis/PodoPathGSet/110podoMotifs.allBckgrnd.ranksum/PodoPathGSet_AMEmotifEnrich.rda")
  
  
  names( PodoPathGSet_new.podo_ame ) <- names( PodoPathGSet_new.podo)
  
  PodoPathGSet_ame.padj <- Reduce( cbind, lapply( seq(PodoPathGSet_ame), function(ii){
    print(ii)
    datt <- PodoPathGSet_ame[[ii]]
    datt.vec <- setNames( datt$adj.pvalue , datt$motif_id)
    datt.vec <- datt.vec[ match( names(mmotifs), names(datt.vec))]
    names(datt.vec) <- names(mmotifs)
    datt.vec[ is.na(datt.vec)]<- max(datt.vec,na.rm = T)
    return(datt.vec)
  }))
  colnames( PodoPathGSet_ame.padj) <- names(PodoPathGSet)
  PodoPathGSet_ame.padj <- t(PodoPathGSet_ame.padj)
  #
  siglbl3 <- ifelse( PodoPathGSet_ame.padj<0.001, "***",
                     ifelse(PodoPathGSet_ame.padj<0.01,"**",
                            ifelse(PodoPathGSet_ame.padj<0.1, "*","")))
  
  better_col_palette <- viridis::mako(30)
  dim(PodoPathGSet_ame.padj)
  gg2 <- ComplexHeatmap::pheatmap( PodoPathGSet_ame.padj , 
                                   row_order= row_order(ht),
                                   column_order= column_order(ht),
                                   color =  rev(better_col_palette)  ,
                                   cluster_cols = F, 
                                   cluster_rows = F,  
                                   number_color = "red", 
                                   # annotation_row = annotation_row ,
                                   # annotation_col = annotation_col,
                                   # annotation_colors = ann_colors ,
                                   annotation_legend = T,
                                   fontsize_number = 18 ,
                                   fontsize = 14,
                                   display_numbers = siglbl3   )
  
  pdf(height = 6, width = 15, file="Supl.Fig4/podoPaths.TFprior.AND.ameMotifEnrich_heatmaps.pdf")
  draw(gg1+gg2)
  dev.off()
  
}



### Suppl. fig TF PWM simmilarity tree
### made by TOBIAS

### Suppl. fig TRN related to Wt1
toPlot <- ATACseq_tgenesM.podo[
  rownames(ATACseq_tgenesM.podo) %in%
    grep("col",  rownames(ATACseq_tgenesM.podo), value = T, ignore.case = T), ]
toPlot <- reshape2::melt(t(toPlot))
net <- igraph::graph_from_data_frame(d= toPlot , directed=T)
plot(net, head.arrow.size=2)

#### 3. show correlation between PDS and circadian (dis)regulation #### 

### Suppl.Fig. plot basics of circadian bulk data
  {
    library( biomaRt)
    # mart_homo <- useMart( "ensembl",dataset="hsapiens_gene_ensembl" , host="www.ensembl.org")
    mart_mouse <- useMart( "ensembl",dataset="mmusculus_gene_ensembl" )
    tx2gene <- getBM( attributes=c( 'ensembl_gene_id', 'external_gene_name',"entrezgene_id"),  mart = mart_mouse)
    
    # load count matrix
   countMat <- read.table( row.names = 1, header = T, "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/counts.tsv")
  countMat <- countMat[,-1]
  # remove one sample with very few counts
    countMat <- countMat[, colnames(countMat)!="SN8580185_21056_806Aligned.sortedByCoord.out.bam"]
    # load normilised by DEseq2 counts
    dds_rlog.norm <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/bulkRNAseq.Kidney.cirk_rlog.rda")
    
    # load annotation
  annot <- gsub(".*_|Aligned.*", "", colnames(countMat))
  annot <-data.frame(  ID= colnames(countMat) ,
                       condition= annot,
                       age=c(rep(10,16), rep(80,15)),
                       time= sub("^..","",annot))
  annot$age <- as.factor( annot$age )
  annot$time <- as.numeric( annot$time )
  
  
  ### load results
  MetaCycle_young <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/MetaCycle_bulk.kidney_young.rda")
  MetaCycle_old <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/MetaCycle_bulk.kidney_old.rda")
  
  datt.y <- MetaCycle_young$meta
  datt.y$gName <- tx2gene$external_gene_name[ match( datt.y$CycID, tx2gene$ensembl_gene_id)]
  datt.y.podo <- datt.y[ datt.y$gName %in% allPodoGenes , ]
  
  datt.o <- MetaCycle_old$meta
  datt.o$gName <- tx2gene$external_gene_name[ match( datt.o$CycID, tx2gene$ensembl_gene_id)]
  datt.o.podo <- datt.o[ datt.o$gName %in% allPodoGenes , ]
  
  
  ### do PCA
  library(factoextra)
  library(viridis)
  
  
  PCAdat_sel <- dds_rlog.norm
  PCAdat_sel <- PCAdat_sel[ rownames(PCAdat_sel) %in% union(
    datt.o.podo$CycID[datt.o.podo$JTK_BH.Q<0.05] ,  datt.y.podo$CycID[datt.y.podo$JTK_BH.Q<0.05]
  ),]
  colnames(PCAdat_sel) <- paste0( annot$age ,"_", annot$time,  "hrs")
  ## run PCA  
  res.pca <- prcomp( t(PCAdat_sel), scale = F)
  
  # plot results
  fviz_eig(res.pca)
  biplot(res.pca, choices = )
  ggl <-  factoextra::fviz_pca_ind(res.pca, 
                                  habillage = ifelse( annot$age==10, "young","old" ))+ 
    scale_color_colorblind()+ theme( text =  element_text(size=20))
  
  pdf( width = 6 , height = 6 , file="Supl.Fig4/bulkRNAseq_circadian_PCAplot.pdf")
  a<-dev.cur()
  png( width =  400 , height = 400 , file="Supl.Fig4/bulkRNAseq_circadian_PCAplot.png")
  dev.control("enable")
  print(ggl)
  dev.copy(which=a)
  dev.off()
  dev.off()
  
  ### 
  toPlot <- dds_rlog.norm[ rownames(dds_rlog.norm)%in% 
                             tx2gene$ensembl_gene_id[tx2gene$external_gene_name%in% 
                                                       c("Bmal1","Nr1d1","Cry1","Per2")],]
  toPlot.melt <- reshape2::melt(toPlot)
  toPlot.melt$Time <- annot$time[ match( toPlot.melt$Var2, annot$ID)]
  toPlot.melt$Gene <- tx2gene$external_gene_name[ match( toPlot.melt$Var1, tx2gene$ensembl_gene_id)]
  toPlot.melt$age <- annot$age[ match( toPlot.melt$Var2, annot$ID)]
  
     # Create the plot
    p <- ggplot(toPlot.melt, aes(x = Time, y = value, color=age)) +
      geom_point( size = 3) +
      geom_smooth( method="loess",span=0.4,aes(fill=age)) +
      labs( x = "Time (hours)",
           y = "Expression Level") + 
      scale_color_colorblind()  + scale_fill_colorblind()  +
      theme_minimal()+ facet_wrap("Gene",scales = "free")+ theme(
        text=element_text( size=20)
      )
    
    pdf( width = 10 , height = 6 , file="Supl.Fig4/bulkRNAseq_circadian_4CoreClockExpr.pdf")
    a<-dev.cur()
    png( width =  800 , height = 400 , file="Supl.Fig4/bulkRNAseq_circadian_4CoreClockExpr.png")
    dev.control("enable")
    print(p)
    dev.copy(which=a)
    dev.off()
    dev.off()
    

  }


### plot summary stat for CRD https://github.com/leihe2021/CRDscore/ 
  {
  
  listSCSN.CRD_table <- readRDS( file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/listSCSN.CRD_table.rda")
  # circ.genes <- circ.genes[ circ.genes%in% allPodoGenes]
  
  listSCSN.CRD_table$sampleDiscr <-  paste( listSCSN.CRD_table$sample , 
                                            listSCSN.CRD_table$group, 
                                            listSCSN.CRD_table$dataSet, sep = "..")
  ssamples <- unique( listSCSN.CRD_table$sampleDiscr)
  
  toPlot <- Reduce( rbind.data.frame , lapply( ssamples, function(ssample){
    
    datt <- listSCSN.CRD_table[ listSCSN.CRD_table$sampleDiscr==ssample, ]
    datt.meta <- unique( datt[  !colnames( datt ) %in%  c( "PDS", "CRD.bulk0.01" ) ] )
    # print(datt.meta)
    
    rho.test <-  tryCatch( cor.test(datt$PDS , 
                                    datt$CRD.bulk0.01, method = "spearman") ,
                           error = function(e) NA)
    
    if( !is.na(rho.test)){
      rho <- rho.test$estimate
      rho.pval <-  rho.test$p.value
    } else rho <- rho.pval <- NA
    
    ccv <- setNames( c(datt.meta,  rho , rho.pval ) ,
                     c( names(datt.meta) , "rho_PDSvsCRD","rho.pval") )
    
    return(ccv)
  }) )
  
  toPlot$sigStat <- ifelse( toPlot$rho.pval < 0.05 , "sig", "non.sig")
  toPlot$sigStat[ is.na(toPlot$sigStat)] <- "non.sig"
  
  toPlot$seq.type <- ifelse( toPlot$dataSet %in% c(
    "Nphs2","Wt1","Pdss2","Lmx1b"), "sn", "sc")
  toPlot$sampleDiscr.sig <- ifelse( toPlot$rho.pval < 0.05 , 
                                    toPlot$sampleDiscr, "")
  
  
  ## plot
  gg1 <- ggplot( data = toPlot, aes( x= gtypeDE, y= rho_PDSvsCRD ,  
                                    shape= seq.type, color=sigStat ))  + 
    # geom_label_repel( aes( label = sampleDiscr.sig ) , 
    #                   color="grey30", size=3,label.size = NA )+
    geom_jitter( width = 0.1, alpha=0.5, size=5)+ theme_bw( )+ 
    ggtitle( "all samples")+
    scale_color_tableau() + theme( text = element_text( size = 20))
  # facet_grid( rows = vars( testRun ) )
  
  
 
}

### Suppl.fig plot CDR results for all individual samples
  {
  
  toPlot <- listSCSN.CRD_table[ listSCSN.CRD_table$sample=="139913",]
  gg2 <- ggplot( data= toPlot,
               aes( y = toPlot[["CRD.bulk0.01"]] , 
                    x = PDS) ) +
    geom_point( alpha=0.3) + ggtitle("Nphs2_139913")+
    # geom_smooth(method = lm) + 
    geom_smooth(method = lm,  se = F, color="black") +
    theme_bw() + theme( text = element_text(size=20), legend.position = "bottom") +
    coord_cartesian(ylim = c( quantile( toPlot[[ "CRD.bulk0.01" ]], 0.025),
                              quantile( toPlot[[ "CRD.bulk0.01" ]], 1)))+
    stat_cor( r.accuracy = 0.01,
              method = "spearman", size=6)+
    scale_color_colorblind()+ labs(y = "CRD") 
  
  pdf(width = 9 , height = 6, file = "CRD_listSCSN_indSampANDsummary.pdf")
  plot_grid(plotlist = list( gg2, gg1))
  
  dev.off()
  
}

### plot summary for Tempo results with Clock as reference  gene
  {    
    tempo.out.list <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/circadian/tempo/tempo.out_probabilities_listSCSN.rda")
    names(tempo.out.list) <- c("Arntl","Clock")

    ### extract max. probability and circ.time from tempo reults
    ## add scsnRNAseq metadata, 
    tempo.out.maxP.list <-lapply( seq(tempo.out.list), function(dd)
      {
      tempo.out <- tempo.out.list[[dd]]
      
      ## extract max. probability and circ. time associated with it from tempo reults
      tempo.out.maxP <- Reduce( rbind , lapply( seq(tempo.out), function(ii)
      {
        
        sampleName <- sub("__.*" , "",  names(tempo.out)[ii])
        dataSet <-  sub(".*__" , "",  names(tempo.out)[ii])
        XX <- data.frame( tempo.maxP= apply( tempo.out[[ii]], 1, max),
                          tempo.maxPtime= apply( tempo.out[[ii]], 1, which.max ),
                          sample = sampleName ,
                          dataSet=  dataSet )
        # XX$tempo.maxPtime[ XX$tempo.maxP< median(XX$tempo.maxP) ] <- NA      
        return(XX)
      }) )
      
      ## combine tempo results with PDS and sn.scRNAseq metadata
      tempo.out.maxP$PDS <- listSCSN.CRD_table[ rownames(tempo.out.maxP),"PDS" ]
      tempo.out.maxP$gtypeDE <- listSCSN.CRD_table[ rownames(tempo.out.maxP),"gtypeDE" ]
      tempo.out.maxP$group <- listSCSN.CRD_table[ rownames(tempo.out.maxP),"group" ]
      tempo.out.maxP$sampleDiscr <- paste0( tempo.out.maxP$sample , "..", tempo.out.maxP$group)
      tempo.out.maxP$testRun <- names(tempo.out.list)[dd]
      
      return(tempo.out.maxP)
    })  
    
    ## calculate likely circadian time of a sample and deviation from it of each cell 
    tempo.out.maxP.diff.list <- lapply( seq(tempo.out.maxP.list), function(dd)
    {
      tempo.out.maxP <- tempo.out.maxP.list[[dd]]
      
      Reduce( rbind, lapply( unique(tempo.out.maxP$sample), function(ss)
      {
         # print(ss)
        # select one sample
        datt <- tempo.out.maxP[tempo.out.maxP$sample==ss,]
        
        ## find most probable sanple time
        # 
        if( table(is.na( datt$tempo.maxPtime))["FALSE"]  <2) {
          datt$timeMost.diff <- NA
        } else {
          dd <- density( datt$tempo.maxPtime , window = "rectangular",na.rm = T)
          
          X<-which.max( dd$y)
          timeMost <- dd$x[X]
          # timeMost <- unique( datt$sample.hrmnTime)
          
          datt$timeMost.diff <-  abs( unlist(
            lapply( datt$tempo.maxPtime, circadian_time_difference , ct1=timeMost)) )
          
        }
       
        
        # peaks <- pracma::findpeaks(dd$y)
        # datt$tempo.maxPtime_scndMax <- max(peaks[,1]) - sort(peaks[,1],partial=1)[1]
        
        return(datt)
      }) )
      
    })
      
    ### correlate circadian time deviation with PDS in each sample
    tempo.out.maxP.diff_vsPDS.rho <-  Reduce( 
        rbind, lapply( seq(tempo.out.maxP.diff.list), function(dd,
                                                               corTO="timeMost.diff")
          {
          
          tempo.out.maxP.diff <- tempo.out.maxP.diff.list[[dd]]
          
      ## calculate most likely circ.time of a sample 
      ## and for each cell calculate deviation from the likely time

      
      # calculate per-sample correlation with PDS
      ssamples <- unique( tempo.out.maxP.diff$sampleDiscr)
      tempo.out.maxP.diff_vsPDS.rho <- Reduce( rbind.data.frame  , lapply( ssamples, function(ssample)
        {
        print( ssample )
        
        datt <- tempo.out.maxP.diff[ tempo.out.maxP.diff$sampleDiscr==ssample, ]
        datt.meta <- unique( datt[  !colnames( datt ) %in%  c( "PDS", "timeMost.diff" ,
                                                               "tempo.maxP" ,"tempo.maxPtime" ) ] )
        # print(datt.meta)
        
        rho.test <-  tryCatch( cor.test(datt$PDS , 
                                        datt[[corTO]], 
                                        method = "spearman" ) ,
                               error = function(e) NA)
     
        if( length(rho.test)>1 ){
          rho <- rho.test$estimate
          rho.pval <-  rho.test$p.value
        } else rho <- rho.pval <- NA
        
        
        ccv <- setNames( c(datt.meta,  rho , rho.pval ) ,
                         c( names(datt.meta) , "rho_PDSvsTEMPO.diff","rho.pval") )
        
        return(ccv)
      }) )
      
      tempo.out.maxP.diff_vsPDS.rho$sigStat <- ifelse( tempo.out.maxP.diff_vsPDS.rho$rho.pval < 0.05 , "sig", "non.sig")
      tempo.out.maxP.diff_vsPDS.rho$sigStat[ is.na(tempo.out.maxP.diff_vsPDS.rho$sigStat)] <- "non.sig"
      
      tempo.out.maxP.diff_vsPDS.rho$seq.type <- ifelse( tempo.out.maxP.diff_vsPDS.rho$dataSet %in% c(
        "Nphs2","Wt1","Pdss2","Lmx1b"), "sn", "sc")
      tempo.out.maxP.diff_vsPDS.rho$sampleDiscr.sig <- ifelse( tempo.out.maxP.diff_vsPDS.rho$rho.pval < 0.05 , 
                                                               tempo.out.maxP.diff_vsPDS.rho$sampleDiscr, "")
    
        return(tempo.out.maxP.diff_vsPDS.rho)
      }) )

    ### plot rho summary for all samples
    {
      toPlot <- tempo.out.maxP.diff_vsPDS.rho[tempo.out.maxP.diff_vsPDS.rho$testRun=="Clock",]
      gg2 <- ggplot( data = toPlot, aes( x= gtypeDE, y= rho_PDSvsTEMPO.diff ,  
                                         shape= seq.type, color=sigStat ))  + 
        ggtitle( "all samples")+
        # geom_label_repel( aes( label = sampleDiscr.sig ) , 
        #                   color="grey30", size=3,label.size = NA )+
        geom_jitter( width = 0.1, alpha=0.5, size=5)+ theme_bw( )+ 
        # facet_wrap( vars( testRun ) )+
        scale_color_tableau() + theme( text = element_text( size = 20))
    }
   
   
    ### plot individual sample scatterplot
    {
      toPlot <- tempo.out.maxP.diff.list[[2]]
      toPlot <- toPlot[ toPlot$sampleDiscr=="CD2AP_KO_2..CD2AP_KO__wk3",]
      
      gg2<- ggplot( data= toPlot,
                    aes( y = toPlot[["timeMost.diff"]] , 
                         x = PDS) ) +
        geom_point( alpha=0.3, aes( color=gtypeDE)) + 
        geom_smooth(method = lm,  se = F) + ggtitle("CD2AP_KO_2")+
        theme_bw() + theme( text = element_text(size=20), legend.position = "none") +
        stat_cor( r.accuracy = 0.01,
                  method = "spearman", size=6)+
        scale_color_colorblind()+ labs(y = "timeMost.diff") 
    }
  
    
    pdf(width = 9 , height = 6, file = "tempo.out_core.Clock_indSampANDsummary.pdf")
      plot_grid(plotlist = list( gg2, gg1))
    dev.off()
  }

### Suppl.fig. summary for Tempo results, Clock and Arntl as references
  {
  toPlot <- tempo.out.maxP.diff_vsPDS.rho
  gg <- ggplot( data = toPlot, aes( x= gtypeDE, y= rho_PDSvsTEMPO.diff ,  
                                     shape= seq.type, color=sigStat ))  + 
    xlab("")+
    geom_label_repel( aes( label = sampleDiscr.sig ) ,
                      color="grey30", size=3,label.size = NA )+
    geom_jitter( width = 0.1, alpha=0.5, size=5)+ theme_bw( )+ 
    facet_wrap( vars( testRun ) )+
    scale_color_tableau() + theme( text = element_text( size = 20))
  
  pdf(width = 8 , height = 6, file = "Supl.Fig4/tempo.out_listSCSN.counts_core.Clock.Arntl.pdf")
    plot_grid(plotlist = list( gg ))
  dev.off()
  
}

### Suppl.fig. scatterplot for all individual samples
  {
  # toPlot <- tempo.out.maxP.diff[ tempo.out.maxP.diff$sample=="Nephritis_d1_2",]
  toPlot <- tempo.out.maxP.diff.list[[2]]

  gg<- ggplot( data= toPlot,
               aes( y = toPlot[["timeMost.diff"]] , 
                    x = PDS) ) +
    geom_point( alpha=0.3, aes( color=gtypeDE)) + 
    geom_smooth(method = lm,  se = F) +
    theme_bw() + theme( text = element_text(size=20), legend.position = "bottom") +
    stat_cor( r.accuracy = 0.01,
              method = "spearman", size=6)+
    facet_wrap(vars(sample ), scales = "free_x") +
    scale_color_colorblind()+ labs(y = "timeMost.diff") 
   
  pdf(width = 20 , height = 20, file = "Supl.Fig4/tempo.out_core.Clock_indSmpls.scatter2.pdf")
    gg
  dev.off()
  
  
}


# #### 4. show deviation of circ.gene correlation in high vs low damage #### 
# listSCSN.1K.sampl <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN_samples.1K.22.12.23.rda")
# listSCSN <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN.PDS_22.12.23.rda")
# 
# circ.genes.bulk0.01 <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/circ.genes.bulk0.01.rda")
# coreClock_genes <- read.table( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/circadian/tempo/core_clock_genes.txt" )$V1 
# 
# SeuDat <- listSCSN
# circGenCor_lowVShighPDS <- Reduce( rbind, lapply(seq(SeuDat), 
#                                   function(ii){
#   print(ii)
# 
#   perStudy <- Reduce( rbind.data.frame, 
#                       lapply( unique( SeuDat[[ii]]$sample), 
#                               function(iid){
#     print(iid)
#     datt <- subset( SeuDat[[ii]] , subset= sample== iid)
#     
#     expr <- datt@assays$RNA@data
#     pds <- datt$PDS[ order( datt$PDS)]
#     
#     if( ncol(datt)>20){
#       lowD <- expr[ , colnames(expr)%in% names(pds)[1:(ncol(expr)/2)] ]
#       highD <-expr[ ,  !colnames(expr) %in% colnames(lowD)]
#       # lowD <- expr[ rownames(expr)%in% allPodoGenes, 
#       #               colnames(expr)%in% names(pds)[1:(ncol(expr)/2)] ]
#       # lowD.r <- apply(lowD, 2, rank)
#       # highD <-expr[  rownames(expr)%in% allPodoGenes,
#       #                !colnames(expr) %in% colnames(lowD)]
#       # highD.r <- apply(highD, 2, rank)
#       # correlation between cycling genes, as defined from bulk
#       # corr_lowD <- mean( abs( cor( t( as.matrix(lowD[  
#       #   rownames(lowD) %in%circ.genes.bulk0.01 , ])), method = "spearman") ), na.rm = T)
#       # corr_highD <- mean( abs( cor( t( as.matrix(highD[
#       #   rownames(highD) %in%circ.genes.bulk0.01 , ])), method = "spearman") ), na.rm = T)
#       # 
#       corr_lowD <- colSums( as.matrix(lowD[  
#         rownames(lowD) %in%  circ.genes.bulk0.01 , ]), na.rm = T)
#       corr_highD <- colSums( as.matrix(highD[
#         rownames(highD) %in% circ.genes.bulk0.01 , ]), na.rm = T)
#       
#       # correlation between core clock genes
#       # coreClock_lowD <- mean( abs( cor( t( as.matrix(lowD[
#       #   rownames(lowD) %in% coreClock_genes , ])), method = "spearman") ), na.rm = T)
#       # coreClock_highD <- mean( abs( cor( t( as.matrix(highD[
#       #   rownames(lowD) %in% coreClock_genes , ])), method = "spearman") ), na.rm = T)
#       coreClock_lowD <- colSums( as.matrix(lowD[
#         rownames(lowD) %in%  coreClock_genes , ]) , na.rm = T)
#       coreClock_highD <- colSums( as.matrix(highD[
#         rownames(highD) %in% coreClock_genes , ]) , na.rm = T)
#       
#       genotypeDE <- unique(datt$gtypeDE)
#       corv <-  c( genotypeDE , cor( corr_lowD,coreClock_lowD,method = "spearman" ),
#                   cor( corr_highD,coreClock_highD,method = "spearman" ) ) 
#       
#       # lowD.circvar <- mean(matrixStats::rowVars( t( lowD.r[ circ.genes.bulk0.01, ])))
#       # highD.circvar <- mean(matrixStats::rowVars( t( highD.r[ circ.genes.bulk0.01, ])))
#       # corv <- setNames( c(lowD.circvar, highD.circvar), nm = c("corr_lowD","corr_highD")  )
#       
#     } else{
#       genotypeDE <- unique(datt$gtypeDE)
#       
#       # corv <-  c(genotypeDE, NA, NA, NA, NA) 
#       corv <-  c(genotypeDE, NA, NA) 
#       
#     }
#     
#   
#     return(corv)
#   }))
#   
#   # colnames(perStudy) <- c("genotypeDE" , "corr_lowD","corr_highD",
#   #                 "coreClock_lowD", "coreClock_highD")
#   colnames(perStudy) <- c("genotypeDE" , "corr_lowD","corr_highD")
#   
#   perStudy$sample <- unique( SeuDat[[ii]]$sample)
#   perStudy$study <- names(SeuDat)[ii]
#   
#   return(perStudy)
# }))
# 
# 
# 
# 
# 
# 
# 
# datt <- reshape2::melt(circGenCor_lowVShighPDS, 
#                        id.vars = c("genotypeDE","sample","study"))
# datt <- datt[datt$value!=1 & !is.na(datt$value),]
# datt$value <- as.numeric(datt$value)
# 
# ggplot( datt , aes(x=study, y=value, color=variable ))+
#   geom_boxplot(outlier.alpha = NULL)+ geom_jitter(width=0.2)+
#   theme_bw()+ scale_color_colorblind()+
#   facet_grid(cols = vars( genotypeDE), scales = "free")
#   
# 
