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
inputdir <-  "/media/tim_nevelsk/WD_tim/PROJECTS/WRITING/PDS_manuscript/Figure_input" 

# load necessary code
source("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/cell-damage-score/AUCell_script.r")
source("/home/tim_nevelsk/PROJECTS/myCode/usefulRfunc.r")

# load damage signatures
DS_all <- read.table( header = T,  file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")
# load 
allPodoGenes <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/SCSN_allPodoGenes.rda")

#### 1. analyse genes that correlate with PDS #### 

### count heatmap, shows N of studies where a gene (row) correlate with PDS in 
### control or experimental setting (columns)
  {
  genecorrPDS.Sprmn_freq <- readRDS(  file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/genecorrPDS.Sprmn_freq.16.05.24.rda")
  
  
  ## decide in how many settings (gtype*study) a gene should correlate with PDS
  ## to be used for the heatmap
  Nn <- 7
  
  toPlot <- genecorrPDS.Sprmn_freq[genecorrPDS.Sprmn_freq$total>=Nn, 1:4]
  cclust <- as.factor( cutree( hclust( dist(as.matrix(toPlot), method = "euclidean"),
                                       method="complete"), k = 4) )
  levels(cclust) <- hue_pal()(4)
  
  # DS42 membership annoatation
  pds42 <-  as.factor( rownames(toPlot)%in%DS_all$gene_symbol[1:42])
  # levels(pds42) <- c("darkgrey" ,"black")
  
  ### annotation
  ha = HeatmapAnnotation(
    df=data.frame( PDS = pds42,
                   clusters = cclust ),
    col = list( 
      PDS = c( "TRUE"="black" , "FALSE"= "darkgrey"),
      clusters =  setNames( cclust , cclust)
    ),
    simple_anno_size = unit(1, "cm")
  )
  
  pdf(height = 4, width = 8, file = "PDScorr_count.heatmap.16.05.24.pdf")
  
  ComplexHeatmap::Heatmap( 
    t(as.matrix(toPlot)),
    column_labels =rep("", nrow(toPlot))  ,
    cluster_rows  = F , 
    # margins= c(1, 8),
    col = viridis(max(toPlot)+1,option="magma"),   #labCol = FALSE,
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "complete" , 
    column_split=cclust ,
    top_annotation= ha, 
    name = "count"
  )
  dev.off()
}


### functional annotation of the heatmap clusters
### Robert GO
{
  ### perform enrichment
  library(DBI)
    source("/home/tim_nevelsk/PROJECTS/myCode/func_analysis.R")  
  

  # perform enrichment
  PDScor.clust_RobertGO.5 <- lapply( unique(cclust),
                                      function( nname ){ 
                                        
                                        # convert to entrez
                                        universe <- unique( entr2gName$entrezgene_id)
                                        universe <- as.character( universe[!is.na(universe)] )
                                        
                                        # prepare gene set, convert ensembleIDs to entrezIDs
                                        gset <- names(cclust)[cclust==nname]
                                        
                                        geneset <- unique(entr2gName$entrezgene_id[match( gset,  entr2gName$external_gene_name)])
                                        geneset <- geneset[!is.na(geneset)]
                                        geneset <- as.character(geneset)
                                        
                                        print(length(intersect(geneset,colnames(gomatrix))))
                                        
                                        # apply Robert's function that given a sparse matrix of GO terms (columns = genes, rows = GO terms)
                                        # a geneset of interest and a background set of genes (universe)
                                        # carry out clustering with members diverging by at most cut_max genes, and do enrichment testing.
                                        # Note, multiplicity adjustment is performed for the representative terms only.
                                        # RobertGO.10 <- sf.clusterGoByGeneset( gomatrix, geneset, universe, 
                                        #                                       min.genes=5, cut_max = 10 )
                                        RobertGO.50 <- tryCatch( sf.clusterGoByGeneset( gomatrix, geneset, universe, 
                                                                                        min.genes= 3 , cut_max = 5  ), 
                                                                 error = function(e) NA )
                                        
                                        # RobertGO.10$results$Term <- goterms[ match( 
                                        #   RobertGO.10$results$GO.ID, names(goterms))]
                                        if( !is.na(RobertGO.50)) RobertGO.50$results$Term <- goterms[ match( 
                                          RobertGO.50$results$GO.ID, names(goterms))]
                                        
                                        return(RobertGO.50)
                                        
                                      })
   names(PDScor.clust_RobertGO.5)<- 
    c("damage signature", "mixed big","mixed small","exprmnt only")
 
   ### top functions per model
  # get a union of top N func in each dataset, select from primary terms only
  thrsh <- 0.1
  goRes <- PDScor.clust_RobertGO.5
  
  topN_GOrb  <- Reduce( union , lapply( seq(goRes) , function( ii,N=5 )
    {
    print(ii)
    datt <- goRes[[ii]]
    # select from primary terms only
    datt <- subset(datt$results,  Primary.Fisher.adj < thrsh )
    datt <- datt[ order(datt$Fisher) ,]
    # return N top primar terms
    return(datt$GO.ID[ 1:min(N, nrow(datt)) ])
  }))
  
  
  # #  order datasets
  # oorder.snNphs2<- c("DE4.Wt1.snRNAseq" , "DE12.Wt1.snRNAseq", "DE25.Wt1.snRNAseq")
  # oorder.snWt1<- c("DE4.Wt1.snRNAseq" , "DE12.Wt1.snRNAseq", "DE25.Wt1.snRNAseq")
  oorder.snAll <- c("damage signature", "mixed big","mixed small","exprmnt only")
  
  gg <- GOrobert_barplot_clustered( iid = seq(goRes) ,
                                    GOtoPlot = topN_GOrb,
                                    datGO = goRes ,
                                    # datLFC= GO_meanLFC
                                     datOrder=oorder.snAll
  )
  
  
  pdf( width =  14 , height = 6 , file="PDScor_robertGO.top8_barplot.pdf")
  a<-dev.cur()
  png( width =  1000 , height = 500 , file="PDScor_robertGO.top8_barplot.png")
  dev.control("enable")
  print(gg)
  dev.copy(which=a)
  dev.off()
  dev.off()
}

### Cluster Profiler
{
  PDScor.clust_GO <- lapply( unique(cclust),
                                     function( nname ){ 
                                       
                                       # convert to entrez
                                       universe <- unique( entr2gName$entrezgene_id)
                                       universe <- as.character( universe[!is.na(universe)] )
                                       
                                       # prepare gene set, convert ensembleIDs to entrezIDs
                                       gset <- names(cclust)[cclust==nname]
                                       
                                       geneset <- unique(entr2gName$entrezgene_id[match( gset,  entr2gName$external_gene_name)])
                                       geneset <- geneset[!is.na(geneset)]
                                       geneset <- as.character(geneset)
                                       
                                       print(length(intersect(geneset,colnames(gomatrix))))
                                       
                                       # apply Robert's function that given a sparse matrix of GO terms (columns = genes, rows = GO terms)
                                       # a geneset of interest and a background set of genes (universe)
                                       # carry out clustering with members diverging by at most cut_max genes, and do enrichment testing.
                                       # Note, multiplicity adjustment is performed for the representative terms only.
                                       # RobertGO.10 <- sf.clusterGoByGeneset( gomatrix, geneset, universe, 
                                       #                                       min.genes=5, cut_max = 10 )
                                       cl.annot <- clusterProfiler::enrichGO( gene = geneset,
                                                                              OrgDb="org.Mm.eg.db",
                                                                                # organism = "mmu",
                                                                                minGSSize = 3)
                                                                                  
                                                                               #                                                                     organism     = 'mmu',
                                                                               #                                                                     pvalueCutoff = 0.05 )@result
                                       
                                       # RobertGO.10$results$Term <- goterms[ match( 
                                       #   RobertGO.10$results$GO.ID, names(goterms))]
                                       # if( !is.na(cl.annot)) RobertGO.50$results$Term <- goterms[ match( 
                                       #   RobertGO.50$results$GO.ID, names(goterms))]
                                       
                                       return(cl.annot)
                                       
                                     })
  names( PDScor.clust_GO )<- 
    c("damage signature", "mixed big","mixed small","exprmnt only")
  
 ggl <- lapply( seq(PDScor.clust_GO), function(ii){
   datt <- PDScor.clust_GO[[ii]]
   gg <- clusterProfiler::dotplot( datt , showCategory=15  ) + 
     ggtitle( names( PDScor.clust_GO )[ii])
   return(gg)
 })
 cowplot::plot_grid( plotlist =  ggl )
  
}

### pathway Fisher test
{
  PopSize <- length(unique( entr2gName$external_gene_name))
  
  sigload_genes_pathAnnot.hyper<- Reduce( cbind,  lapply( unique(cclust) , function(nname)
    {
    
    print(nname)
    
    # prepare gene set, convert ensembleIDs to entrezIDs
    geneset <- names(cclust)[cclust==nname]
    list21 <- length(geneset )
    # geneset <- unique(datt$symbol[ datt$fromOverlappingOrNearest == "Overlapping"])
    # test all genes, cis and trans regulated genes
    
    XX <- sapply( seq(pathDB), function(jj){
      list1 <-  length( pathDB[[jj]] )
      overlap <- length( intersect( geneset ,
                                    pathDB[[jj]]) )

      if(overlap>1) phyper( overlap -1 , list21 , PopSize , list1  , lower.tail= FALSE) else NA
    }) 
    names(XX) <- names(pathDB)

    return(XX)
  } ) )
  
  
  # adjust p-values
  sigload_genes_pathAnnot.hyper.adj <- apply( sigload_genes_pathAnnot.hyper, 2, p.adjust , "fdr")
  colnames(sigload_genes_pathAnnot.hyper.adj) <- colnames(sigload_genes_pathAnnot.hyper) <- 
    c("damage signature", "mixed big","mixed small","exprmnt only")
  
  
  topN_fisher  <- Reduce( union , sapply( seq(unique(cclust)) , function( ii,N=10 )
    {
    print(ii)
    datt <- sigload_genes_pathAnnot.hyper[,ii]
    # select from primary terms only
    datt <- datt[ datt<0.05 &  !is.na(datt)]
    datt <- sigload_genes_pathAnnot.hyper[names(datt),ii]
    datt <- datt[order(datt)]
    # return N top primar terms
    return( names( datt[ 1:min(N, length(datt)) ]) )
  }))
  topN_fisher <- topN_fisher[ ! topN_fisher %in% c(
    "mediumpurple1","turquoise3","salmon", "chartreuse3","SARS-CoV-2 Infection__REACT")]
 
  # cluster 
  toPlot <- sigload_genes_pathAnnot.hyper[ topN_fisher , ]
  # toPlot[ toPlot >=0.05 ] <- NA
  toPlot <-  apply( apply(toPlot, 2, rank,na.last='keep'), 2, rank_normalize)
  
  rownames(toPlot) <- sapply( sub( "__.*","",rownames(toPlot)) , str_trunc, 50)
  toPlot[is.na(toPlot) ]<- 1.1
  better_col_palette <- viridis::viridis(9,direction = -1)
  ggl <- pheatmap::pheatmap( toPlot , cluster_cols = F , 
                             color = better_col_palette , fontsize = 16)
  
  
  
  pdf( width =  8 , height = 8 , file="PDScor_path.Fisher_barplot.pdf")
  a<-dev.cur()
  png( width =  600 , height = 600 , file="PDScor_path.Fisher_barplot.png")
  dev.control("enable")
  print(ggl)
  dev.copy(which=a)
  dev.off()
  dev.off()
  
}

#### 2. combining TRN with PDS to find interesting TFs #### 
# load podo Paths/gsets
PodoPathGSet <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PodoPathGSet.05.01.24.7plus.rda")
circ.genes.bulk0.01 <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/circ.genes.bulk0.01.rda")
PodoPathGSet_plus <- c(PodoPathGSet , circ.genes.bulk0.01=list(circ.genes.bulk0.01 ))
### load the prior
# ## Wt1 and Tead1 Chipseq targets
WT1wt2014_TFtc <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/Chipseq/Chipseq.WT1wt2014_TFtc.rda")
WT1wt2014_TFtc.sig <- WT1wt2014_TFtc$score * ( WT1wt2014_TFtc$qvalue <0.3)
Tead1wt2018_TFtc <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/Chipseq/Chipseq.Tead1wt2018_TFtc.rda")
Tead1wt2018_TFtc.sig <-  Tead1wt2018_TFtc$score * ( Tead1wt2018_TFtc$qvalue <0.3)
## ATACseq prior
# ATACseq_tgenesM.TFtc <- readRDS( file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/ATACseq_TFtc.tgenesM_77.TF.rda")
# ATACseq_tgenesM <- readRDS( file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/ATACseq_tgenesM.rda")
ATACseq_tgenes <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/ATACseq/Podocytes/TOBIAS/FIMO/atac.podo_tobias.fimo110TF.cut1.p5e4_TFtc.rda")
ATACseq_tgenes.S <- Reduce( cbind , lapply( ATACseq_tgenes , "[[", "score" ))
ATACseq_tgenes.Q <- Reduce( cbind , lapply( ATACseq_tgenes , "[[", "qvalue" ))
colnames(ATACseq_tgenes.S) <- colnames(ATACseq_tgenes.Q) <- names(ATACseq_tgenes)
ATACseq_tgenesM <- ATACseq_tgenes.S * ( ATACseq_tgenes.Q <0.1)
rownames(ATACseq_tgenesM) <- rownames(ATACseq_tgenes.S) <- rownames(ATACseq_tgenes.Q) <- rownames(ATACseq_tgenes[[1]])
# combine ATAC and Chip
identical(rownames(WT1wt2014_TFtc), rownames(ATACseq_tgenesM))
ATAC.Chip <- ATACseq_tgenesM
ATAC.Chip[ , "Cluster_1__Wt1"] <- (ATAC.Chip[ , "Cluster_1__Wt1"] + WT1wt2014_TFtc.sig)/2
ATAC.Chip[ , "Cluster_5__Tead1"] <- (ATAC.Chip[ , "Cluster_5__Tead1"] + Tead1wt2018_TFtc.sig)/2

# ATACseq_tgenesM <- ATACseq_tgenesM[ rowSums(ATACseq_tgenesM)>0 , colSums(ATACseq_tgenesM)>0  ]
ATACseq_tgenesM <- ATAC.Chip[ rowSums(ATAC.Chip)>0 , colSums(ATAC.Chip)>0  ]
ATACseq_tgenesM.podo <- ATACseq_tgenesM[ rownames( ATACseq_tgenesM ) %in% allPodoGenes ,  ]

### visualising TRN
  {
    ### black-white heatmap
    toPlot <- ATACseq_tgenesM.podo
    # toPlot <- toPlot[ rowSums(toPlot)>0 , colSums(toPlot)>0  ]
    # toPlot <- toPlot[ rownames( toPlot ) %in% allPodoGenes ,  ]
    toPlot <- apply(toPlot , 2, rank,ties.method= "min")
    gg1 <- pheatmap::pheatmap(toPlot, treeheight_row = 100, show_rownames = F,
                       clustering_method = "ward.D2",
                       col= brewer.pal(9, "Greys"), treeheight_col = 0)
    

    pdf( width =  5 , height = 6 , file="TRN.14.TFtcCut0.1_blwh.heatmap.pdf")
    a<-dev.cur()
    png( width =  400 , height = 400 , file="TRN.14.TFtcCut0.1_blwh.heatmap.png")
    dev.control("enable")
    print( gg1 )
    dev.copy(which=a)
    dev.off()
    dev.off()
    
    
    ### network vew
    library(igraph)
    TF_MeanMed <- readRDS(  file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/TFA/TFs_SNSC.podo.MeanMed.08.01.24.rda")
    
    PodoTFnet <-ATACseq_tgenesM.podo[ (rownames(ATACseq_tgenesM.podo) %in% TF_MeanMed) & 
                                        rownames(ATACseq_tgenesM.podo) %in%  rownames(genecorrPDS.Sprmn_freq)[genecorrPDS.Sprmn_freq$total>=4], ]
    
    PodoTFnet <- PodoTFnet[ , colSums(PodoTFnet!=0)!=0 ]
    PodoTFnet <- PodoTFnet[ rowSums(PodoTFnet!=0)!=0,]

    edge_list <- reshape2::melt(t(PodoTFnet))
    edge_list <- edge_list[edge_list$value!=0,]
    toPlot <- igraph::graph_from_edgelist( as.matrix( edge_list[,1:2] ), directed = T )
    
    # 
    
  plot( toPlot, vertex.size=6, directed=T,
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

### heatmap of intersect of gSets and TFtargets
  {
    
### test enrichment of TF targets in gSets 
    PodoPathGSet_PDS.test <- sapply( PodoPathGSet , function(X) 
      {
      # print(head(X))
      gset.TFqval.xprmnt <-  gset.TFqval( geneSet =  X ,
                                          GRN = ATACseq_tgenesM.podo )
      
    } )
    colnames(PodoPathGSet_PDS.test) <-names( PodoPathGSet ) 
    
    
    # PodoPathGSet_PDS.test[ PodoPathGSet_PDS.test >  0.05] <- NA
    PodoPathGSet_PDS.test <- PodoPathGSet_PDS.test[ rowSums(!is.na(
      PodoPathGSet_PDS.test))>1,
      colSums(!is.na(
        PodoPathGSet_PDS.test))>1]
    PodoPathGSet_PDS.test <- t(PodoPathGSet_PDS.test)
    siglbl <- round( PodoPathGSet_PDS.test , 4)
    siglbl <- ifelse( siglbl < 0.001, "***", 
                      ifelse( siglbl< 0.01, "**",
                              ifelse( siglbl< 0.1,"*","")))

### calculate intersect between gSets and TFtargets  
  datt <- ATACseq_tgenesM.podo
  toPlot <- Reduce( rbind, lapply( PodoPathGSet , function(gSet){ 
    apply( datt, 2, function(TF ){
      TF <- TF[ TF>0 & names(TF) %in% gSet] 
      print( length(TF))
    })
  }))
  rownames(toPlot) <-  names( PodoPathGSet ) 
  
  # toPlot <- toPlot[rowSums(toPlot)>10,]
  # # scale
  toPlot <-t(scale(t(toPlot),  center = F ))
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
   tgene_N = log( colSums( ATACseq_tgenesM.podo[ ,colnames(toPlot)]>0 )) 
 )
 annotation_row = data.frame(
   row.names = rownames(toPlot) ,
   gSet_size =  log( Reduce( c, lapply( PodoPathGSet_plus[ rownames(toPlot) ], length)) )
 )
 # colors for annotations
 ann_colors = list(
   
   gSet_size =  ( setNames( map2color( annotation_row$gSet_size ,
                                       rev(grey.colors(nrow(toPlot)+1))), 
                                       annotation_row$gSet_size ) ),
   tgene_N =  (  setNames(  map2color( annotation_col$tgene_N  ,
                                       rev(grey.colors(ncol(toPlot)+1))), 
                           annotation_col$tgene_N ) )
 )
 
 # specify a color pallete
 better_col_palette <- viridis::magma(30)
 ### make a plot
  gg1 <- ComplexHeatmap::pheatmap( toPlot , 
                             color =  better_col_palette ,
                             cluster_cols = T, 
                             cluster_rows = T,  
                             number_color = "white", 
                              annotation_row = annotation_row ,
                              annotation_col = annotation_col,
                             annotation_colors = ann_colors ,
                             annotation_legend = F,
                             fontsize_number = 18 ,
                             fontsize = 14,
                             display_numbers = siglbl2 
                             
                             )
  
  
  pdf(height = 6, width = 8, file="podoPaths.TFprior_overlap.test_heatmap.rank.pdf")
    ht <- draw(gg1) 
  dev.off()
  
  }


### Suppl. fig TF PWM simmilarity tree
### made by TOBIAS

### Suppl. fig motif enrichment of geneSet promoters
  { 
  mmotifs <- universalmotif::read_meme( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/ATACseq/Podocytes/TOBIAS/motif_clusters/motif_comparison_seqcor_Ward.cut.1_consensus_motifs.meme" )
  # mmotifs <- universalmotif::read_meme( "/media/tim_nevelsk/WD_tim/ANNOTATIONS/CISBP/CISBP2mouse_podoTF_08.01.24.meme" )
    
    names(mmotifs) <- sapply(seq(mmotifs) , function(ii) mmotifs[[ii]]@name )
    mmotifs <- mmotifs[colnames(ATACseq_tgenesM.podo)]
  # PodoPathGSet_ame <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/AME_motifAnalysis/PodoPathGSet/15meta.podoMotifs.allBckgrnd.ranksum/PodoPathGSet_AMEmotifEnrich.rda")
  PodoPathGSet_ame <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/AME_motifAnalysis/PodoPathGSet/15meta.podoMotifs.podoBckgrnd.ranksum/PodoPathGSet_AMEmotifEnrich.rda")
  #  PodoPathGSet_ame <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/AME_motifAnalysis/PodoPathGSet/110podoMotifs.allBckgrnd.ranksum/PodoPathGSet_AMEmotifEnrich.rda")
  
  
  names(PodoPathGSet_ame) <- names(PodoPathGSet)
  
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

# ### Suppl. fig TRN related to Wt1
# toPlot <- ATACseq_tgenesM.podo[
#   rownames(ATACseq_tgenesM.podo) %in% 
#     grep("col",  rownames(ATACseq_tgenesM.podo), value = T, ignore.case = T), ]
# toPlot <- reshape2::melt(t(toPlot))
# net <- igraph::graph_from_data_frame(d= toPlot , directed=T) 
# plot(net, head.arrow.size=2)

#### 3. show correlation between PDS and circadian (dis)regulation #### 

### Suppl.Fig. plot basics of circadian bulk data
  {
  # load normilised by DEseq2 counts
  dds_rlog.norm <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/bulkRNAseq.Kidney.cirk_rlog.rda")
  # remove one sample with very few counts
  countMat <- countMat[, colnames(countMat)!="SN8580185_21056_806Aligned.sortedByCoord.out.bam"]
  annot <- gsub(".*_|Aligned.*", "", colnames(countMat))
  annot <-data.frame(  ID= colnames(countMat) ,
                       condition= annot,
                       age=c(rep(10,16), rep(80,15)),
                       time= sub("^..","",annot))
  annot$age <- as.factor( annot$age )
  annot$time <- as.numeric( annot$time )
  
  
  ### load results
  MetaCycle_young <- readRDS( "MetaCycle_bulk.kidney_young.rda")
  MetaCycle_old <- readRDS( "MetaCycle_bulk.kidney_old.rda")
  
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
                                                       c("Bmal1","Nr1d1","Cry2","Per1")],]
  toPlot.melt <- reshape2::melt(toPlot)
  toPlot.melt$Time <- annot$time[ match( toPlot.melt$Var2, annot$ID)]
  toPlot.melt$Gene <- tx2gene$external_gene_name[ match( toPlot.melt$Var1, tx2gene$ensembl_gene_id)]
  
     # Create the plot
    p <- ggplot(toPlot.melt, aes_string(x = "Time", y = "value")) +
      geom_smooth(color = "blue", size = 1,span=0.3) +
      geom_point(color = "red", size = 2) +
      labs( x = "Time (hours)",
           y = "Expression Level") +
      theme_minimal()+ facet_wrap("Gene",scales = "free")+ theme(
        text=element_text( size=20)
      )
    
    pdf( width = 8 , height = 6 , file="Supl.Fig4/bulkRNAseq_circadian_4CoreClockExpr.pdf")
    a<-dev.cur()
    png( width =  600 , height = 400 , file="Supl.Fig4/bulkRNAseq_circadian_4CoreClockExpr.png")
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
    
    ## calculate likely circadian time of a sample deviation from it of each cell 
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

