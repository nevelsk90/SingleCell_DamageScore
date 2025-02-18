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
library(cowplot)
library( reshape2 )
library( plyr )
library( Seurat )
library( viridis )
library( ggthemes )
library( AUCell )
library( ggpubr )
library( biomaRt )
library( RColorBrewer )


# setwd
setwd("/home/tim_nevelsk/PROJECTS/PODOCYTE/WRITING/PDS_manuscript/Figures/Figure3")
inputdir <-  "/media/tim_nevelsk/WD_tim/PROJECTS/WRITING/PDS_manuscript/Figure_input" 

# load necessary code
source("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/cell-damage-score/AUCell_script.r")
source("/home/tim_nevelsk/PROJECTS/myCode/usefulRfunc.r")

# load damage signatures
DS_all <- read.table( header = T,  file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")
# podocyte expr.genes
allPodoGenes <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/SCSN_allPodoGenes.rda")

## load expression data
listSCSN.1K.sampl <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN_samples.1K.22.12.23.rda")
listSCSN <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN.PDS_22.12.23.rda")

## pathway info
keggPath_pathlist <- readRDS("/media/tim_nevelsk/WD_tim/ANNOTATIONS/Pathways/kegg_pathsList_gName.07.09.23.rda")
reactPath_pathlist <- readRDS("/media/tim_nevelsk/WD_tim/ANNOTATIONS/Pathways/reactome_pathsList_gName.07.09.23.rda")

### load the data and functions
### load AUCell score
# listSCSN.1K.sampl_PAaucell <- readRDS( file= "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PathwayActivity/listSCSN.PDS_1K_PAaucell.29.09.2023rda")
listSCSN.1K.sampl_PAaucell <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN.PDS_samples.1K_PAaucell.13.05.24.rda")

# separate pathways by database pathways
listSCSN.PDS_1K_PAaucell.kegg <- lapply(listSCSN.1K.sampl_PAaucell, "[[", 1)
listSCSN.PDS_1K_PAaucell.react <- lapply(listSCSN.1K.sampl_PAaucell, "[[", 2)
listSCSN.PDS_1K_PAaucell.podoPath <- lapply(listSCSN.1K.sampl_PAaucell, "[[", 3)
# set smoothing factors
WsizeF=10
WstepF=30

### a function to plot pathway fingerprint
fingerprintPlot <-function( PAtoPlot , 
                            collAnnot , 
                            podoPaths=F ,
                            rowNames = T, 
                            globalAnnotLimit=T , # adjust limits for AlbCr and PDS to see difference in spans across studies
                            rowOrder= rownames(MMpval)[rev(gsea.order$rowInd)] )
  {
  # print(ii)
  require(ComplexHeatmap)
  
  
  ### prepare annotation rows that show percentage 
  ### of cells of a certain gtype in a window
  annot <- (  collAnnot )
  annot <- annot[ ,!is.na(colSums(annot)), drop=F ]
  
  if( min(annot) == max(annot)) { 
    minn <- 0 
    maxx = max(annot)
    
  } else if( isTRUE(globalAnnotLimit)) { # adjust limits for AlbCr and PDS to see difference in spans across studies
    minn1 <- PDSmin
    maxx1 <- PDSmax
    col_fun1 <- circlize::colorRamp2( 
      breaks = c( minn1 , (minn1+maxx1)/2, maxx1 ) ,
      colors= c(  "white","orange2", "brown4") )
    # breaks = c( minn1 , maxx1 ) ,
    # 
    # colors= c(  "white","orange2") ) 
    
    minn2 <- AlbCr.min
    maxx2 <- AlbCr.max
    col_fun2 <- circlize::colorRamp2( 
      breaks = c( minn2 , maxx2 ) ,
      colors= c(  "white","black") ) 
    
    ### annotation
    ha = HeatmapAnnotation( df=as.data.frame( annot ),
                            col = list( AlbCr = col_fun2 ,
                                        PDS = col_fun1 ),
                            simple_anno_size = unit(1.5/ncol(annot), "cm"),
                            show_legend =  FALSE
                            # show_heatmap_legend==F
                            
    )
  } else {
    maxx = max(annot)
    minn = min(annot)
    
    col_fun <- circlize::colorRamp2( 
      breaks = c( minn , (minn+maxx)/2, maxx ) ,
      colors= c(  "white","orange2", "brown4") ) 
    
    ### annotation
    ha = HeatmapAnnotation( foo=annot ,
                            col = list(all = col_fun ),
                            simple_anno_size = unit(1.5/ncol(annot), "cm"),
                            show_legend =  FALSE
                            # show_heatmap_legend==F
                            
    )
  }
  
  
  
  
  
  # treat differently if podopaths are to plot
  if( isTRUE(podoPaths)) {
    paths <- names(PodoPathGSet_pathlist)
    rowOrder<- seq(PodoPathGSet_pathlist)
  }
  
  toPlot <- t( scale(as.matrix(PAtoPlot)))
  toPlot <- toPlot[ rowOrder , ]
  if(isFALSE(rowNames)) rowNames =rep("", nrow(toPlot))  else if(isTRUE(rowNames)) {
    rowNames <- rowOrder }
  
  gg <-  ComplexHeatmap::Heatmap(  
    toPlot ,
    row_names_gp = gpar(fontsize = 15),
    row_order = rowOrder ,
    column_labels =rep("", ncol(toPlot))  ,
    row_names_side = "left",
    row_labels = rowNames  ,
    cluster_columns =  F , 
    cluster_rows = F , 
    col = viridis(100),
    # col = viridis(max(PAtoPlot, na.rm=T)+1),   #labCol = FALSE,
    # column_split=gtype ,
    top_annotation= ha, 
    name = "AUCell PAscores"
  )
  
  return(gg)
}

# a function to plot dot heatmap of changes in activity 
# of an individual pathway, across FSGS models
PAdotHeatPlot <- function( cells="all", 
                           single_path=sepLath[[ii]],
                           dataSetIndex = c(1:3,5:9) )
  {
  
  ## cells - a charachter vector which specifies whether to use all cells,
  # when set to "all", or a specific genotype, for the analysis.
  
  require(scales)
  require(colorspace)
  
  # path to plot
  toPlot <- Reduce( rbind, lapply( dataSetIndex , function( ii ) 
  {
    
    print(ii)
    
    if( single_path %in% names(reactPath_pathlist) ) {
      datt <- listSCSN.PDS_1K_PAaucell.react[[ii]]
    } else datt <- listSCSN.PDS_1K_PAaucell.kegg[[ii]]
    
    dataSet <- listSCSN.1K.sampl[[ii]]
    # downsample
    Idents(dataSet) <- "group"
    dataSet <- subset( dataSet, downsample=1000 )
    
    if( cells=="all" )  {
      Idents(dataSet) <- "gtypeDE"
      ## check if GSE data has control cells 
      # and if not (nephritis day 5 and doxo) then add controls  
      if( length(unique(dataSet$gtypeDE))==1) dataSet <- merge( 
        dataSet,   subset( listSCSN.1K.sampl[[8]], gtypeDE=="control") )
          # dataSet <- subset( dataSet, downsample=1000 )
          
          ## add controls to PA data
          if( single_path %in% names(reactPath_pathlist) ) {
            XX <- listSCSN.PDS_1K_PAaucell.react[[8]]
          } else XX <- listSCSN.PDS_1K_PAaucell.kegg[[8]]
          datt <- cbind( datt, XX[ match(rownames(datt) , rownames(XX)),] )
        
      
      ### normalise PA by control means
      # toPlot <- sweep(t( datt[ rownames(datt) == single_path, colnames(dataSet), drop=F] ), 
      #                 2, rowMeans(datt[ rownames( datt ) == single_path,
      #                                   colnames(dataSet)[dataSet$gtypeDE=="control"], drop=F]), 
      #                 FUN = '/')            
      # toPlot <- log2(toPlot) #  get log2FC of PA
          toPlot <-  scale( t( datt[rownames(datt) == single_path, 
                                    colnames(datt) %in% colnames(dataSet), drop=F]) , center =T)
      colnames(toPlot) <- single_path
      
     

    } else  {  
      #select cells of spec. genotype(s) 
      dataSet <- subset( dataSet, gtypeDE==cells )
      
      ### normalise PA by scaling  
      toPlot <-  scale( t( datt[rownames(datt) == single_path, 
                                colnames(datt) %in% colnames(dataSet), drop=F]) , center =T)
      colnames(toPlot) <- single_path
      
    } 
    
    toPlot <- cbind( toPlot, PDS=dataSet$PDS ) # add PDS
    toPlot <- toPlot[ order(toPlot[,"PDS"]),]  # order by PDS
    toPlot <- t(toPlot)
    toPlot[is.infinite(toPlot)] <- NA
    ## smooth both PA and PDS
    toPlot_smooth <-  apply( toPlot, 1 , function(x) slideFunct(
      x, ncol(toPlot)/WsizeF,  ncol(toPlot)/WstepF ))
    toPlot_smooth <- as.data.frame(toPlot_smooth)
    toPlot_smooth$dataSet <- names(listSCSN.1K.sampl)[[ii]]
    toPlot_smooth$PDS.index <- 1:nrow(toPlot_smooth)
    return(toPlot_smooth)
    
  } ) )
  
  
  if( cells=="experimental"){
    ### add controll cells
    XX <-  lapply(  listSCSN.1K.sampl[c(1:6,8)], subset, gtypeDE=="control")
    dataSet <- Reduce( rbind, lapply( seq(XX), function(ii){
      XX[[ii]][["PDS"]]
    }))
    # dataSet <- subset(listSCSN.1K.sampl[[1]], gtypeDE=="control")[["PDS"]]
    
    if( single_path %in% names(reactPath_pathlist))  {
      datt.wt <- listSCSN.PDS_1K_PAaucell.react[c(1:6,8)]
    } else datt.wt <- listSCSN.PDS_1K_PAaucell.kegg[c(1:6,8)]
    
    datt.wt <- Reduce( cbind, lapply( seq(datt.wt) , function(ii){
      datt.wt[[ii]][ rownames(datt.wt[[ii]]) == single_path ,, drop=F]
    }))
    
    datt.wt <- datt.wt[ , rownames(dataSet),  drop=F]
    datt.wt <- t( datt.wt )
    datt.wt <- scale(datt.wt)
    toPlot.wt <- cbind( datt.wt, PDS=dataSet[,"PDS"] ) # add PDS
    toPlot.wt <- toPlot.wt[ order(toPlot.wt[,"PDS"]),] # order by PDS
    toPlot.wt <- t(toPlot.wt)
    toPlot.wt[is.infinite(toPlot.wt)] <- NA
    ## smooth both PA and PDS
    toPlot_smooth.wt <-  apply( toPlot.wt, 1 , function(x) slideFunct(
      x, ncol(toPlot.wt)/WsizeF,  ncol(toPlot.wt)/WstepF ))
    toPlot_smooth.wt <- as.data.frame(toPlot_smooth.wt)
    toPlot_smooth.wt$dataSet <- "controls"
    toPlot_smooth.wt$PDS.index <- 1:nrow(toPlot_smooth.wt)
    
    toPlot <- rbind( toPlot, toPlot_smooth.wt )
    nnames <- c("controls", names(listSCSN.1K.sampl )[dataSetIndex] )
  } else nnames <- names(listSCSN.1K.sampl)[dataSetIndex]
  
  toPlot.m <- reshape2::melt(toPlot,id.vars=c("dataSet","PDS","PDS.index") )
  toPlot.m$dataSet <- factor( toPlot.m$dataSet , 
                              levels = rev(  nnames ))
  
  toPlot.m <- droplevels.data.frame( toPlot.m )
  
  
  
  gg<- ggplot( toPlot.m , aes(x=PDS, y= dataSet ))+
    # geom_boxplot(aes(color=experiment))+
    geom_point(aes( color=as.numeric(value)), size=8) +
    # scale_color_viridis( midpoint=0) +
    labs(x = "PDS",color = "PA Zscore") +
    facet_wrap( facets  = vars(variable)) +
    theme_bw() + 
    theme(text=element_text(size=24 ))
  
  # decide on a color pallete depending on wether only one or 2 gtypes
  # if(  cells=="all" ) {
  #   gg <- gg + scale_color_continuous_divergingx(palette = 'RdBu',
  #                                                mid = 0.0 , rev = T,
  #                                                p1=0.8, p2=0.8, p3=0.8, p4=0.8 ) 
  # } else  gg <- gg + scale_color_viridis() 
  
  gg <- gg + scale_color_viridis() 
  return(gg)
}


#### 1. Pathway fingerprint heamaps #### 


### KFO studies, annotated by PDS and AlbCr
  {

  ### load selected paths
    pthsPlt <- read.table( sep = "\t","/home/tim_nevelsk/PROJECTS/PODOCYTE/WRITING/PDS_manuscript/Figures/Figure3/KFO_select.paths.txt")

    pthsPlt <- pthsPlt[!(pthsPlt$V1 %in% c( "Muscarinic acetylcholine receptors *REACT",
                                            "Metabolism *REACT")),]
  
  ### cluster pathways
    {

    # table of similarity
    library(GeneOverlap)
    library(dendextend)
    
      RR <- reactPath_pathlist[ paste0( names(reactPath_pathlist),  " *REACT") %in% pthsPlt]
      names(RR) <- paste( names(RR), "*REACT")
      KK <-  keggPath_pathlist[paste0( names(keggPath_pathlist), " *KEGG") %in% pthsPlt ]
      names(KK) <- paste( names(KK), "*KEGG")
    patLlist <- c( RR, KK )
    
    MM<- GeneOverlap::newGOM( patLlist , patLlist ,genome.size=20000)
    MMpval <- getMatrix(MM, name="odds.ratio")
    MMpval[ !is.finite(MMpval)] <- max( MMpval[ is.finite(MMpval)] )
    # MMpval[MMpval==0] <- NA
    # diag(MMpval) <- NA
    MMpval <- log10(MMpval+1)
    
    ## plot
    # pdf( width =  12 , height = 9, file="gseaKEGG.REACT.top3_pvalHeat.clust.pdf")
    # a<-dev.cur()
    # png( width =  1200 , height = 800 , file="gseaKEGG.REACT.top3_pvalHeat.clust.png")
    # dev.control("enable")
    gsea.order <- gplots::heatmap.2( MMpval,
                                     col = RColorBrewer::brewer.pal( name = "Greens", n=9), 
                                     trace = "none", symkey=T,
                                     na.color = "#00441B",
                                     hclustfun = function(X) hclust(X, "ward.D2"),
                                     dendrogram = "row",margins = c(2,25) ,
                                     labCol = FALSE, cexRow= 1.2,
                                     ColSideColors = map2color(log10(sapply( patLlist , length)),
                                                               pal = gray.colors(9, rev = T)))
    
    
    
  }

  #### plot KFO samples, annotate by albuminuria
    {
    ### compute number of cells within low/mid/high PDS per window
    KFOannot <- read.table( sep="\t", header = T,
                            "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/Sample_Names_KFO.csv")
    
    samplTag <- c("BDP","AKK","BWV|CDC")
    
    # make plots for 3 KFO datasets
    ggData<-  lapply( 1:3, function(ii){
      
      # print(ii)
      dataSet <- listSCSN.1K.sampl[[ii]]
      dataSet$sample <- sub( "SID", "", dataSet$sample )
      # select metadata for samples on one model
      Albuminuria <- KFOannot[ grep( samplTag[ii] , 
                                     KFOannot$Sample_Name, ignore.case = T),]
      # Albuminuria <- KFOannot[ KFOannot$Genotype =="Nphs2",]
      Albuminuria <-  Albuminuria[!is.na(Albuminuria$AlbCrRatio_v2),]
      # get cells of a specific genotype
      datt <- subset( dataSet , 
                      subset = sample %in% as.character( Albuminuria$CCG_Sample_ID ))
      # downsample
      Idents(datt) <- "group"
      datt <- subset( datt, downsample=1000 )
      Idents(datt) <- "gtypeDE"
      # datt <- subset( datt, downsample=1000 )
      # datt <- subset( datt, gtypeDE=="experimental" )
      
      datt$AlbCr <- Albuminuria$AlbCrRatio_v2[ match( datt$sample, Albuminuria$CCG_Sample_ID)]
      
      # order according to PDS
      expPDSorder <- colnames( datt )[ order( datt$PDS) ]
      
      # create column annotations
      Alb.binFr <- datt@meta.data[  expPDSorder , c("AlbCr","PDS")]
      Alb.binFr_smooth <- apply( Alb.binFr, 2, function(XX) slideFunct(XX, nrow(Alb.binFr)/WsizeF, 
                                                                       nrow(Alb.binFr)/WstepF ))
      # Alb.binFr_smooth <- scale(Alb.binFr_smooth)
      
      
      # extract PA info, select pathways
      # order PA by PDS
      bothDB <-  rbind( listSCSN.1K.sampl_PAaucell[[ ii ]]$cells_AUC.kegg[ , expPDSorder ],
                        listSCSN.1K.sampl_PAaucell[[ ii ]]$cells_AUC.react[ , expPDSorder ] )
      bothDB <- bothDB[ match( sub(" \\*KEGG| \\*REACT","",pthsPlt), 
                               rownames( bothDB) ),  ]
      rownames( bothDB ) <-  pthsPlt
      
      
      datTOplot <- bothDB
      
      ## smooth pathway activity
      path_smooth <- apply( datTOplot , 1 , function(x) slideFunct(x, ncol(datTOplot)/WsizeF,
                                                                  ncol(datTOplot)/WstepF ))

      ll <- list(path_smooth, Alb.binFr_smooth)
      names( ll) <- c("path_smooth", "Alb.binFr_smooth")
      return(ll )

    })
    names(ggData) <- names(listSCSN.1K.sampl)[1:3]
    
    
    ### plot
    PDSmin <- min( unlist( lapply( seq(ggData), function(ii){
      ggData[[ii]]$Alb.binFr_smooth[,2]
    })), na.rm = T)
    PDSmax <- max( unlist( lapply( seq(ggData), function(ii){
      ggData[[ii]]$Alb.binFr_smooth[,2]
    })), na.rm = T)
    AlbCr.min <-  min( unlist( lapply( seq(ggData), function(ii){
      ggData[[ii]]$Alb.binFr_smooth[,1]
    })))
    AlbCr.max <- max( unlist( lapply( seq(ggData), function(ii){
      ggData[[ii]]$Alb.binFr_smooth[,1]
    })))
      
   gglist <-  lapply( seq(ggData), function(ii)
     {
      
      path_smooth <- ggData[[ii]]$path_smooth
      Alb.binFr_smooth <- ggData[[ii]]$Alb.binFr_smooth
        
      fingerprintPlot( PAtoPlot = path_smooth ,
                              collAnnot= Alb.binFr_smooth ,
                              rowNames = F,
                              globalAnnotLimit = T,
                              rowOrder= rownames(MMpval)[rev(gsea.order$rowInd)]
      )
    })
   

    gg.ctrl <- lapply( c(1:3) , function(ii)
      { 
      
       print(ii)
      datt <- listSCSN.1K.sampl[[ii]]
      datt$sample <- sub( "SID", "", datt$sample )
     

      # downsample
      Idents(datt) <- "group"
      datt <- subset( datt, downsample=1000 )
      # select controls
      Idents(datt) <- "gtypeDE"
      datt <- subset( datt, gtypeDE=="control" )
      

      # order according to PDS
      expPDSorder <- colnames( datt )[ order( datt$PDS) ]
      
      # create column annotations, use only PDS
      Alb.binFr <- datt@meta.data[  expPDSorder ,"PDS"]
      Alb.binFr_smooth <- slideFunct(Alb.binFr, length(Alb.binFr)/WsizeF, 
                                     length(Alb.binFr)/WstepF )
      # Alb.binFr_smooth <- scale(Alb.binFr_smooth)
      
      
      # extract PA info, select pathways
      # order PA by PDS
      bothDB <-  rbind( listSCSN.1K.sampl_PAaucell[[ii]]$cells_AUC.kegg[ , expPDSorder ],
                        listSCSN.1K.sampl_PAaucell[[ii]]$cells_AUC.react[ , expPDSorder ])
      bothDB <- bothDB[ match( sub(" \\*KEGG| \\*REACT","",pthsPlt), 
                               rownames( bothDB) ),  ]
      rownames( bothDB ) <-  pthsPlt
      
      
      datTOplot <- bothDB
      
      ## smooth pathway activity
      path_smooth <- apply( datTOplot, 1 , function(x) slideFunct(x, ncol(datTOplot)/WsizeF, 
                                                                  ncol(datTOplot)/WstepF ))
      
      ll <- list(path_smooth, Alb.binFr_smooth )
      names(ll) <-  c( "path_smooth", "Alb.binFr_smooth" )
      return( ll )
    
      })
    
    path_smooth.ctrl <- Reduce( "+", lapply(gg.ctrl, "[[", 1 ) )/3
    Alb.binFr_smooth.ctrl <- Reduce(  "+", lapply(gg.ctrl, "[[", 2 ) )/3
    gg.ctrl.plot <-  fingerprintPlot( PAtoPlot = path_smooth.ctrl , 
                            collAnnot= cbind( AlbCr=Alb.binFr_smooth.ctrl ,
                                              PDS= Alb.binFr_smooth.ctrl),
                            globalAnnotLimit = T,
                            rowNames = T,
                            rowOrder= rownames(MMpval)[rev(gsea.order$rowInd)]
    )
    
    gg.list.all <- c( gg.ctrl.plot, gglist )
    
    pdf(height = 8, width = 18, "podoPaths.fingerprints_KFO_both.select.14.11.24.pdf")
    
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nr = 1, nc = 4)))
    
    lapply( seq(gg.list.all), function(ii)
      {
      pushViewport(viewport(layout.pos.row = 1,
                            layout.pos.col = c(1:4)[ii]))
      draw( gg.list.all[[ ii ]], newpage = FALSE, 
            column_title= c( "controls_mean", names(listSCSN.1K.sampl)[1:3])[ii],
            heatmap_legend_side = "bottom"      )
      upViewport()
    })
    
    dev.off()
    
  }
  
}

### Suppl.fig. all studies, annotated by PDS only
  {
    
  ### select top or common Paths
  {
    PAaucell.react_vs_PDS <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PathwayActivity/Multi.pathways/PAaucell.react_vs_PDSrank.gtypeAgg.rda")
    PAaucell.kegg_vs_PDS <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PathwayActivity/Multi.pathways/PAaucell.kegg_vs_PDSrank.gtypeAgg.rda")
    
    N <- 5
    PAaucell.react_topN.KFO <- Reduce( union, apply(
      PAaucell.react_vs_PDS, 2, function(X) names(X[order(X)])[1:N]))
    PAaucell.kegg_topN.KFO <- Reduce( union, apply(
      PAaucell.kegg_vs_PDS, 2, function(X) names(X[order(X)])[1:N]))
    
    
    
  }
  
  ### cluster pathways
  {
    
    # table of similarity
    library(GeneOverlap)
    library(dendextend)
    
    RR <- reactPath_pathlist[names(reactPath_pathlist) %in% PAaucell.react_topN.KFO]
    names(RR) <- paste( names(RR), " *REACT")
    KK <-  keggPath_pathlist[names(keggPath_pathlist) %in% PAaucell.kegg_topN.KFO]
    names(KK) <- paste( names(KK), " *KEGG")
    patLlist <- c( RR, KK )
    
    MM<- GeneOverlap::newGOM( patLlist , patLlist ,genome.size=20000)
    MMpval <- getMatrix(MM, name="odds.ratio")
    MMpval[ !is.finite(MMpval)] <- max( MMpval[ is.finite(MMpval)] )
    # MMpval[MMpval==0] <- NA
    # diag(MMpval) <- NA
    MMpval <- log10(MMpval+1)
    
    ## plot
    # pdf( width =  12 , height = 9, file="gseaKEGG.REACT.top3_pvalHeat.clust.pdf")
    # a<-dev.cur()
    # png( width =  1200 , height = 800 , file="gseaKEGG.REACT.top3_pvalHeat.clust.png")
    # dev.control("enable")
    gsea.order <- gplots::heatmap.2( MMpval,
                                     col = RColorBrewer::brewer.pal( name = "Greens", n=9), 
                                     trace = "none", symkey=T,
                                     na.color = "#00441B",
                                     hclustfun = function(X) hclust(X, "ward.D2"),
                                     dendrogram = "row",margins = c(2,25) ,
                                     labCol = FALSE, cexRow= 1.2,
                                     ColSideColors = map2color(log10(sapply( patLlist , length)),
                                                               pal = gray.colors(9, rev = T)))
    
    
    write.table(rownames(MMpval)[rev(gsea.order$rowInd)], sep = "\t",quote = F, row.names = F, 
                file="Supl.Fig3/allData_top5.paths_clust.tsv")
  }
  
  ### make plots for all datasets
  gglist <-  lapply( seq(listSCSN.1K.sampl), function(ii)
    {
    
    
    print(ii)
    datt <- listSCSN.1K.sampl[[ii]]
    
    # downsample
    Idents(datt) <- "group"
    datt <- subset( datt, downsample=1000 )
    Idents(datt) <- "gtypeDE"
    datt <- subset( datt, downsample=1000 )
    
    
    # order according to PDS
    expPDSorder <- colnames( datt )[ order( datt$PDS) ]
    
    # create column annotations
    Alb.binFr <- datt@meta.data[  expPDSorder , "PDS", drop=F]
    Alb.binFr_smooth <- apply( Alb.binFr, 2, function(XX) slideFunct(XX, nrow(Alb.binFr)/WsizeF, 
                                                                     nrow(Alb.binFr)/WstepF ))
    Alb.binFr_smooth <- scale(Alb.binFr_smooth)
    
    
    # extract PA info, select pathways
    # order PA by PDS
    if( !is.na(PAaucell.kegg_topN.KFO) & !is.na(PAaucell.react_topN.KFO)){
      KEGG <-  listSCSN.1K.sampl_PAaucell[[ii]]$cells_AUC.kegg[ 
        match( PAaucell.kegg_topN.KFO, rownames( listSCSN.1K.sampl_PAaucell[[ii]]$cells_AUC.kegg) ), expPDSorder ]
      rownames( KEGG ) <-  paste(  PAaucell.kegg_topN.KFO, " *KEGG")
      REACT <- listSCSN.1K.sampl_PAaucell[[ii]]$cells_AUC.react[ 
        match( PAaucell.react_topN.KFO , rownames(listSCSN.1K.sampl_PAaucell[[ii]]$cells_AUC.react)), expPDSorder ]
      rownames( REACT ) <-  paste(  PAaucell.react_topN.KFO , " *REACT")
      
      datTOplot <- rbind( KEGG ,REACT )
    } else if( is.na(PAaucell.react_topN.KFO) ) {
      datTOplot <-  listSCSN.1K.sampl_PAaucell[[ii]]$cells_AUC.kegg[ 
        match( PAaucell.kegg_topN.KFO, rownames( 
          listSCSN.1K.sampl_PAaucell[[ii]]$cells_AUC.kegg) ), expPDSorder ]
      rownames( datTOplot ) <-  paste(  PAaucell.kegg_topN.KFO, " *KEGG")
    } else if( is.na(PAaucell.kegg_topN.KFO) ) {
      datTOplot <-  listSCSN.1K.sampl_PAaucell[[ii]]$cells_AUC.react[ 
        match( PAaucell.react_topN.KFO, rownames( 
          listSCSN.1K.sampl_PAaucell[[ii]]$cells_AUC.react) ), expPDSorder ]
      rownames( datTOplot ) <-  paste(  PAaucell.react_topN.KFO , " *REACT")
    }
    
    
    ## smooth pathway activity
    path_smooth <- apply( datTOplot, 1 , function(x) slideFunct(x, ncol(datTOplot)/WsizeF, 
                                                                ncol(datTOplot)/WstepF ))
    
    
    # toPlot <- t( scale(as.matrix(PAtoPlot)))
    # toPlot <- toPlot[ rowOrder , ]
    
    ### plot
    gg <-  fingerprintPlot( PAtoPlot = path_smooth , 
                            collAnnot= Alb.binFr_smooth ,
                            # podoPaths= T ,
                            rowNames = F
                            # rowOrder= NULL
    )
    return(gg)
  })
  
  
  
  
  
  
  pdf(height = 18, width = 18, "Supl.Fig3/podoPaths.fingerprints_allData_REACT.pdf")
  
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nr = 3, nc = 3)))
  
  lapply( seq(gglist), function(ii){
    pushViewport(viewport(layout.pos.row = rep(1:3,each=3)[ii],
                          layout.pos.col = rep(1:3,3)[ii]))
    draw( gglist[[ii]], newpage = FALSE, column_title= names(listSCSN.1K.sampl)[ii])
    upViewport()
  })
  
  dev.off()
  
  } 

### Suppl.Fig. expression of Oxphos and Glycolisis genes
  {
  ggenes1 <- keggPath_pathlist[["Glycolysis / Gluconeogenesis"]][ keggPath_pathlist[["Glycolysis / Gluconeogenesis"]] %in%  allPodoGenes ]
  ggenes2 <- keggPath_pathlist[["Oxidative phosphorylation"]][ keggPath_pathlist[["Oxidative phosphorylation"]] %in%  allPodoGenes ]
  
  ### transcript lvls 
  {
    GlyOxph.pathMean <- lapply( seq(listSCSN.1K.sampl), function(ii)
    {
      print( ii )
     
      datt <-  Seurat::ScaleData( listSCSN.1K.sampl[[ii]], features = c( ggenes1, ggenes2))
      
      datt1 <- cbind( rowMeans( datt@assays$RNA@data[ ggenes1 ,  datt$gtypeDE=="control"] ), 
                      rowMeans(  datt@assays$RNA@data[ ggenes1 ,  datt$gtypeDE!="control"] ) )
      
      datt2 <- cbind( rowMeans( datt@assays$RNA@data[ ggenes2 ,  datt$gtypeDE=="control"] ),
                      rowMeans(  datt@assays$RNA@data[ ggenes2 ,  datt$gtypeDE!="control"] ) )
      
      ll <- list( datt1 , datt2)
      names(ll) <- c("glyc","oxphos")
      return(ll)
    })
    
    Gly.pathMean <- apply(simplify2array(lapply( GlyOxph.pathMean, "[[",1) ), 
                          c(1, 2), FUN = mean, na.rm = TRUE)
    
    colnames(Gly.pathMean) <-  paste( c("ctrl","xprmnt"),"Glyc" ,  sep = "_")
    
    Oxph.pathMean <-apply(simplify2array(lapply( GlyOxph.pathMean, "[[",2) ), 
                          c(1, 2), FUN = mean, na.rm = TRUE)
    colnames(Oxph.pathMean) <-  paste(c("ctrl","xprmnt"),  "Oxphos" ,  sep = "_")
    
    
    gg1 <- ComplexHeatmap::Heatmap(( Gly.pathMean), col=viridis(10), 
                                   cluster_columns  = F, name = "Glycolisis KEGG")
    gg2 <- ComplexHeatmap::Heatmap(( Oxph.pathMean), col=viridis(10), 
                                   cluster_columns = F, name="Oxphos KEGG")
    
    pdf( width = 4, height = 24, file = "XX.pdf")
    plot_grid(grid.grabExpr(draw(gg1)), grid.grabExpr(draw(gg2)), nrow = 2, rel_heights = c(1,3))
    dev.off()
  }
 
   
   
   ### protein expression of markers
   {
     # convert gene to protID
     mart_mouse <- useMart( "ensembl",dataset="mmusculus_gene_ensembl" )
     gene2prot <- biomaRt::getBM( attributes=c( 'external_gene_name',"uniprotsptrembl","uniprotswissprot"),  mart = mart_mouse)
     gene1.prot <-  c( toupper( ggenes1 ), 
                       gene2prot$uniprotsptrembl[ ( gene2prot$external_gene_name ) %in%  ggenes1 ] , 
                       gene2prot$uniprotswissprot[ ( gene2prot$external_gene_name ) %in%  ggenes1 ] )
     gene1.prot <- gene1.prot[gene1.prot!=""]
     gene2.prot <-  c( toupper( ggenes2 ), 
                       gene2prot$uniprotsptrembl[ ( gene2prot$external_gene_name ) %in%  ggenes2 ] , 
                       gene2prot$uniprotswissprot[ ( gene2prot$external_gene_name ) %in%  ggenes2 ] )
     gene2.prot <- gene2.prot[gene2.prot!=""]
     
     library(DEP)
     ll <- list.files( pattern = "proteinGroups", full.names = T , 
                       path = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_charachterise/ProteinExpr/data/")
     
     # read
     podoGlom_proteome.raw <- lapply( ll , read.csv, header = TRUE,
                                      stringsAsFactors = FALSE, sep = "\t")
     # add gene names to first dataset
     X <- (strsplit(podoGlom_proteome.raw[[1]]$Majority.protein.IDs, split = "\\||;"))
     X1 <- unlist( lapply(X, function(x) gsub( "_MOUSE", "",paste( unique(x[grep("MOUSE",x)]), collapse = ";") )))
     X2 <- unlist( lapply(X, function(x) gsub( ".*CON__", "",paste( unique(x[grep("P\\d.*|Q\\d.*",x)]), collapse = ";") )))
     podoGlom_proteome.raw[[1]]$Gene.names <- ifelse( X1=="", X2 , X1)  
     
     # read the data, make unique protein names and filter out  contaminants
     podoGlom_proteome <- lapply( podoGlom_proteome.raw , function(proteome) {
       # remove decoy matches and matches to contaminant
       proteome = proteome[!proteome$Reverse=="+",]
       proteome = proteome[!proteome$Potential.contaminant=="+",]
       
       # Make unique names using the annotation in the "Gene.names" column as primary names 
       # and the annotation in "Protein.IDs" as name for those that do not have an gene name.
       proteome_unique <- DEP::make_unique( proteome , names="Gene.names", 
                                            ids="Protein.IDs", delim = ";")
       proteome_unique$name <- toupper( proteome_unique$name )
       return(proteome_unique)
       
     })
     
     # Generate a SummarizedExperiment object using an experimental design
     ll.meta <- list.files( pattern = "experimental_design", full.names = T , 
                            path = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_charachterise/ProteinExpr/data/")
     
     # creating se objects
     podoGlom_proteome_se <- lapply(seq(podoGlom_proteome), function( ii ){
       proteome <- podoGlom_proteome[[ii]]
       print(ii)
       MRpodo_LFQ_columns <- grep("LFQ.", colnames( proteome )) # get LFQ column numbers
       experimental_design <- read.table(header=T, sep = "\t",ll.meta[[ii]])
       MRpodo_proteome_se <- make_se( proteome , MRpodo_LFQ_columns, experimental_design )
       return(MRpodo_proteome_se)
     })
     # filtering
     podoGlom_proteome_filt <- lapply(seq(podoGlom_proteome_se), function( ii ){
       print(ii)
       proteome <- podoGlom_proteome_se[[ii]]
       experimental_design <- read.table(header=T, sep = "\t",ll.meta[[ii]])
       
       # Filter for proteins that are identified in all replicates of at least one condition
       proteome <- filter_missval(proteome, thr = round(min(summary( as.factor(experimental_design$condition)))*0.49))
       
       # Plot a barplot of the protein identification overlap between samples
       print( plot_coverage(proteome) )
       
       # Normalize the data
       proteome <- normalize_vsn(proteome)
       meanSdPlot( proteome )
       
       # # Visualize normalization by boxplots for all samples before and after normalization
       # plot_normalization(Nphs2_proteome_filt, Nphs2_proteome_Norm)
       
       # imputation
       plot_detect( proteome )
       proteome  <- impute( proteome, fun = "MLE")
       
     })

     

     
     ### custom jitterplot
     ggl <- lapply( seq(podoGlom_proteome_filt) , function( ii )
      {
       
       ### prepare the data
       XX <- podoGlom_proteome_filt[[ii]]@assays@data[[1]]
       XX <- XX[  rownames(XX) %in% c( gene1.prot, gene2.prot) , ]
       # remove genes with less than 3 measurments 
       XX <- as.data.frame(  XX[rowSums( !is.na(XX) )>=3, ] )
       
       XX$name <- rownames(XX)
       
       XX.m <- reshape2::melt( XX )
       ## add metadata
       XX.m$condition <- sub( "_.*", "", XX.m$variable)
       XX.m$pathN <- ifelse(XX.m$name %in%  gene1.prot , "Glycolisis","Oxphos" )
       
       # plot 
       gg <- ggplot2::ggplot( data = XX.m ,  
                              aes( x = value , 
                                   y = reorder(name, value, FUN = mean, na.rm=T) ) ) + 
         geom_jitter(shape=16, position=position_jitter(width =0.1, height = 0.1),
                     size = 5, aes( color=condition)) + 
         facet_grid(rows = vars(pathN) ,scales = "free_y", space = "free_y") + 
         scale_colour_colorblind(name = "proteomics\ndataset") +
         theme_bw() + labs(title = sub("_.*","",basename(ll)[ii]) ) +
         labs(y = "Gene symbol", x = "LFQ itensity") +
         theme(text=element_text(size=20 )) +
         guides(color = guide_legend(override.aes = list(size=5)))
       return(gg)
       # save plot
       # pdf(width = 10, height =  nrow( XX)*0.4, file = paste( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/ProteinValidation/figures/",
       #                                                        sub("_.*","",basename(ll)[ii]),"_PDS.50_jitterplotNOTscaled.pdf", sep = ""))
       # print(gg)
       # dev.off()
       # png(width = 800, height =  nrow( XX)*40, file = paste( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/ProteinValidation/figures/",
       #                                                        sub("_.*","",basename(ll)[ii]),"_PDS.50_jitterplotNOTscaled.png", sep = ""))
       # print(gg)
       # dev.off()
     })
     
     pdf( height = 24, width = 20, file = "XX.pdf")
   cowplot::plot_grid( plotlist = ggl[4:5])
   dev.off()
   
  
}
}

### Suppl.Fig. granger causality
  {
  lmtest::grangertest
  iids <- seq(listSCSN.1K.sampl)
  
  Grangertest_list <-  lapply( iids, function(ii)
    {
    
    print(ii)
    dataSet <- listSCSN.1K.sampl[[ii]]

    
    # # downsample
    # Idents(datt) <- "group"
    # datt <- subset( datt, downsample=1000 )
    # Idents(datt) <- "gtypeDE"
    # datt <- subset( datt, downsample=1000 )
    # datt <- subset( datt, gtypeDE=="experimental" )
    

    # order according to PDS
    expPDSorder <- colnames( dataSet )[ order( dataSet$PDS) ]
    

    # extract PA info, select pathways
    # order PA by PDS
    bothDB <-  rbind( listSCSN.1K.sampl_PAaucell[[ ii ]]$cells_AUC.kegg[ , expPDSorder ],
                      listSCSN.1K.sampl_PAaucell[[ ii ]]$cells_AUC.react[ , expPDSorder ] )
    bothDB <- bothDB[ match( sub(" \\*KEGG| \\*REACT","",pthsPlt), 
                             rownames( bothDB) ),  ]
    rownames( bothDB ) <-  pthsPlt
    
    
    datTOplot <- bothDB
    
    # ## smooth pathway activity
    # path_smooth <- apply( datTOplot , 1 , function(x) slideFunct(x, ncol(datTOplot)/WsizeF,
    #                                                              ncol(datTOplot)/WstepF ))
    
   #  # calculate LFC of the path change
   # LFC_vec <- path_smooth[1,]/path_smooth[ nrow(path_smooth),]
   # # swap the direction of curve to make all pathways changing in the same direction
   # path_smooth[ , LFC_vec >1 & !is.na(LFC_vec)] <- path_smooth[, LFC_vec >1 & !is.na(LFC_vec) ]*-1
   # 
    ll <- sapply( 1:(nrow(bothDB)), function(jj){
       print(jj)
      sapply( 1:(nrow(bothDB)) , function( ff){
        # XX <- tryCatch( lmtest::grangertest( x=path_smooth[ , jj], y=path_smooth[ , ff] )$`Pr(>F)`[2], 
        XX <- tryCatch( lmtest::grangertest( x=bothDB[ jj, ], y=bothDB[ ff, ] )$`Pr(>F)`[2], 
                                          error=function(e) NA )
      })
    })
    
    rownames(ll) <- colnames(ll) <- rownames(datTOplot)
    return( ll )
    
  })
  
  
  names( Grangertest_list ) <- names( listSCSN.1K.sampl )[iids]
  
  Grangertest_combineP <- sapply( seq(Grangertest_list[[1]]), 
                                  function(ii){
  datt <- sapply(seq(Grangertest_list), function(jj) Grangertest_list[[jj]][ii])
  datt[ is.na(datt)] <- 1
   metap::sumlog( datt)$p
  })
  MM <- matrix(Grangertest_combineP, nrow = 26, ncol = 26 )
  colnames(MM) <- rownames(MM) <- colnames(Grangertest_list[[1]])
  
  
  MM_adj <- apply( MM, 2, p.adjust,"fdr")
  toPlot <-  -log10( MM_adj )
  hist(toPlot)
  toPlot[ toPlot< 90] <- 0
 
  ComplexHeatmap::Heatmap( (toPlot) , cluster_columns = F, cluster_rows = F, 
                           col=viridis(20) )
 
  library(igraph)
  ggraph <- reshape2::melt(toPlot)
  ggraph <- ggraph[ggraph$value!=0,]
 
  pdf(width = 5,height = 5, file="selected_path_Grngr.graph.pdf")
    plot( igraph::graph_from_edgelist(as.matrix(ggraph[,1:2])),edge.arrow.size=0.5 )
  dev.off()
  
  ### calculate direction
  LFCs <-  sapply( iids, function(ii)
    {
    
    print(ii)
    dataSet <- listSCSN.1K.sampl[[ii]]
    
    
    # # downsample
    # Idents(datt) <- "group"
    # datt <- subset( datt, downsample=1000 )
    # Idents(datt) <- "gtypeDE"
    # datt <- subset( datt, downsample=1000 )
    # datt <- subset( datt, gtypeDE=="experimental" )
    
    
    # order according to PDS
    expPDSorder <- colnames( dataSet )[ order( dataSet$PDS) ]
    
    
    # extract PA info, select pathways
    # order PA by PDS
    bothDB <-  rbind( listSCSN.1K.sampl_PAaucell[[ ii ]]$cells_AUC.kegg[ , expPDSorder ],
                      listSCSN.1K.sampl_PAaucell[[ ii ]]$cells_AUC.react[ , expPDSorder ] )
    bothDB <- bothDB[ match( sub(" \\*KEGG| \\*REACT","",pthsPlt), 
                             rownames( bothDB) ),  ]
    rownames( bothDB ) <-  pthsPlt
    
    
    datTOplot <- bothDB
    
    ## smooth pathway activity
    path_smooth <- apply( datTOplot , 1 , function(x) slideFunct(x, ncol(datTOplot)/WsizeF,
                                                                 ncol(datTOplot)/WstepF ))

     # calculate LFC of the path change
    LFC_vec <- log10(path_smooth[1,]/path_smooth[ nrow(path_smooth),])
  })
  
  rowMeans(LFCs)
}



#### 2. show changes in individual pathway(s) across all studies #### 
  {
  
  # pathway activity
  sepLath <- c( "Neutrophil degranulation" ,
                "Ephrin signaling",
                "Focal adhesion" ,
                "Oxidative phosphorylation",
                "Glycolisis")
  
 
  
  
  ggl<-  lapply(seq(sepLath) , function(pp){
    print(sepLath[pp])
    PAdotHeatPlot( cells="experimental",
                   single_path=sepLath[pp],
                   dataSetIndex=c(1:3,5:9))
   })
  
  pdf( height = 8 , width = 20, file="Reactome_singlePath_heatDots.pdf")
  cowplot::plot_grid( plotlist = ggl, nrow = 2)
  dev.off()
  }




#### 3. show corr. of individual genes with PDS on a pathway diagram #### 
### combine corr tables with networks and TF info for cytoscape 
{ 
  ### data to plot on Net
  toPlot.all <- genecorrPDS.Sprmn_r.cntrd
  # combine all controls and all exprmnt in indiv. columns
  toPlot.all <- cbind(  ctrl_mean=rowMeans(toPlot.all[,1:7], na.rm = T), 
                        toPlot.all[,8:14])
  
  colnames(toPlot.all) <- c( "ctrl_mean","Nphs2","Wt1","Pdss2" ,      
                             "Btbr","Cd2ap","doxo","nephr.day5")
  
  ### load network information
  Net <- read.table(sep = "\t", header = T, file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/cytoscapeVIZ/Podo_areas/C.shell/FAnet.tsv")
  Net.nodes <- unique( union(Net$Symbol, Net$Interactors))
  Net.nodes <- Net.nodes[Net.nodes!="" & !is.na(Net.nodes)]
  PDScorr <- toPlot.all[ match( Net.nodes , rownames(toPlot.all) ) ,  ]
  rownames(PDScorr) <- Net.nodes
  
 
  # aggegate TF reg information in ove vector
  TFregVec <- Reduce( c, apply( TFreg, 1, function(X){
    stringr::str_c(colnames(TFreg)[ X==1 & !is.na(X)],collapse = ",")
  }))
  Net.nodeTab <- cbind.data.frame( gName=Net.nodes,
                                   DSgenes = ifelse( rownames(PDScorr) %in% DS_all$gene_symbol[1:42], "DS.42", 
                                                     ifelse( rownames(PDScorr) %in% DS_all$gene_symbol , 
                                                             "DS.rest", "")),
                                   Nstd.PDScor0.01.xprmnt =rowSums( genecorrPDS.Sprmn_qval[match( Net.nodes , 
                                                                                                  rownames(genecorrPDS.Sprmn_qval)),
                                                                                           8:14] < 0.01, na.rm = T ) ,
                                   gNameUP = toupper(Net.nodes) ,
                                   SCSNmeasured = Net.nodes%in%rownames(toPlot.all) ,
                                   PDScorr ,
                                   Wt1chipseq=Net.nodes %in% rownames(WT1chipTgenes)[WT1chipTgenes$x<0.3],
                                   TFregSum=TFregVec,
                                   TFreg)
  
  write.table( Net.nodeTab , sep = "\t", quote = F, col.names = NA, 
               file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/crosscorr/ECMnet.PDS.TFs3plus.tsv")
  
}


