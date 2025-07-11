options( connectionObserver = NULL )
.libPaths(c("/home/tim_nevelsk/R/x86_64-pc-linux-gnu-library/4.0", 
            "/media/tim_nevelsk/WD_tim/SOFT/R"))
library( GSEABase )
library( biomaRt )
library( ggplot2 )
library( ggthemes)

#### load Data
keggPath_pathlist <- readRDS( "/media/tim_nevelsk/WD_tim/ANNOTATIONS/Pathways/kegg_pathsList_gName.07.09.23.rda")
reactPath_pathlist <- readRDS( "/media/tim_nevelsk/WD_tim/ANNOTATIONS/Pathways/reactome_pathsList_gName.07.09.23.rda")
PodoPathGSet_pathlist <- readRDS( file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PodoPathGSet.05.01.24.7plus.rda")

# podocyte genes
allPodoGenes <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/SCSN_allPodoGenes.rda")


#### estimate pathway activities in KFO snRNAseq with AUCell ####

  # prepare reactorme and kegg gsets
  pathCollect__list <- lapply( list( keggPath_pathlist, 
                                     reactPath_pathlist, 
                                     PodoPathGSet_pathlist ) , function( pathlist ) {
                                       print( head( names(pathlist)))
                                       GSEABase::GeneSetCollection(
                                         lapply(seq(pathlist), function(jj){
                                           ppath <- unique(pathlist [[jj]])
                                           # if( length(ppath)>2 ){}
                                           GSEABase::GeneSet( ppath , 
                                                              setName= names(pathlist)[jj] )    
                                         }) )
                                     } )
  names(pathCollect__list) <- c("kegg", "react","podoPaths")
  
  
  
#### calculate AUCell score
  {
    
    
    listSCSN.1K.sampl_AUCell.ranks <- lapply( seq(listSCSN.1K.sampl) , 
                                              function(ii , datt=listSCSN.1K.sampl)
                                              {
                                                print(ii)
                                                
                                                # load expression data 
                                                exprMatrices <- datt[[ii]]@assays$RNA@data
                                                exprMatrices <- exprMatrices[ 
                                                  rowSums( round(exprMatrices) > 0 ) > ncol(exprMatrices)*0 , ]
                                                
                                                # 1. Build gene-expression rankings for each cell  
                                                AUCell::AUCell_buildRankings( exprMatrices,  nCores=1, plotStats=F)
                                              } )
    
    
    ## calculate for pathways
    listSCSN.1K.sampl_PAaucell <- lapply( seq( listSCSN.1K.sampl_AUCell.ranks ) , 
                                          function( ii )
                                          {
                                            
                                            cat( "experiment N ", names(listSCSN.1K.sampl)[ii] , "\n" , "*****************************", "\n" )
                                            
                                            # 1. Build gene-expression rankings for each cell  
                                            cells_rankings <-  listSCSN.1K.sampl_AUCell.ranks[[ii]]
                                            
                                            # 2. Calculate enrichment for the gene signatures (AUC)
                                            # store NAs if less than 20% of geneSet genes are expressed in snRNAseq
                                            cells_AUC.kegg <- AUCell::AUCell_calcAUC( pathCollect__list$kegg ,
                                                                                      cells_rankings ,
                                                                                      verbose = T)
                                            cells_AUC.kegg <- AUCell::getAUC( cells_AUC.kegg)
                                            
                                            
                                            cells_AUC.react <- AUCell::AUCell_calcAUC( pathCollect__list$react ,
                                                                                       cells_rankings ,
                                                                                       verbose = T)
                                            cells_AUC.react <- AUCell::getAUC( cells_AUC.react)
                                            
                                            cells_AUC.podoPth <- AUCell::AUCell_calcAUC( pathCollect__list$podoPaths , 
                                                                                         cells_rankings , 
                                                                                         verbose = T)
                                            cells_AUC.podoPth <- AUCell::getAUC( cells_AUC.podoPth)
                                            
                                            
                                            ll <- list( cells_AUC.kegg, cells_AUC.react ,cells_AUC.podoPth)
                                            
                                            names(ll) <- c("cells_AUC.kegg", "cells_AUC.react","cells_AUC.podoPaths")
                                            return( ll )
                                          } )
    names(listSCSN.1K.sampl_PAaucell) <- names(listSCSN.1K)
    saveRDS( listSCSN.1K.sampl_PAaucell , file= "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN.PDS_samples.1K_PAaucell.13.05.24.rda")
    
  }
  


#### select pathways to plot #### 

#### prepare Podo Paths and PDScorr gene.sets
## cluster from analysis of PDS correlations
toPlot <- genecorrPDS.Sprmn_freq[genecorrPDS.Sprmn_freq$total>=Nn, 1:4]
cclust <- as.factor( cutree( hclust( dist(as.matrix(toPlot), method = "euclidean"),
                                     method="complete"), k = 4) )
levels(cclust) <- hue_pal()(4)


toTestGSET <- lapply(seq( levels(cclust)), function(ii){
  names(cclust)[cclust==levels(cclust)[ii]]
}) 
names(toTestGSET) <- sapply(levels(cclust), plotrix::color.id)

## load podocyte gene sets
# PMID:36307401 tab 1 and 2
slitDiaphr <-  read.table(sep = "\t", header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/cytoscapeVIZ/Podo_areas/PPIsource/PMID:36307401_tab1tab2.slitDiaphr.csv")$GeneName
# Schell et al. PMID:33514561 tab 3 and 4
actCtsklt <-  read.table(sep = "\t", header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/cytoscapeVIZ/Podo_areas/PPIsource/PMID:33514561_SupplTab3.4_ActinCytoSkeleton.csv")$name
# Schell et al. PMID:33514561 fig1f
actCtsklt_fig1f <-read.table(sep = "\t", header = T, fill = T,  "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/cytoscapeVIZ/Podo_areas/C.shell/ActnCtskl_node.csv")
actCtsklt_fig1f <- unique( actCtsklt_fig1f$gene_symbol)
# PMID:28536193 , figure 1E in Schell et al. 2017
fcladhsn <-  read.table(sep = "\t", header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/cytoscapeVIZ/Podo_areas/PPIsource/PMID:28536193_SupplTab4_focalAdhesion.csv")$Gene.names...primary.

# adhesome
cnsrvAdhesome <-  read.table(sep = "\t", header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/cytoscapeVIZ/Podo_areas/PPIsource/PMID:32147508_SupplTab3.ConsAdhesome.csv")$Mouse.gene.name
fltrdAdhesome <- read.table(sep = "\t", header = T,"/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/cytoscapeVIZ/Podo_areas/PPIsource/PMID:33761352_SupplTab10.fltrdAdhesome.csv")$Mouse.gene.name
Adhesome <- union(cnsrvAdhesome , fltrdAdhesome)
# matrisome
Mtrxsome <- read.table(sep = "\t", fill = T , header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/cytoscapeVIZ/Podo_areas/PPIsource/PMID:33761352_SupplTab7.fltrdMtxsome.csv")$Gene_names
Mtrxsome <- unlist( strsplit(Mtrxsome, split = ";") )
Mtrxsome <- unique(unlist( fun_homoTO.FROMmouse(Mtrxsome)$MUS) )
# ECM from MArtin
ECM <- readxl::read_xlsx("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/cytoscapeVIZ/Podo_areas/PPIsource/Network_structure_ECM-updated.xlsx",
                         sheet = 1)
ECM <- unique( ECM$Symbol )
ECM <- ECM[ECM!=""]

# minimal gene set
PodoPathGSet <- c( toTestGSET , list( slitDiaphr, actCtsklt, actCtsklt_fig1f, 
                                      fcladhsn, ECM, Adhesome , Mtrxsome ) )
names(PodoPathGSet) <- c(names(toTestGSET),"Slit.Diaphr", "Actin.Ctsklt", "Actin.Ctsklt.fig1f",
                         "Focal.Adhesion", "ECM","Adhesome","Matrisome" )
PodoPathGSet <- lapply(PodoPathGSet, unique)

Reduce( cbind, lapply(seq(PodoPathGSet), function(ii){
  Reduce( c, lapply(seq(PodoPathGSet), function(jj) length( intersect(
    PodoPathGSet[[ii]],PodoPathGSet[[jj]] ) ) ) )
}) )
# saveRDS(PodoPathGSet , file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PodoPathGSet.29.09.23.rda")
saveRDS(PodoPathGSet , file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PodoPathGSet.05.01.24.7plus.rda")

## create gene set collection
PodoPathGSet_genesets <- lapply( seq(PodoPathGSet) , function( jj ) {
  GeneSet( PodoPathGSet[[jj]], 
           setName= names(PodoPathGSet)[jj] ) } )
PodoPathGSet_genesets <-  GSEABase::GeneSetCollection( PodoPathGSet_genesets )







### select top or common Paths
{
  
  # aggregate over samples, groups and genotypes
  PAaucell.react_vs_PDS <- Reduce( union, lapply( seq( listSCSN.1K.sampl), 
                                                  function(ii, N=10){
                                                    print(ii)
                                                    
                                                    # split by gtype
                                                    Reduce( union , lapply( list("control","experimental"), 
                                                                            function(ggtype){
                                                                              datt <- subset(  listSCSN.1K.sampl[[ii]] , 
                                                                                               gtypeDE==ggtype)
                                                                              ggroup<- as.character( unique( datt$group) )
                                                                              # split by group
                                                                              Reduce( intersect , lapply( seq(ggroup), 
                                                                                                          function(zz){
                                                                                                            # print(zz)
                                                                                                            datGroup <- subset(  datt , group==ggroup[zz] )
                                                                                                            ssample <- as.character( unique( datGroup$sample ) )
                                                                                                            # split by sample
                                                                                                            Reduce( union, lapply( seq(ssample), 
                                                                                                                                   function(vv){
                                                                                                                                     # print(vv)
                                                                                                                                     datTest <- subset(  datGroup , sample==ssample[vv] )
                                                                                                                                     dattPDS <- datTest$PDS
                                                                                                                                     dattPath <- listSCSN.PDS_1K_PAaucell.react[[ii]][ 
                                                                                                                                       , colnames(datTest)]
                                                                                                                                     PAvsPDScorr <- psych::corr.test( t( dattPath ) ,dattPDS , method="spearman") 
                                                                                                                                     ress <- PAvsPDScorr$r[ PAvsPDScorr$p[,1]<0.05 & 
                                                                                                                                                              !is.na(PAvsPDScorr$p[,1]), ]
                                                                                                                                     names( ress[ order(ress)][1:N] )
                                                                                                                                   }))
                                                                                                          }))
                                                                            })
                                                    )
                                                    
                                                    
                                                    
                                                  }))
  
  
  # aggregate over samples, groups and genotypes
  PAaucell.KEGG_vs_PDS <- Reduce( union, lapply( seq( listSCSN.1K.sampl), 
                                                 function(ii, N=10){
                                                   print(ii)
                                                   
                                                   # split by gtype
                                                   Reduce( union , lapply( list("control","experimental"), 
                                                                           function(ggtype){
                                                                             datt <- subset(  listSCSN.1K.sampl[[ii]] , 
                                                                                              gtypeDE==ggtype)
                                                                             ggroup<- as.character( unique( datt$group) )
                                                                             # split by group
                                                                             Reduce( union , lapply( seq(ggroup), 
                                                                                                     function(zz){
                                                                                                       # print(zz)
                                                                                                       datGroup <- subset(  datt , group==ggroup[zz] )
                                                                                                       ssample <- as.character( unique( datGroup$sample ) )
                                                                                                       # split by sample
                                                                                                       Reduce( intersect, lapply( seq(ssample), 
                                                                                                                                  function(vv){
                                                                                                                                    # print(vv)
                                                                                                                                    datTest <- subset(  datGroup , sample==ssample[vv] )
                                                                                                                                    dattPDS <- datTest$PDS
                                                                                                                                    dattPath <- listSCSN.PDS_1K_PAaucell.kegg[[ii]][ 
                                                                                                                                      , colnames(datTest)]
                                                                                                                                    PAvsPDScorr <- psych::corr.test( t( dattPath ) ,dattPDS , method="spearman") 
                                                                                                                                    ress <- PAvsPDScorr$r[ PAvsPDScorr$p[,1]<0.05 & 
                                                                                                                                                             !is.na(PAvsPDScorr$p[,1]), ]
                                                                                                                                    names( ress[ order(ress)][1:N] )
                                                                                                                                  }))
                                                                                                     }))
                                                                           })
                                                   )
                                                   
                                                   
                                                   
                                                 }))
  
  
}


### histogram of pathway sizes
hist(sapply(wikiPath_gName, function(x) log10(length(x))), 
     main = "wikiPath pathway size distribution", 
     xlab = "log10 pathway size") 
hist(sapply(keggPath_gName, function(x) log10(length(x))), 
     main = "KEGG pathway size distribution", 
     xlab = "log10 pathway size") 
hist(sapply(reactPath_gName, function(x) log10(length(x))), 
     main = "Reactome pathway size distribution", 
     xlab = "log10 pathway size") 

#### PPS https://github.com/paulvanderlaken/ppsr 


#### compute and plot smoothed signal ####

### load AUCell score
# listSCSN.1K.sampl_PAaucell <- readRDS( file= "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PathwayActivity/listSCSN.PDS_1K_PAaucell.29.09.2023rda")
listSCSN.1K.sampl_PAaucell <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN.PDS_samples.1K_PAaucell.13.05.24.rda")
# separate pathways by database pathways
listSCSN.PDS_1K_PAaucell.kegg <- lapply(listSCSN.1K.sampl_PAaucell, "[[", 1)
listSCSN.PDS_1K_PAaucell.react <- lapply(listSCSN.1K.sampl_PAaucell, "[[", 2)
listSCSN.PDS_1K_PAaucell.podoPath <- lapply(listSCSN.1K.sampl_PAaucell, "[[", 3)

# load PDS
listSCSN.1K.sampl <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN_samples.1K.22.12.23.rda")
names(listSCSN.1K.sampl) <- names(listSCSN.1K.sampl_PAaucell)
# set smoothing factors
WsizeF=10
WstepF=30

### load results of PDS correlation with PA
PAaucell.react_vs_PDS <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PathwayActivity/Multi.pathways/PAaucell.react_vs_PDSrank.gtypeAgg.rda")
PAaucell.kegg_vs_PDS <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PathwayActivity/Multi.pathways/PAaucell.kegg_vs_PDSrank.gtypeAgg.rda")

### a function to plot pathway fingerprint
fingerprintPlot <-function( PAtoPlot , 
                            collAnnot= gtypeFr_smooth[[ii]], 
                            podoPaths=F ,
                            rowNames = T, 
                            rowOrder= rownames(MMpval)[rev(gsea.order$rowInd)] )
{
  print(ii)
  require(ComplexHeatmap)
  
  
  ### prepare annotation rows that show percentage 
  ### of cells of a certain gtype in a window
  annot <- (  collAnnot )
  annot <- annot[ ,!is.na(colSums(annot)), drop=F ]
  if( min(annot) == max(annot)) { minn =0 
  } else  minn = min(annot)
  col_fun <- circlize::colorRamp2( 
    breaks = c( minn ,  max(annot)) ,
    colors=c( "white", "black"))                                     
  
  ### annotation
  ha = HeatmapAnnotation( foo=annot ,
                          col = list(foo = col_fun),
                          simple_anno_size = unit(1.5/ncol(annot), "cm"),
                          show_legend =  FALSE
                          # show_heatmap_legend==F
                          
  )
  
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
PAdotHeatPlot <- function(cells="all", 
                          single_path=sepLath[[ii]],
                          dataSetIndex = c(1:3,5:9) ){
  
  ## cells - a charachter vector which specifies whether to use all cells,
  # when set to "all", or a specific genotype, for the analysis.
  
  require(scales)
  require(colorspace)
  
  print( sepLath[[ii]])
  
  
  # path to plot
  toPlot <- Reduce( rbind, lapply( dataSetIndex , function( ii ) 
  {
    
    # print(ii)
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
      if( length(unique(dataSet$gtypeDE))>1){
        dataSet <- subset( dataSet, downsample=1000 ) } else {
          ## add controls to expression data
          dataSet <- merge( dataSet,   subset( listSCSN.1K.sampl[[8]], gtypeDE=="control") )
          dataSet <- subset( dataSet, downsample=1000 )
          
          ## add controls to PA data
          if( single_path %in% names(reactPath_pathlist) ) {
            XX <- listSCSN.PDS_1K_PAaucell.react[[8]]
          } else XX <- listSCSN.PDS_1K_PAaucell.kegg[[8]]
          datt <- cbind( datt, XX[ match(rownames(datt) , rownames(XX)),] )
        }
      
      ### normalise PA by control means
      toPlot <- sweep(t( datt[ rownames(datt) == single_path, colnames(dataSet)] ), 
                      2, rowMeans(datt[ rownames( datt ) == single_path,
                                        colnames(dataSet)[dataSet$gtypeDE=="control"]]), 
                      FUN = '/')            
      toPlot <- log2(toPlot) #  get log2FC of PA
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
    toPlot <- toPlot[ order(toPlot[,"PDS"]),] # order by PDS
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
    nnames <- c("controls", names(listSCSN.1K.sampl ) )
  } else nnames <- names(listSCSN.1K.sampl)
  
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
  if(  cells=="all" ) {
    gg <- gg + scale_color_continuous_divergingx(palette = 'RdBu',
                                                 mid = 0.0 , rev = T,
                                                 p1=0.8, p2=0.8, p3=0.8, p4=0.8 ) 
  } else  gg <- gg + scale_color_viridis() 
  return(gg)
}


### KFO studies, annotated by PDS and AlbCr
{
  
  ### select top or common Paths
  {
    
    N <- 10
    PAaucell.react_topN.KFO <- Reduce( union, apply(
      PAaucell.react_vs_PDS[,1:6], 2, function(X) names(X[order(X)])[1:N]))
    PAaucell.kegg_topN.KFO <- Reduce( union, apply(
      PAaucell.kegg_vs_PDS[,1:6], 2, function(X) names(X[order(X)])[1:N]))
    
    
    
  }
  
  ### cluster pathways
  {
    
    # table of similarity
    library(GeneOverlap)
    library(dendextend)
    
    RR <- reactPath_pathlist[ paste0( names(reactPath_pathlist),  " *REACT") %in% pthsPlt$V1]
    names(RR) <- paste( names(RR), "*REACT")
    KK <-  keggPath_pathlist[paste0( names(keggPath_pathlist), " *KEGG") %in% pthsPlt$V1 ]
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
  
  #### plot KFO samples, annotate by albuninuria
  {
    ### compute number of cells within low/mid/high PDS per window
    KFOannot <- read.table( sep="\t", header = T,
                            "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/Sample_Names_KFO.csv")
    
    samplTag <- c("BDP","AKK","BWV|CDC")
    
    # make plots for 3 KFO datasets
    gglist <-  lapply( 1:3, function(ii,
                                     pathType= "podoPaths"){
      
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
      datt <- subset( datt, downsample=1000 )
      
      datt$AlbCr <- Albuminuria$AlbCrRatio_v2[ match( datt$sample, Albuminuria$CCG_Sample_ID)]
      
      # order according to PDS
      expPDSorder <- colnames( datt )[ order( datt$PDS) ]
      
      # create column annotations
      Alb.binFr <- datt@meta.data[  expPDSorder , c("AlbCr","PDS")]
      Alb.binFr_smooth <- apply( Alb.binFr, 2, function(XX) slideFunct(XX, nrow(Alb.binFr)/WsizeF, 
                                                                       nrow(Alb.binFr)/WstepF ))
      Alb.binFr_smooth <- scale(Alb.binFr_smooth)
      
      
      # extract PA info, select pathways
      # order PA by PDS
      bothDB <-  rbind( listSCSN.1K.sampl_PAaucell[[ii]]$cells_AUC.kegg[ , expPDSorder ],
                        listSCSN.1K.sampl_PAaucell[[ii]]$cells_AUC.react[ , expPDSorder ])
      bothDB <- bothDB[ match( sub(" \\*KEGG| \\*REACT","",pthsPlt$V1), 
                               rownames( bothDB) ),  ]
      rownames( bothDB ) <-  pthsPlt$V1
      
      
      datTOplot <- bothDB
      
      ## smooth pathway activity
      path_smooth <- apply( datTOplot, 1 , function(x) slideFunct(x, ncol(datTOplot)/WsizeF, 
                                                                  ncol(datTOplot)/WstepF ))
      
      
      ### plot
      gg <-  fingerprintPlot( PAtoPlot = path_smooth , 
                              collAnnot= Alb.binFr_smooth ,
                              # podoPaths= T ,
                              rowNames = F,
                              rowOrder= rownames(MMpval)[rev(gsea.order$rowInd)]
      )
      return(gg)
    })
    
    
    
    
    pdf(height = 8, width = 18, "podoPaths.fingerprints_KFO_both.select.pdf")
    
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nr = 1, nc = 3)))
    
    lapply( seq(gglist), function(ii){
      pushViewport(viewport(layout.pos.row = 1,
                            layout.pos.col = c(1:3)[ii]))
      draw( gglist[[ii]], newpage = FALSE, column_title= names(listSCSN.1K.sampl)[ii])
      upViewport()
    })
    
    dev.off()
    
  }
  
}

### all studies, annotated by PDS only
{
  
  ### select top or common Paths
  {
    
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

### one gene in all studies
{
  
  # ppath <- "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/PA/single_path/KEGG/"
  # pNames <- names(keggPath_pathlist)
  
  ppath <- "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/PA/single_path/REACT/"
  pNames <- names(reactPath_pathlist)
  
  
  lapply( seq(pNames) , function(pp){
    # print(pNames[pp])
    gg <- tryCatch( PAdotHeatPlot(cells="experimental",
                                  single_path=pNames[pp],
                                  dataSetIndex=c(1:3,5:9)), error = function(e) NA)
    if(!is.na(gg)) {
      nnames <- gsub(" |-|/|:|>|<","_",pNames[pp])
      pdf( height = 4 , width = 10, file= paste0(ppath,nnames,".pdf"))
      print(gg)
      dev.off()
    }
  })
}



#### granger causality for selected pathways ####
  
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
  


