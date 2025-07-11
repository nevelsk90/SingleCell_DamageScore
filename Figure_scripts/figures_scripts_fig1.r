###======# Single-Cell Resolution of Cellular Damage Illuminates Disease Progression #=======###
###======# publication figure 1 PDS code #=======###

# release memory
mallinfo::malloc.trim()

#### load DS and code ####
.libPaths(c("/home/tim_nevelsk/R/x86_64-pc-linux-gnu-library/4.0", "/media/tim_nevelsk/WD_tim/SOFT/R"))
options( connectionObserver = NULL )
library( org.Mm.eg.db )
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
library(rstatix)


# setwd
setwd( "/home/tim_nevelsk/PROJECTS/PODOCYTE/WRITING/PDS_manuscript/Figures/Figure1" )
inputdir <-  "/media/tim_nevelsk/WD_tim/PROJECTS/WRITING/PDS_manuscript/Figure_input" 

# load necessary code
source("https://raw.githubusercontent.com/PauUng/HepatocyteDamageScore/master/SharedFunctions.R")
source("https://raw.githubusercontent.com/nevelsk90/R_scripts/master/usefulRfunc.r")  


# load expression data
listSCSN.1K <- readRDS( paste0( inputdir, "/listSCSN_1K.22.12.23.rda") )


#### 1. panel A. Scheme of the damage score pipeline, made in inkscape ####

#### 2. panel B. Dim.Red. plots showing PDS gradient in the population of damaged podocytes ####
### UMAP of podocyte colored by PDS, from Wt1het.del snRNAseq study
  {
    toPlot <- listSCSN.1K$Wt1
    
   gg1 <- FeaturePlot( toPlot , features = "PDS", 
                        cells = sample(colnames(toPlot)),
                        cols = brewer.pal( n = 9 , name = "YlOrBr"),
                       pt.size =1.5)+ 
      theme(legend.position = "bottom", 
            legend.key.height=unit(5,"mm"),
            text = element_text(size=20),
            axis.title=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank()) 
      # coord_cartesian(xlim = c(-10,-6), ylim = c(-2,3))+
      # theme(legend.position = "none")
    
    gg2 <- DimPlot( toPlot , group.by = "gtype",
                   pt.size = 1.5 , shuffle = T)+
      # coord_cartesian(xlim = c(-10,-6), ylim = c(-2,3))+
      scale_colour_colorblind()+
      theme(legend.position = "bottom",
            legend.key.height=unit(15,"mm"),
            text = element_text(size=20),
            axis.title=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank())
   
    
    ## densito plot for PDS
    gg3 <- ggplot(data=toPlot@meta.data,
                  aes( x=PDS, color=gtype))+
      geom_density(lwd=3)+ theme_bw()+
      scale_colour_colorblind()+
      theme(legend.position = "none",
            axis.title.y =element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            text = element_text(size=20))
    
    ggpnl <-cowplot::plot_grid(plotlist = list(gg2, gg1, gg3), 
                               rel_heights = c(1,0.25), nrow = 2)
    
    
    pdf(height = 6, width = 8, file = "Wt1hd.snRNAseq_podo_PDS.umap.NOdim.pdf")
    ggpnl
    dev.off()
    png(height = 600, width = 600, file = "Wt1hd.snRNAseq_podo_PDS.umap.NOdim.png")
    ggpnl
    dev.off() 
  }


#### 3. panel C. model cross validation plots ####

CV_model.list_PDS42.agg <- readRDS( file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/Dev.Valid/CV_model.list_PDS42.sampAgg.rda")

# scale PDS
toPlot.list <- lapply(seq(CV_model.list_PDS42.agg), function(ii){
  toPlot <- CV_model.list_PDS42.agg[[ii]]
  Reduce( rbind, lapply( seq(toPlot), function(jj){
    datt <- toPlot[[jj]]
    datt$PDS.42 <- scale(datt$PDS.42)
    return(datt)
  }))
}) 
names(toPlot.list) <- names(CV_model.list_PDS42.agg)

### make one plot for all models
  {

  toPlot <- Reduce( rbind, lapply(seq(toPlot.list), 
                                  function(ii){
    datt <- toPlot.list[[ii]]
    datt <- aggregate( PDS.42 ~ gtypeDE + dataset + CV_model, data=datt , FUN=mean )
    names(datt)[4]<- "PDS.42"
    return(datt)
  }))
  
  gg<-  ggplot2::ggplot( toPlot , 
                         aes( y= PDS.42, x=gtypeDE, color=gtypeDE)) +
    scale_color_colorblind() +  geom_boxplot(lwd=1.5) + 
    geom_jitter(size=5)+
    theme_bw() + theme( text = element_text(size = 24) , axis.title.x=element_blank(),
                        axis.text.x = element_text(angle = 45, hjust = 1)) +
    stat_compare_means( size = 6   ) +
    facet_grid(cols = vars(CV_model),scales = "free")
  
  pdf( height = 6, width = 12, file = "Figure1/PDS_cv.models_boxplot_agg.pdf" )
  print( gg )
  dev.off()
  
  # png( height = 600, width = 1000, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/WRITING/PDS_manuscript/Figures/Figure1/PDS_cv.models_boxplot_agg.png" )
  # print( gg )
  # dev.off()
  }

### Suppl.Fig 2B plot individual models
  {
  ppglist <- lapply(seq(toPlot.list), function(ii)
  {
    datt <- toPlot.list[[ii]]
    ggplot2::ggplot( datt , 
                     aes( y= PDS.42, x=gtypeDE, color=gtypeDE)) +
      scale_color_colorblind() +  geom_boxplot(lwd=1.5) + 
      geom_jitter(size=3)+
      theme_bw() + theme( text = element_text(size = 24) , axis.title.x=element_blank(),
                          axis.text.x = element_text(angle = 45, hjust = 1)) +
      # stat_compare_means( size = 6   )+
      facet_grid(cols = vars(dataset),scales = "free")
  })
  ppglist <- cowplot::plot_grid( plotlist = ppglist, ncol = 1)
  
  pdf( height = 20, width = 20, file = "Supl.Fig1/PDS_cv.models_boxplot_Indv.pdf" )
  print( ppglist )
  dev.off()
  png( height = 1800, width = 1500, file = "Supl.Fig1/PDS_cv.models_boxplot_Indv.png" )
  print( ppglist )
  dev.off()
}

#### Signature size test ####
{
  require(GSEABase)
  require( AUCell)
  require( ggplot2)
  require( reshape2 )
  
  size_vec <- c(  1:nrow(DSignature) )
  size_Gsets <- Reduce( c , lapply(size_vec, function(ssize){
    
    genesUP = DSignature$gene_symbol[1:ssize][ 
      which( DSignature$direction_foldchange[1:ssize]==1)]
    genesDOWN =  DSignature$gene_symbol[1:ssize][ 
      which( DSignature$direction_foldchange[1:ssize]== -1)]
    
    list( GeneSet( genesUP , 
                   setName= paste0("UP.",ssize )) , 
          GeneSet( genesDOWN , 
                   setName= paste0("DOWN.", ssize)))
  }) )
  size_Gsets <- GSEABase::GeneSetCollection( size_Gsets  )
  
  ### calculate PDS for a list of studies (listSCSN.1K.sampl) 
  # using a range of damage signature sizes (size_vec)
  {
    ceilThrsh <- 0.05
    
    ## inject controls in doxo and NTS.d5
    dattM <- listSCSN.1K.sampl
    dattM$doxo <- merge( dattM$doxo, 
                         subset(dattM$nephr.D1, subset=gtypeDE=="control"))
    dattM$nephr.D5 <- merge( dattM$nephr.D5, 
                             subset(dattM$nephr.D1, subset=gtypeDE=="control"))
    
    # iterate over sc studies and calculate PDS
    
    set.seed(42)
    listSCSN.PDS_DSsizeTest <- lapply( seq(dattM),
                                       function(ii)
                                       {
                                         ### uppercase row. names
                                         print(ii)
                                         newSeu <- dattM[[ii]]
                                         # subset 
                                         Idents(newSeu) <- newSeu$sample
                                         
                                         # get counts
                                         expr <- newSeu@assays$RNA@counts
                                         expr <- expr[ rowSums(round(expr) > 0) > 0.01*ncol(expr) , ]
                                         
                                         ###  1. Build gene-expression rankings for each cell  
                                         cellRanks <- AUCell_buildRankings( expr, nCores=4, plotStats=F, 
                                                                            verbose = T )
                                         
                                         ###  2. Calculate enrichment for the gene signatures (AUC)
                                         # gene sets should include UP and DOWN regulated genes seperately
                                         cells_AUC <- AUCell_calcAUC( size_Gsets, cellRanks , verbose = T,
                                                                      # aucMaxRank parameter is a tricky one,
                                                                      # for single cell 0.05 is a reasonable 
                                                                      # but for bulk a bigger, e.g. 0.2 can be used
                                                                      aucMaxRank = ceiling( ceilThrsh * nrow(cellRanks))
                                         )
                                         cells_AUC <- getAUC( cells_AUC) #this will extract results 
                                         
                                         ## make sure even empty gSets are in the matrix
                                         cells_AUC <- cells_AUC[ match(
                                           names(size_Gsets), rownames(cells_AUC)),]
                                         cells_AUC[is.na(cells_AUC)] <- 0
                                         
                                         # calculate coefficients to balance by N of UP and DOWN
                                         coeff <- sapply( seq(size_Gsets), function(jj) length(
                                           size_Gsets[[jj]]@geneIds)/rep(size_vec, each=2)[jj])
                                         ## adjust 
                                         cells_AUC.wghtd <- cells_AUC*coeff
                                         
                                         # subtract DOWN from UP AUCell scores
                                         DSsizeTest = cells_AUC.wghtd[ seq( 1,length(size_Gsets), by=2) ,] - 
                                           cells_AUC.wghtd[seq( 1,length(size_Gsets), by=2)+1,]
                                         rownames(DSsizeTest) <- paste0("PDS.size",size_vec)
                                         
                                         # add PDS to the metadata
                                         newSeu@meta.data <- cbind(newSeu@meta.data, t(DSsizeTest))
                                         return(newSeu)
                                       })
    names(listSCSN.PDS_DSsizeTest) <- names(listSCSN.1K.sampl)
    # saveRDS(listSCSN.PDS_DSsizeTest, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/Dev.Valid/listSCSN.PDS42_DSsizeTest.22.05.25.rda")
    
  }
  
  ### plot PDS distributions
  {
    # listSCSN.PDS_DSsizeTest<- readRDS(file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/Dev.Valid/listSCSN.PDS42_DSsizeTest.22.05.25.rda")
    
    ## prepare the table
    toPlot <- Reduce( rbind,  lapply( seq(listSCSN.PDS_DSsizeTest) , function(ii){
      datt <- listSCSN.PDS_DSsizeTest[[ii]]@meta.data
      # datt <-   newSeu@meta.data
      datt <- datt[,c( "gtypeDE","gtype","group","sample", 
                       grep("PDS",names(datt), value = T)) ]
      
      datt<- datt[ , ! colnames(datt) %in% 
                     c("PDS.42.005",  "PDS.42.005.2" ,"PDS") ]
      datt[,grep("PDS.*",names(datt), value = T)] <- scale(
        datt[,grep("PDS.*",names(datt), value = T)], scale = T )
      datt$dataSet <- names( listSCSN.1K.sampl)[ii]
      print( dim(datt))
      return(datt)
    }) )
    
    ### plot mean of the differences between means of ctrl and exprmntal
    toPlot4 <- Reduce( rbind,  lapply(size_vec, function(ii){
      data.frame( mean.of.means = mean( toPlot3$PDS.diff.Mean[toPlot3$size==ii] ) , 
                  var.of.means = var(toPlot3$PDS.diff.Mean[toPlot3$size==ii] ))
    }))
    toPlot4$size <- size_vec 
    toPlot4.m <- reshape2::melt(toPlot4, id.var="size")
    # toPlot4.m.sel <- toPlot4.m[toPlot4.m$size %in% 20:50,]
    # toPlot4.m.sel$size <- as.factor(toPlot4.m.sel$size)
    toPlot4.m.sel <- toPlot4.m[ toPlot4.m$size %in% 
                                  c(2,12,22,32,42,52,62,72,82,92,102,152,381) , ]
    toPlot4.m.sel$size <- as.factor(toPlot4.m.sel$size)
    
    ppg4 <-  ggplot(toPlot4.m.sel, aes( x=size , y=value , group=1)) +
      geom_line()+
      geom_point( alpha=0.4) +
      geom_vline(xintercept  = 20 , color="red" )+
      geom_vline(xintercept  = 50 , color="red"  )+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      facet_wrap(vars(variable), scales = "free")
    # ggbreak::scale_y_break(c(0.2, 0.65)) +
    # ggbreak::scale_y_break(c(0.675, 0.925))
    # save plots for 2 groups, used for DE
    
    
    pdf(height = 6, width = 12, file =  paste( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/SetSize/",
                                               "PDS_sizeTest_ctrlVSxprmnt.medianDiff.plot.BY10.pdf", sep = "") )
    print(ppg4)
    dev.off()
    
    
  }
  
  
  
}




#### 4. panel D. validate PDS in mouse spatial transcriptomics #### 
### read glom coordinates
TimGlom_coord <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/podoClust/TimGlom_coord.rda" )

### GSM5713367 Spatial plot, with glom morphology and PDS
  {
    
    ### read cell-type annotation
    ll <-list.files( pattern = "glom_annot" , full.names = T,
                     path = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/podoClust/BeadLocationsPodoTim/" )
    TimGlom_ctypes <- lapply( ll,  read.csv , row.names=1)
    names(TimGlom_ctypes) <- sub( "_.*","",basename(ll))
    
    
    ### Spatial plot: highlight gloms with circles colored by PDS
    library(ggforce)
    toPlotID <- "GSM5713367"
    toPlot <- TimGlom_coord[TimGlom_coord$geoID== "GSM5713367",] 
    
    ### read seurat object
    seur <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Spatial.Transcr/GSE190094_RAW/BTBR/podoClust/filtSeur/GSM5713367_filtSeur.rda")
    # annotate all glomerular cells
    seur$glom.KNN <-  TimGlom_ctypes$GSM5713367$cell_type[ match( colnames(seur),
                                                                  TimGlom_ctypes$GSM5713367$barcode)]
    seur$glom.KNN[!is.na(seur$glom.KNN)] <- "glom" 
    
   
    # spatial plot of Slide'seq.v2 data from kidney-tissue slice
    pp<- SpatialDimPlot(seur, group.by = "glom.KNN",stroke = 0.1)+
      scale_fill_colorblind()+ theme( text = element_text( size=16))

    # plot circles coloured by PDS
    cc <- ggplot(toPlot, aes(x0=y, y0=x, r=radii, color=PDS, fill=PDS)) + 
      geom_circle(alpha=0.4) + coord_equal() + 
      ggeasy::easy_remove_axes()+theme_void()+  
      theme( text = element_text( size=16))+
      theme( legend.position = "bottom") +  scale_y_reverse() +
      scale_color_gradientn( colours = brewer.pal( n = 9 , name = "YlOrBr"))+
      scale_fill_gradientn( colours = brewer.pal( n = 9 , name = "YlOrBr"))
    
      # combine on 2 panels, merge panels in the illustrator
   ppl  <- cowplot::plot_grid( cc, pp , rel_heights = c(1, 1), rows = 2)
      
   
   # save plots
   pdf(height = 15, width = 15, file = "Dimplot_glom.PDS_GSM5713367.pdf")
    ppl
   dev.off()
   
  }

### regression line of PDS and podo N (normalised by the glom area)
  {
  
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
  

  pdf( width = 8, height = 6 , file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Spatial.trnscrptm/GSE190094/BTBR/podoClust/Glom.metrics_scatter.bxplts.test.pdf")
  print(gg2)
  dev.off()
  
  }

### Suppl.Fig 3A boxplots of Boxplots showing glomerular size estimates, 
# podocyte bead number, and normalized podocyte density for diabetic Btbrob/ob 
# and control mice as calculated on Slide-seqV2 spatial transcriptomics data (GSE190094)
  {
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
  
}








