###=====# Single-Cell Resolution of Cellular Damage Illuminates Disease Progression #=======###
###=====# supplementary figure 1 PDS code #=======###

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

# setwd
setwd( "/home/tim_nevelsk/PROJECTS/PODOCYTE/WRITING/PDS_manuscript/Figures/Figure1" )
inputdir <-  "/media/tim_nevelsk/WD_tim/PROJECTS/WRITING/PDS_manuscript/Figure_input" 

# load necessary code
source("https://raw.githubusercontent.com/PauUng/HepatocyteDamageScore/master/SharedFunctions.R")
source("https://raw.githubusercontent.com/nevelsk90/R_scripts/master/usefulRfunc.r")  


# load damage signatures
DS_all <- read.table( header = T,  file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")
# load expression data
listSCSN.1K <- readRDS( paste0( inputdir, "/listSCSN_1K.22.12.23.rda") )

#### Suppl.Fig 1A Signature size test ####

  library(GSEABase)
  DSignature <- DS_all
  # DSignature <- DSignature[order(DSignature$mean_rankSC),]
  
  
  size_vec <- c( 5, 10, 15, 20, 30, 42 , 50, 100, 200, nrow(DSignature) )
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
  
  ### analyse effect of a set size on individual studies
  ceilThrsh <- 0.05
  
  # iterate over sc studies and calculate PDS
  set.seed(42)
  listSCSN.PDS_DSsizeTest <- lapply( seq(listSCSN.1K.sampl),
                                     function(ii)
                                     {
                                       ### uppercase row. names
                                       print(ii)
                                       newSeu <- listSCSN.1K.sampl[[ii]]
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
  # saveRDS(listSCSN.PDS_DSsizeTest, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/Dev.Valid/listSCSN.PDS42_DSsizeTest.rda")
  # saveRDS(listSCSN.PDS_DSsizeTest, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/Dev.Valid/listSCSN.PDS_DSsizeTest.rda")
  
  ### plot PDS distributions
  {
    ## prepare the table
    toPlot <- lapply( seq(listSCSN.PDS_DSsizeTest) , function(ii){
      datt <- listSCSN.PDS_DSsizeTest[[ii]]@meta.data
      datt <-   newSeu@meta.data
      datt <- datt[,c( "gtypeDE","gtype","group","sample", 
                       grep("PDS",names(datt), value = T))]
      datt[,grep("PDS.*",names(datt), value = T)] <- scale(
        datt[,grep("PDS.*",names(datt), value = T)], scale = T )
      return(datt)
    })
    
    toPlot.melt <- reshape2::melt( toPlot , id.vars=c("gtypeDE","gtype","group","sample"))
    toPlot.melt$L1 <- factor( toPlot.melt$L1, labels =(names(listSCSN.1K.sampl)))
    # remove pdss2 6 w as having no effect, and old PDS
    toPlot.melt <- toPlot.melt[toPlot.melt$group!="Pdss2_6" &
                                 toPlot.melt$variable!="PDS" ,]
    toPlot.melt$variable <- droplevels(toPlot.melt$variable)
    
    ### make plots
    ppg1 <- ggplot( toPlot.melt , aes( y=value , 
                                       x=gtypeDE , 
                                       color=gtypeDE)) +
      geom_violin(trim = T,scale="width") + theme_bw() + scale_color_colorblind() + 
      stat_summary(
        fun.data = "mean_sdl",  fun.args = list(mult = 1), 
        geom = "pointrange", color = "black"
      )+ theme( text = element_text(size=20), legend.position = "bottom", 
                axis.text.x = element_blank()
                # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
      ) +
      facet_grid( cols = vars(variable))
    # # separate each sc.sn study
    # ppg2 <- ppg1+facet_grid(rows = vars(L1), scales = "free")+
    #   stat_summary(fun=median, geom="line",
    #                aes(group=gtypeDE ,color=gtypeDE))
    
    # ##  add t stats
    # t_stat<- sapply( unique(toPlot.melt$variable), function(ii){
    #   tres <- t.test(x= toPlot.melt$value[toPlot.melt$variable==unique(toPlot.melt$variable)[ii] & 
    #                                         toPlot.melt$gtypeDE=="experimental" ],
    #                  y= toPlot.melt$value[toPlot.melt$variable==unique(toPlot.melt$variable)[ii] &
    #                                         toPlot.melt$gtypeDE=="control"] )
    #   tres$statistic
    # })
    # names(t_stat) <- unique(toPlot.melt$variable)
    # # create df of tests
    # test.means <- compare_means(
    #   value ~ gtypeDE, data = toPlot.melt, group.by = "variable",
    #   method = "t.test", ref.group = "control",p.adjust.method = "bonferroni"
    # )
    # # add tstat to the test df
    # test.means$tstat <- round(t_stat,2)
    # # plot
    # ppg3 <- ppg1 +
    #   stat_pvalue_manual(data = test.means, label="tstat",
    #                      x = "variable", 
    #                      y.position = 4 )+ 
    #   scale_x_discrete(limits = rev(levels(toPlot.melt$variable)))+
    #   stat_summary(fun=mean, geom="line",
    #                aes(group=gtypeDE ,color=gtypeDE))+ 
    #   coord_flip()+ facet_grid( cols = vars(variable))
    
    
    
    # save plots for 2 groups, used for DE
    pdf(height = 6, width = 12, file =  paste( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/SetSize/",
                                               "PDS_sizeTest_scsn.VlnPlot.pdf", sep = "") )
    print(ppg1)
    dev.off()
    
    
  }
  
  
  ### plot correlation with PDS
  {
    # combine metadata from 3 experiments 
    datt <- Reduce( rbind , lapply( 1:3, function(ii) {
      XX <- listSCSN.PDS_DSsizeTest[[ii]]@meta.data 
      return(XX[,c( "gtypeDE","group","sample","gtype", 
                    grep("PDS",names(XX), value = T)
      )])
    }))
    datt <- datt[, colnames(datt)!="PDS.42"]
    toPlot.melt <- reshape2::melt( datt , id.vars=c("gtypeDE","gtype","group","sample"))
    # toPlot.melt$L1 <- factor( toPlot.melt$L1, labels =(names(listSCSN.PDS)))       
    
    # treat carefully 21 week Pdss2 samples since they have only per group measurements
    datt1 <- datt[ datt$sample %in% c( "146985", "146986", "143485" , "143486") ,]
    # aggregate
    aggPDS1 <- aggregate( .~group , FUN = mean , 
                          data= datt1[ , c("group", 
                                           grep("PDS",names(datt), value = T))] )
    aggPDS1 <- aggPDS1[rep(seq_len(nrow(aggPDS1)), each = 2), ]
    aggPDS1$group <-  c( "146985", "146986", "143485" , "143486")
    colnames(aggPDS1)[1] <- "sample"
    # the rest of samples
    datt2 <- datt[ !(datt$sample %in%  c("146985", "146986", "143485" , "143486")), ]
    aggPDS2 <- aggregate( .~sample , FUN=mean,
                          data= datt2[ , c( "sample",
                                            grep("PDS",names(datt), value = T))] ) 
    
    aggPDS <- rbind(aggPDS2, aggPDS1)
    
    aggPDS$AlbCrRatio <- annot_tab$AlbCrRatio[ match( sub("SID","" ,aggPDS$sample) , 
                                                      annot_tab$CCG_Sample_ID)]
    aggPDS$group <- annot_tab$group[ match( sub("SID","" ,aggPDS$sample ), 
                                            annot_tab$CCG_Sample_ID)]
    aggPDS$gtype <- as.factor(annot_tab$Genotype[ match( sub("SID","" ,aggPDS$sample ), 
                                                         annot_tab$CCG_Sample_ID)])
    # melt 
    aggPDS.melt <- reshape2::melt( aggPDS ,
                                   id.vars= names(aggPDS)[!names(aggPDS) %in%
                                                            grep("PDS",names(aggPDS), value = T)])
    gg <- ggplot2::ggplot( data = aggPDS.melt, aes( x=value, y=log(AlbCrRatio))) +
      geom_point(  size=6, aes(col=gtype)) +
      theme_bw() +  theme( text = element_text(size = 22)) + 
      geom_smooth(method='lm', se = FALSE) + ggtitle( "Nphs2mut podocytes")+ 
      stat_cor( size=5, method = "spearman" )  +  
      geom_text(aes(label = sample  ), size=6, position = "dodge")
    # facet_grid(cols = vars(variable), scales = "free")
    gglist <-  gg+ ggforce::facet_wrap_paginate( ~variable, scales = "free", 
                                                 ncol = 2 ) 
    
    
    # distribution plot for NPHS2
    dattX <- datt[datt$sample%in% listSCSN.PDS$Nphs2$sample,]
    dattX.melt <- reshape2::melt( dattX ,
                                  id.vars= names(dattX)[!names(dattX) %in%
                                                          grep("PDS",names(dattX), value = T)])
    ggplot2::ggplot( data = dattX.melt, aes( x=value, color=group)) +
      geom_density(size=1.5) +
      theme_bw() +  theme( text = element_text(size = 22))+
      ggforce::facet_wrap_paginate( ~variable, scales = "free", 
                                    ncol = 2 )
    
    XX <- aggPDS[aggPDS$sample%in% listSCSN.PDS$Nphs2$sample,]
    XX.melt <- reshape2::melt(XX,  id.vars= names(XX)[!names(XX) %in%
                                                        grep("PDS",names(XX), value = T)])
    ggplot(XX.melt, aes(x=group,y=value, color=group))+geom_jitter()+theme_bw()+
      geom_label(aes(label=sample))+
      ggforce::facet_wrap_paginate( ~variable, scales = "free", 
                                    ncol = 2 )
  }
  
#### Suppl Fig 1C randomisation test ####
#### Suppl Fig 1E Cell-type specificity test #### 
