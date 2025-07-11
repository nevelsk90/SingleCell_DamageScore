# ### Signature size test ### #


mallinfo::malloc.trim()
gc()

options( connectionObserver = NULL )
.libPaths(c("/home/tim_nevelsk/R/x86_64-pc-linux-gnu-library/4.0", "/media/tim_nevelsk/WD_tim/SOFT/R"))


library( Seurat)
require( GSEABase)
require( AUCell)
require( ggplot2)
library( ggpubr)
library( ggrepel)
library( ggthemes)
require( reshape2 )

#### load functions
# generate damage signature
source("https://raw.githubusercontent.com/PauUng/HepatocyteDamageScore/master/SharedFunctions.R")
# calculate damage score
source("https://raw.githubusercontent.com/nevelsk90/R_scripts/master/usefulRfunc.r")

DS_all <- read.table( header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")

allPodoGenes <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/SCSN_allPodoGenes.rda")

### load subsampled sc data
listSCSN.1K.sampl <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN_samples.1K.22.12.23.rda")


### define a vector of siyes to test 
  size_vec <- c(  1:nrow(DS_all) )
### create the respective gene set collection
  size_Gsets <- Reduce( c , lapply(size_vec, function(ssize)
    {
    
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
   listSCSN.PDS_DSsizeTest<- readRDS(file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/Dev.Valid/listSCSN.PDS42_DSsizeTest.22.05.25.rda")
    
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
      datt$dataSet <- names( listSCSN.PDS_DSsizeTest)[ii]
      print( dim(datt))
      return(datt)
    }) )
    
    
    ### plot mean of the differences between means of ctrl and exprmntal
    toPlot3 <-  Reduce( rbind,  lapply( seq(listSCSN.PDS_DSsizeTest) , function(ii){
      datt <- listSCSN.PDS_DSsizeTest[[ii]]@meta.data
      # datt <-   newSeu@meta.data
      
      datt<- datt[ , c("gtypeDE", grep("PDS.size.*",names(datt), value = T)) ]
    # scale
      datt[,grep("PDS.size.*",names(datt), value = T)] <- scale(
        datt[,grep("PDS.size.*",names(datt), value = T)], scale = T )
      # print( dim(datt))
      
      datt.ag <- aggregate( . ~ gtypeDE , dat= datt,  FUN = mean)
      datt.diff <- datt.ag[2,-1] - datt.ag[1,-1]
      
      XX <- data.frame( PDS.diff.Mean = t( datt.diff )[,1], dataSet= names(listSCSN.PDS_DSsizeTest)[ii],
                  size=size_vec)

    }) ) 
    
    
    
    toPlot4 <- Reduce( rbind,  lapply(size_vec, function(ii){
      data.frame( mean.of.means = mean( toPlot3$PDS.diff.Mean[toPlot3$size==ii] ) , 
                  var.of.means = var(toPlot3$PDS.diff.Mean[toPlot3$size==ii] ))
    }))
    toPlot4$size <- size_vec 
    toPlot4.m <- reshape2::melt(toPlot4, id.var="size")
    # toPlot4.m.sel <- toPlot4.m[toPlot4.m$size %in% 20:50,]
    # toPlot4.m.sel$size <- as.factor(toPlot4.m.sel$size)
    
    sel.size <- c(2,12,22,32,42,52,62,72,82,92,102,152,202)
    sel.size <- c(5,10,15,20,30,42,50,100,150,200)
    
    toPlot4.m.sel <- toPlot4.m[ toPlot4.m$size %in% sel.size ,   ]
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
                                               "PDS_sizeTest_ctrlVSxprmnt.medianDiff.plot.v3.pdf", sep = "") )
    print(ppg4)
    dev.off()
    
    
  }
  
### test
  {
    listSCSN.PDS_DSsizeTest<- readRDS(file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/Dev.Valid/listSCSN.PDS42_DSsizeTest.22.05.25.rda")



    ### plot mean of the differences between means of ctrl and exprmntal
    sel.size <- c(2,12,22,32,42,52,62,72,82,92,102,152,202)

    toPlot5 <- sapply( seq(size_vec), function(jj){

      print(jj)
      ssize = sel.size[jj]

      toTest <- Reduce( rbind, lapply( seq(listSCSN.PDS_DSsizeTest) , function(ii){
        datt <- listSCSN.PDS_DSsizeTest[[ii]]@meta.data
        # datt <-   newSeu@meta.data
        datt <- datt[,c( "gtypeDE","gtype","group","sample", paste0("PDS.size", ssize)) ]
        # exclude PDS columns not related to the size test
        datt<- datt[ , ! colnames(datt) %in%
                       c("PDS.42.005",  "PDS.42.005.2" ,"PDS") ]
        # scale
        datt[,paste0("PDS.size", ssize)] <- scale(
          datt[,paste0("PDS.size", ssize)], scale = T )

        return(datt)
      }))

      toTest$gtypeDE <- as.factor(toTest$gtypeDE)

            t.test( toTest[,ncol(toTest)] ~gtypeDE, toTest )$statistic

    })

    toPlot <- data.frame( size= ( sel.size), value = abs(toPlot5))
    ppg5 <-  ggplot(toPlot, aes( x=size , y= value , group=1)) +
      geom_line()+
      geom_point( alpha=0.4) +
      # geom_vline(xintercept  = 20 , color="red" )+
      # geom_vline(xintercept  = 50 , color="red"  )+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.text = element_text(size = 18),
            text = element_text(size = 18)) +
      scale_x_continuous(
        breaks = toPlot$size,
        labels = toPlot$size
      ) + ylab("T-statistic") + xlab("Gene signature size")
    ppg5
    
    pdf(height = 5, width = 10, file =  paste( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/SetSize/",
                                               "PDS_sizeTest_ctrlVSxprmnt.tstat.plot.pdf", sep = "") )
    print(ppg5)
    dev.off()
    
  }
  
  

  



