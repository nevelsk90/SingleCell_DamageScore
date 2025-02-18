# ###################################################### #
# release memory
mallinfo::malloc.trim()
gc()

options( connectionObserver = NULL )
.libPaths(c("/home/tim_nevelsk/R/x86_64-pc-linux-gnu-library/4.0", 
            "/media/tim_nevelsk/WD_tim/SOFT/R"))
library( AnnotationDbi )
library( org.Mm.eg.db )
library( BSgenome.Mmusculus.UCSC.mm10  )
library(devtools)
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

# mart_homo <- useMart( "ensembl",dataset="hsapiens_gene_ensembl" , host="www.ensembl.org")
 mart_mouse <- useMart( "ensembl",dataset="mmusculus_gene_ensembl" )
 tx2gene <- getBM( attributes=c( 'ensembl_gene_id', 'external_gene_name',"entrezgene_id"),  mart = mart_mouse)

#### load functions
# calculate damage score
 source("https://raw.githubusercontent.com/PauUng/HepatocyteDamageScore/master/SharedFunctions.R")
 source("https://raw.githubusercontent.com/nevelsk90/R_scripts/master/func_analysis.R")  

#### load data #### 
# ### load KFO and GSE146912 data
# listSCSN.1Kcells.PDS <- readRDS( file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN.1Kcells.PDSlist_23.12.2022.rda")
# Wt1Nphs2.4Kcells.PDS <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/Wt1Nphs2.4Kcells.PDSlist_23.12.2022.rda")
listSCSN.1K <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN_1K.22.12.23.rda")
listSCSN <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN.PDS_22.12.23.rda")
listSCSN.1K.sampl <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN_samples.1K.22.12.23.rda")

# # load PDS subsampled to 1K cells per experiment and limited to -0.5 ; 0 PDS range
# listSCSN.PDSlimit.1K <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN.PDSlimit.1K.rda")
annot_tab <- read.table(sep = "\t",header = T, "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/Sample_Names_KFO.csv")
annot_tab$group <- paste( annot_tab$Genotype,annot_tab$Age_weeks,sep = "_" )




##  load Damage signature
# DS_all <- read.table( header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.02.06.2022.tsv")
DS_all <- read.table( header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")

### load all genes expressed in SCSN podocytes
allPodoGenes <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/SCSN_allPodoGenes.rda")
allPodoGenes_mean <-  readRDS(file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/SCSN_allPodoGenes_mean.rda")

### load TFs expressed in SCSN podocytes
TF_MeanMed <- readRDS(  file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/TFA/TFs_SNSC.podo.MeanMed.rda")

### load the prior
  # ATACseq_tgenesM.TFtc <- readRDS( file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/ATACseq_TFtc.tgenesM_77.TF.rda")
  # ATACseq_tgenesM <- readRDS( file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/ATACseq_tgenesM.rda")
  ATACseq_tgenes <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/ATACseq/Podocytes/TOBIAS/FIMO/atac.podo_tobias.fimo110TF.p5e4_TFtc.rda")
  ATACseq_tgenes.S <- Reduce( cbind , lapply( ATACseq_tgenes , "[[", "score" ))
  ATACseq_tgenes.Q <- Reduce( cbind , lapply( ATACseq_tgenes , "[[", "qvalue" ))
  colnames(ATACseq_tgenes.S) <- colnames(ATACseq_tgenes.Q) <- names(ATACseq_tgenes)
  ATACseq_tgenesM <- ATACseq_tgenes.S * ( ATACseq_tgenes.Q <0.1)
  rownames(ATACseq_tgenesM) <- rownames(ATACseq_tgenes[[1]])
  ATACseq_tgenesM <- ATACseq_tgenesM[ rowSums(ATACseq_tgenesM)>0 , colSums(ATACseq_tgenesM)>0  ]
  # ATACseq_tgenesFiltNames <- read.table( sep = "\t", file="/media/tim_nevelsk/WD_tim/PROJECTS/MiraldiGRN/MiraldiStyle_Tim/inputs/geneLists/mouseTFprior_podoExprsd.tsv")
  ## filter TFprior based on expression
  # ATACseq_tgenesM.podo <- ATACseq_tgenesM[ ,  colnames(ATACseq_tgenesM)%in% TF_MeanMed ]
  # ATACseq_tgenesM.podo <- ATACseq_tgenesM.podo[ rowSums(ATACseq_tgenesM.podo)>0 , colSums(ATACseq_tgenesM.podo)>0  ]
  
### list of all mouse TFs
  TFmouse_MAT <- read.table( header = T, sep = "\t", fill = T , "/media/tim_nevelsk/WD_tim/ANNOTATIONS/CISBP/Mus_musculus_2021_01_19_5 55_pm/TF_Information.txt")
  TFmouse_all <- unique( TFmouse_MAT$TF_Name)
  
### Wt1 targets
  WT1chipTgenes <- read.table( sep="\t","/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/GRN/ATACseq_GRNprior/Tgenes/Podo_chipWt1_TFtc.csv")

  
#### calculate and plot correlation of genes with PDS ####
 
  ### calculate correlation 
    {
      # calculate correlation for full data, with sampling
      # to get comparable rho's for experiments with different N of cells
      # turns out this is almost equivalent to using the full data
    
      # full data without sampling
      genecorrPDS_Sprmn <- lapply( seq( listSCSN ), function(ii)
        {
        require(metap)
        print(ii)
        newSeu <- listSCSN[[ii]]


          # get cell IDs persample type
          if( "control" %in% newSeu$gtypeDE ) {
            cell.wt <- colnames(newSeu)[newSeu$gtypeDE == "control" ]
            expr <- subset(newSeu, cells= cell.wt)
            expr <- expr@assays$RNA@data[ rowSums( round(expr@assays$RNA@counts) > 0 ) > ncol(expr)*0.005 , ]
            cor.wt <- psych::corr.test( as.matrix(t(expr )) , newSeu$PDS[ cell.wt ] ,
                                        method = "spearman" )
            gc()
          } else cor.wt <- NA
        
        if("experimental" %in% newSeu$gtypeDE) {
          cell.mut <- colnames(newSeu)[newSeu$gtypeDE != "control" ]
          
          expr <- subset(newSeu, cells= cell.mut)
          expr <- expr@assays$RNA@data[ rowSums( round( expr@assays$RNA@counts) > 0 ) > ncol(expr)*0.01 , ]
          
          cor.mut <- psych::corr.test( as.matrix(t(expr ) ), newSeu$PDS[ cell.mut ] ,
                                       method = "spearman" )
          gc()
        } else  cor.mut <- NA
         
        corrList <-  list( cor.wt ,
                            cor.mut )

        return(corrList)
      })
      
      names(genecorrPDS_Sprmn) <-names( listSCSN )
      saveRDS( genecorrPDS_Sprmn , file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/SCSN_genecorrPDS_Sprmn.22.12.23.rda")
      
      # calculate sample-wise correlations
      genecorrPDS_smpls.Sprmn <- lapply( seq( listSCSN.1K.sampl ), 
                                         function(ii )
        {
          print(ii)
  
           datt <- listSCSN.1K.sampl[[ii]]
          ## fit curves for each sample
          snames <- names(table(datt$sample))[ table(datt$sample)> 50 ]
          
          persample <- lapply( seq(snames),
                                          function(jj){
                                            print(snames[jj])
                                            # extract normalised data to avoid spurios corr, 
                                            # due to FSGS related change the lib.size
                                            datt.smpl <-datt@assays$RNA@data[
                                              , datt$sample==snames[jj] ]
                                            datt.smpl <- datt.smpl[ rowSums(datt.smpl>0) > sqrt(ncol(datt.smpl)) , ]
                                            datt.smpl <- t( as.matrix( datt.smpl ))
                                            
                                            
                                            # datt.smpl <- datt.smpl[ datt.smpl$PDS.42>  mean(datt.smpl$PDS.42)-3*sd(datt.smpl$PDS.42) & 
                                            #                           datt.smpl$PDS.42<  mean(datt.smpl$PDS.42)+3*sd(datt.smpl$PDS.42), ]
                                            PDS.42 <- datt$PDS[ datt$sample == snames[jj] ]
                                            
                                            cor.smpl <- psych::corr.test(   y= PDS.42 ,
                                                                            x = datt.smpl,
                                                                            method = "spearman" )
                                            
                                           
                                            # calculate q-values
                                            cor.qval.smpl <- qvalue::qvalue( cor.smpl$p )$qvalues

                                            return( cbind.data.frame( cor.r=cor.smpl$r , 
                                                                      cor.p=cor.smpl$p  , 
                                                                      cor.qval = cor.qval.smpl ,
                                                                      sample=snames[jj] ,
                                                                      gtypeDE=unique( datt$gtypeDE[
                                                                        datt$sample == snames[jj] ]
                                                                      )) )
          })
          
          names(persample) <- snames
          return( persample )

        } )
      saveRDS( genecorrPDS_smpls.Sprmn , file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/SCSNsmpls_genecorrPDS_Sprmn.13.05.24.rda")
      

     
    }
  

  ### extract q and r value tables
    {
      
      # genecorrPDS_Sprmn <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PDScorr/SCSN_genecorrPDS_Sprmn.29.09.23.rda")
      genecorrPDS_Sprmn <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/SCSN_genecorrPDS_Sprmn.22.12.23.rda")

   
    ### extract r
      {
        
        X1 <-  lapply( lapply( genecorrPDS_Sprmn,"[[", 1),"[[", 1 )
        X2 <-  lapply( lapply( genecorrPDS_Sprmn,"[[", 2),"[[", 1 )
        
        allgenes <- Reduce( union, lapply( seq(X2) , function(ii) union( 
          rownames(X1[[ii]]) , rownames(X2[[ii]]) )))

        names(X1) <- names(X2) <- names(genecorrPDS_Sprmn)
        
        X1[is.na(X1)] <- list( X1$nephr.D1 , X1$nephr.D1)

        cnames <- names(X1)
        X1 <- Reduce( cbind.data.frame,  lapply( seq(X1) , function(ii) X1[[ii]][
            match( allgenes, rownames( X1[[ii]])),] ) )
        colnames(X1) <-paste0( cnames, "_ctrl")
        
        cnames2 <- names(X2)
        X2 <- Reduce( cbind.data.frame,  lapply( seq(X2) , function(ii) X2[[ii]][
          match( allgenes, rownames( X2[[ii]])),] ) )
        colnames(X2) <-paste0( cnames2, "_xprmnt")
        rownames(X1) <- rownames(X2) <- allgenes
        
        ### r-value table
        genecorrPDS.Sprmn_r <- cbind(X1,X2)
        
        genecorrPDS.Sprmn_r.cntrd <- apply( genecorrPDS.Sprmn_r, 2, scale, scale = F)
        rownames(genecorrPDS.Sprmn_r.cntrd) <- rownames(genecorrPDS.Sprmn_r)
    
        saveRDS( genecorrPDS.Sprmn_r.cntrd , file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/SCSN_genecorrPDS_Sprmn.r.cntrd.rda")
        # # group all controls and all exprmntl smpls
        # genecorrPDS.Sprmn_qval <- genecorrPDS.Sprmn_qval[,c(seq(from=1,to=14,by=2), seq(from=2,to=14,by=2)) ]
        # genecorrPDS.Sprmn_p <- genecorrPDS.Sprmn_p[,c(seq(from=1,to=14,by=2), seq(from=2,to=14,by=2)) ]
        # genecorrPDS.Sprmn_r.cntrd <- genecorrPDS.Sprmn_r.cntrd[,c(seq(from=1,to=14,by=2), seq(from=2,to=14,by=2)) ]
        # genecorrPDS.Sprmn_r <- genecorrPDS.Sprmn_r[,c(seq(from=1,to=14,by=2), seq(from=2,to=14,by=2)) ]
      }
      
    ### extract p
      {
        X1 <- lapply( genecorrPDS_Sprmn,"[[", 1)
        X2 <- lapply( genecorrPDS_Sprmn,"[[", 2)
        names(X1) <- names(X2) <- names(genecorrPDS_Sprmn)
       
         X1[is.na(X1)] <- list( X1$nephr.D1 , X1$nephr.D1)
        
        
        X1 <- lapply( X1 ,"[[", 4)
        X2 <- lapply( X2 ,"[[", 4)
        
  
        cnames <- names(X1)
        X1 <- Reduce( cbind.data.frame,  lapply( seq(X1) , function(ii) X1[[ii]][
          match( allgenes, rownames( X1[[ii]])),] ) )
        colnames(X1) <-paste0( cnames, "_ctrl")
        
        cnames2 <- names(X2)
        X2 <- Reduce( cbind.data.frame,  lapply( seq(X2) , function(ii) X2[[ii]][
          match( allgenes, rownames( X2[[ii]])),] ) )
        colnames(X2) <-paste0( cnames2, "_xprmnt")
        rownames(X1) <- rownames(X2) <- allgenes
        
        ### r-value table
        genecorrPDS.Sprmn_p <- cbind(X1,X2)
        genecorrPDS.Sprmn_q <- apply(genecorrPDS.Sprmn_p, 2, 
                                     function(X) qvalue::qvalue(X)$qvalues)
        
        # # group all controls and all exprmntl smpls
        # genecorrPDS.Sprmn_qval <- genecorrPDS.Sprmn_qval[,c(seq(from=1,to=14,by=2), seq(from=2,to=14,by=2)) ]
        # genecorrPDS.Sprmn_p <- genecorrPDS.Sprmn_p[,c(seq(from=1,to=14,by=2), seq(from=2,to=14,by=2)) ]
        # genecorrPDS.Sprmn_r.cntrd <- genecorrPDS.Sprmn_r.cntrd[,c(seq(from=1,to=14,by=2), seq(from=2,to=14,by=2)) ]
        # genecorrPDS.Sprmn_r <- genecorrPDS.Sprmn_r[,c(seq(from=1,to=14,by=2), seq(from=2,to=14,by=2)) ]
      }
      
  }
  
  ### convert to binary matrix
    {
    genecorrPDS.Sprmn_bin <- genecorrPDS.Sprmn_qval
    sigL <- 0.01
    genecorrPDS.Sprmn_bin[genecorrPDS.Sprmn_bin< sigL ] <- 0
    genecorrPDS.Sprmn_bin[genecorrPDS.Sprmn_bin> sigL ] <- 1
    # modify column names
    colnames(genecorrPDS.Sprmn_bin ) <- c( paste( names(listSCSN.PDS), "ct", sep = "_") ,
                                           paste( names(listSCSN.PDS), "expr", sep = "_")     )
    
    
  }

  ### make a matrix that counts FSGS studies where a gene significantly correlaties with PDS 
    {
      sigL <- 0.01
      
     
    genecorrPDS.Sprmn_freq <- t( sapply( seq(rownames(genecorrPDS.Sprmn_p)), 
                                         function( ii )
      {
      rowSums( sapply( seq( length(genecorrPDS_Sprmn)), function(jj){
        c( genecorrPDS.Sprmn_p[ii,jj] < sigL , 
           genecorrPDS.Sprmn_p[ii,jj+ length(genecorrPDS_Sprmn)] < sigL , 
           genecorrPDS.Sprmn_p[ii,jj] < sigL &  
             genecorrPDS.Sprmn_p[ii,jj+  length(genecorrPDS_Sprmn)] > sigL ,
           genecorrPDS.Sprmn_p[ii,jj] > sigL &  
             genecorrPDS.Sprmn_p[ii,jj+  length(genecorrPDS_Sprmn)] < sigL,
           genecorrPDS.Sprmn_p[ii,jj] < sigL &  
             genecorrPDS.Sprmn_p[ii,jj+  length(genecorrPDS_Sprmn)] < sigL)
      }), na.rm = T)
    }))
    
    rownames(genecorrPDS.Sprmn_freq) <- rownames(genecorrPDS.Sprmn_p)
    genecorrPDS.Sprmn_freq <- as.data.frame(genecorrPDS.Sprmn_freq)
    genecorrPDS.Sprmn_freq$total <- genecorrPDS.Sprmn_freq$V1 + genecorrPDS.Sprmn_freq$V2
    colnames(genecorrPDS.Sprmn_freq) <- c( "control", "experiment", "ctr.only","exp.only", "both", "total" ) 
    
    }
  
  ### prepare gene data, all datasets
    {
      PDS42.SpCor_barcode <- genecorrPDS.Sprmn_r.cntrd
      # combine all controls and all exprmnt in indiv. columns
      PDS42.SpCor_barcode <- cbind(  ctrl_mean=rowMeans(
        PDS42.SpCor_barcode[ , grepl( "ctrl" , colnames(PDS42.SpCor_barcode) ) ], na.rm = T), 
        PDS42.SpCor_barcode[, grepl( "xprmnt" , colnames(PDS42.SpCor_barcode) )])
    
    colnames(PDS42.SpCor_barcode) <- c( "ctrl_mean", names(genecorrPDS_Sprmn ))
    saveRDS( PDS42.SpCor_barcode , file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PDS42.SpCor_barcode.03.05.24.rda")
    
    ### save a table with extra info
    # readTF clustering info
    mouseTFmm.podo_clust <- yaml::read_yaml("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/ATACseq/Podocytes/TOBIAS/motif_clusters/motif_comparison_seqcor0.3_clusters.yml")
    names(mouseTFmm.podo_clust) <-sapply(seq(mouseTFmm.podo_clust), 
                                         function(ii)  paste( gsub( ".*__| ", "", mouseTFmm.podo_clust[[ii]]), collapse=":") )
    mouseTFmm.podo_clust <- setNames( rep(names(mouseTFmm.podo_clust), lengths(mouseTFmm.podo_clust)),
              unlist(mouseTFmm.podo_clust, use.names=F))
    
     # add TFtgt info 
    TFreg <- ATACseq_tgenesM[ rownames(ATACseq_tgenesM) %in% rownames(PDS42.SpCor_barcode),
                              colnames(ATACseq_tgenesM) %in% sub( " ", "",names(mouseTFmm.podo_clust)) ] 
    # combine targets of all TFs belonging to the same PWM cluster
    cclust_names <-  mouseTFmm.podo_clust[ match( 
      colnames(TFreg) , sub( " ", "",names(mouseTFmm.podo_clust)) ) ] 
    TFreg_clust <- t( rowsum( t(TFreg), group=cclust_names) )
    saveRDS( TFreg_clust , file="/home/tim_nevelsk/PROJECTS/PODOCYTE/ATACseq/atac.podo_tobias_fimo_tftc.0.01.GRN.PWMclust.rda")
    
    # vectorise TF-target gene list
    TFregVec <- Reduce( c, apply( TFreg, 1, function(X)
      {
      stringr::str_c( unique( colnames(TFreg)[ X>0 & !is.na(X)]),collapse = ",")
    }))
    names(TFregVec) <- rownames(TFreg)
    # vectorise clustered TF-target list
    TFreg_clust.Vec <- Reduce( c, apply( TFreg_clust, 1, function(X)
    {
      stringr::str_c( unique( colnames(TFreg_clust)[ X>0 & !is.na(X)]),collapse = ",")
    }))
    names(TFreg_clust.Vec) <- rownames(TFreg_clust)
    
    # table
    PDS42.SpCor_barcode_tab <- cbind.data.frame( 
      DSgenes = ifelse( rownames(PDS42.SpCor_barcode) %in% DS_all$gene_symbol[1:42], "DS.42",
                        ifelse( rownames(PDS42.SpCor_barcode) %in% DS_all$gene_symbol , 
                                                                    "DS.rest", "")),
                                          Nstd.PDScor0.01.xprmnt =rowSums( genecorrPDS.Sprmn_p[,10:18] < sigL, na.rm = T ) ,
                                          TFtrgt= TFregVec[match( rownames(PDS42.SpCor_barcode), names(TFregVec) )], 
                                          gNameUP = toupper(rownames(PDS42.SpCor_barcode)) ,
                                          PDS42.SpCor_barcode ,
      gName= rownames(PDS42.SpCor_barcode))
    
    # calculate average cor for TF clust
    PDS42.SpCor_barcode.TFclust <- as.data.frame( t( sapply( seq( ncol( TFreg_clust) ), 
                                           function(ii){
                                             print(ii)
                                             ggenes <- colnames(TFreg_clust)[ii]
                                             ggenes <- unlist( strsplit(ggenes,split = ":") )
                                             datt <- PDS42.SpCor_barcode[ 
                                               rownames(PDS42.SpCor_barcode) %in% ggenes,]
                                             if( length(ggenes)>1) {
                                                expr <- colSums( datt , na.rm = T)
                                             } else  expr <- datt

                                             c( 
                                               DSgenes=ifelse( sum(ggenes %in% DS_all$gene_symbol[1:42])>0 , "DS.42", 
                                                                    ifelse( sum(ggenes %in% DS_all$gene_symbol)>0 , 
                                                                            "DS.rest", "")),
                                                 Nstd.PDScor0.01.xprmnt=sum( colSums(genecorrPDS.Sprmn_p[ggenes,10:18] < 
                                                                                        sigL), na.rm = T ),
                                                 TFtrgt="",
                                                 gNameUP= "", 
                                                 expr )
                                           } ) ) )
    rownames(PDS42.SpCor_barcode.TFclust) <- paste0( sub(":.*","",colnames( TFreg_clust )), "_clust")
    PDS42.SpCor_barcode.TFclust$gName <-  paste0( sub(":.*","",colnames( TFreg_clust )), "_clust")

    PDS42.SpCor_barcode.combo <- rbind(  PDS42.SpCor_barcode_tab, PDS42.SpCor_barcode.TFclust)
      
      write.table( PDS42.SpCor_barcode.combo , sep = "\t", quote = F, col.names = NA, 
                   file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/PDS42.sprmCorr.barcode_centered.05.01.24.tsv")
  

      
    ### plot corelation of PDS with specific genes
      {
        genes <- grep("nup", rownames(toPlot.all_tab), ignore.case = T, value = T)
        
        # extract significance info 
        qvalTab <- genecorrPDS.Sprmn_qval[ genes , ]
        qvalTab[is.na(qvalTab)]<- 1
        qvalTab[qvalTab < 0.01] <- "**"
        qvalTab[qvalTab < 0.1 & qvalTab != "**"] <- "*"
        qvalTab[ qvalTab!="*" & qvalTab!="**"] <- ""
        # plot
        toPlot <-   genecorrPDS.Sprmn_r.cntrd[ genes , ]
        toPlot[is.na(toPlot)]<- 0
        colnames( toPlot ) <- c( paste( "ctrl", names(listSCSN.PDS), sep = "_") ,
                                 names(listSCSN.PDS))
        gplots::heatmap.2(as.matrix(toPlot ), notecol="red", notecex = 2, 
                          cellnote =   qvalTab , 
                          col = rev(brewer.pal(11,"RdBu")) , margins = c(6,6), Colv = F, trace = "none")
        
        ### plot expr. lvls
        genes_mean <- Reduce( cbind, lapply( seq( listSCSN.PDS ) , function(ii ){
          datt <-  listSCSN.PDS[[ii]]@assays$RNA@data[genes, ] 
          datt <- rowMeans( datt )
        }) )
        colnames(genes_mean) <- names(listSCSN.PDS)
        
        library( viridisLite )
        
        gplots::heatmap.2( as.matrix( genes_mean ), scale = "row", 
                           col = viridis , margins = c(6,6), Colv = F, trace = "none", 
                           cellnote =  qvalTab )
      }
     
      
     ### what's the heck, why PDS correlation distributions are shifted from zero
      {
        nnames <-  c("Nphs2","Wt1","Pdss2","btbr","cd2ap","doxo","nephr.D5")
        # nnames <-  names(Wt1Nphs2.4Kcells.PDS)
        
        # make a function for repetative plotting
        plotit<-function(ii,datt){
          p1 <- hist( datt[,ii] , breaks=50)                     # centered at 4
          p2 <- hist( datt[,ii+ncol(datt)/2], breaks=100 )                     # centered at 6
          maxY <- max(p2$counts, max(p1$counts))
          plot( p1, col=rgb(0,0,1,1/4) , ylim=c(0,maxY),
                      xlim= c(-0.3, 0.3), xlab="Spearman rho",
                      main=paste( "Distribution of correlations between PDS and gene lvls.\nlibnorm counts,", 
                                  nnames[ii],"data" , sep = " "))  # first histogram
          plot( p2, col=rgb(1,0,0,1/4) ,  xlim= c(-0.3, 0.3), ylim=c(0,2000) , add=T)  # second
          abline( v = 0 , lwd=5, col="red")
        }
       
        pdf( width = 12 , height = 15 , file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/scRNAseq_GSE146912/raw_counts/PDScorr.pdf")
        par(mfrow=c(7,3))
        pplot <-  lapply( seq(nnames), function(jj){
          ZZ <- function(){plotit( ii=jj, datt=genecorrPDS.Sprmn_r )}
          print(ZZ())
        })
       
        dev.off()
        # plot rho distribution for nephritis
       
        
        # read the expr data
        listSCSN.PDS.sct <- readRDS( file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN.")
        listSCSN.PDS.sct <- readRDS( file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN.PDSlist.SCT_23.12.2022.rda")
       
        # select genes to plot(
        ggenes <- c( "Wt1", "Nphs1")
        
         datt_list <-lapply( seq(listSCSN.PDS.sct) , function(ii){
          datt <- listSCSN.PDS.sct[[ii]]@meta.data[ , c( "gtypeDE" , "PDS.42.005") ]
          datt$countMean <- colMeans( listSCSN.PDS.sct[[ii]]@assays$RNA@counts[ allPodoGenes, ] )
          datt$libnormCountMean <- colMeans( listSCSN.PDS.sct[[ii]]@assays$RNA@data[ allPodoGenes, ] )
          datt$nonzeroN <- colSums( listSCSN.PDS.sct[[ii]]@assays$RNA@counts[ allPodoGenes, ] > 0 )
          datt$sctMean <- colMeans( listSCSN.PDS.sct[[ii]]@assays$SCT@counts )
          datt$sctNormMean <- colMeans( listSCSN.PDS.sct[[ii]]@assays$SCT@data)
          datt$sctScaleMean <- colMeans( listSCSN.PDS.sct[[ii]]@assays$SCT@scale.data)
        
          datt$libnormCountSum <- colSums( listSCSN.PDS.sct[[ii]]@assays$RNA@data[ allPodoGenes, ] )
          datt$libnormCountSum <- colSums( listSCSN.PDS.sct[[ii]]@assays$RNA@data[ allPodoGenes, ] )
          gdat <- t(rbind(  listSCSN.PDS.sct[[ii]]@assays$RNA@counts[ ggenes, ], 
                          listSCSN.PDS.sct[[ii]]@assays$RNA@data[ ggenes, ], 
            listSCSN.PDS.sct[[ii]]@assays$SCT@counts[ ggenes, ], 
                         listSCSN.PDS.sct[[ii]]@assays$SCT@data[ ggenes, ], 
                         listSCSN.PDS.sct[[ii]]@assays$SCT@scale.data[ ggenes, ]) )
          colnames(gdat) <-  paste(colnames(gdat) , rep( c("RNA.counts","RNA.data",
                              "SCT.counts","SCT.data","SCT.scaledData"), each=2) ,
                                     sep = "_")
          datt$dataset <- names(listSCSN.PDS.sct)[ii]
          datt$index <- order( datt$PDS.42.005 )
          datt <- cbind( datt ,gdat  )
          return(datt)
        }) 
        
        library(ggpubr)
        gglist <- lapply( seq( datt_list ) , function(ii){
          toPlot <- datt_list[[ii]]
          ggl <- lapply( seq(2, 10, by = 2), function(jj){
            gg <- ggplot2::ggplot( data = toPlot , aes( x = countMean  , y = toPlot[, ncol(toPlot)-10+jj] )) + 
              geom_jitter(alpha = 0.05)+
              geom_smooth( method = "lm") + 
              ggtitle(paste( names(listSCSN.PDS.sct)[ii], 
                             colnames(toPlot)[ncol(toPlot)-10+jj], sep = " " ) )+
              facet_grid( cols = vars( gtypeDE ) , 
                          scales = "free" , shrink=T ) + 
              stat_cor( size=4 , col="red", method = "spearman") 
            return(gg)
          })
          ggll <- cowplot::plot_grid( plotlist = ggl , ncol = 5)

          return(ggll)
        })
        
        
        cowplot::plot_grid( plotlist = gglist , ncol = 1)
        
      }
    
     
  
    }

#### find features commonly perturbed in FSGS #### 
    Nn <- 7
    
    ### select genes that correlate in N and more cases
    {
      
      
      # all samples
      gene.topCorPDS <- rownames(genecorrPDS.Sprmn_p)[
        rowSums(genecorrPDS.Sprmn_p < sigL, na.rm = T) >= Nn  ]
      # # only experimental samples
      # gene.topCorPDS.xprmnt <- rownames(genecorrPDS.Sprmn_qval)[
      #   rowSums(genecorrPDS.Sprmn_qval[,8:14] < sigL) > Nn ]
      # # only controls
      # gene.topCorPDS.ctrl <- rownames(genecorrPDS.Sprmn_qval)[
      #   rowSums(genecorrPDS.Sprmn_qval[,1:7] < sigL) > Nn ]
      
      
      
      
    }
    
    ### explore relation between the N of identified PDS correlates and the dataset metrics 
    {
      genecorrPDS.Sprmn_sigN <- t( sapply( seq(listSCSN), 
                                           function(ii){
                                             q.wt <- genecorrPDS.Sprmn_q[ , ii ]
                                             q.mut <- genecorrPDS.Sprmn_q[ , ii+9]
                                             as.numeric( c( NsigGenes.ctr=summary( q.wt < sigL )["TRUE"] , 
                                                            NsigGenes.expr=summary( q.mut < sigL )["TRUE"] ,
                                                            Ncells.ctrl= sum( listSCSN[[ii]]$gtypeDE=="control"), 
                                                            Ncells.expr = sum( listSCSN[[ii]]$gtypeDE!="control"),
                                                            PDSspan.ctrl= ( max( listSCSN[[ii]]$PDS[
                                                              listSCSN[[ii]]$gtypeDE=="control" ] ) - min( listSCSN[[ii]]$PDS[
                                                                listSCSN[[ii]]$gtypeDE=="control" ] ) ),
                                                            PDSspan.exp=  max( listSCSN[[ii]]$PDS[
                                                              listSCSN[[ii]]$gtypeDE!="control" ] ) - min( listSCSN[[ii]]$PDS[
                                                                listSCSN[[ii]]$gtypeDE!="control" ] ),
                                                            PDSvar.ctrl= var( listSCSN[[ii]]$PDS[
                                                              listSCSN[[ii]]$gtypeDE=="control" ] ) ,
                                                            PDSvar.exp=  var( listSCSN[[ii]]$PDS[
                                                              listSCSN[[ii]]$gtypeDE!="control" ] ) ) ) 
                                             
                                           }) ) 
      colnames(genecorrPDS.Sprmn_sigN) <- c("NsigGenes.ctr", "NsigGenes.expr", "Ncells.ctrl", 
                                            "Ncells.expr","PDSspan.ctrl", "PDSspan.exp",
                                            "PDSvar.ctrl","PDSvar.exp")
      rownames(genecorrPDS.Sprmn_sigN) <- names(listSCSN)
      
      genecorrPDS.Sprmn_sigN <- genecorrPDS.Sprmn_sigN[ order(genecorrPDS.Sprmn_sigN[,2]),]
   
      # correlate N of correlates and N of cells
      # PDS score span 
      library( ggrepel)
      toPlot <- data.frame( Npds.cor=c(genecorrPDS.Sprmn_sigN[,1], genecorrPDS.Sprmn_sigN[,2]),
                            gtype=c( rep("control",9), rep("experimental",9) ),
                            Ncells_log10=log10( c(genecorrPDS.Sprmn_sigN[,3], genecorrPDS.Sprmn_sigN[,4]) ),
                            Ncells=( c(genecorrPDS.Sprmn_sigN[,3], genecorrPDS.Sprmn_sigN[,4]) ),
                            
                            dataSet = rep( rownames(genecorrPDS.Sprmn_sigN) , 2)  , 
                            PDSspan=c(genecorrPDS.Sprmn_sigN[,5], genecorrPDS.Sprmn_sigN[,6]),
                            PDSvar=c(genecorrPDS.Sprmn_sigN[,7], genecorrPDS.Sprmn_sigN[,8]))

      toPlot[sapply(toPlot, is.infinite)] <- NA
      # plot
      gg0 <- ggplot(toPlot , aes(x=Ncells_log10 , y=log10( Npds.cor), color=gtype)) + 
        geom_jitter( size=3)+         
        coord_cartesian(ylim = c(1,3.25))+
        geom_smooth( method = "lm" )+ 
        stat_cor( size=6 , method = "spearman" ) +
        theme_bw() + geom_label_repel(aes(label = dataSet), size = 5) +
        # ggtitle("How N of cells affects N of PDS correlates") + 
        theme( text = element_text(size=20), legend.position = "none")+
        scale_color_colorblind()
      
      gg2 <- ggplot(toPlot , aes(x= PDSspan, y=log10( Npds.cor), color=gtype)) + 
        geom_jitter( size=3)+
        geom_smooth( method = "lm" )+ 
        stat_cor( size=6 , method = "spearman" ) +
        coord_cartesian(ylim = c(1,3.25))+
        theme_bw() + geom_label_repel(aes(label = dataSet), size = 5) +
        # ggtitle("How PDS span affects N of PDS correlates") + 
        theme( text = element_text(size=20), legend.position = "none")+ 
        scale_color_colorblind()
      
      gg3 <- ggplot(toPlot , aes(x= PDSvar, y=log10( Npds.cor), color=gtype)) + 
        geom_jitter( size=3)+
        geom_smooth( method = "lm" )+ 
        stat_cor( size=6 , method = "spearman" ) +
        theme_bw() + geom_label_repel(aes(label = dataSet), size = 5) +
        coord_cartesian(ylim = c(1,3.25))+
        # ggtitle("How PDS span affects N of PDS correlates") + 
        theme( text = element_text(size=20), legend.position = "none")+ 
        scale_color_colorblind()
      
      ####
      ggl <- cowplot::plot_grid( plotlist = list(gg0, gg2, gg3), nrow = 1) 
      pdf( height = 5, width = 16, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PDScorr/NgnsPDScor.VS.dataMeta_scatterPlot.pdf")
        ggl
      dev.off()
      
     
      ### corr heatmap
      ttt <- toPlot[ !colnames(toPlot) %in% c("dataSet","Ncells_log10")]
      ttt$gtype <- as.numeric(as.factor(ttt$gtype))
      ttt <- ttt[!is.na(ttt$PDSvar),]
      
      gg1 <- corrplot::corrplot( cor(ttt , method = "spearman",
                                     use = "pairwise.complete.obs"),
                                 tl.col =  "black", 
                                 title = "spearman rho")
      
      ## lm
      lms <- glm( Npds.cor~ gtype + Ncells* PDSspan* PDSvar  , data=ttt)
      
      # parcor
      gg2 <- corrplot::corrplot( ppcor::pcor(ttt , method = "spearman")$estimate, 
                          tl.col =  "black", title = "partial spearman rho")
      
      pdf( height = 4, width = 8, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PDScorr/NgnsPDScor.VS.dataMeta_corrHeatmap.pdf")
      par(mfrow=c(1,2))
      corrplot::corrplot( cor(ttt , method = "spearman",
                              use = "pairwise.complete.obs"),
                          tl.col =  "black", 
                          title = "spearman rho")
      corrplot::corrplot( ppcor::pcor(ttt , method = "spearman")$estimate, 
                          tl.col =  "black", title = "partial spearman rho")
      dev.off()
      
    }
    
    ### plot count matrix of sig.correlations in Nn and more studies 
    {
      library(scales)
      library(ComplexHeatmap)
      
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
      
      # save 
      pdf(height = 4, width = 12, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PDScorr_count.heatmap.05.01.24.7plus.pdf")
      
      ComplexHeatmap::Heatmap( 
        t(as.matrix(toPlot)),
        column_labels =rep("", nrow(toPlot))  ,
        cluster_rows  = F , 
        # margins= c(1, 8),
        col = viridis(max(toPlot)+1),   #labCol = FALSE,
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "complete" , 
        column_split=cclust ,
        top_annotation= ha, 
        name = "count heatmap of sc/sn FSGS models\nwhere a gene correlates with PDS"
      )
      dev.off()
      
      ## save as a table
      write.table( cbind( toPlot ,clustID=cclust, clustName=sapply(cclust, plotrix::color.id), 
                          DS= ifelse(rownames(toPlot)%in%DS_all$gene_symbol[1:42],
                                     "DS.top42", ifelse( rownames(toPlot) %in% 
                                                           DS_all$gene_symbol[43:nrow(DS_all)],
                                                         "DS","")) ), sep = "\t",
                   file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PDScorr_count.heatmap.05.01.24.7plus.tsv")
    }
    
    ### sic! use Roberts function for nonredundant GO annotations 
    
    ### functional annotation with clusterprofile
    {    
      
      # convert background to entrez IDs
      gene.bckgrnd_eID <- tx2gene$entrezgene_id[ match( allPodoGenes, tx2gene$external_gene_name)]
      gene.bckgrnd_eID <- as.character( gene.bckgrnd_eID[!is.na(gene.bckgrnd_eID)] )
      
      ## to test
      library(plotrix)
      library( "reactome.db" , lib.loc = "/media/tim_nevelsk/WD_tim/SOFT/R")
      library( "ReactomePA" , lib.loc = "/media/tim_nevelsk/WD_tim/SOFT/R")
      
      toTest <-  as.data.frame( cclust)
      toTest$eID <- tx2gene$entrezgene_id[ match( rownames(toTest), tx2gene$external_gene_name)]
      toTest$gName  <- rownames(toTest)
      toTestGSET.eID <- lapply(seq( length(unique(cclust))), function(ii){
        NN <-  toTest$eID[toTest$cclust==levels(cclust)[ii]]
        NN <- unique( NN[!is.na(NN)] )
      }) 
      toTestGSET <- lapply(seq( length(unique(cclust))), function(ii){
        NN <- toTest$gName[toTest$cclust==levels(cclust)[ii]]
        NN <- unique( NN[!is.na(NN)] )
      }) 
      names(toTestGSET) <- names(toTestGSET.eID) <- sapply(levels(cclust), plotrix::color.id)
      
      
      # annotate
      gene.topCorPDS_clustAnnot <- lapply( seq( toTestGSET.eID ), 
                                           function(ii){
                                             print(ii)
                                             ggenes.eID <- toTestGSET.eID[[ii]]
                                             ### GO annotations 
                                             cProfiler.GKR( ggenes.eID=ggenes.eID , 
                                                            gene.bckgrnd_eID=gene.bckgrnd_eID)
                                             
                                           })
      
      saveRDS(gene.topCorPDS_clustAnnot, file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PDScorr/gene.topCorPDS_clustAnnot.rda")
      
      
      # make a df for visualising
      toPlot <- Reduce( rbind, lapply( 1:4, function(ii){
        toPlot <- gene.topCorPDS_clustAnnot[[ii]]
        toPlot <- Reduce( rbind , lapply(toPlot, function(X) X@result))
        toPlot$cluster <-   names(toTestGSET)[ii]
        return(toPlot)
        
      }))
      # toPlotGO <- toPlot[ toPlot$ID %in% toPlot$ID[toPlot$qvalue<0.01] &
      #                       toPlot$ONTOLOGY%in%c("BP","MF","CC"), ]
      # toPlotGO$Description <- sapply( toPlotGO$Description , wrap_text, 40)
      toPlotPTH <- unique( toPlot[(toPlot$pvalue< 0.01 & toPlot$ONTOLOGY%in%c("KEGG","REACT")) | 
                                    (toPlot$qvalue< 0.01 & !toPlot$ONTOLOGY%in%c("KEGG","REACT")), ] )
      
      toPlotPTH$Description <- sapply( toPlotPTH$Description , substring, 1, 50)
      toPlotPTH <- toPlotPTH[toPlotPTH$Count>1,]
      # toPlotPTH$Description <- sapply( toPlotPTH$Description , wrap_text, 40)
      
      # visualise
      ggplot2::ggplot(data=toPlotPTH, aes(
        x=-log10(pvalue), 
        y=reorder(Description, -log10(pvalue))))+
        geom_bar(stat="identity") + 
        # geom_vline(xintercept = 1, linetype="dashed", color="red")+
        theme_bw() + facet_grid( rows=vars(ONTOLOGY), 
                                 cols = vars(factor(toPlotPTH$cluster)), 
                                 scales = "free_y", space = "free" ) + theme( text = element_text(size=20)) + 
        ylab("pathway names")
      
      
      ## cluster terms
      pplist <-  lapply( seq(gene.topCorPDS_clustAnnot), function(ii){
        print(ii)
        datt <- gene.topCorPDS_clustAnnot[[ii]]
        
        p1 <- enrichplot::emapplot_cluster(enrichplot::pairwise_termsim(
          datt$GO.enrich), color = "pvalue", cex_label_group=1.2)
        # p3 <-  enrichplot::emapplot_cluster(enrichplot::pairwise_termsim(
        #   datt$KEGG.enrich), color = "pvalue")
        p2  <- enrichplot::emapplot_cluster(enrichplot::pairwise_termsim(
          datt$REACT.enrich), color = "pvalue", cex_label_group=1.2)
        cowplot::plot_grid( plotlist =  list(p1,p2), nrow = 1, labels = c("GO","REACT"))
        
      })
      cowplot::plot_grid( plotlist =  pplist , ncol = 1, labels =   names(toTestGSET) )
      
      
    }
    

    
    # #### GSEA plot
    #   library( fgsea)
    #   allPodoGenes_mean <- readRDS(  file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/SCSN_allPodoGenes_mean.rda")
    #   
    #   gene.topCorPDS <- rownames(genecorrPDS.Sprmn_qval)[rowSums(
    #     genecorrPDS.Sprmn_qval < sigL) > 2 ]
    #   
    #   plotEnrichmentPDS( ticksSize=0.2, pathway= (gene.topCorPDS),  stats = rowMeans(allPodoGenes_mean)) +
    #     ggtitle("enrichemtn of podocyte genes, ranked by the expr.lvls.,\nin 61 genes correlating with PDS (qval<0.01) in >= 3 models")
    #  
    
    ### study features unique to specific FSGS models
    {
      ### plot selected gene correlations with PDS (rows) over models (columns)
      sigL <- 0.01
      
      
      ### all unique
      FSGS.unqPDS.01.list <- lapply( seq(listSCSN.PDS), function(ii){
        rownames(genecorrPDS.Sprmn_qval)[ genecorrPDS.Sprmn_qval[,7+ii] < 0.01 & 
                                            ( rowSums(genecorrPDS.Sprmn_qval[,8:14] < sigL ) ==1) ]
      })
      names(FSGS.unqPDS.01.list)  <- names(listSCSN.PDS)
      
      ### select unique by ranks
      genecorrPDS.Sprmn_r.rank <- apply( abs( genecorrPDS.Sprmn_r[,8:14]) , 2, rank )
      # difference between the gene rank in the model and the average gene rank across other models
      genecorrPDS.Sprmn_r.rank.diff <- sapply( 1:7 , function(ii){
        genecorrPDS.Sprmn_r.rank[,ii] - rowMeans( genecorrPDS.Sprmn_r.rank[,-ii] )
      } )
      # top 100 with the largest positive rank differnce
      genecorrPDS.Sprmn_r.rank.5K <- lapply(1:7, function(ii){
        X <- genecorrPDS.Sprmn_r.rank.diff[,ii]
        X <- names(X)[X>5000]
      })
      names(genecorrPDS.Sprmn_r.rank.5K) <- names(listSCSN.PDS)
      
      ### functional annotation by topGO
      {

        ### functional annotation by Robert
        # run function on list of DE results
        
        # use Robert's function to create sparse binary matrix indicating membership of genes (columns) in GO terms (rows)
        gomatrix=sf.createGoMatrix()
        # load biomart
        mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
        tx2gene_GO <- biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id"), mart = mart)
        
        uniqueBYmodel.UPSTRM_robertGO <- lapply(seq(uniqueBYmodel.UPSTRM), 
                                                function(ii){
                                                  print( ii )
                                                  
                                                  # define background gene set
                                                  universe <- unique( tx2gene_GO$entrezgene_id[ 
                                                    tx2gene_GO$external_gene_name %in% allPodoGenes ] ) 
                                                  universe <- as.character( universe[!is.na(universe)] )
                                                  
                                                  
                                                  # prepare gene set, convert ensembleIDs to entrezIDs
                                                  geneset <- uniqueBYmodel.UPSTRM[[ii]]
                                                  geneset <- unique(tx2gene_GO$entrezgene_id[
                                                    tx2gene_GO$external_gene_name %in% geneset ])
                                                  geneset <- geneset[!is.na(geneset)]
                                                  geneset <- as.character(geneset)
                                                  print(length(intersect(geneset,colnames(gomatrix))))
                                                  if( length(geneset)==0) return(NA)
                                                  
                                                  # apply Robert's function that given a sparse matrix of GO terms (columns = genes, rows = GO terms)
                                                  # a geneset of interest and a background set of genes (universe)
                                                  # carry out clustering with members diverging by at most cut_max genes, and do enrichment testing.
                                                  # Note, multiplicity adjustment is performed for the representative terms only.
                                                  RobertGO=sf.clusterGoByGeneset( gomatrix, geneset, universe, 
                                                                                  min.genes=5, cut_max = 5 )
                                                  
                                                  return(RobertGO)
                                                })
        saveRDS(uniqueBYmodel.UPSTRM_robertGO, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/crosscorr/uniqueBYmodel.UPSTRM_robertGO.v2.rda")
        
        
        
        # make plots
        
      }
      
      ### functional annotation by clusterprofile
      {
        # convert background to entrez IDs
        gene.bckgrnd_eID <- tx2gene$entrezgene_id[ match( allPodoGenes, tx2gene$external_gene_name)]
        gene.bckgrnd_eID <- as.character( gene.bckgrnd_eID[!is.na(gene.bckgrnd_eID)] )
        
        ## to test
        toTestGSET <- genecorrPDS.Sprmn_r.rank.5K
        # remove gene sets with less than 50 genes
        toTestGSET <- toTestGSET[ lengths(toTestGSET)>50]
        # convert to entrez IDs
        toTestGSET <- lapply(toTestGSET, function(X){
          
          X<- tx2gene$entrezgene_id[ match( X , tx2gene$external_gene_name)]
          X <- X[!is.na(X)]
        })
        
        
        ## annotate
        gene.UniCorPDS_clustAnnot <- lapply( seq( toTestGSET ), 
                                             function(ii){
                                               print(ii)
                                               ggenes.eID <- toTestGSET[[ii]]
                                               ###  annotations 
                                               cProfiler.GKR( ggenes.eID=ggenes.eID , 
                                                              gene.bckgrnd_eID=gene.bckgrnd_eID)
                                             })
        
        
        # make a df for visualizing
        toPlot <- Reduce( rbind, lapply( seq(gene.UniCorPDS_clustAnnot), function(ii){
          toPlot <- gene.UniCorPDS_clustAnnot[[ii]]
          toPlot <- Reduce( rbind , lapply(toPlot, function(X) X@result))
          toPlot$cluster <-   names(toTestGSET)[ii]
          return(toPlot)
          
        }))
        toPlotGO <- toPlot[ toPlot$ID %in% toPlot$ID[toPlot$qvalue<0.05] &
                              toPlot$ONTOLOGY%in%c("BP","MF","CC"), ]
        toPlotGO$Description <- sapply( toPlotGO$Description , substring, 1, 50)
        toPlotPTH <-  toPlot[ toPlot$ID %in% toPlot$ID[toPlot$qvalue<0.05 &
                                                         toPlot$ONTOLOGY%in%c("KEGG","REACT")], ]
        toPlotPTH$Description <- sapply( toPlotPTH$Description , substring, 1, 50)
        
        # toPlotPTH$Description <- sapply( toPlotPTH$Description , wrap_text, 40)
        
        # visualise
        ggplot2::ggplot(data=toPlotGO, aes(x=-log10(qvalue), y=reorder(Description, -log10(pvalue))))+
          geom_bar(stat="identity") + 
          geom_vline(xintercept = 1, linetype="dashed", color="red")+
          theme_bw() + facet_grid( rows=vars(ONTOLOGY), 
                                   cols = vars(cluster), 
                                   scales = "free", space = "free_y" ) + theme( text = element_text(size=20)) + 
          ylab("pathway names")
        
      }
      
      
    }
  
  
#### plot correlation of PA/genes with PDS on pathway diagrams ####
  {
    # select genes that are significantly expressed in control or 
    # experimental samples of at least 2 models
    sigL <- 0.01
    
    gene.topCorPDS_eID <- tx2gene$entrezgene_id[ match( gene.topCorPDS, tx2gene$external_gene_name)]
    gene.bckgrnd_eID <- tx2gene$entrezgene_id[ match( allPodoGenes, tx2gene$external_gene_name)]
    
    ### data to plot on Net
    toPlot.all <- genecorrPDS.Sprmn_r.cntrd
    # combine all controls and all exprmnt in indiv. columns
    toPlot.all <- cbind(  ctrl_mean=rowMeans(toPlot.all[,1:7], na.rm = T), 
                          toPlot.all[,8:14])
    
    colnames(toPlot.all) <- c( "ctrl_mean","Nphs2","Wt1","Pdss2" ,      
                               "Btbr","Cd2ap","doxo","nephr.day5")
    

    ### visualise Reactome
    {
      
      
      library(SBGNview)
      data("sbgn.xmls")
      setwd("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/reactome/13.10.23/")
      ensembl.pathway <- sbgn.gsets( id.type = "ENSEMBL",
                                     species = "mmu",
                                     mol.type = "gene",
                                     truncate.name.length = 500,
                                     output.pathway.name = TRUE )
      
      ## prepare gene info 
      toPlotR <- toPlot.all
      rownames(toPlotR) <- tx2gene$ensembl_gene_id[ match(   rownames(toPlotR)  ,  tx2gene$external_gene_name)]
      # toPlot<- t(c(-0.99,-0.75,-0.5,0.5,0.75,0.99) * t(toPlot))
      
      # # selec paths
      ppname <- "RHO GTPases activate IQGAPs"
      pName <- pathways.info[ pathways.info$pathway.name %in% ppname, ]
      # pName <-  findPathways(c("cytoskeletal"))
      
      lapply(seq( pName$pathway.id), function(ii){
        require(SBGNview)
        
        nname <- gsub(" |\\(|\\)|/|:", "_", pName$pathway.name[ii])
        SBGNview.obj <- SBGNview( gene.data = toPlotR ,
                                  # input.sbgn =  pathways ,
                                  input.sbgn = pName$pathway.id[ii] ,
                                  gene.id.type = "ENSEMBL",
                                  output.formats =  c("png", "pdf"),
                                  org = "mmu" , 
                                  output.file =paste( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/reactome/13.10.23/", 
                                                      nname, sep = ""),
                                  
                                  node.sum = "sum",
                                  min.gene.value = -0.2,
                                  max.gene.value = 0.2 ,
                                  
                                  font.size = 2,
                                  text.length.factor.complex = 3,
                                  if.scale.compartment.font.size = TRUE,
                                  node.width.adjust.factor.compartment = 0.04 )+
          highlightNodes(select.glyph.class = "macromolecule",
                         stroke.width = 4, 
                         stroke.color = "green")
        print(SBGNview.obj)
      })
      
      
    }
    
    ### visualise KEGG
    {
      
      
      library( pathview)
      library(KEGG.db)
      keggid2keggname <- as.list(KEGGPATHID2NAME)
      # mmu04510 , mmu04512, 04810, 04370, 
      # mmu04520 - adherens junction
      # mmu04060
      # toPlot <- KEGG.fsgs_tab$Term[KEGG.fsgs_tab$Freq>1]
      setwd("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/kegg/05.01.24/")
      # pathsPlot <- "Calcium signaling pathway"
      pathsPlot <- c( "ECM-receptor interaction" ,
                      "Focal adhesion",
                      "Regulation of actin cytoskeleton",
                      "Axon guidance")
      
      # ### annotate
      # nphs2wt1.spec <- rownames(toPlot)[ ( rowSums(genecorrPDS.Sprmn.Max_qval[, c(1,2,10)]> sigL)==3 & 
      #                                        genecorrPDS.Sprmn.Max_qval[, 9 ]< sigL) | 
      #                                      (rowSums(genecorrPDS.Sprmn.Max_qval[, c(1,2,9)]> sigL)==3 & 
      #                                         genecorrPDS.Sprmn.Max_qval[, 10 ]< sigL)  ]
      # nphs2wt1.spec.ENTRZ <- tx2gene$entrezgene_id[ match( nphs2wt1.spec, tx2gene$external_gene_name)]
      # 
      # PDS.0.01.nphs2wt1.spec_KEGG.enrich <-  clusterProfiler::enrichKEGG( gene = nphs2wt1.spec.ENTRZ,
      #                                                                     organism     = 'mmu',
      #                                                                     pvalueCutoff = 0.05 )@result
      
      # plot!
      KEGGgraph_fsgs <- lapply( seq(pathsPlot), function(ii , expr=toPlot.all )
      {
        require(pathview)
        IDD <- pathsPlot[[ii]] 
        IDD <- names(keggid2keggname)[which(keggid2keggname==IDD)]
        nname1 <- paste( "_Mctrl.Mexpr.9expr_All.sum",gsub( " ", "_", pathsPlot[ii]) , sep ="_")
        pathview( gene.data = expr , pathway.id = IDD, node.sum="sum",
                  gene.idtype="SYMBOL",kegg.dir="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/kegg/13.10.23/",
                  species = "mmu", out.suffix = nname1 , kegg.native = F,
                  expand.node=F, limit = list(gene = 0.2, cpd = 1) )
        pathview(gene.data = expr ,pathway.id = IDD ,node.sum="sum",
                 gene.idtype="SYMBOL",kegg.dir="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/kegg/13.10.23/", 
                 species = "mmu", out.suffix = nname1 , kegg.native = T ,
                 limit = list(gene = 0.2, cpd = 1))
        
        # nname2 <- paste( "exprmnt.500",sigL, gsub( " ", "_", pathsPlot[ii]) , sep ="_")
        # pathview(gene.data = toPlot[, 9:16] , pathway.id = IDD,
        #          gene.idtype="SYMBOL",kegg.dir="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/crosscorr/pathway_viz/KEGG/500_0.01/",
        #          species = "mmu", out.suffix = nname2 , kegg.native = F,
        #          expand.node=F, limit = list(gene = 0.2, cpd = 1))
        # pathview(gene.data = toPlot[, 9:16] ,pathway.id = IDD ,
        #          gene.idtype="SYMBOL",kegg.dir="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/crosscorr/pathway_viz/KEGG/500_0.01/",
        #          species = "mmu", out.suffix = nname2 , kegg.native = T ,
        #          limit = list(gene = 0.2, cpd = 1))
      })
      
      
      
    }
    
    
    ### combine corr tables with networks and TF info for cytoscape
    {
      Net <- read.table(sep = "\t", header = T, file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/crosscorr/pathway_viz/cytoscapeVIZ/C.shell/ECMnet.tsv")
      Net.nodes <- unique( union(Net$Symbol, Net$Interactors))
      Net.nodes <- Net.nodes[Net.nodes!="" & !is.na(Net.nodes)]
      PDScorr <- toPlot.all[ match( Net.nodes , rownames(toPlot.all) ) ,  ]
      rownames(PDScorr) <- Net.nodes
      # include TFtget ATACseq net
      TF.01 <- readRDS(file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/TFA/TFs_SNSC.podo.0.01.rda")
      
      TFreg <- ATACseq_tgenesM[ match(Net.nodes , rownames(ATACseq_tgenesM) ),
                                colnames(ATACseq_tgenesM)%in%TF.01 ] 
      rownames(TFreg) <- Net.nodes
      TFreg <- TFreg[ , colSums(TFreg, na.rm = T)>0 ]
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
    
  }
  
#### cluster features using "barcode" #### 
  
  PDS42.SpCor_barcode <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PDS42.SpCor_barcode.03.05.24.rda")
  
    # select "pathway genes"
    # acyskl <- read.table( sep = "\t", header = T,"/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/crosscorr/pathway_viz/cytoscapeVIZ/PPIsource/PMID:33514561_SupplTab3.4_ActinCytoSkeleton.csv")
    # ggenes <- acyskl$name
    
    # # select genes that correlate with PDS
    # ggenes <- (rownames(genecorrPDS.Sprmn_freq)[(genecorrPDS.Sprmn_freq$total > 4)])
    # # only wt1 and nphs2
    
    # ### select pathway genes
    # # KEGG
    # ggenes <- All_kegg$Gene.symbol[All_kegg$pName=="Regulation of actin cytoskeleton"]
    # # reactome
    # pt <- read.table("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/crosscorr/pathway_viz/Participating_Molecules_[R-MMU-2022090].tsv", fill = T, header = T)
    # ggenes <- tx2prot$mgi_symbol[ match( pt$Identifier, tx2prot$uniprotswissprot)]
    
    ### tab
    sigL <- 0.1
    
    ggenes <- rownames(genecorrPDS.Sprmn_q)[rowSums(genecorrPDS.Sprmn_q < sigL, na.rm = T)>3 ]
    # # ggenes <- ggenes[!ggenes%in% DS_all$gene_symbol[1:42]]
    # toPlot <- toPlot.all[ rownames(toPlot.all)%in% ggenes, 3:4]
    # toPlot[ is.na(toPlot)]<-0
    
    # plot vars
    myColor <- rev( RColorBrewer::brewer.pal(n = 11, name = "RdBu") )
    mypal <- viridis( 2 )
    # rowcol <- map2color( rowMeans(allPodoGenes_mean[ rownames(toPlot),]),  mypal) 
    colcol <- ifelse( genecorrPDS.annot$gtypeDE[ 
      genecorrPDS.annot$sample%in% colnames(toPlot)] == "control", mypal[1] ,mypal[2]) 
    mybreaks <- seq(-0.3,0.3,length=12)
    # make meaningful columnames
    colnames(toPlot) <- sub("SID","",colnames(toPlot))
    colnames(toPlot)[colnames(toPlot)%in% annot_tab$CCG_Sample_ID]<- paste( 
      "KFO",sep = "_", annot_tab$group[match(  colnames(toPlot)[colnames(toPlot)%in% annot_tab$CCG_Sample_ID],  annot_tab$CCG_Sample_ID)])
    
    
    cclust <-   as.factor( cutree( hclust(dist((as.matrix(toPlot))),
                                                    method="ward.D2"), k = 7) )
    levels(cclust) <- hue_pal()(7)
    
    
    # plot
    
    toPlot <- as.matrix( genecorrPDS.Sprmn_r.cntrd[ 
      ggenes, !colnames(genecorrPDS.Sprmn_r.cntrd) %in% 
        c("nephr.D1_ctrl","nephr.D5_ctrl")] )
    toPlot[is.na(toPlot)]<- 0
    gg<- gplots::heatmap.2( toPlot,  
                            plot.col.key=T,
                            Colv = F,
                       na.color = "darkgrey", 
                       breaks = mybreaks,
                       # ColSideColors  = colcol ,
                       # RowSideColors  = as.character( cclust ),
                       hclustfun =function(x) hclust(x, method="ward.D2"),
                       col = myColor,  trace = "none",
                       # labRow = FALSE,
                       cexCol=1 , 
                       margins=c(10,5) , dendrogram = "both"
                       # colRow = ifelse( rownames(toPlot)%in% DS_all$gene_symbol[1:42],"red" ,"black") 
    )
    
    # legend("topright", title = "5ean.expr",legend=c(min( rowMeans(allPodoGenes_mean[ rownames(toPlot),])),max( rowMeans(allPodoGenes_mean[ rownames(toPlot),]))), 
    #        fill= viridis( 2 ), cex=0.8, box.lty=0)
    
    # for paths
    hplots <- lapply(seq(ppname), function(ii){
      library(gplots)
      
      # select gene names
      pthw <- ppname[ii]
      ggenes <- reactPath_gName[[pthw]]
      
      # tab
      toPlot <- cbind(rowMeans(genecorrPDS.Sprmn.Max_r[,1:8], na.rm = T), genecorrPDS.Sprmn.Max_r[,8:14])
      toPlot <- toPlot[ ( rownames(toPlot) ) %in% ggenes,]
      colnames(toPlot) <- c("mean wild-type", gsub("podo_|sn_|cem.","" , names(SCSNdata_list_sub)))
      llables <- genecorrPDS.Sprmn_qval[( rownames(genecorrPDS.Sprmn_qval)) %in% ggenes, 8:14]
      llables <- ifelse( llables < 0.01 , "**", ifelse(llables<0.1, "*",""))
      llables <- cbind( "", llables)
      
      
      qq <- tryCatch( gplots::heatmap.2(toPlot, col = redblue(31),Colv = F, trace = "none", 
                                        margins=c(10,5), 
                                        cellnote= llables,
                                        notecex=1.2,notecol="black",
                                        main = pthw), error = function(e) NA)
      
      
      return(qq)
      
    })
    
    
  
#### estimate pathway activities with AUCell ####
  {
    # prepare reactorme and kegg gsets
    source( "/home/tim_nevelsk/PROJECTS/myCode/func_analysis.R")
    library(GSEABase)
    keggPath_pathlist <- readRDS("/media/tim_nevelsk/WD_tim/ANNOTATIONS/Pathways/kegg_pathsList_gName.07.09.23.rda")
    reactPath_pathlist <- readRDS("/media/tim_nevelsk/WD_tim/ANNOTATIONS/Pathways/reactome_pathsList_gName.07.09.23.rda")
    PodoPathGSet_pathlist <- readRDS(file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PodoPathGSet.05.01.24.7plus.rda")
      
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
      
    

#### visualise activity curves ####
  # PAaucell_reactome <- readRDS(file= "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PathwayActivity/PAaucell_list.reactome.rda") 

  
### smooth expr./activity curves
  {
    TFlist <- c( "Klf4", "Notch1", "Foxc1", "Foxc2", "Wt1", "Mafb", "Tcf21", "Lmx1b" )
    
    PAaucell_TFs<- readRDS(file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PathwayActivity/PAaucell_list.TFs.rda")
    # listSCSN.PDS_1K <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN.PDSlimit.1K.Test.rda")
    listSCSN.1K.sampl <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN_samples.1K.22.12.23.rda")
    TF_MeanMed <- readRDS(  file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/TFA/TFs_SNSC.podo.MeanMed.08.01.24.rda")
    
    ## calculate smoothed activity profile of mRNA lvls and make a DF for plotting
    iids <- c(1:6,8)
    geneExpr_loess.1K <- Reduce( rbind , lapply( iids , 
                                                 function(ii, datt= listSCSN.1K.sampl , 
                                                          ggenes = TFlist[ TFlist%in%TF_MeanMed], 
                                                          sspan=0.3, scale.arg=F )
                                                 {
                                                   print(ii)
                                                   ## select TFs and order cells by PDS
                                                     expr <- datt[[ii]]@assays$RNA@data
                                                   expr <- expr [ match( ggenes , rownames( expr ) ) ,  ]
                                                   
                                                   expr[is.na(expr)]<- 0
                                                   rownames( expr) <- ggenes
                                                   ## separate wt and mut cells and  limit to span -0.5 to 0 PDS
                                                   # balance groups
                                                   col.wt <- intersect( names(  which(datt[[ii]]$gtypeDE=="control") )  ,
                                                                        names( which( datt[[ii]]$PDS < - 0.2  &
                                                                                        datt[[ii]]$PDS > -0.4 ) ) )
                                                   
                                                   col.mut <- intersect( names( which(datt[[ii]]$gtypeDE!="control") ) ,
                                                                         names( which( datt[[ii]]$PDS <  0.2 &
                                                                                         datt[[ii]]$PDS > -0.40 ) ) )
                                                   #                  
                                                   print( length( col.wt))
                                                   print( length( col.mut))
                                                   
                                                   if( isTRUE(scale.arg)){
                                                     # scale
                                                     datt.wt <- scale( t( as.matrix(expr[, col.wt] )), center = F )
                                                     datt.wt[is.na(datt.wt)]<- 0
                                                     
                                                     datt.mut <- scale( t( as.matrix(expr[, col.mut] )) , center = F)
                                                     datt.mut[is.na(datt.mut)]<- 0
                                                   } else{
                                                     datt.wt <- t( as.matrix(expr[, col.wt] ))
                                                     
                                                     datt.mut <- t( as.matrix(expr[, col.mut] ))
                                                   }
                                                   
                                                   
                                                   
                                                   # add damage scores
                                                   datt.wt <- cbind.data.frame(  
                                                     PDS= datt[[ii]]$PDS[col.wt],
                                                     (datt.wt) )
                                                   datt.mut <- cbind.data.frame( 
                                                     PDS= datt[[ii]]$PDS[col.mut],
                                                     (datt.mut)  )
                                                   
                                                   
                                                   # # add index
                                                   datt.wt <- cbind.data.frame( index=seq(nrow(datt.wt)) , datt.wt)
                                                   datt.mut <- cbind.data.frame( index=seq(nrow(datt.mut)) , datt.mut)
                                                   
                                                   # smoothen
                                                   loessRes.wt <- Reduce( cbind, lapply( (ncol(datt.wt)-length(ggenes)+1):ncol(datt.wt) ,
                                                                                         function(jj){
                                                                                           # print(jj)
                                                                                           
                                                                                           loessMod1 <- loess(  datt.wt[,jj] ~ datt.wt[,"index"], data = datt.wt, span = sspan)
                                                                                           loessMod2 <- loess( datt.wt[,jj] ~ datt.wt[,"PDS"], data = datt.wt, span = sspan)

                                                                                           smoothed <- c(  range01(predict(loessMod1)),
                                                                                                           range01(predict(loessMod2))  )
                                                                                         }))
                                                   
                                                   loessRes.mut <- Reduce( cbind , lapply( (ncol(datt.mut)-length(ggenes)+1):ncol(datt.mut),
                                                                                           function(jj){
                                                                                             # print(jj)
                                                                                             loessMod1 <- loess( datt.mut[,jj] ~ datt.mut[,"index"], data = datt.mut, span = sspan)
                                                                                             loessMod2 <- loess( datt.mut[,jj] ~ datt.mut[,"PDS"], data = datt.mut, span = sspan)

                                                                                             smoothed <- c( range01( predict(loessMod1)),
                                                                                                            range01( predict(loessMod2)) )
                                                                                           }) )
                                                   
                                                   # ad feature names
                                                   colnames(loessRes.wt ) <-   colnames(loessRes.mut ) <- ggenes
                                                   
                                                   # combine with metadata
                                                   loessRes.wt <- cbind.data.frame( 
                                                     sample= datt[[ii]]$sample[col.wt],
                                                     response= c( datt.wt[,"index"] ,
                                                                  datt.wt[,"PDS"] ) ,
                                                     rType = c( rep("index",nrow(datt.wt)),
                                                                rep("PDS",nrow(datt.wt))) ,
                                                     loessRes.wt)
                                                   
                                                   loessRes.mut <- cbind.data.frame( 
                                                     sample= datt[[ii]]$sample[col.mut],
                                                     response= c( datt.mut[,"index"] ,
                                                                  datt.mut[,"PDS"] ) ,
                                                     rType = c( rep("index",nrow(datt.mut)),
                                                                rep("PDS",nrow(datt.mut))) ,
                                                     loessRes.mut)
                                                   
                                                   # combine wt and mut
                                                   toPlot<- cbind( gtypeDE= c( rep( "control", nrow(loessRes.wt) ), 
                                                                               rep( "experiment", nrow(loessRes.mut) )),
                                                                   dataSet=names(datt)[ii] , 
                                                                   rbind( loessRes.wt , loessRes.mut) )
                                                   return( toPlot )
                                                 }) )
    
    # saveRDS(TFmRNA_loess.1K, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/TFA/TFmRNA_PDS42_loess1.RDS")
    
    # melt for plotting
    geneExpr_loess.1K.melt <- reshape2::melt( geneExpr_loess.1K ,
                                              id=c("gtypeDE", "dataSet","sample", "response", "rType"))
    colnames( geneExpr_loess.1K.melt )[ colnames(geneExpr_loess.1K.melt)=="value" ] <- "predicted"
    
      }

### TF activity profiles
  {
    
  ## plot selected genes for many studies
    TFset <- c( "Klf4", "Notch1", "Foxc1", "Foxc2", "Wt1", "Mafb", "Tcf21", "Lmx1b" )
    
  datTOplot <- geneExpr_loess.1K.melt
  datTOplot <- datTOplot[datTOplot$variable%in%TFset &
                           datTOplot$rType=="PDS",]
  ggplot2::ggplot( data = datTOplot , 
                   aes(x=response , 
                       color=variable ,
                       y=predicted) )+
      geom_line( lwd=1  )+ 
    geom_jitter( ) +
    facet_grid( rows = vars(dataSet), cols = vars(gtypeDE),
                scales = "free", space="free") + theme_bw() 
  
  ## indiv samples
  datTOplot <- TFmRNA_loess.smpls.melt
  datTOplot <- datTOplot[datTOplot$variable%in%TFset &
                           datTOplot$rType=="PDS.42" ,]
  ## select samples where TFs sig cor with PDS
  # TFsigsamp<- colnames(genecorrPDS.Sprmn_qval)[colSums(genecorrPDS.Sprmn_qval[TFset,]<0.05)==2]
  # datTOplot <- datTOplot[datTOplot$sample %in% TFsigsamp,]
  # make descriptive sample names
  datTOplot$sample <- sub("SID","",datTOplot$sample)
  
  datTOplot$sample[datTOplot$sample%in% annot_tab$CCG_Sample_ID]<- paste( 
    datTOplot$sample[datTOplot$sample%in% annot_tab$CCG_Sample_ID],sep = "_", annot_tab$group[match(  datTOplot$sample[datTOplot$sample%in% annot_tab$CCG_Sample_ID],  annot_tab$CCG_Sample_ID)])
  
  
  # plot
  gglist<- lapply( seq(unique(datTOplot$sample)), function(ii){
    toPlot <- datTOplot[ datTOplot$sample== unique(datTOplot$sample)[ii] , ] 
    gg <- ggplot2::ggplot( data = toPlot , aes(x=response , y=predicted, 
                                            color=variable) )+
      geom_line( lwd=2  )+ ggtitle( unique(datTOplot$sample)[ii] )+
      # geom_jitter( ) +
      #  facet_grid( rows = vars(variable), cols = vars(dataSet),scales = "free", space="free") + 
      theme_bw() 
    if(ii!=length(unique(datTOplot$sample))) gg<- gg+ theme(legend.position = "none")
    return(gg)
  })
  
  pdf(width = 28, height = 24, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFmRNA_loess.KFOsmpl.pdf")
      cowplot::plot_grid( plotlist = gglist , nrow = 7)
  dev.off()
  
  

  
  
}

### test for Granger causality
  {
    library(tseries)
    library(bruceR)
    
    TFmRNA_loess.smpls <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/TFA/TFmRNA_loess1.smpl.RDS")
    ### do for one sample and a subset of TFs 
    ### to test if granger is appropriate and choose parameters
    {
      ggenes <- c("Wt1","Zbtb20","Lmx1b")
      expr <- t(listSCSN.PDS[[1]]@assays$RNA@data[ggenes, 
                                                  listSCSN.PDS[[1]]$sample=="143489"]  )
      expr <- expr[ order(listSCSN.PDS[[1]]$PDS.42[ listSCSN.PDS[[1]]$sample=="143489"]), ]
      # loess
      expr.loess <- TFmRNA_loess.smpls[TFmRNA_loess.smpls$rType=="PDS.42" & 
                                         TFmRNA_loess.smpls$sample=="143489",]
      expr.loess <- expr.loess[ order(expr.loess$response), ggenes]
      
      # maw 
      WsizeF=100
      WstepF=100
      expr.maw <- apply( expr , 2, slideFunct, window =nrow(expr)/WsizeF ,
                         step =nrow(expr)/WstepF )
      expr.maw.diff <- apply(expr.maw, 2, diff)
      ### test for if time series are stationary
      Reduce( c, lapply( 1:ncol(expr) , 
                         function(ii) adf.test(expr[,ii])$p.value))
      Reduce( c, lapply( 1:ncol(expr.loess) , 
                         function(ii) adf.test(expr.loess[,ii])$p.value))
      Reduce( c, lapply( 1:ncol(expr.maw) , 
                         function(ii) adf.test(expr.maw[,ii])$p.value))
  
      
      ## test if trends are cointegrated
      library(dynlm)
      tseries::po.test(expr)
      tseries::po.test(expr.loess)
      tseries::po.test(expr.maw)

      ### define number of lags
      ss1 <- select.lags(x=expr[,"Wt1"],y=expr[,"Zbtb20"],max.lag=100)
      ss2 <- select.lags(x=expr.loess[,"Wt1"],y=expr.loess[,"Zbtb20"],max.lag=100)
      ss3 <- select.lags(x=expr.maw[,"Wt1"],y=expr.maw[,"Zbtb20"],max.lag=20)
      
      Reduce( c , lapply( seq(ggenes), function(ii){
         Reduce( c, lapply(seq(ggenes), function(jj){
          if(ggenes[ii]==ggenes[jj]) return(NA) else{
            X <- select.lags(x=expr[,ii],y=expr[,jj],max.lag=20)$selection$bic
            names(X) <- paste(ggenes[ii],ggenes[jj] ,sep = "_")
            return(X)
          }
          
        }))
      }) )
      
      # plot.ts(ss1$ic)
      # plot.ts(ss2$ic)
      plot.ts(ss3$ic)
     
      ### multivariate granger
      #  type="both" include "constant" and "trend" deterministic regressors,
      # making VAR suitable for non-stationary data
      varmodel <- vars::VAR(expr.maw,lag.max = 50 ,ic="AIC",  type="both")
      grc <- granger_causality(varmodel)
        
      
    }
    
   

    ggenes <- colnames(TFmRNA_loess.smpls)[ colnames(TFmRNA_loess.smpls)%in% colnames(ATACseq_tgenesM.podo)]
    ggenes <- c("Wt1","Zbtb20","Lmx1b")
    ggenes <- TFset.down
    # ssample <- unique(TFmRNA_loess.smpls$sample)
    
    # do for maw data for individual samples
    # maw 
    WsizeF=50
    WstepF=WsizeF*3
    WsizeF=200
    WstepF=WsizeF
   
 
    TFmRNA_maw.grangerList <- lapply( 1 , function(ss){
      print( ss)
      WsizeF=50
      WstepF=WsizeF*ss
      Reduce( c , lapply( seq(listSCSN.PDS) ,
                          function(ii, datt=listSCSN.PDS)
                          {
                            require(SimDesign)
                            ssname <- table(datt[[ii]]$sample)
                            ssname <- names(ssname)[ssname>100]
                            # print(ssname)
                            grngList <- lapply( seq( ssname), function(jj){
                              # print(jj)
                              # print(ssname[jj])
                              expr <- t(datt[[ii]]@assays$RNA@data[ggenes, 
                                                                   datt[[ii]]$sample==ssname[jj]]  )
                              expr <- expr[ order(datt[[ii]]$PDS.42[ datt[[ii]]$sample==ssname[jj]]), ]
                              
                              if( nrow(expr) < (WstepF*2) ) expr.maw <- expr else{
                                # maw
                                expr.maw <- apply( expr , 2, slideFunct, 
                                                   window = nrow(expr)/WsizeF ,
                                                   step =nrow(expr)/WstepF )
                                
                              }
                              
                              Reduce( rbind, lapply( seq(ggenes), function(gg){
                                
                                Reduce( rbind, lapply(seq(ggenes)[seq(ggenes)!=gg], function(qq){
                                  # do granger 
                                  # print( paste(ggenes[gg], ggenes[qq] ))
                                  varmodel <- vars::VAR(expr.maw[, c(gg, qq)], 
                                                        lag.max = nrow(expr.maw)/2 ,ic="AIC",  type="both")
                                  grc <- quiet(  granger_causality(varmodel))
                                  # grc$result$p.F.qvalue <- qvalue::qvalue(grc$result$p.F)$qvalue
                                  # grc$result$p.Chisq.qvalue <- qvalue::qvalue(grc$result$p.Chisq)$qvalue
                                  
                                  return( grc$result )
                                }))
                                
                              }))
                              
                            })
                            
                            names(grngList)<- ssname
                            
                            
                            return(grngList)
                          }))
    })
    
       
    TFmRNA_maw.grangerList.qF <- lapply( seq(TFmRNA_maw.grangerList), 
                                         function(ii){
      TFmRNA_maw.granger <- TFmRNA_maw.grangerList[[ii]]
      
      TFmRNA_maw.granger.pF <- Reduce( cbind.data.frame, 
                                       lapply(TFmRNA_maw.granger, function(X) X["p.F"]))
      TFmRNA_maw.granger.pF$causality <-   TFmRNA_maw.granger[[1]]$Causality
      TFmRNA_maw.granger.pF <- TFmRNA_maw.granger.pF[!grepl("ALL", TFmRNA_maw.granger.pF$causality),]
      colnames(TFmRNA_maw.granger.pF) <-  c(  make.unique( names(TFmRNA_maw.granger)),"causality")
      
      TFmRNA_maw.granger.pF <- aggregate( . ~ causality , data=TFmRNA_maw.granger.pF, FUN=mean, na.rm=T)
      # TFmRNA_maw.granger.p.Chisq <- Reduce( cbind.data.frame, 
      #                                  lapply(TFmRNA_maw.granger, function(X) X["p.Chisq"]))
      
      
      rownames(TFmRNA_maw.granger.pF) <- TFmRNA_maw.granger.pF$causality
      TFmRNA_maw.granger.qF <- apply(TFmRNA_maw.granger.pF[, -1], 1, p.adjust, method="fdr")
      XX <- colSums(TFmRNA_maw.granger.qF<0.1)
      XX <- XX[order(XX)]
      print( head( XX ))
      return(TFmRNA_maw.granger.qF)
    })
   
    lapply(TFmRNA_maw.grangerList.qF, function(X){
      XX<- colSums(X<0.1)
      XX <- XX[order( XX )]
      print( tail( XX ) )
    })

  }
  
### combine expression and AUcell activity tables, + metadata
  {

   
  ### combine expression and AUcell activity tables, + metadata
  Alist <- lapply( seq(SCSNdata_list_sub) , function(ii){
    PDS <- listSCSN.PDSlimit.1K[[ii]]$PDS
    # react <-  PAaucell_reactome[[ii]][ , match( names(PDS), colnames(PAaucell_reactome[[ii]]))]
    expr <-  if( !is.null(SCSNdata_list_sub[[ii]]@assays$RNA) ) {
      as.matrix( SCSNdata_list_sub[[ii]]@assays$RNA@data )
    }  else {   as.matrix( SCSNdata_list_sub[[ii]]@assays$data@data ) }
    expr <- expr[ , match( names(PDS), colnames(expr))]
    # limit to genes expressed higher than a certan thrshld
    # expr <- expr[log10(rowMeans(expr))>=-1.5,]
    # combine expression and activities
    XX<- cbind.data.frame( t(react), t(expr ))
    # # add PDS
    XX$PDS <- PDS
    # add vector of genotypes
    XX$gtype <- SCSNdata_list_sub[[ii]]$gtypeDE[
      match( names(PDS), names( SCSNdata_list_sub[[ii]]$gtypeDE))]
    return(XX)
    
  })
    
    # all features
    allFeatures <- Reduce(intersect, lapply(Alist, colnames))
    names( Alist ) <- names(SCSNdata_list_sub)
    # featureSelect <- c( React_common.pName, TF_MeanMed)
    # PDSmarkers42 <- DS_all$gene_symbol[1:42]
    
    # save the object for shiny app
    # save( Alist , featureSelect , PDSmarkers42,
    # file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDStrajectories_app/PDS_shiny.dataBig.RData" )
    
}

### plot activity curves
  {

    ## make a table for plotting
    datTOplot <- Reduce( rbind, lapply(1:8, function(ii, 
                                                     features = input$featureList )
    {
      
      # select only cells in a desired PDS range
      Act <- Alist[[ii]]
      # # select features,  make sure we don't drop metadata
      Act <- as.data.frame(  Act[ , match( c(features , "PDS", "gtype") , 
                                           colnames(Act))] )
      # colnames(Act) <- features
      # scale all features except metadata
      if( isTRUE( input$sscale) ) Act[,  seq(features) ] <- 
        as.data.frame( scale(Act[ ,  seq(features) ]))
      
      ### compute smooth PDS
      # create bins  with equal PDS for all experiment
      spots <- seq( -0.4, 0, 0.01)
      # compute smooth PDS
      Act.melt <- Reduce( rbind,
                          lapply( 1:(length(spots)-1),  function( i )
                          {
                            # select cells within range of PDS
                            data <- Act[ Act$PDS > spots[i] & Act$PDS < spots[i+1],  ]
                            
                            # calculate average within bin for all, experimental and control cells only
                            if( nrow(data) > 3 ) {
                              all <- colMeans( data[ , seq(features) , drop=F], na.rm = T)
                            } else all <- setNames( rep(NA, length(features)), (features)) # if no cells within the bin - return NA
                            if( nrow(data[data$gtype=="control", ]) > 3 ) {
                              control <- colMeans( data[ data$gtype=="control", 
                                                         seq(features), drop=F ], na.rm = T)
                            } else control <- setNames( rep(NA, length(features)), (features))
                            if(nrow(data[data$gtype!="control", ]) > 3 ) {
                              experiment <- colMeans( data[ data$gtype!="control",
                                                            seq(features), drop=F ], na.rm = T)
                            } else experiment <- setNames( rep(NA, length(features)), (features))
                            
                            RR <- as.data.frame( melt( rbind( all, control , experiment) ) )
                            RR$PDS <- spots[i+1]
                            return(RR)
                          }))
      
      
      Act.melt$dataSet <- names(Alist)[ii]
      
      # assign proper colnames
      colnames(Act.melt)[1:2] <- c( "gtype" , "gene.sets" )
      
      return(Act.melt)
    }))
    
    # 
    datTOplot$gene.sets <- sapply( as.character(datTOplot$gene.sets) , wrap_text, n = 40)
    
    # draw the plot
   gg <-  ggplot( data = datTOplot,
            aes(x=PDS, y=value, color= dataSet ))  +
      geom_smooth(method = "loess", show.legend = T, se = T, 
                  lwd=1.2, span=1.2) +
      facet_grid( rows = vars(gtype),
                  scales ="free_y") + theme_bw()  +
      theme( text = element_text(size = 16))  
  
 print(gg)
 

}

### study cross-correlation of genes and pathways between studies, ordered by PDS
  {
    PDS.crosscorr.func  <- function( features , allFeatures,
                                     corrM="spearman")
      {
      require(reshape2)
      require(psych)
      
      # features - either a reactome pathway or a list of genes
      
      # if a reactome pathway name is provided - select all gene members of this pathway
      if( features %in% names(reactPath_gName) ) features <- c( features , reactPath_gName[[features]])
      features <- features[ features %in% allFeatures]
      
      ## make a table for plotting
      Alist2 <- Reduce( rbind, lapply(1:8, function(ii)
        {

        # select only cells in a desired PDS range
        Act <- Alist[[ii]]
        # # select features,  make sure we don't drop metadata
        Act <- as.data.frame(  Act[ , match( c(features , "PDS", "gtype") , 
                                             colnames(Act))] )
        
        ### compute smooth PDS
        # create bins  with equal PDS for all experiment
        spots <- seq( -0.4, 0, 0.01)
        # compute smooth PDS
        Act.melt <- Reduce( rbind,
                            lapply( 1:(length(spots)-1),  function( i )
                            {
                              # select cells within range of PDS
                              data <- Act[ Act$PDS > spots[i] & Act$PDS < spots[i+1],  ]
                              
                              # calculate average within bin for all, experimental and control cells only
                              if( nrow(data) > 3 ) {
                                all <- colMeans( data[ , seq(features) , drop=F], na.rm = T)
                              } else all <- setNames( rep(NA, length(features)), (features)) # if no cells within the bin - return NA
                              if( nrow(data[data$gtype=="control", ]) > 3 ) {
                                control <- colMeans( data[ data$gtype=="control", 
                                                           seq(features), drop=F ], na.rm = T)
                              } else control <- setNames( rep(NA, length(features)), (features))
                              if(nrow(data[data$gtype!="control", ]) > 3 ) {
                                experiment <- colMeans( data[ data$gtype!="control",
                                                              seq(features), drop=F ], na.rm = T)
                              } else experiment <- setNames( rep(NA, length(features)), (features))
                              
                              RR <- as.data.frame( melt( rbind( all, control , experiment) ) )
                              RR$PDS <- spots[i+1]
                              return(RR)
                            }))
        
        
        Act.melt$dataSet <- names(Alist)[ii]
        
        # assign proper colnames
        colnames(Act.melt)[1:2] <- c( "gtype" , "gene.sets" )
        
        
        
        return(Act.melt)
      }))
      
      # cast from long to wide
      Alist2.wide <- reshape(Alist2, idvar = c("PDS","gtype","gene.sets"), 
                             timevar = "dataSet", direction = "wide")
      
      corrMats <- lapply( seq(features), function(jj){
        datt <- Alist2.wide[Alist2.wide$gene.sets==features[jj] , ]
        
       
        ct <- corr.test( datt[datt$gtype=="control", -c(1:3)], adjust= "fdr",
                     method = corrM , use =  "pairwise.complete.obs")
        

         ex <- corr.test( datt[datt$gtype=="experiment", -c(1:3)], adjust= "fdr",
                           method = corrM , use =  "pairwise.complete.obs")
        return( list( control.r =ct$r ,experimental.r=ex$r ,
                      control.p =ct$p ,experimental.p=ex$p ) )
        
      })
      names(corrMats) <- features
      return(corrMats)
        
    }
   
    ### calculate cross-correlation for all podocyte genes
    # collect all podocyte genes
    allPodoGenes <- Reduce(intersect, lapply(1:8, function(ii){
      expr <-  if( !is.null(SCSNdata_list_sub[[ii]]@assays$RNA) ) {
        as.matrix( SCSNdata_list_sub[[ii]]@assays$RNA@data )
      }  else as.matrix( SCSNdata_list_sub[[ii]]@assays$data@data )
      rownames(expr)
    }  ))
    # # calculate  correlations
    # crosscorr.mats <- PDS.crosscorr.func(features = allPodoGenes ,
    #                                      allFeatures =allPodoGenes )
    crosscor.Prs.genes <- PDS.crosscorr.func(features = allPodoGenes ,
                                         allFeatures = allPodoGenes ,  corrM="pearson")
    saveRDS(crosscor.Prs.genes, file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/crosscorr/FSGScrosscor.Prs.genes")
    
    # compute average cross-crrelation between models
    crosscorr.ave <- Reduce( rbind, lapply( crosscor.Prs.genes, function(X) {
      
      sapply( X, function(XX) {
        XX <- XX[lower.tri(XX)]
        XX[is.na(XX)] <- 0 
        mean(XX, na.rm=T ) 
      } )
    }))
    
    rownames( crosscorr.ave )<- names( crosscorr.mats )
    crosscorr.ave <- as.data.frame(crosscorr.ave)
    crosscorr.ave$dirct <-  ifelse( crosscorr.ave$control.r > 
                                      crosscorr.ave$experimental.r , 
                                       "down", "up")
    crosscorr.ave$gene.set <- rownames(crosscorr.ave)
    crosscorr.ave$diff <- crosscorr.ave$experimental.r -  crosscorr.ave$control.r
    crosscorr.ave$ave <- (crosscorr.ave$experimental.r +  crosscorr.ave$control.r)/2
    crosscorr.ave$PDmarker <- ifelse( crosscorr.ave$gene.set%in%DS_all$gene_symbol, "yes","no" )
    crosscorr.ave$PDmarker42 <- ifelse( crosscorr.ave$gene.set%in%DS_all$gene_symbol[1:42] , "yes","no" )
    crosscorr.ave <- crosscorr.ave[ order( crosscorr.ave$ave , na.last = F ) , ]
    crosscorr.ave$p.diff <- -log10(crosscorr.ave$experimental.p/ crosscorr.ave$control.p )
    
    write.table(crosscorr.ave , sep = "\t", quote = F, row.names = F, 
                file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/crosscorr/crossCor.allgenes_tab.tsv")
    
    ### plot distribution of values
    crosscorr.ave <- read.table(header = T,  "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/crosscorr/crossCor.allgenes_tab.tsv" )
    datTOplot <- crosscorr.ave
    datTOplot <- melt( as.data.frame(datTOplot), 
                       measure.var = c("control.r","experimental.r"))
    library( reshape2)
    ggplot( datTOplot , aes( x= value, color= variable)) + 
      geom_density(lwd=2, adjust = 2)+ theme_bw()+ 
      ggtitle("crosscorrelation of genes\nbetween FSGS models")+
      theme(text=element_text(size=24))
    
    ### plot individual genes
    {
      # gain some ideas about possible candidate distributions
      library(fitdistrplus)
      descdist( as.numeric(na.omit(crosscorr.ave$ave)), discrete = FALSE)
      
      ### select genes more than 2sd from the mean
      crosscorr_2.5sdAway <- crosscorr.ave$gene.set[abs(crosscorr.ave$ave) > (
        mean(crosscorr.ave$ave, na.rm=T)+ 2.58*sd(crosscorr.ave$ave, na.rm=T)
      ) & !is.na(crosscorr.ave$ave )]
      # how many genes are described PDS markers
      summary(as.factor(
        crosscorr.ave$PDmarker[crosscorr.ave$gene.set%in%crosscorr_2.5sdAway])   )
      ### select genes with diffeerence between control and exp. more than 2sd from the mean
      crosscorr.diffUP_2.5sdAway <- crosscorr.ave$gene.set[ (crosscorr.ave$diff) > (
        mean(crosscorr.ave$diff, na.rm=T)+ 2.58*sd(crosscorr.ave$diff, na.rm=T)
      ) & !is.na(crosscorr.ave$diff ) & 
        ( abs(crosscorr.ave$control.r)  > 0.1 | abs(crosscorr.ave$experimental.r) > 0.1)]
      crosscorr.diffDOWN_2.5sdAway <- crosscorr.ave$gene.set[ (crosscorr.ave$diff) < (
        mean(crosscorr.ave$diff, na.rm=T) - 2.58*sd(crosscorr.ave$diff, na.rm=T)
      ) & !is.na(crosscorr.ave$diff )  & 
        ( abs(crosscorr.ave$control.r) > 0.1 | abs(crosscorr.ave$experimental.r) > 0.1)]
      
      # plot the genes  
      library(ggrepel)
      datTOplot <- tail(crosscorr.ave ,20) 
      datTOplot <- crosscorr.ave[crosscorr.ave$gene.set %in% crosscorr.diffUP_2.5sdAway,]
      
      datTOplot <- melt( datTOplot, 
                         measure.var = c("control.r","experimental.r"))
      datTOplot$label <- ifelse( datTOplot$variable=="experimental.r",
                                 as.character( datTOplot$gene.set), "" )
      
      ggplot( datTOplot, aes(x=variable, y=value, group = gene.set))+
        geom_line(aes(col=dirct), lwd=2)  +
        scale_color_manual(values = c("down" = "cornflowerblue",
                                      "up"="salmon",
                                      "NA"="grey")) +
        geom_label_repel(aes(label = label), nudge_x = 0.35, size = 5) +
        geom_point( size=5, col="red") + ylab("ρ") + ggtitle("average spearman correlation across models")+
        theme_bw()+ theme( text = element_text( size=18))
      
      ### ephrin signalling
      {
        ggenes <- reactPath_gName[["Ephrin signaling"]]
        ggenes <- crosscorr.diffDOWN_2.5sdAway
        
        
        plotExpr <- function(ggenes=crosscorr.diffDOWN_2.5sdAway , 
                 ggtitle="genes which crosscorr. goes DOWN in FSGS", foldchange=F) 
          {
        
        ll <- Reduce( cbind, lapply( 1:8, function(ii){
          
          datt <- SCSNdata_list_sub[[ii]]
          expr <-  if( !is.null(datt@assays$RNA) ) { as.matrix( datt@assays$RNA@data )
          }  else {   as.matrix( datt@assays$data@data ) }
          # compute fold change or mean
          if(isTRUE(foldchange)){
            wt <- rowMeans( expr[ match(  ggenes, rownames(expr)),
                                                        datt$gtypeDE=="control"]) 
            mut <-  rowMeans( expr[ match(  ggenes, rownames(expr)),
                                                          datt$gtypeDE!="control"])
            log10(mut/wt)
           } else  rowMeans( expr[ match(  ggenes, rownames(expr)),] )
        }) )
        
        
        rownames(ll) <- ggenes
        colnames(ll) <- names(SCSNdata_list_sub)
        
        ll.melt <- melt(ll)
        
        if( isTRUE(foldchange) ) xxlab="log10FC in FSGS relative to controls" else xxlab="mean expression in datasets"
        ggplot( ll.melt, aes( x=value, y = reorder(Var1, value,na.rm = TRUE))) + 
          geom_boxplot( )+ 
          geom_jitter(height = 0.2, aes( color=Var2), size=3 )+ 
          theme_bw()+ theme(text=element_text(size=20))+ ggtitle(ggtitle)+
          ylab("genes") + xlab( xxlab )
        }
        
        gglist <-  list(crosscorr.diffUP_2.5sdAway,
                        crosscorr.diffDOWN_2.5sdAway
                        )
        namelist <- list( "genes which crosscorr. goes UP in FSGS", 
                          "genes which crosscorr. goes DOWN in FSGS")
        ll <- lapply( 1:2, function(ii) plotExpr( ggenes=gglist[[ii]],
                                                 ggtitle=namelist[[ii]], 
                                                 foldchange=T) )
        ll2 <- lapply( 1:2, function(ii) plotExpr( ggenes=gglist[[ii]],
                                                  ggtitle=namelist[[ii]], 
                                                  foldchange=F) )
        
        
        gridExtra::grid.arrange( grobs = c(ll, ll2), nrow=2)
        
        }
    }
    
    
    


    
    }

### order using AUC
  {
    
  # TFs
    TFmRNA_smooth.common_AUC <- Reduce( cbind, lapply( 1:8 , function(ii){
      XX <-  data.frame((TFmRNA_smooth.norm[[ii]] ))
      XX$ID <- (1:nrow(XX))
      
      AUC <- Reduce( c, lapply( 1:(ncol(XX)-1), function(jj){    sum(XX[,jj])
        # datt <- ( XX[ , c(jj,ncol(XX))] )
        # loessMod70 <- loess( as.formula(paste0(colnames(XX)[1], "~ ID")) , data= XX, span=0.70) # 70% smoothing span
        # smoothed70 <- predict(loessMod70) 
      }))
      AUC <- AUC[ match( TF_MeanMed , colnames(TFmRNA_smooth.norm[[ii]]))]
    }))
    rownames(TFmRNA_smooth.common_AUC) <- TF_MeanMed
      
  # reactomr  
  PAaucell.react_smooth.common_AUC <- Reduce( cbind, lapply( 1:8 , function(ii){
    XX <-  data.frame((PAaucell.react_smooth.common.norm[[ii]] ))
    colnames(XX) <- paste( "var" , 1:42, sep = "_")
    XX$ID <- (1:nrow(XX))
    
    AUC <- Reduce( c, lapply( 1:42, function(jj){    sum(XX[,jj])
      # datt <- ( XX[ , c(jj,ncol(XX))] )
      # loessMod70 <- loess( as.formula(paste0(colnames(XX)[1], "~ ID")) , data= XX, span=0.70) # 70% smoothing span
      # smoothed70 <- predict(loessMod70) 
    }))
  }))
  
    colnames(PAaucell.react_smooth.common_AUC) <- names(SCSNdata_list_sub)[1:8]
  rownames(PAaucell.react_smooth.common_AUC) <- React_common.pName
  
  # mutated cells
  PAaucell.react_smooth.mut.common_AUC <- Reduce( cbind, lapply( 1:8 , function(ii){
    XX <-  data.frame((PAaucell.react_smooth.common.norm[[ii]] ))
    AUC <- Reduce( c , lapply( 1:42, function(jj) sum(XX[,jj])))
  }))
  
  # select pathways that go down
  PAaucell.react_smooth.sameDirDown_AUC <- PAaucell.react_smooth.common_AUC[React_common.pName.sameDirUp,]
  # make a plot of pathways that go down
  pheatmap::pheatmap(apply(PAaucell.react_smooth.sameDirDown_AUC,2,rank), 
                     color = rev(inferno(27)))
  }

#### relate pathways/gSets and TFs via GRNs ####

  # load podo Paths/gsets
  PodoPathGSet <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PodoPathGSet.05.01.24.7plus.rda")
  circ.genes.bulk0.01 <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_circadian/circ.genes.bulk0.01.rda")
  
  
  ### ATACseq based GRN
    {
    # ATACseq_hint.fimo_tftc.GRN <-   readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/ATACseq/atac.podo_hint_fimo.77TF_GRN.rda")
      ATACseq_hint.fimo_tftc.GRN <-   readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/ATACseq/atac.podo_tobias_fimo_tftc.0.01.GRN.PWMclust.rda")
      ATACseq_hint.fimo_tftc.GRN <-   readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/ATACseq/atac.podo_tobias_fimo_tftc.0.1.GRN.PWMclust.rda")
      
      ATACseq_tgenesM.podo <- ATACseq_hint.fimo_tftc.GRN[ 
      rownames(ATACseq_hint.fimo_tftc.GRN) %in% allPodoGenes, ]
    


  }
  
  ### create REGnets for specific gene sets
    {
     lapply( seq(PodoPathGSet), function(ii){
        nett <- ATACseq_tgenesM.podo[rownames(ATACseq_tgenesM.podo)%in% 
                                       PodoPathGSet[[ii]],]
        nett <- nett[ rowSums(nett)>0, colSums(nett)>0]
        nett <- reshape2::melt(nett)[,c(2,1,3)]
        nett <- nett[nett$value>0,]
        colnames(nett) <- c("TFreg", "tgene","TFtc.SUMscore")
        nett$edgeID <- paste0("edge_",seq(nrow(nett)))
        nett$TFname <- paste0( sub(":.*","",nett$TFreg), "_clust")
        # print(dim(nett))
        write.table(nett , row.names = F, sep = "\t", quote = F,
                  file=paste0( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/pathway_viz/cytoscapeVIZ/Podo_areas/PodoPathGSet_GRNs/",
                                      names(PodoPathGSet)[ii] ,"_tftc.0.1.GRN.clust0.3.csv"))
      })
      

  }
  
  ### test TF regltrs of gene sets
    {
  
    
  ### Hypergeometric test of enrichment 
    { 
        ## function to test TF target enrichment
        gset.TFqval  <- function( geneSet , GRN=ATACseq_tgenesM.podo) 
          {
          p.adjust( apply( GRN, 2, function(TF){
            TFtg <- rownames(GRN)[TF>2]
            print( length(TFtg))
            # phyper=(overlap,list1,PopSize,list2,lower.tail = FALSE, log.p = FALSE)
            overlap <- length( intersect( geneSet,  TFtg ) )
            list1 <- length( geneSet )
            list2 <- length( TFtg )
            popSize <- nrow( GRN )
            ph <- phyper(overlap-1, list2, popSize - list2, list1, lower.tail=FALSE )
            
            
           
            return(ph )
          }), method="fdr" ) }
        
        # apply 
        # gset.TFqval(  geneSet = gene.topCorPDS.xprmnt )
        PodoPathGSet_PDS.test <- sapply( c(list(circ.genes.bulk0.01), PodoPathGSet), function(X) {
          # print(head(X))
          gset.TFqval.xprmnt <-  gset.TFqval( geneSet =  X ,
                                              GRN = ATACseq_tgenesM.podo )
         
        } )
       colnames(PodoPathGSet_PDS.test) <- c( "circ.genes.bulk0.01", 
                                           names( PodoPathGSet ) )
          
          
        # PodoPathGSet_PDS.test[ PodoPathGSet_PDS.test >  0.05] <- NA
        PodoPathGSet_PDS.test <- PodoPathGSet_PDS.test[ rowSums(!is.na(
          PodoPathGSet_PDS.test))>1,
          colSums(!is.na(
            PodoPathGSet_PDS.test))>1]
        siglbl <- round( PodoPathGSet_PDS.test , 4)
        pheatmap::pheatmap(-log10(PodoPathGSet_PDS.test), cluster_cols = T,
                           cluster_rows = T,  display_numbers = siglbl, 
                           number_color = "red" )
      }
   
      ### test for TF motif enrichment
      {
        ddir <-"/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/AME_motifAnalysis/PodoPathGSet/"
        
        # run the test
        PodoPathGSet_AMEmotifEnrich <- lapply( seq(PodoPathGSet), function(ii){
          print(ii)
          TF_motifEnrichTest( gset =  PodoPathGSet[[ii]] , 
                              bckgrGenes = "all", 
                              motifDir = "/media/tim_nevelsk/WD_tim/ANNOTATIONS/CISBP/CISBP2mouse_podoTF_08.01.24.meme",  
                              outDir=  paste0( ddir,names(PodoPathGSet)[ii],"/") )
        })
        
        saveRDS(PodoPathGSet_AMEmotifEnrich , "PodoPathGSet_AMEmotifEnrich.rda" ) 
        
        
        
      }
      
    
    }
  
  ### heatmap of intersect of gSets and TFtargets
    {
    datt <- ATACseq_tgenesM.podo
    XX <- Reduce( rbind, lapply( c(list(circ.genes.bulk0.01), PodoPathGSet), function(gSet){ 
      apply( datt, 2, function(TF ){
        TF <- TF[ TF>0 & names(TF) %in% gSet] 
        print( length(TF))
      })
    }))
    rownames(XX) <- c( "circ.genes.bulk0.01", 
                       names( PodoPathGSet ) )
    
    toPlot <- XX[rowSums(XX)>10,]
    # # scale
    toPlot <-t(scale(t(toPlot),  center = F ))
    toPlot <-  toPlot[,colSums(toPlot)>5]
    
    # # or rank
    # toPlot <- t( apply( toPlot ,1, rank, ties.method = "min" ) )
    # # col1 <-RColorBrewer::brewer.pal(10,"Paired")
    
    # plot 
    colCol <- map2color( colSums( ATACseq_tgenesM.podo[ ,colnames(toPlot)] ), grey.colors(ncol(toPlot)))
    rowCol <- map2color( Reduce( c, lapply( PodoPathGSet[rownames(toPlot)], length)), grey.colors(20))
    bbreaks <- seq( ncol(toPlot) )
    
    pdf(height = 8, width = 24, file="circ.podoPaths.TFprior.overlap_heatmap.rank")
     # png(height = 400, width = 1000, file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/PDS.TFprior.overlap_heatmap.rank.png")
    
    gplots::heatmap.2( toPlot , col = viridis_pal(option="inferno" ) , 
                         Rowv = T , trace="none",
                          # breaks = bbreaks, 
                         RowSideColors= rowCol , 
                         ColSideColors= colCol , 
                         margins = c(20,14) , key=T, 
                         cexRow=2.0 , cexCol=1.8 )
    dev.off()
  }
  
  ### Heatmaps that relates TF activity to PDS
    {
    ### Heatmap of correlation between TF mRNA lvls ad PDS
      {
        # load correlations of individual genes 
        toPlot.all <- genecorrPDS.Sprmn_r
        
        # combine all controls and all exprmnt in indiv. columns
        toPlot.all <- cbind(  rowMeans(toPlot.all[,1:7], na.rm = T), 
                              rowMeans(toPlot.all[,8:14], na.rm = T),
                              toPlot.all[,8:14])
        
        colnames(toPlot.all) <- c("ctrl_mean","exprmnt_mean", names( listSCSN.PDS))
        toPlot <- scale(toPlot.all, scale = F)
        toPlot <- t( toPlot[ rownames(toPlot) %in% TF_MeanMed , ] )
        
      # plot 
      # colCol <- map2color( colSums( ATACseq_tgenesM[ , colnames(toPlot)] ), grey.colors(ncol(toPlot)))
      mybreaks <- seq(-0.2,0.2,length=22)
      gplots::heatmap.2(toPlot, col = colorspace::diverge_hcl(21), 
                        Rowv = F, trace="none",
                        breaks=mybreaks, 
                        # ColSideColors= colCol, 
                        margins = c(5,12), key=T, 
                        cexRow=2, cexCol=1.0 )
    }
    
    ###  Heatmap of correlation between AUCell of TF genes and PDS
      {
        PAaucell_TFs <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PathwayActivity/PAaucell_list.TFs.unique.rda")
        ### calculate corr with PDS
        PAaucell_TFs.PDScorr <- lapply(1:7, function(ii){
          cell.wt <- colnames(listSCSN.PDS[[ii]])[listSCSN.PDS[[ii]]$gtypeDE=="control"]
          cells <- colnames(listSCSN.PDS[[ii]])[listSCSN.PDS[[ii]]$gtypeDE!="control"]
          
          wt <- psych::corr.test( t(PAaucell_TFs[[ii]][,cell.wt]) , 
                                  y =listSCSN.PDS[[ii]]$PDS.42[cell.wt] ,
                                  method = "spearman")
          mut <- psych::corr.test( t(PAaucell_TFs[[ii]][,cells]) , 
                                   y =listSCSN.PDS[[ii]]$PDS.42[cells] ,
                                   method = "spearman")
          XX <- cbind.data.frame(wt.r=wt$r,
                                 wt.p=wt$p,
                                 mut.r=mut$r,
                                 mut.p=mut$p)
          # XX <- XX[rownames(XX)!="Hoxa7",]
          print( dim(XX) )
          return( XX) 
        })
        
        ggenes <- Reduce( union , lapply(PAaucell_TFs.PDScorr, rownames))
        toPlot <- Reduce(cbind, lapply(PAaucell_TFs.PDScorr, function(X){
          X <- X[match( ggenes, rownames(X)), c("wt.r","mut.r"), drop=F]
          rownames(X) <- ggenes
          return(X)
        }))
        toPlot <- toPlot[ , c(seq(from=1,to=14,by=2), seq(from=2,to=14,by=2)) ]
        toPlot <- cbind(  rowMeans(toPlot[,1:7], na.rm = T), 
                          rowMeans(toPlot[,8:14], na.rm = T),
                          toPlot[,8:14])
        
        colnames(toPlot) <- c("ctrl_mean","exprmnt_mean", names( listSCSN.PDS))
        # colnames(toPlot)<- paste( rep( names(listSCSN.PDS) ,each=2),
        #                            rep(c("ctrl","xprmnt"),7) , sep ="_" )
        ## plot
        # colCol <- map2color( colSums( ATACseq_tgenesM.podo[ ,colnames(toPlot)] ), grey.colors(ncol(toPlot)))
        mybreaks <- seq(-0.2,0.2,length=22)
        gplots::heatmap.2( toPlot, col = colorspace::diverge_hcl(21), 
                           Rowv = F, trace="none",
                           breaks=mybreaks, 
                           # ColSideColors= colCol, 
                           margins = c(5,12), key=T, 
                           cexRow=2, cexCol=1.0 )
      }
      
    ###  Heatmap of correlation between AUCell of unique TFtargets and PDS
      {
        PAaucell_TFs.unique <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PathwayActivity/PAaucell_list.TFs.unique.rda")
        ### calculate corr with PDS
        PAaucell_TFs.unique.PDScorr <- lapply(1:7, function(ii){
          cell.wt <- colnames(listSCSN.PDS[[ii]])[listSCSN.PDS[[ii]]$gtypeDE=="control"]
          cell.mut <- colnames(listSCSN.PDS[[ii]])[listSCSN.PDS[[ii]]$gtypeDE!="control"]
          
          wt <- psych::corr.test( t(PAaucell_TFs.unique[[ii]][,cell.wt]) , 
                                  y =listSCSN.PDS[[ii]]$PDS.42[cell.wt] ,
                                  method = "spearman")
          mut <- psych::corr.test( t(PAaucell_TFs.unique[[ii]][,cell.mut]) , 
                                   y =listSCSN.PDS[[ii]]$PDS.42[cell.mut] ,
                                   method = "spearman")
          XX <- cbind.data.frame(wt.r=wt$r,
                                 wt.p=wt$p,
                                 mut.r=mut$r,
                                 mut.p=mut$p)
          # XX <- XX[rownames(XX)!="Hoxa7",]
          print( dim(XX) )
          return( XX) 
        })
        
        toPlot <- Reduce(cbind, lapply(PAaucell_TFs.unique.PDScorr, function(X){
          X[,c("wt.r","mut.r"), drop=F]
        }))
        toPlot <- toPlot[ , c(seq(from=1,to=14,by=2), seq(from=2,to=14,by=2)) ]
        toPlot <- cbind(  rowMeans(toPlot[,1:7], na.rm = T), 
                          rowMeans(toPlot[,8:14], na.rm = T),
                          toPlot[,8:14])
        
        colnames(toPlot) <- c("ctrl_mean","exprmnt_mean", names( listSCSN.PDS))
        # colnames(toPlot)<- paste( rep( names(listSCSN.PDS) ,each=2),
        #                            rep(c("ctrl","xprmnt"),7) , sep ="_" )
        ## plot
        # colCol <- map2color( colSums( ATACseq_tgenesM.podo[ ,colnames(toPlot)] ), grey.colors(ncol(toPlot)))
        mybreaks <- seq(-0.2,0.2,length=22)
        # toPlot <- t( scale(toPlot, scale = F))
        gplots::heatmap.2( t(toPlot), col = colorspace::diverge_hcl(21), 
                          Rowv = F, trace="none",
                          breaks=mybreaks, 
                          # ColSideColors= colCol, 
                          margins = c(5,12), key=T, 
                          cexRow=2, cexCol=1.0 )
      }
      
    ###  Heatmap of correlation between AUCell of all TFtargets and PDS
      {
        PAaucell_TFs.GRN <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PathwayActivity/PAaucell_list.TFs.unique.rda")
        ### calculate corr with PDS
        PAaucell_TFs.GRN.PDScorr <- lapply(1:7, function(ii){
          cell.wt <- colnames(listSCSN.PDS[[ii]])[listSCSN.PDS[[ii]]$gtypeDE=="control"]
          cell.mut <- colnames(listSCSN.PDS[[ii]])[listSCSN.PDS[[ii]]$gtypeDE!="control"]
          
          wt <- psych::corr.test( t(PAaucell_TFs.GRN[[ii]][,cell.wt]) , 
                                  y =listSCSN.PDS[[ii]]$PDS.42[cell.wt] ,
                                  method = "spearman")
          mut <- psych::corr.test( t(PAaucell_TFs.GRN[[ii]][,cell.mut]) , 
                                   y =listSCSN.PDS[[ii]]$PDS.42[cell.mut] ,
                                   method = "spearman")
          XX <- cbind.data.frame(wt.r=wt$r,
                                 wt.p=wt$p,
                                 mut.r=mut$r,
                                 mut.p=mut$p)
          # XX <- XX[rownames(XX)!="Hoxa7",]
          print( dim(XX) )
          return( XX) 
        })
        
        toPlot <- Reduce(cbind, lapply(PAaucell_TFs.GRN.PDScorr, function(X){
          X[,c("wt.r","mut.r"), drop=F]
        }))
        toPlot <- toPlot[ , c(seq(from=1,to=14,by=2), seq(from=2,to=14,by=2)) ]
        toPlot <- cbind(  rowMeans(toPlot[,1:7], na.rm = T), 
                          rowMeans(toPlot[,8:14], na.rm = T),
                          toPlot[,8:14])
        
        colnames(toPlot) <- c("ctrl_mean","exprmnt_mean", names( listSCSN.PDS))
        # colnames(toPlot)<- paste( rep( names(listSCSN.PDS) ,each=2),
        #                            rep(c("ctrl","xprmnt"),7) , sep ="_" )
        ## plot
        # colCol <- map2color( colSums( ATACseq_tgenesM.podo[ ,colnames(toPlot)] ), grey.colors(ncol(toPlot)))
        mybreaks <- seq(-0.2,0.2,length=22)
        toPlot <- t( scale(toPlot, scale = F))
        gplots::heatmap.2( (toPlot), col = colorspace::diverge_hcl(21), 
                           Rowv = F, trace="none",
                           breaks=mybreaks, 
                           # ColSideColors= colCol, 
                           margins = c(5,12), key=T, 
                           cexRow=2, cexCol=1.0 )
      }
    
      ###  Heatmap of correlation between Wt1 measures of activity and PDS
      {
        PAaucell_Wt1 <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PathwayActivity/PAaucell_Wt1.rda")
        ### calculate corr with PDS
        PAaucell_Wt1.PDScorr <- lapply(1:7, function(ii){
          cell.wt <- colnames(listSCSN.PDS[[ii]])[listSCSN.PDS[[ii]]$gtypeDE=="control"]
          cell.mut <- colnames(listSCSN.PDS[[ii]])[listSCSN.PDS[[ii]]$gtypeDE!="control"]
          
          wt <- psych::corr.test( t(PAaucell_Wt1[[ii]][,cell.wt]) , 
                                  y =listSCSN.PDS[[ii]]$PDS.42[cell.wt] ,
                                  method = "spearman")
          mut <- psych::corr.test( t(PAaucell_Wt1[[ii]][,cell.mut]) , 
                                   y =listSCSN.PDS[[ii]]$PDS.42[cell.mut] ,
                                   method = "spearman")
          XX <- cbind.data.frame(wt.r=wt$r,
                                 wt.p=wt$p,
                                 mut.r=mut$r,
                                 mut.p=mut$p)
          # XX <- XX[rownames(XX)!="Hoxa7",]
          print( dim(XX) )
          return( XX) 
        })
        
        toPlot <- Reduce(cbind, lapply(PAaucell_Wt1.PDScorr, function(X){
          X[,c("wt.r","mut.r"), drop=F]
        }))
        toPlot <- toPlot[ , c(seq(from=1,to=14,by=2), seq(from=2,to=14,by=2)) ]
        toPlot <- cbind(  rowMeans(toPlot[,1:7], na.rm = T), 
                          rowMeans(toPlot[,8:14], na.rm = T),
                          toPlot[,8:14])
        
        colnames(toPlot) <- c("ctrl_mean","exprmnt_mean", names( listSCSN.PDS))
        # colnames(toPlot)<- paste( rep( names(listSCSN.PDS) ,each=2),
        #                            rep(c("ctrl","xprmnt"),7) , sep ="_" )
        ## plot
        # colCol <- map2color( colSums( ATACseq_tgenesM.podo[ ,colnames(toPlot)] ), grey.colors(ncol(toPlot)))
        mybreaks <- seq(-0.4,0.4,length=22)
        pdf(width = 8, height = 8, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/TFA/PDScorr.Wt1_various_heatmap.pdf")
        # toPlot <- t( scale(toPlot, scale = F))
        gplots::heatmap.2( t(toPlot), col = colorspace::diverge_hcl(21), 
                           Rowv = F, trace="none",
                           breaks=mybreaks, 
                           # ColSideColors= colCol, 
                           srtCol=45,
                           margins = c(12,10), key=T, 
                           cexRow=2,cexCol=2 )
        dev.off()
      }

    }
  
  ### heatmap that relates TF activity, pathways and PDS
    {
    PAaucell_TFs.Pthws <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PathwayActivity/PAaucell_TFs.Pthws.rda")
   
      
       ### calculate corr with PDS
   
       PAaucell_TFs.Pthws.PDScorr <- lapply(1:7, function(ii){
      cell.wt <- colnames(listSCSN.PDS[[ii]])[listSCSN.PDS[[ii]]$gtypeDE=="control"]
      cell.mut <- colnames(listSCSN.PDS[[ii]])[listSCSN.PDS[[ii]]$gtypeDE!="control"]
      
      wt <- psych::corr.test( t(PAaucell_TFs.Pthws[[ii]][,cell.wt]) , 
                              y =listSCSN.PDS[[ii]]$PDS.42[cell.wt] ,
                              method = "spearman")
      mut <- psych::corr.test( t(PAaucell_TFs.Pthws[[ii]][,cell.mut]) , 
                               y =listSCSN.PDS[[ii]]$PDS.42[cell.mut] ,
                               method = "spearman")
      XX <- cbind.data.frame(wt.r=wt$r,
                             wt.p=wt$p,
                             mut.r=mut$r,
                             mut.p=mut$p)
      # XX <- XX[rownames(XX)!="Hoxa7",]
      print( dim(XX) )
      return( XX) 
    })
    
    allgsets<- Reduce( union, lapply(PAaucell_TFs.Pthws.PDScorr, rownames))
   
    TFs.Pthws_tab <- Reduce(cbind, lapply(PAaucell_TFs.Pthws.PDScorr, function(X){
      X <-  X[match( allgsets  , rownames(X)) ,c("wt.r","mut.r"), drop=F]
      rownames(X) <- allgsets
      return(X)
    }))
    
    
    TFs.Pthws_tab <- TFs.Pthws_tab[ , c(seq(from=1,to=14,by=2), 
                                        seq(from=2,to=14,by=2)) ]
    TFs.Pthws_tab <- cbind(  rowMeans(TFs.Pthws_tab[,1:7], na.rm = T), 
                      rowMeans(TFs.Pthws_tab[,8:14], na.rm = T),
                      TFs.Pthws_tab[,8:14])
    
    colnames(TFs.Pthws_tab) <- c("ctrl_mean","exprmnt_mean", names( listSCSN.PDS))
    # select a gene set 
    toPlotList  <- lapply( seq((PodoPathGSet)), function(ii){
      print(ii)
      X <- TFs.Pthws_tab[ grepl( paste( names(PodoPathGSet)[ii],"_" , sep = "") ,
                                        rownames(TFs.Pthws_tab)),]
      rownames(X)<- sub(  "_", "",sub( names(PodoPathGSet)[ii], "", rownames(X)))
      # number of tgenes per TF
      X$annotCl <- sapply( ATACseq_tgenesM.Pths.List[  grepl( paste( names(PodoPathGSet)[ii],"_" , sep = "") ,
                                                              names(ATACseq_tgenesM.Pths.List))] ,
                           length)
      return(X)
      
    })
    

# plot
    lapply( seq(toPlotList), function(ii){
      
      toPlot <- toPlotList[[ii]]
      # select TFs with 3 and more  targets
      toPlot <- toPlot[toPlot$annotCl>2,]
      # select expressed in more than one sample
      toPlot <- toPlot[rowSums( is.na(toPlot[,1:9]))<8,]
      # 
      colCol <- map2color( toPlot$annotCl, grey.colors(ncol(toPlot)))
      # mybreaks <- seq(-0.2,0.2,length=22)
      
      pdf( height = 8, width = 20, file=paste("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/TFA/TFAaucell/Pthws.Targets/PDScorr.TFs_aucell.Pthw.tg_heatmap_",
      names(PodoPathGSet)[ii]  ,"_small.pdf", sep = ""))
      # png( height = 800, width = 2000, file=paste("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/TFA/TFAaucell/Pthws.Targets/PDScorr.TFs_aucell.Pthw.tg_heatmap_",
      #                                         names(PodoPathGSet)[ii]  ,"_small.png", sep = ""))
      gplots::heatmap.2( (t(toPlot[,-10])), col = colorspace::diverge_hcl(21), 
                         Rowv = F, trace="none",symbreaks = T, 
                          # breaks=mybreaks,
                         main = paste( "correlation of PDS and TFA (TFtg AUCell score)\nin ", 
                                       names(PodoPathGSet)[ii],sep = ""),
                         ColSideColors= colCol, 
                         margins = c(8,12), key=T, 
                         cexRow=2, cexCol=2 )
      dev.off()
    })

    }

#### visualise expression changes maped on GRN  ####
  
### make a movie!
  {
      ### prepare a network
      library(igraph)
      
      ATACseq_deTF <- ATACseq_tgenesM[ (rownames(ATACseq_tgenesM) %in% (TF_MeanMed)) , 
                                       colnames(ATACseq_tgenesM) %in% names(TFcore.42) ]
      # ATACseq_deTF <- ATACseq_deTF[, !colnames(ATACseq_deTF)%in%c("Max", "Nr1d1")]
      # remove TF2 without regulators
      ATACseq_deTF <- ATACseq_deTF[rowSums(ATACseq_deTF>0)>=1, colSums(ATACseq_deTF>0)>=1]
      # make an edge list from non-square adjacency matrix
      ATACseq_deTF_melt <- reshape2::melt(t(ATACseq_deTF))
      ATACseq_deTF_melt <- ATACseq_deTF_melt[ATACseq_deTF_melt$value>0,]
      # create a graph
      gg <- graph_from_edgelist( as.matrix(ATACseq_deTF_melt[,-3]) )
      
      V(gg)$size <- sqrt(degree(gg, mode = "out")+1)*5
      
      pdf(width=8, height = 8, file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFA/19TF_ATACnet.pdf")
      plot(gg,   edge.arrow.size=.2,
           vertex.label.cex=1.2,
           vertex.label.dist=0.1 ,
           # layout=layout.fruchterman.reingold(gg, niter=50) ,
           vertex.label = ifelse( degree(gg, mode = "all") > 0, V(gg)$name, NA))
      dev.off()
      
      #### color edges
      exprMat <- as.matrix( SCSNdata_list[[1]]@assays$RNA@data )
      exprMat <- exprMat[match( V(gg)$name , rownames(exprMat)),]
      corMat <- psych::corr.test( t(exprMat))
      corMat.melt <- reshape2::melt((corMat$r) )
      corMat.melt$pval <- reshape2::melt((corMat$p))$value
      corMat.melt$padj <- p.adjust( corMat.melt$pval  )
      
      #### prepare expression
      exprMat <- as.matrix( SCSNdata_list[[1]]@assays$RNA@data )
      
      exprMat <-  as.data.frame ( t(apply( exprMat[match( V(gg)$name , rownames(exprMat)),], 1 , function(x) slideFunct(x,500, 25) ) ) )
      exprMat <- t(scale(t(exprMat)))
      exprMat_mod <- exprMat
      exprMat_mod[exprMat_mod>2] <- 2
      exprMat_mod[exprMat_mod< -2] <- - 2
      
      ### make a movie!
      library('animation')
      library(viridis)
      library(ggraph)
      
      animation::saveVideo( {  
        for( i in 1:ncol(exprMat_mod)){
          pp <- ggraph(gg,"stress",bbox = 15) +
            geom_edge_link0(edge_colour = "grey66",edge_width = 0.2)+
            geom_node_point(aes(fill = exprMat_mod[,i]),shape = 21,size =  sqrt(degree(gg, mode = "all"))*3 )+
            geom_node_text(aes(label = ifelse( degree(gg, mode = "all") > 12, V(gg)$name, NA),size = degree(gg)),
                           family = "serif",repel = F,
                           colour = "darkblue" , size = 8 )+
            scale_fill_viridis( limits=c(-2.01,2.01))+
            # scale_size(range=c(2,5),guide = FALSE)+
            theme_graph()+theme(legend.position = "bottom")
          print( pp )
        }
      },
      interval = 0.04166667, 
      ani.height = 600 ,
      ani.width = 600 ,
      video.name ="/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/snRNAseq_podo/PDS/disease.score/animation/networkScaled_movie3.mp4" 
      )
      
      animation::saveHTML( 
        {  
          for( i in 1:ncol(exprMat)){
            pp <- ggraph(gg,"stress",bbox = 15) +
              geom_edge_link0(edge_colour = "grey66",edge_width = 0.2,
                              arrow = grid::arrow(angle = 15,type = "closed"))+
              geom_node_point(aes(fill = exprMat[,i]),shape = 21,size =  sqrt(degree(gg, mode = "all"))*3 )+
              geom_node_text(aes(label = ifelse( degree(gg, mode = "all") > 10, V(gg)$name, NA)),
                             family = "serif",repel = T,
                             colour = "darkblue" , size = 8 )+
              scale_fill_viridis()+
              # scale_size(range=c(2,5),guide = FALSE)+
              theme_graph()+theme(legend.position = "bottom")
            print( pp )
          }
        } ,
        interval = 0.1, 
        ani.height = 600 ,
        ani.width = 600 ,
        # htmlfile ="/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/snRNAseq_podo/PDS/disease.score/animation/networkScaled_animated.html" 
      )
      
      
    }

### combine network view nets with GRN
  {
    
     colnames(ATACseq_tgenesM.podo)
    ### compare PWM to get nonredundant set
      {
      library(universalmotif)
        
      motifs <-read_meme("/media/tim_nevelsk/WD_tim/ANNOTATIONS/CISBP/CISBP2core_mouse28.03.20.meme")
      nname <- sapply( motifs, function(x) x["name"])
      names(motifs) <-sub( ".*__", "", nname)
      # extract Tfs expressed well in podocytes
      motifs.podo <- motifs[ grep(paste(TF_MeanMed,collapse="|"), 
                                  names(motifs))]
      # remove artifact
      mmethods <- c("PCC", "EUCL", "SW", "KL", "ALLR", "BHAT", "HELL", "SEUCL", 
                    "MAN", "ALLR_LL", "WEUCL", "WPCC")
      motifs.podo <- motifs.podo[ !names(motifs.podo)%in% c("Fosl1","Fosl2","Creb3l1")]
      motif.compare <- lapply( seq(mmethods) , function(ii){
       X <- compare_motifs( motifs.podo, method=mmethods[ii], min.mean.ic =0)
        colnames(X) <- rownames(X) <- sub(".*__","",colnames(X))
        return(X)
      })
      names(motif.compare) <- mmethods

      
      # heatmap

      gplots::heatmap.2( motif.compare , trace="none", 
                         hclustfun = function(x) hclust(x, method="ward.D2"), 
                         col = viridis(20, option = "B",direction = 1),
                         # breaks = bbreaks , 
                         margins = c( 10,10)      )
      
      # cluster tree
      lapply(seq(motif.compare), function(ii){
        pdf(height = 7, width = 14, file=paste("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/motifCompare/TFmotifs.77_treePlot_",
                                               names(motif.compare)[ii],".pdf",sep = ""))
        # png(height = 8, width = 21, file=paste("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/motifCompare/TFmotifs.77_treePlot_",
        #                                        names(motif.compare)[ii],".png",sep = ""))
        # 
         motif.c <- motif.compare[[ii]]
        tree <- hclust( dist(motif.c) , method = ("ward.D2"))
        cl_members <- cutree(tree = tree, h = 1)
        plot(x = tree, labels =  row.names(tree), cex = 0.5,
             xlab="", main=paste("clustegram of 78 podocyte TFs\n",
                                 names(motif.compare)[ii],sep = ""))
        # rect.hclust(tree = tree, k = max(cl_members), which = 1:max(cl_members), 
        #             border = 1:max(cl_members), cluster = cl_members)
        # abline(h=1, col="red")
        dev.off()
      })
     
      # make a table with simmilar motifs
      cl_table <- as.data.frame(cl_members )
      
      cl_names <- sapply(seq( table(cl_members)) , function(ii){
        cl_members[cl_members==ii][1] } ) 
      cl_table$cl_name <- names(cl_names)[match( cl_table$cl_members , cl_names )] 
      # rownames(cl_table) <- sub(":.*","", rownames(cl_table))
      
     }
   
    ### load nets
      {
      # PMID:36307401 tab 1 and 2
      slitDiaphr <-  read.table(sep = "\t", header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PDScorr/pathway_viz/cytoscapeVIZ/Podo_areas/PPIsource/PMID:36307401_tab1tab2.slitDiaphr.csv")
      names(slitDiaphr)[4] <- "Symbol"
      # Schell et al. PMID:33514561 fig1f
      actCtsklt_fig1f <-read.table(sep = "\t", header = T, fill = T,  "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PDScorr/pathway_viz/cytoscapeVIZ/Podo_areas/C.shell/ActnCtskl_node.csv")
      names(actCtsklt_fig1f)[  names(actCtsklt_fig1f) == "gene_symbol"] <- "Symbol"
      actCtsklt_fig1f <- actCtsklt_fig1f[ names(actCtsklt_fig1f)%in% c("shared.name", "Symbol")]
      
      # PMID:28536193 , figure 1E in Schell et al. 2017
      fcladhsn <-  read.table(sep = "\t", header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PDScorr/pathway_viz/cytoscapeVIZ/Podo_areas/C.shell/FAnet.tsv")
      
      # # adhesome
      # cnsrvAdhesome <-  read.table(sep = "\t", header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PDScorr/pathway_viz/cytoscapeVIZ/Podo_areas/PPIsource/PMID:32147508_SupplTab3.ConsAdhesome.csv")[1:3]
      # fltrdAdhesome <- read.table(sep = "\t", header = T,"/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PDScorr/pathway_viz/cytoscapeVIZ/Podo_areas/PPIsource/PMID:33761352_SupplTab10.fltrdAdhesome.csv")
      # 
      # Adhesome <- union(cnsrvAdhesome , fltrdAdhesome)
      
      # matrisome
      Mtrxsome <- read.table(sep = "\t", fill = T , header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PDScorr/pathway_viz/cytoscapeVIZ/Podo_areas/PPIsource/PMID:33761352_SupplTab7.fltrdMtxsome.csv")
      ggenes <-  unlist(sapply( Mtrxsome$Gene_names , function(X) {
        strsplit( X , split = ";")[[1]][1] }))
      Mtrxsome$Symbol <- fun_homoTO.FROMmouse(ggenes)$MUS
      Mtrxsome <- Mtrxsome[c(1,2,10:13)]
      Mtrxsome$Symbol <- as.character(Mtrxsome$Symbol)
      # ECM updated
      ECM <- read.table(sep="\t", header = T, file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PDScorr/pathway_viz/cytoscapeVIZ/Podo_areas/PPIsource/Network_structure_ECM-updated.csv")
      
      NNets <- c( lapply(PodoPathGSet[1:3], function(xx) data.frame(Symbol=xx)),
                  list( slitDiaphr, actCtsklt_fig1f, fcladhsn, ECM, Mtrxsome ))
      names(NNets) <- c( names(PodoPathGSet)[1:3], "Slit.Diaphr", "Actin.Ctsklt.fig1f",
                         "Focal.Adhesion", "ECM","Matrisome")
    }
   
    
    ### ad TFreg GRN info
   NNets.TFs<- lapply(seq(NNets), function(ii)
     {
     print(ii)
     nett <- NNets[[ii]]
     nett.TFs <- Reduce( rbind , lapply( seq(nrow(nett)), function(jj){
       # print(jj)
       datt <- nett[jj,,drop=F]
       XX <- data.frame( datt,
                         TFreg= if( datt$Symbol %in% rownames(ATACseq_tgenesM.podo)){
                           names( ATACseq_tgenesM.podo[datt$Symbol,])[
                             ATACseq_tgenesM.podo[datt$Symbol,]>0 ] } else "" )
     }))
  
    
    # add TF fcount
    nett.TFs$TFtgCount <- table( nett.TFs$TFreg)[ match( nett.TFs$TFreg ,
                                                           names( table( nett.TFs$TFreg)))] 
    ## add TF cluster/group name
    iid <- sapply(  nett.TFs$TFreg , function(xx){
      X <- which( lapply(rownames(cl_table), function(yy) xx %in% unlist( 
        strsplit(yy, split = ":"))) =="TRUE" )
      if( identical(X, integer(0))) X <- NA
      
      return(X)
    }) 
    nett.TFs$TFgroup <- cl_table$cl_name[ iid ]
    nett.TFs$TFgroup[is.na(nett.TFs$TFgroup)] <- ""
    
    nett.TFs$TFgroup <- sub(":.*","",nett.TFs$TFgroup)
    nett.TFs$TFgroup[nett.TFs$TFgroup!=""] <- paste(  nett.TFs$TFgroup[nett.TFs$TFgroup!=""],
                                                      "clst",sep = "_")
    # add group count name
    nett.TFs$TFgroupCount <- table( nett.TFs$TFgroup)[ match( nett.TFs$TFgroup ,
                                                              names( table( nett.TFs$TFgroup)))]
    return(nett.TFs)
    })
   names(NNets.TFs) <- names(NNets)
   
   saveRDS(NNets.TFs, file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/TFreg/PodoPathGSet.TFs.rda")
    # # remove TFs with less than 3 targets
    # nett.TFs$TFreg <- ifelse( nett.TFs$TFreg %in% names(table(nett.TFs$TFreg))[
    #   table(nett.TFs$TFreg)<3 ]  , "",
    #       nett.TFs$TFreg )
    write.table(nett.TFs , sep = "\t", quote = F, col.names = NA ,
                file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PDScorr/pathway_viz/cytoscapeVIZ/Podo_areas/PPIsource/Network_structure_ECM-updated.TFs.tsv")
    # write tables
    lapply( seq(NNets.TFs), function(ii){
      write.table(NNets.TFs[[ii]], sep = "\t", file= paste( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/PDScorr/pathway_viz/cytoscapeVIZ/Podo_areas/PPIsource/",names(NNets.TFs)[ii], "_TFreg.podo.tsv",sep = ""), quote = F, col.names = NA)
    })
    }

#### which TFs regulate pathways of interest
  {
   
  ATACseq_tgenesM.TFtc.P <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/GRN/ATACseq_GRNprior/ATACseq_TFtc.pvalue.tgenesM_77tf.rda" )
  
  pathG <- as.character( unlist(path5) )
  
  wiki_TFstat <-  Reduce( cbind , lapply( seq(wikiPath_gName), function(ii){
    pathG <- as.character( wikiPath_gName[[ii]] )
    p.adjust( apply( ATACseq_tgenesM.TFtc.P, 2, function(tgenes){
      tgenes <- names(tgenes)[tgenes<0.5]
      phyper( (length( intersect(tgenes, pathG))-1) , length(tgenes) , nrow(ATACseq_tgenesM.TFtc.P), length(pathG ) ,
              lower.tail= FALSE)
  }))
  }) )
  colnames(wiki_TFstat) <- names(wikiPath_gName)
  
  wiki_TFstat[,grep("primary " , colnames(wiki_TFstat),ignore.case = T)]
  
  # kegg
  kegg_TFstat <-  Reduce( cbind , lapply( seq(keggPath_gNames), function(ii){
    pathG <- as.character( keggPath_gNames[[ii]] )
    p.adjust( apply( ATACseq_tgenesM.TFtc.P, 2, function(tgenes){
      tgenes <- names(tgenes)[tgenes<0.5]
      phyper( (length( intersect(tgenes, pathG))-1) , length(tgenes) , nrow(ATACseq_tgenesM.TFtc.P), length(pathG ) ,
              lower.tail= FALSE)
    }))
  }) )
  colnames(kegg_TFstat) <- names(keggPath_gNames)
  
  kegg_TFstat["Wt1",][order(kegg_TFstat["Wt1",])]
 
  # reactome
  react_TFstat <-  Reduce( cbind , lapply( seq(reactPath_gName), function(ii){
    pathG <- as.character( reactPath_gName[[ii]] )
    p.adjust( apply( ATACseq_tgenesM.TFtc.P, 2, function(tgenes){
      tgenes <- names(tgenes)[tgenes<0.5]
      phyper( (length( intersect(tgenes, pathG))-1) , length(tgenes) , nrow(ATACseq_tgenesM.TFtc.P), length(pathG ) ,
              lower.tail= FALSE)
    }))
  }) )
  colnames(react_TFstat) <- names(reactPath_gName)
  
  head(react_TFstat["Wt1",][order(react_TFstat["Wt1",])])

  }

#### use dispersion estimate
  {
  
}
    
  