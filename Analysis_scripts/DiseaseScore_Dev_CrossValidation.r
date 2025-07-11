# ================================================================= #
# ### script for validation and description of the damage score ### # 
# ================================================================= #
mallinfo::malloc.trim()
gc()

options( connectionObserver = NULL )
.libPaths(c("/home/tim_nevelsk/R/x86_64-pc-linux-gnu-library/4.0", "/media/tim_nevelsk/WD_tim/SOFT/R"))

library( AUCell)
library( ggplot2)
library( ggrepel)
library( ggthemes)
library( plyr)
library( biomaRt)
library( Seurat)
library( ggpubr)

mart_mouse <- useMart( "ensembl",dataset="mmusculus_gene_ensembl" )

tx2gene <- biomaRt::getBM(attributes=c( "ensembl_gene_id", "external_gene_name", 
                              "uniprotsptrembl",  "uniprotswissprot"),  mart = mart_mouse)

#### load functions
# generate damage signature
source("https://raw.githubusercontent.com/PauUng/HepatocyteDamageScore/master/SharedFunctions.R")
# calculate damage score
source("https://raw.githubusercontent.com/nevelsk90/R_scripts/master/usefulRfunc.r")

DS_all <- read.table( header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")

allPodoGenes <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/SCSN_allPodoGenes.rda")

### load subsampled sc data
listSCSN.1K.sampl <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN_samples.1K.22.12.23.rda")
## triangle_vec_to_matrix
triangle_vec_to_matrix <- function(vec,
                                   diag_val = 0,
                                   side = c("lower", "upper"),
                                   order = c("column", "row")) {
  side  <- match.arg(side)   # which half is encoded?
  order <- match.arg(order)  # in what order are the elements stored?
  
  k <- length(vec)           # number of elements supplied
  n <- (1 + sqrt(1 + 8 * k)) / 2   # solve k = n · (n-1) / 2
  if (n != floor(n)) stop("Vector length must be n*(n-1)/2 for some integer n.")
  n <- as.integer(n)
  
  mat <- matrix(diag_val, n, n)    # initialise with chosen diagonal value
  
  ## Decide the index sequence that fills the triangle -----------------------
  if (side == "lower") {
    idx <- lower.tri(mat, diag = FALSE)
    if (order == "row") {                 # convert to row-wise ordering
      idx <- idx[order(row(mat)[idx], col(mat)[idx])]
    }
    mat[idx] <- vec
    mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]   # mirror to upper part
  } else {                                # side == "upper"
    idx <- upper.tri(mat, diag = FALSE)
    if (order == "row") {
      idx <- idx[order(col(mat)[idx], row(mat)[idx])]
    }
    mat[idx] <- vec
    mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]   # mirror to lower part
  }
  
  mat
}
#### load the data ####
    
      # genes expressed in sc.sn
      scsnRNAseq_mean.rank <- readRDS( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/scsnRNAseq.podo_gene.inter_mean.rank.rda")
      
      ### load DE results
      FSGS_MA_DE <- readRDS(  file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DE/MA_DElist.04.04.22.rda")
      FSGS_bulk_DE <- readRDS(  file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DE/bulk_DElist.04.04.22.rda")
      FSGS_sc_DE <- readRDS(  file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DE/sc_DElist.04.04.22.rda")
      
      names(FSGS_sc_DE) <- paste0(c("Coq2","Pdss2","GSE127235" ,"GSE139506",
                                 "btbr", "cd2ap" ,"doxo", "nephr.D1", "nephr.D5", 
                                 "GSE164273","GSE174013","GSE174102", "Nphs2","Wt1" ),"_sc")
      DE_list <- c( FSGS_MA_DE , FSGS_bulk_DE , FSGS_sc_DE)
  
      
    ### load MA expression matrices
    ll <- list.files( full.names = T, path= "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/PDS/MA/expr_tabs")
    MA_explist <- lapply( ll, readRDS )
    names( MA_explist ) <- sub( "_.*","",basename(ll))
  
    ### load bulk expression matrices
    ll <- list.files( full.names = T, path= "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/PDS/bulkRNAseq/expr_tabs")
    XX <- lapply( ll, readRDS )
    bulk_explist <- sapply( XX,"[[",1)
    bulk_annlist <- sapply( XX,"[[",2)
    names(bulk_explist) <-names(bulk_annlist) <- sub( "_.*","",basename(ll))
  
    ### load expression matrices
    # names(SCSNdata_list_sub) <- names(SCSNdata_list12)
    # saveRDS( SCSNdata_list_sub,  file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/PDS/sc/expr_tabs/SCSNdata_list12_sub.rda")
    SCSNdata_list_sub <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/PDS/sc/expr_tabs/SCSNdata_list12_sub.rda")
    
   
    # extract counts, downsample for faster results
    sc_explist <- lapply( seq( SCSNdata_list_sub ) , function( ii )
      {
      print(ii)
      newSeu <- SCSNdata_list_sub[[ii]]
  
      if( !is.null(newSeu@assays$RNA)) {
        newSeu@assays$RNA@counts} else {
          newSeu@assays$data@counts
        }
    } )
    names(sc_explist) <- names(SCSNdata_list_sub)
    
  #### combine all
    # expression Matrices
    explist <- c( MA_explist, bulk_explist, sc_explist )
    saveRDS( explist , file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/PDS/explist_44.ma.bulk.sc.rda")

    # annotations
    annotlist <- c( lapply( MA_explist , function(X) sub(".*__", "",colnames(X))),
                    lapply( bulk_explist , function(X) sub(".*__", "",colnames(X))) ,
                    lapply( seq(SCSNdata_list_sub) , function(ii){
                      SCSNdata_list_sub[[ii]]$groups
                    })  )

    annotlist_gtype <-c( lapply( MA_explist , function(X) ifelse( sub(".*__", "",colnames(X))=="mutant",
                                                "control", sub(".*__", "",colnames(X)))),
                              lapply( bulk_explist , function(X) sub(".*__", "",colnames(X))) ,
                              lapply( seq(SCSNdata_list_sub) , function(ii){
                                SCSNdata_list_sub[[ii]]$gtype
                              })  )
  
    annotlist_gtype <- lapply( seq(annotlist_gtype), function(ii){
      X <- annotlist_gtype[[ii]]
      X <- ifelse( X %in% c("control","wt"), "control","experiment")
    })
    # equate names of experiments
    names(annotlist) <-names(annotlist_gtype) <- names(explist)
    # ll <- list( annotlist, annotlist_gtype , annotlist_scSample)
    # names( ll ) <- c( "annotlist", "annotlist_gtype" , "annotlist_scSample")
    saveRDS( ll, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/PDS/explist_44.ma.bulk.sc_annot.rda" )
  
#### cross validate platforms base DS + pseudo-time based DS ####  

    ### make signatures combining results of DE
      {
        # get list of genes expressed in sc|sn RNAseq
        SCSNpodogenes <- readRDS(  file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DE/podoSNgenes.rda")
        
        ### generated damage signatures
        MA_PDS <- damage_signature.func( FSGS_MA_DE, thresh = 0.75 , 
                                         geneFilter = SCSNpodogenes)
        bulk_PDS <- damage_signature.func( FSGS_bulk_DE, thresh = 0.75, 
                                           geneFilter = SCSNpodogenes)
        sc_PDS <- damage_signature.func( FSGS_sc_DE , thresh = 0.75, 
                                         geneFilter = SCSNpodogenes)
        DS_all <- damage_signature.func( DE_list , thresh = 0.75, 
                                         geneFilter = SCSNpodogenes )
        DS_all.2 <- damage_signature.func( DE_list , thresh = 0.75, 
                                         geneFilter = names(scsnRNAseq_mean.rank) )
        bulkMA_PDS <- damage_signature.func( c(FSGS_MA_DE , FSGS_bulk_DE ), 
                                             thresh = 0.75, geneFilter = SCSNpodogenes )
        
        # use genes within top 2000 mean rank of expr. in sc.sn podocytes
        scsnRNAseq_mean.rank <- readRDS( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/scsnRNAseq.podo_gene.inter_mean.rank.rda")
        DS_all.2K <- damage_signature.func( DE_list , thresh = 0.75, 
                                            geneFilter = names(scsnRNAseq_mean.rank[1:2000]  ))
        
        ### signature based on ctype markers
        podo_sc_DE <- readRDS(  file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DE/podoSCSN_DElist.04.04.22.rda")
        ctype_PDS <- damage_signature.func( podo_sc_DE, thresh = 0.75, 
                                            geneFilter = SCSNpodogenes )
        # reverse fold change so podocytes are "control"
        ctype_PDS$direction_foldchange <- -1*ctype_PDS$direction_foldchange
        
        ### sc/sn ptime signature
        pseudo <- readRDS( file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/SCSNpseudoSig.rds")
        
        ptime_PDS <-  damage_signature.func( list(pseudo,pseudo), thresh = 0 , 
                                             geneFilter = SCSNpodogenes )
          
        ###  previously used PDS
        old_PDS <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/FSGS_markers0.75_36stud.rda")[[1]]
        old_PDS$gene_symbol <- rownames(old_PDS)
        # convert to lower case gene names
        allGenes <- Reduce( union , lapply(explist, rownames))
        old_PDS$gene_symbol <- allGenes[ match( old_PDS$gene_symbol , 
                                                toupper( allGenes) )]
        colnames(old_PDS)[1] <- "mean_rank"
        old_PDS$direction_foldchange <- ifelse(old_PDS$meanLFC<0, -1 , 1)
        old_PDS <- old_PDS[order(old_PDS$mean_rank),]
       
        
        # list all signatures
        DS_cv.list <- list( MA_PDS , bulk_PDS , sc_PDS ,bulkMA_PDS , 
                             DS_all , old_PDS , ctype_PDS , ptime_PDS)
        names(DS_cv.list) <- c( "MA_DS" , "bulk_DS" , "sc_DS" , "bulkMA_DS" , 
                                 "all_DS" , "old_DS" , "ctype_DS","ptime_PDS" )
        saveRDS(DS_cv.list , file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/platforms/DS_platform.list.03.06.22.tsv")
        
        # save tables
        write.table( MA_PDS , sep = "\t", row.names = F, 
                    file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/platforms/MA_PDS.tsv")
        write.table( bulk_PDS , sep = "\t", row.names = F, 
                    file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/platforms/bulk_PDS.tsv")
        write.table( sc_PDS , sep = "\t", row.names = F, 
                    file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/platforms/sc_PDS.tsv")
        # write.table( DS_all , sep = "\t", row.names = F, 
        #             file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/platforms/DamageSignatures/DS_all.02.06.2022.tsv")
        write.table( bulkMA_PDS , sep = "\t", row.names = F, 
                    file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/platforms/bulkMA_PDS.tsv")
        write.table( ctype_PDS , sep = "\t", row.names = F, 
                     file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/platforms/ctype_PDS.tsv")
        write.table( old_PDS , sep = "\t", row.names = F, 
                     file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/platforms/old2020_PDS.tsv")
        write.table( ptime_PDS , sep = "\t", row.names = F, 
                     file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/platforms/ptime_PDS.tsv")
        write.table( DS_all.2K , sep = "\t", row.names = F, 
                     file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.19.09.2023.tsv")
        write.table( DS_all.2 , sep = "\t", row.names = F, 
                     file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")
      }
      
    ### make a corr plot of disease signatures
      {
        DS_cv.list <- readRDS(file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/platforms/DS_platform.list.rda")
        
        ### combine scores in one df
        gg <-  Reduce( union, lapply( DS_cv.list , function(X) X$gene_symbol))
       
        DS_union <- Reduce( cbind, lapply( DS_cv.list, function(X){
          X$mean_rank[ match( gg , (X$gene_symbol))] } ))
        DS_unionLFC <- Reduce( cbind, lapply( DS_cv.list, function(X){
          X$direction_foldchange[ match( gg , (X$gene_symbol))] } ))
        
        colnames(DS_union) <- colnames(DS_unionLFC) <- names(DS_cv.list)
        
        ### calculate spearman correlation for p-values
        spCor_rank <- cor(DS_union , method = "spearman", 
                    use = "pairwise.complete.obs")
        ### calculate kendal correlation for binary lfc vectors
        kCor_lfc <- cor(DS_unionLFC , method = "kendall", 
                    use = "pairwise.complete.obs")
        ### plot 
        corrplot::corrplot( spCor_rank , p.mat = spCor_rank ,tl.col="black",
                            insig =  "p-value" ,  sig.level=0 , 
                            addCoef.col="red")
        corrplot::corrplot( kCor_lfc , p.mat = kCor_lfc ,tl.col="black",
                            insig =  "p-value" ,  sig.level=0 , 
                            addCoef.col="red")
      }

    ### calculate the damage score, using various signatures
      {
        DS_cv.list <- readRDS(file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/platforms/DS_platform.list.tsv")
        
        ### calculate DS using platform specific signatures
        # iterate over platforms
        PDS.platforms <- lapply( seq(DS_cv.list), 
                                    function( id )
          {
            print(id)

            # iterate over datasets
            lapply( seq(explist) , function(ii )
              {
            DS_calc.func( exprMatrices = explist[[ii]], ceilThrsh = 0.05 , 
                          DSignature= DS_cv.list[[id]] , ntop=42  )
          })
          
        } )
        
        # name platforms 
        names(PDS.platforms) <- sub( "_",".",names(DS_cv.list) )
      
      ### make individual DS plots 
      lapply( seq(PDS.platforms), function( jj )
        {
        print(jj)
        DSlist <- PDS.platforms[[jj]]
        
        PDS_plot <-lapply( seq(DSlist), function(ii )
          {
          GSE <-  names(DSlist)[ii]  # get GSE id of a current experiment
          toPlot <- data.frame( score= DSlist[[ii]] , gt =annotlist[[ii]])
          if ( length(DSlist)<50) {
            gg <- ggplot2::ggplot( toPlot , aes( y=score, x=gt, color=gt)) +
              scale_color_colorblind() +  geom_jitter(width = 0.1, size=3) + 
              theme_bw() + theme( text = element_text(size = 18) , legend.position = "none") +
              ggtitle(GSE) +
              stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean,
                           geom = "crossbar", width = .2, color = "red")  
          } else {
            # make a density plot
            gg <-  ggplot2::ggplot( toPlot , aes( x=score, color=gt)) + 
              scale_color_colorblind() +  geom_density(size=1.5) + 
              theme_bw() +  ggtitle( GSE ) +
              theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                     text = element_text(size = 18)) +
              geom_vline(data=ddply( toPlot, "gt",
                                     summarise, grp.mean=mean(score)), 
                         aes(xintercept=grp.mean, color=gt), linetype="dashed")
          }
          return(gg)
          
        })
        
        ppg1 <- cowplot::plot_grid( PDS_plot[1:14], nrow = 4)
        ppg2 <- cowplot::plot_grid( PDS_plot[15:32],nrow = 5)
        ppg3  <- cowplot::plot_grid( PDS_plot[32:44], nrow = 4)

      pdf(height = 12, width = 12, file = paste( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/platforms/",
                          names(PDS.platforms)[jj],"_plot.pdf", sep = ""))
        print(ppg , ppg2, ppg3)
      dev.off()

      png(height = 1000, width = 1000, file = paste( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/platforms/",
                                                     names(PDS.platforms)[jj],"_plot.png", sep = ""))
        print(ppg , ppg2, ppg3)
      dev.off()
        } )
      
      saveRDS(PDS.platforms , file = " PDS.platforms.ctype.ptime.rda")
     
      ### combine all studies by platform DS plots 
        {
        platCV.onePlot <- Reduce( rbind, 
                                  lapply( which(!names(PDS.platforms)%in%
                                                  # exclude not relevant signatures
                                                  c("old_DS","bulkMA_DS")), 
                                          function( jj )
          {
          print(jj)
          DSlist <- PDS.platforms[[jj]]
          
          PDS_plot <- Reduce( rbind , lapply( seq(DSlist), function(ii ){
            nname <- names(explist)[ii]
            
            # agregate sc to pseudo-bulk if needed and create ggplot-ready long df
            if( grepl( "_sc", nname )){
              XX<- data.frame( score= DSlist[[ii]] ,
                               centered_score= scale( DSlist[[ii]] , scale = F),
                               sample = as.factor(annotlist_scSample[[nname]]) ,
                               gt =annotlist_gtype[[ii]] , 
                               EXPsource= "sc" )  
              
              XX <- aggregate( c )
            }  else  XX <- data.frame( sample = sub( "__.*", "", names(DSlist[[ii]])) ,
                                   gt =annotlist_gtype[[ii]] , 
                                  EXPsource= ifelse( ii%in% 1:14, "MA", "bulk"),
                                  score= DSlist[[ii]] ,
                                  centered_score= scale( DSlist[[ii]] , scale = F) )  
            return( XX )
          } ) )
          
          PDS_plot$DSsource=names(PDS.platforms)[jj]

        
                                   
          return(PDS_plot)
          
        } ) )
        
          # reorder factor
          platCV.onePlot$DSsource <- factor( as.factor(platCV.onePlot$DSsource), 
                               levels = c("all_DS", "MA_DS", "bulk_DS", "sc_DS", "ctype_DS", "ptime_PDS" ))

          
        ppg <- ggplot2::ggplot( platCV.onePlot , 
                                aes( y=centered_score, x=EXPsource , color=gt)) +
          scale_color_colorblind() +  geom_boxplot(lwd=1.2) + 
          coord_cartesian(ylim=c(-0.3, 0.3)) +
          theme_bw() + theme( text = element_text(size = 24) , axis.title.x=element_blank(),
                              axis.text.x = element_text(angle = 45, hjust = 1)) +
          facet_grid(cols = vars(DSsource)) + stat_compare_means( size = 6, method = "wilcox.test",
                                                                  aes(color=gt)  ,
                                                                  label.y.npc = c(0.76,0.72,0.68) ,
                                                                  label.x.npc= c('center','center','center') )
        
        pdf( height = 6, width = 12, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/platforms/platformsCV_onePlot_boxplot.pdf" )
        print( ppg )
        dev.off()
        
      }
      
    }
  
    ### correlate PDS calculated with different signatures 
      {
      PDS.platforms_corr <-  lapply( 1:length(explist), function(jj) 
        {
        datt <-  Reduce( cbind, lapply( PDS.platforms , function(X) {
          X[[jj]] 
        } ) )
        colnames(datt) <- names(PDS.platforms)
        pCor <- cor(datt , method = "spearman", 
                    use = "pairwise.complete.obs")
        return(pCor)
      })
      names(PDS.platforms_corr) <- names(explist)
      pCor <- Reduce("+", PDS.platforms_corr) / length(PDS.platforms_corr)

      ### plot average 
      corrplot::corrplot( pCor , p.mat = pCor ,tl.col="black",
                            insig =  "p-value" ,  sig.level=0 , 
                            title = "average correlation pf PDS across sc-sn datasets" ,
                            addCoef.col="red", mar=c(0,0,1,0) )
      ### scatter plots for sc data 
      # pdf( width = 8, height = 6, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PseudoTimeTrajectory/DSall.VS.ptime_SCSN.pCorr.pdf")
      png( width = 400, height = 300, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PseudoTimeTrajectory/DSall.VS.ptime_SCSN.pCorr.png")
      
      lapply( 1:12, function(ii){
        
        plot( x = PDS.platforms[[8]][[32+ii]], y=PDS.platforms[[5]][[32+ii]] ,
                xlab= "pseudo-time DS" , 
                ylab= "supervised DS" , main=names(explist)[32+ii])
        abline( lm(PDS.platforms[[5]][[32+ii]] ~ PDS.platforms[[8]][[32+ii]]) , 
                col="red" , lwd=3)
          legend( "bottomright" ,
                  legend = paste( "rho = ", round( cor(PDS.platforms[[8]][[32+ii]], 
                                      y=PDS.platforms[[5]][[32+ii]] ), digits = 2 ), sep = ""  ) )
      })
      dev.off()
      
    }
    
    ### test questionable datasets
      {
      ### 117571 MA
      {
        DS <- MA_DSall[[which(names( MA_explist )=="GSE117571")]]
        
        toPlot <- data.frame( score= DS , gt =sub(".*__","",names(DS)) , 
                              groups=paste( sub(".*__","",names(DS)),
                                            c(rep("gloms",4),rep("podo.cult",8)) ,sep = "_"))
        
       ggplot2::ggplot( toPlot , aes( y=score, x=groups, color=groups)) +
          scale_color_colorblind() +  geom_jitter(width = 0.1, size=5) + 
          theme_bw() + theme( text = element_text(size = 18) , 
                              axis.text.x = element_text(angle = 45, hjust = +1) ,
                              legend.position = "none") +
          ggtitle("GSE117571") +
          stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean,
                       geom = "crossbar", width = .2, color = "red") 
      }
    
      ### 168676 bulk
      {
        DS <- PDS.all[[which(names( PDS.all )=="GSE168676")]]
        annot <- bulk_annlist[[which(names( bulk_annlist )=="GSE168676")]]
        
        toPlot <- data.frame( score= DS , gt =sub(".*__","",names(DS)) , 
                              groups= annot$groups )
        
        ggplot2::ggplot( toPlot , aes( y=score, x=groups, color=groups)) +
          scale_color_colorblind() +  geom_jitter(width = 0.1, size=5) + 
          theme_bw() + theme( text = element_text(size = 18) , 
                              axis.text.x = element_text(angle = 45, hjust = +1) ,
                              legend.position = "none") +
          ggtitle("GSE168676") +
          stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean,
                       geom = "crossbar", width = .2, color = "red")
      }
      
      ### 79291 bulk 
      {
        DS <- bulk_DSall[[which(names( bulk_explist )=="GSE79291")]]
        annot <- bulk_annlist[[which(names( bulk_annlist )=="GSE79291")]]
        
        toPlot <- data.frame( score= DS , gt =sub(".*__","",names(DS)) , 
                              groups= annot$groups )
        
        ggplot2::ggplot( toPlot , aes( y=score, x=groups, color=groups)) +
          scale_color_colorblind() +  geom_jitter(width = 0.1, size=5) + 
          theme_bw() + theme( text = element_text(size = 18) , 
                              axis.text.x = element_text(angle = 45, hjust = +1) ,
                              legend.position = "none") +
          ggtitle("GSE79291") +
          stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean,
                       geom = "crossbar", width = .2, color = "red")
      }
    }
    

#### use ML to predict Alb/Cre or Class labels #### 
#  see how well each gene predicts FSGS
### correlate mean g.expr lvl with FSGS
  {
  
  # gene sets
  gSets.DS_all <- lapply( seq(DS_all$gene_symbol), function(ii) DS_all$gene_symbol[ii])
  names(gSets.DS_all) <- DS_all$gene_symbol
  
  ## sc/sn datasets
    {

    # laod annotation
    library(ggpubr)
    annot_tab <- read.table(sep = "\t",header = T, "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/Sample_Names_KFO.csv")
    annot_tab$group <- paste(annot_tab$Genotype,annot_tab$Age_weeks,sep = "_")
    
    ### load KFO data 
    listSCSN.PDS <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN.PDSlist_23.12.2022.rda")
    listSCSN.PDS[[2]]$sample <- sub( "SID", "", listSCSN.PDS[[2]]$sample)
    listKFO_sub.smpl <- lapply( 1:3 , function(ii){
      datt <- listSCSN.PDS[[ii]]
      Idents(datt)<- datt$sample
      
      newDatt <- subset( datt, downsample=250 )
      
      # Idents(newDatt)<- newDatt$gtypeDE
      # 
      # newDatt <- subset( newDatt, downsample=1000 )
      return(newDatt)
    })
      
    # calculate AUCell for each gene
    # split computation in 2 batches due to the size
    KFO_AUCell.indGene.b1 <- lapply( seq(listKFO_sub.smpl), function(ii){
        print(ii)
        
      # remove genes that expressed in less than  0.005 Ncells
      expr <- listKFO_sub.smpl[[ii]]@assays$RNA@counts
      # expr <- expr[ rowSums( expr > 0 ) > ncol(expr)*0.005 , ]
      
       calculateAUCcell( gSets.DS_all[1:200] , ddata= expr , AUCthrsh = 0.05 )
      })
    KFO_AUCell.indGene.b2 <- lapply( seq(listKFO_sub.smpl), function(ii){
      print(ii)
      
      # remove genes that expressed in less than  0.005 Ncells
      expr <- listKFO_sub.smpl[[ii]]@assays$RNA@counts
      # expr <- expr[ rowSums( expr > 0 ) > ncol(expr)*0.005 , ]
      
      
      calculateAUCcell( gSets.DS_all[201:length(gSets.DS_all)] , ddata= expr , AUCthrsh = 0.05 )
    })

    # combine metadata from KFO experiments 
    KFO_AUCell.indGene <- Reduce( rbind, lapply( seq(listKFO_sub.smpl), function(ii) {
     print(ii)
      XX <- cbind( t( KFO_AUCell.indGene.b1[[ii]]) ,
                   t( KFO_AUCell.indGene.b2[[ii]]) )
      XX <- XX[ ,match( DS_all$gene_symbol, colnames(XX))]
      XX[is.na(XX)] <- 0
      colnames(XX) <- DS_all$gene_symbol
      # add meta
      x <- listKFO_sub.smpl[[ii]]@meta.data 
      x$group <- paste( x$gtypeDE , x$group, sep = "."  )
      XX <- cbind( x[,c( "group","sample")] , XX )

            return(XX )
    }))
    saveRDS(KFO_AUCell.indGene, file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/snKFO.pergeneAUcell.05.rds")
    
    
    
    # treat carefully 21 week Pdss2 samples since they have only per group measurements
    datt1 <- KFO_AUCell.indGene[ KFO_AUCell.indGene$sample %in% c( "146985", "146986", "143485" , "143486") ,]
    # aggregate
    aggPDS1 <- aggregate( .~group , data= datt1[ , -2], 
                          mean )
    aggPDS1 <- aggPDS1[rep(seq_len(nrow(aggPDS1)), each = 2), ]
    aggPDS1$group <-  c( "146985", "146986", "143485" , "143486")
    colnames(aggPDS1)[1] <- "sample"
    # the rest of samples
    datt2 <- KFO_AUCell.indGene[ !(KFO_AUCell.indGene$sample %in%  c("146985", "146986", "143485" , "143486")), ]
    aggPDS2 <- aggregate( .~sample , data= datt2[ , -1], 
                          mean )
    
    aggPDS <- rbind(aggPDS2, aggPDS1)
    
    aggPDS.sc <- cbind( AlbCrRatio = annot_tab$AlbCrRatio[ match( sub("SID","" ,aggPDS$sample) , 
                                                      annot_tab$CCG_Sample_ID)],
                    group = annot_tab$group[ match( sub("SID","" ,aggPDS$sample ), 
                                            annot_tab$CCG_Sample_ID)] ,
                    gtype = as.factor(annot_tab$Genotype[ match( sub("SID","" ,aggPDS$sample ), 
                                                         annot_tab$CCG_Sample_ID)] ) , 
                    aggPDS)
    

   #  DS_all.AUcell <-  cor(aggPDS.sc$AlbCrRatio, aggPDS.sc[ , 5:ncol(aggPDS.sc)] , 
   #                        method = "spearman", use = "pairwise.complete")
   #  DS_all.AUcell <-   setNames( as.vector(DS_all.AUcell) , colnames(DS_all.AUcell) )
   # DS_all.AUcell <- DS_all.AUcell[ order(- abs(DS_all.AUcell))]
   # saveRDS(DS_all.AUcell, file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.AUcell.FSGScorr.rds")
   # 

    
  }
  
  ## bulk datasets
    {
    stud <- c( "GSE117571" ,  "GSE108629" , "GSE17709" , "GSE117987" , 
               "GSE126217", "GSE154955" , "GSE112116" , "GSE131266", 
               "GSE134327", "GSE110092", "GSE77717", "KFO.Wt1" )
    # annotation
    annot_bulkStage <- readRDS(file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Disease_stages/annot_bulkStage.rda")
    # expression
    listBULK <- c( MA_explist , bulk_explist)
    
    # calculate PDS
    bulk_AUCell.indGene <- lapply( seq( listBULK ) , function(ii  )
                                {
                              calculateAUCcell( geneSets = gSets.DS_all , 
                                                ddata= listBULK[[ii]] ,
                                                AUCthrsh=0.2)
                                } )
    names(bulk_AUCell.indGene) <- names(listBULK)
    # GSE117571 use only gloms
    bulk_AUCell.indGene[["GSE117571"]] <- bulk_AUCell.indGene[["GSE117571"]][,1:4]
    
 ### studies with avalable proteinuria
    stud <- c( "GSE117571" ,  "GSE108629" , "GSE17709" , "GSE117987" , 
               "GSE126217", "GSE154955" , "GSE112116" , "GSE131266", 
               "GSE134327", "GSE110092", "GSE77717", "KFO.Wt1" )
    # aggregate data
    datMean <- Reduce( rbind , lapply( seq(stud) , 
                                      function( ii ){
                                        print(ii)
                                        id <- stud[ii]
                                        score <- bulk_AUCell.indGene[[id]]
                                        score <- score[ match(names(gSets.DS_all), rownames(score)),]
                                        rownames(score) <- names(gSets.DS_all)
                                        score[is.na(score)] <- 0
                                        # if a multistaged study - use the prepared annotation list,
                                        # otherwise extract annotation from sample names
                                        if( id %in% names(annot_bulkStage)) {
                                          annot <- annot_bulkStage[[id]]$groups[ 
                                            match( sub( "__.*", "", colnames(score)), 
                                                   rownames( annot_bulkStage[[id]] ))]
                                        } else {
                                          annot <- sub( ".*__", "", colnames(score))
                                        }
                                        
                                        XX <- data.frame(groups= as.factor(annot)  , t(score))
                                        XX <- aggregate( .~groups, data=XX , FUN=mean)
                                        XX <-  data.frame( study=id , XX)
                                        
                                        return(XX)
                                      }))
    
    ### create extra-annotation table
    XX <-  data.frame( 
      # add proteinuria data
      AlbCrRatio = c(0.1, 80, 
                       0.1, 0.2, 100,
                       0.04, 0.4,
                       0.01, 10,
                       0.01, 20,
                       0.005, 0.02, 0.065,
                       0.01, 1.5, # GSE112116
                       0.01, 5.5 , # GSE131266
                       0.1, 3,
                       0.04, 0.125,
                       0.1, 0.45, 
                       0.033 , 0.486 , 0.019 , 0.519 ) ,  # group  Wt1h.d. KFO
    
    # add platform type annotation
    platform = c( rep( "MA" , 2) , rep( "MA" ,3 ) ,  rep( "MA" , 2) , 
                           rep(  "bulk" , 2) , rep( "bulk", 2) ,  rep( "bulk" , 3),
                           rep(  "MA" , 2) , rep( "MA", 2) , rep( "bulk", 2) ,  
                           rep( "bulk", 2) , rep(  "bulk", 2), rep(  "bulk", 4)), 
    # add stage
    gtype = ifelse( datMean$groups %in% c("experiment","mutant", "ko_4w","D9 after ADR injection",
                                          "LMB2day4", "ko_12w", "D14 after ADR injection","LMB2day7"),
                             "experiment", "control"))
  
    ##combine annotation and AUCell
    aggPDS.bulk <- cbind( XX, datMean )
    
  }
 
  ### make lm to find predictors
    {
      # datt <- aggPDS.sc
      datt.sc <- aggPDS.sc[ , !colnames(aggPDS.sc) %in% 
                          c("group","gtype","sample","study","PDS.42.005",
                            "platform", "groups")] 
      datt.bulk <- aggPDS.bulk[ , !colnames(aggPDS.bulk) %in% 
                             c("group","gtype","sample","study","PDS.42.005",
                               "platform", "groups")]
      datT <- rbind( datt.bulk, datt.sc )

      # cv.glmnet
      require(caret)
      ll <- list(  datt.bulk, datt.sc, datT)
      cv.lasso_list <- lapply( 1:3, function(ii){
        datt <- ll[[ii]]
        datt$AlbCrRatio <- log(datt$AlbCrRatio)
        datt <- datt[!is.na(datt$AlbCrRatio),]
        datt <- as.matrix(datt)
        cv.lasso_res <- cv.glmnet( y=datt[,1], x=datt[,-1], 
                                   type.measure="mse", alpha=1 , 
                                   family="gaussian", nfolds=5)
        outcome <- coef(cv.lasso_res, s=cv.lasso_res$lambda.min)
        cv.lasso_lm.lmin <-  outcome[outcome[,1]!=0,][-1]
        return(cv.lasso_lm.lmin)

      })
      names(cv.lasso_list) <- c(  "bulk", "sc.KFO", "bulkANDscKFO")

   
      # kinda stability
      stabcv.lasso_list <- lapply( 1:3, function(ii){
        datt <- ll[[ii]]
        datt$AlbCrRatio <- log(datt$AlbCrRatio)
        datt <- datt[!is.na(datt$AlbCrRatio),]
        datt <- as.matrix(datt)
        freqq <- lapply( 1:100, function(jj){
          cv.lasso_res <- cv.glmnet( y=datt[,1], x=datt[,-1], 
                                     type.measure="mse", alpha=1 , 
                                     family="gaussian", nfolds=4)
          outcome <- coef(cv.lasso_res, s=cv.lasso_res$lambda.min)
          cv.lasso_lm.lmin <-  outcome[outcome[,1]!=0,][-1]
        })
        freqq <- table(unlist(lapply(freqq, names)))
        
        return(freqq)
        
      })
      
      
      # coeff
      X <-  Reduce( c, lapply( seq( stabcv.lasso_list[[2]]), function(ii ){
        datt<-  datt.sc
        datt$AlbCrRatio <- log(datt$AlbCrRatio)
        datt <- datt[!is.na(datt$AlbCrRatio),]
        datt <- datt[,colnames(datt)%in% c("AlbCrRatio", names(stabcv.lasso_list[[2]])[ii])]
        coefficients( lm( AlbCrRatio ~ ., data =datt ) )[2]
        }))
      
      
      names(cv.lasso_list) <- c(  "bulk", "sc.KFO", "bulkANDscKFO")
      nn<- intersect( names(freqq.5[freqq.5>90]),names(freqq[freqq.5>90]))
      datt.29 <- as.data.frame( datt[,c("AlbCrRatio",nn)] )
      lasso.sc.29 <- glm(AlbCrRatio ~ . , data = datt.29, 
                             
                             family="gaussian")
       summary(lasso.sc.29)
      X1 <-  coefficients(lasso.sc.29)[-1]
     
      
    }

}

###  predict control VS experiment 
  { 
    SCSN_PDSlist_limit <- readRDS( file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN.PDS.42limit.rda")
    
    ### sc
    {
      listSCSN.PDS_sub.gt <- lapply( seq(SCSN_PDSlist_limit), function(ii){
        datt <- listSCSN.PDS[[ii]]
        Idents( datt ) <- datt$gtypeDE
        datt <- subset( datt, downsample=500)
        return(datt)
      })
     
      # calculate AUCell
      sc_AUCell.indGene.b1 <- lapply( seq(listSCSN.PDS_sub.gt), function(ii){
        print(ii)
        
        # remove genes that expressed in less than  0.005 Ncells
        expr <- listSCSN.PDS_sub.gt[[ii]]@assays$RNA@counts
        expr <- expr[ rowSums( expr > 0 ) > ncol(expr)*0.005 , ]
        
        calculateAUCcell( gSets.DS_all[1:150] , ddata= expr )
      })
      sc_AUCell.indGene.b2 <- lapply( seq(listSCSN.PDS_sub.gt), function(ii){
        print(ii)
        
        # remove genes that expressed in less than  0.005 Ncells
        expr <- listSCSN.PDS_sub.gt[[ii]]@assays$RNA@counts
        expr <- expr[ rowSums( expr > 0 ) > ncol(expr)*0.005 , ]
        
        calculateAUCcell( gSets.DS_all[151:300] , ddata= expr )
      })
      sc_AUCell.indGene.b3 <- lapply( seq(listSCSN.PDS_sub.gt), function(ii){
        print(ii)
      
        # remove genes that expressed in less than  0.005 Ncells
        expr <- listSCSN.PDS_sub.gt[[ii]]@assays$RNA@counts
        expr <- expr[ rowSums( expr > 0 ) > ncol(expr)*0.005 , ]
        
        calculateAUCcell( gSets.DS_all[301:length(gSets.DS_all)] , ddata= expr )
      })
      
      # combine AUCell calculation
      sc_AUCell.indGene <- Reduce( rbind , lapply( seq(listSCSN.PDS_sub.gt), function(ii) {
        print(ii)
        XX <- cbind( t( sc_AUCell.indGene.b1[[ii]]) ,
                     t( sc_AUCell.indGene.b2[[ii]]) , 
                     t (sc_AUCell.indGene.b3[[ii]]) )
        XX <- XX[ , match( names(gSets.DS_all), colnames(XX))]
        colnames(XX) <- names(gSets.DS_all)
        
        # combine sample and genotype
        XX <- cbind.data.frame( study= names(SCSN_PDSlist_limit)[ii],
                                sample=listSCSN.PDS_sub.gt[[ii]]@meta.data$sample ,
                                gtypeDE = listSCSN.PDS_sub.gt[[ii]]@meta.data$gtypeDE,
                                group = listSCSN.PDS_sub.gt[[ii]]@meta.data$group,
                                PDS.42 = listSCSN.PDS_sub.gt[[ii]]@meta.data$PDS.42 ,
                                lasso.scKFO = listSCSN.PDS_sub.gt[[ii]]@meta.data$lasso.scKFO , 
                                XX)
        print( dim(XX) )
        
        return( XX)
      }))
      sc_AUCell.indGene.list <- lapply( seq(listSCSN.PDS_sub.gt), function(ii) {
        print(ii)
        XX <- cbind( t( sc_AUCell.indGene.b1[[ii]]) ,
                     t( sc_AUCell.indGene.b2[[ii]]) , 
                     t (sc_AUCell.indGene.b3[[ii]]) )
        XX <- XX[ , match( names(gSets.DS_all), colnames(XX))]
        colnames(XX) <- names(gSets.DS_all)
        
        # combine sample and genotype
        XX <- cbind.data.frame( study= names(SCSN_PDSlist_limit)[ii],
                                sample=listSCSN.PDS_sub.gt[[ii]]@meta.data$sample ,
                                gtypeDE = listSCSN.PDS_sub.gt[[ii]]@meta.data$gtypeDE,
                                group = listSCSN.PDS_sub.gt[[ii]]@meta.data$group,
                                PDS.42 = listSCSN.PDS_sub.gt[[ii]]@meta.data$PDS.42 ,
                                lasso.scKFO = listSCSN.PDS_sub.gt[[ii]]@meta.data$lasso.scKFO , 
                                XX)
        print( dim(XX) )
        
        return( XX)
      })
      saveRDS(sc_AUCell.indGene.list, file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/sc_AUCell.indGene.list.rds")
      
    }
    
  ### make lm to find predictors for single cell
    {
      # select only samples that were used for training models that predict AlbCre
      snKFO_AUCell <- sc_AUCell.indGene[ sc_AUCell.indGene$sample %in% annot_tab$CCG_Sample_ID ,
                                 !colnames(sc_AUCell.indGene) %in% c("sample","group","study",
                                                                     "PDS.42","lasso.scKFO")] 
      scGSE_AUCell <- sc_AUCell.indGene[ !(sc_AUCell.indGene$sample %in% annot_tab$CCG_Sample_ID) ,
                                         !colnames(sc_AUCell.indGene) %in% c("sample","group","study",
                                                                             "PDS.42","lasso.scKFO" )]
      snsc_AUCell <- rbind( snKFO_AUCell , scGSE_AUCell )
      
      
      require(glmnet)
      llscsnCL <- list( snKFO_AUCell, scGSE_AUCell  )
      lasso_scsnCL <- lapply( seq(llscsnCL), function(ii){
        print(ii)
        datt <- llscsnCL[[ii]]
        datt$gtypeDE <- ifelse(datt$gtypeDE=="control", 0,1 )
        datt[is.na(datt)] <- 0
        datt <- as.matrix(datt)
        cv.lasso_res <- cv.glmnet( y=datt[,1], x=datt[,-1], 
                                   type.measure="class", alpha=1 , 
                                   family="binomial", nfolds= 30 )
        outcome <- coef(cv.lasso_res, s=cv.lasso_res$lambda.1se)
        cv.lasso <-  outcome[outcome[,1]!=0,][-1]
        
        sspath<- c060::stabpath(y=datt[,1], x=datt[,-1] ,size=0.6,steps= 50,
                 weakness=1, mc.cores=4, family="binomial", type.measure="class")
        sstab<- c060::stabsel(sspath , error=0.05, type="pfer",pi_thr=0.6)$stable
        
        return(list(cv.lasso, sstab))
        
      } )
      lasso_scsnCL <- Reduce( c , lasso_scsnCL)
      names( lasso_scsnCL )  <- c(  "cvlasso.CL.snKFO","sslasso.CL.snKFO", 
                                    "cvlasso.CL.scGSE", "sslasso.CL.scGSE" )
      
   
      cv.lasso_listALLcv.lasso_listALL# saveRDS(cv.lasso_listALL, file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/ML/cv.lasso_listALL.rda")
      
      
    }
    
  ### sc aggregated
    {

     # combine AUCell calculation
     datMean.scCL <- sc_AUCell.indGene 
     datMean.scCL[is.na(datMean.scCL)] <- 0
     datMean.scCL <- aggregate( .~sample , data=datMean.scCL[,-c(2,3)] ,   mean)
     gg <- unique(sc_AUCell.indGene[, c("gtypeDE","sample")])
     datMean.scCL <- cbind( gtypeDE = gg$gtypeDE[ match( datMean.scCL$sample, gg$sample)] , 
                            datMean.scCL )
     datMean.scCL$gtypeDE <- ifelse(datMean.scCL$gtypeDE=="control", "control", "experiment")
    
  }

  ### bulk
    # reaggregate based on sample, not group!
    {
        # aggregate data
        datMean.bulkCL <- Reduce( rbind , lapply( seq(listBULK) , 
                                                function( ii ){
                                                  print(ii)
                                                  id <- names(listBULK)[ii]
                                                  score <- bulk_AUCell.indGene[[id]]
                                                  score <- score[ match(names(gSets.DS_all), rownames(score)),]
                                                  rownames(score) <- names(gSets.DS_all)
                                                  score[is.na(score)] <- 0
                                                  # if a multistaged study - use the prepared annotation list,
                                                  # otherwise extract annotation from sample names
                                                  if( id %in% names(annot_bulkStage)) {
                                                    annot <- annot_bulkStage[[id]]$groups[ 
                                                      match( sub( "__.*", "", colnames(score)), 
                                                             rownames( annot_bulkStage[[id]] ))]
                                                  } else {
                                                    annot <- sub( ".*__", "", colnames(score))
                                                  }
                                                  
                                                  XX <- data.frame( groups= as.factor(annot)  , t(score))
                                                  # XX <- aggregate( .~groups, data=XX , FUN=mean)
                                                  XX <-  data.frame( study=id , XX)
                                                  
                                                  return(XX)
                                                }))
        
      # categorise all samples in  experimental and control samples
        datMean.bulkCL <- cbind(  gtypeDE = ifelse( datMean.bulkCL$groups %in% 
                                                  c("experiment","mutant", "ko_4w","D9 after ADR injection",
                                                    "LMB2day4", "ko_12w", "D14 after ADR injection","LMB2day7",
                                                    "ob 4 weeks","ob 8 weeks", "ob 16 weeks","ob 24 weeks",
                                                    "mut_4w", "mut_12w", "Ercc1[-/Δ], 4 weeks",
                                                    "Ercc1[-/Δ], 14 weeks", "wild type, 96 weeks"),
                                                "experiment", "control") , datMean.bulkCL)
      }
   
  
  ### make lm to find predictors for pseudo(bulk)
    {
      # select only samples that were used for training models that predict AlbCre
      datt.scCL <- datMean.scCL[, !colnames(datMean.scCL) %in% 
                                      c("group","gtype","sample","study","PDS.42",
                                        "platform", "groups","sgt")] 
      datt.scCL.sn <- datMean.scCL[ datMean.scCL$sample %in% annot_tab$CCG_Sample_ID ,
                                 !colnames(datMean.scCL) %in% 
                             c("group","gtype","sample","study","PDS.42",
                               "platform", "groups","sgt")] 
      datt.scCL.sc <- datMean.scCL[ !(datMean.scCL$sample %in% annot_tab$CCG_Sample_ID ),
                                 !colnames(datMean.scCL) %in% 
                                   c("group","gtype","sample","study","PDS.42",
                                     "platform", "groups","sgt")] 
      datt.bulkCL <- datMean.bulkCL[ (datMean.bulkCL$study) %in% stud , 
                                    !colnames(datMean.bulkCL) %in% 
                                 c("group","gtype","sample","study","PDS.42.005",
                                   "platform", "groups","sgt")] 
      datt.bulkCLall <- datMean.bulkCL[, 
                                     !colnames(datMean.bulkCL) %in% 
                                       c("group","gtype","sample","study","PDS.42.005",
                                         "platform", "groups","sgt")] 
      dattCL <- rbind(datt.scCL, datt.scCL.sn )
      dattCLall <- rbind(datt.bulkCLall, datt.scCL )
      
      
      require(glmnet)
      llCL <- list( datt.bulkCL , datt.scCL.sn ,datt.scCL.sc, datt.scCL , dattCL, dattCLall )
      cv.lasso_listCL <- lapply( seq(llCL) , function(ii){
        datt <- llCL[[ii]]
        datt$gtypeDE <- as.numeric( as.factor( datt$gtypeDE) )
        datt <- as.matrix(datt)
        cv.lasso_res <- cv.glmnet( y=datt[,1], x=datt[,-1], 
                                   type.measure="mse", alpha=1 , 
                                   family="gaussian", nfolds= 5 )
        outcome <- coef(cv.lasso_res, s=cv.lasso_res$lambda.min)
        cv.lasso_lm.lmin <-  outcome[outcome[,1]!=0,][-1]
        return(cv.lasso_lm.lmin)
        
      })
      names(cv.lasso_listCL)  <- c(   "datt.bulkCL" , "datt.scCL.sn" ,
                                      "datt.scCL.sc", "datt.scCL" , "dattCL", "dattCLall")
      
    
      # kinda stability
      stabcv.lassoCL_list <- lapply( seq(llCL), function(ii){
        print(ii)
        datt <- llCL[[ii]]
        datt$gtypeDE <- ifelse(datt$gtypeDE=="control", 0,1 )
        datt <- as.matrix(datt)
        freqq <- lapply( 1:100, function(jj){
          cv.lasso_res <- cv.glmnet( y=datt[,1], x=datt[,-1], 
                                     type.measure="class", alpha=1 , 
                                     family="binomial", nfolds=round( sqrt(nrow(datt))))
          outcome <- coef(cv.lasso_res, s=cv.lasso_res$lambda.min)
          cv.lasso_lm.lmin <-  outcome[outcome[,1]!=0,][-1]
        })
        freqq <- table(unlist(lapply(freqq, names)))
        freqq <- freqq[order(freqq)]
        return(freqq)
        
      })
      saveRDS( stabcv.lassoCL_list , "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/ML/tabcv.lassoCL_list.rda" )
      # find coefficients
      Y <-  stabcv.lassoCL4f_list[[6]][stabcv.lassoCL4f_list[[6]]==100]
      Y<- stabcv.lasso_list[[3]][stabcv.lasso_list[[3]]>30]
      Y <-  Reduce( c, lapply( seq( Y), function(ii, datt=datT ){
        # datt$gtypeDE <- ifelse(datt$gtypeDE=="control", 0,1 )
        datt <- as.data.frame(datt)
        datt$AlbCrRatio <- log( datt$AlbCrRatio )
        datt <- datt[,colnames(datt)%in% c("AlbCrRatio", names(Y)[ii])]
        coefficients( lm( AlbCrRatio ~ ., data =datt ) )[2]
      }))
      
      # combine models predicting Alb/Cre and class labels and save
      cv.lasso_listALL <- c( cv.lasso_list , cv.lasso_listCL)
      cv.lasso_listALL
      names(cv.lasso_listALL) <- c("lasso.bulk" ,  "lasso.scKFO", "lasso.bulkANDsc" , 
                                   "lassoCL.bulk", "lassoCL.scKFO" , "lassoCL.bulkANDsc","PDS.42")
      # saveRDS(cv.lasso_listALL, file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/ML/cv.lasso_listALL.rda")
      
      
    }

    

}

### plot score VS alb/cre
  { 
  #   # combine models predicting FSGS lvls and class labels
   cv.lasso_listALL <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/ML/cv.lasso_listALL.rda")
    modelsPlot <- list( X[stabcv.lasso_list[[2]]>30], Y[stabcv.lassoCL_list[[2]]>30], cv.lasso_listALL[[7]])
    modelsPlot <- list(cv.lasso_list[[2]], cv.lasso_listCL[[2]],  cv.lasso_listALL[[7]])
    modelsPlot <- c(cv.lasso_listALL[c(2,3)] ,list(Y ), list(X ) )
    names(modelsPlot) <- c("lasso.sc", "lasso.both","lasso.all","stabselCLall")

### scatterplot
    {
    # prepare data to plot
    toPlot.sc <- aggPDS.sc
    toPlot.sc$gtype <- ifelse( toPlot.sc$gtype=="wt", "control", 
                               as.character( toPlot.sc$gtype) )
    toPlot.bulk <- aggPDS.bulk
    toPlot <- rbind( toPlot.bulk[, -c(2:5)], toPlot.sc[, -c(2:5)])
    toPlot$platform <- c( aggPDS.bulk$platform , rep("sc", nrow(aggPDS.sc)) )
    toPlot$gtype <- c( aggPDS.bulk$gtype , ifelse( 
      aggPDS.sc$gtype=="wt", "control","experiment" ) )
    
    # best fit 
    # [1] "Snx31"   "Itgb8"   "Myo1d"   "Cxcl1"   "Celf2"   "Soat1"   "Anxa3"   "Tnfaip2"
    # [9] "Cd59a"   "Rhoq"    "Snx13"   "Nsmce2"  "Kdm6a"   "Tsr2"    "Dnajc3" 
    
  
    ### make plots
    datt.list <- list( toPlot.bulk, toPlot.sc )
    names(datt.list) <- c( "bulk.data", "scKFO.data" )
    
      gglist <-  Reduce( c , lapply( seq(datt.list), 
                                     function(jj, modelList=modelsPlot){ 
        toPlot <- datt.list[[jj]]
        
        glist <- lapply( seq(modelList), function( ii ){
          
         tt1 <- rowSums( as.data.frame(
            toPlot[,names(modelList[[ii]])[modelList[[ii]]>0]]) ) 
           tt2 <-   rowSums( as.data.frame(
              toPlot[,names(modelList[[ii]])[modelList[[ii]]<0]]) )
           
         toPlot$gset <- tt1-tt2 
          gg <-  ggplot2::ggplot( data = toPlot , aes( x=gset , y=log(AlbCrRatio))) +
            geom_point(  size=6, aes(  color=gtype )) +
            theme_bw() +  theme( text = element_text(size = 22)) + 
            geom_smooth( method='lm', se = FALSE) + 
            xlab(names(modelList)[ii]) +
            ggtitle( names(datt.list)[jj])+ 
            stat_cor( size=6, method = "pearson" )  +
            stat_cor( size=8, method = "spearman",  cor.coef.name ="rho", label.y.npc = 0.8  )
          
            
          return(gg)
        })
        return(glist)
        
      }))
      names(gglist)<- apply( expand.grid( names(modelsPlot) , names(datt.list) ), 
               1, paste, collapse="_")
      
      cowplot::plot_grid( plotlist = gglist , ncol = 2 )
   
  }
  
### boxplots
    {
    ### get the data
    toPlot.bulkCL <- datMean.bulkCL
    # train
    toPlot.bulkCL.train <- toPlot.bulkCL[ (toPlot.bulkCL$study %in%  aggPDS.bulk$study), ]
    # test
    toPlot.bulkCL.test <- toPlot.bulkCL[ !(toPlot.bulkCL$study %in%  aggPDS.bulk$study), ]
    #======================
    toPlot.scCL <- datMean.scCL
    # train
     toPlot.scCL.train <- toPlot.scCL[ toPlot.scCL$sample   %in%  aggPDS.sc$sample, ]
     # test
     toPlot.scCL.test <- toPlot.scCL[ !(toPlot.scCL$sample   %in%  aggPDS.sc$sample ), ]
     #==========
     toPlotCL.train <- rbind(toPlot.bulkCL.train[,-c(2:3)] , toPlot.scCL.train[,-2] )
     toPlotCL.train$type <- c( rep("bulk/MA", nrow(toPlot.bulkCL.train)), rep("sc", nrow(toPlot.scCL.train)))
     toPlotCL.test <- rbind(toPlot.bulkCL.test[,-c(2:3)] , toPlot.scCL.test[,-2] )
     toPlotCL.test$type <- c( rep("bulk/MA", nrow(toPlot.bulkCL.test)), rep("sc", nrow(toPlot.scCL.test)))
     
    # dat list
     datt.listCL <- list( toPlot.bulkCL.train , toPlot.bulkCL.test ,
                          toPlot.scCL.train, toPlot.scCL.test   )
     names(datt.listCL) <- c( "bulk.train" , "bulk.test" , "sc.train" , "sc.test"  )
     
     # make plots
     gglist.CL <-  Reduce( c , lapply( seq(datt.listCL), 
                                       function(jj,modelList=modelsPlot){ 
       toPlot <- datt.listCL[[jj]]

             glist <- lapply( seq(modelList), function(ii ){
               
               tt1 <- rowSums( as.data.frame(
                 toPlot[,names(modelList[[ii]])[modelList[[ii]]>0]]) ) 
               tt2 <-   rowSums( as.data.frame(
                 toPlot[,names(modelList[[ii]])[modelList[[ii]]<0]]) )
               
               toPlot$gset <- tt1-tt2 
 
         
         gg <-  ggplot2::ggplot( data = toPlot, aes( 
           y=gset , x=gtypeDE)) +
           geom_boxplot(  ) + geom_jitter ( width = 0.25, alpha= 0.2 ) + 
           theme_bw() +  theme( text = element_text(size = 16)) + 
           ggtitle(names(datt.listCL)[jj]) + ylab( names(modelList)[ii]) +
           stat_compare_means(  method = "t.test",  size=5, label.y.npc = 0.95)
        
          return(gg)
       })
      return(glist)
     
      }))
    names(gglist.CL) <- apply(expand.grid( names(modelsPlot), names(datt.listCL)  ), 
          1, paste, collapse="_")
    # print plots
    cowplot::plot_grid( plotlist = gglist.CL , ncol = 8 )
    
    
    # DS_all$stage.pred.SC <- ifelse(DS_all$gene_symbol %in% names( cv.lasso_listALL[[2]]) , "YES", "NO")
    # DS_all$label.pred.SC <- ifelse(DS_all$gene_symbol %in%  names(cv.lasso_listALL[[5]]), "YES", "NO")
    # # DS_all$stage.pred.SC42 <- ifelse(DS_all$gene_symbol %in% names( cv.lasso_listALL[[2]]) , "YES", "NO")
    # DS_all$label.pred.SC42 <- ifelse(DS_all$gene_symbol %in%  names(cv.lasso_listALL[[5]]), "YES", "NO")
    # 
    write.table(DS_all , sep = "\t", quote = F, col.names = NA,
                file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/ML/DS_all.12.04.2023.tsv")
 
    
    ## combine plots for scKFO data
    g1 <-  cowplot::plot_grid( plotlist = gglist[
      c( "lasso.scKFO_scKFO.data","lassoCL.scKFO_scKFO.data","PDS.42_scKFO.data" ) ],
      ncol = 3)
    g2 <- cowplot::plot_grid( plotlist = gglist.CL[ 
      c( "lasso.scKFO_sc.train", "lasso.scKFO_sc.test",
         "lassoCL.scKFO_sc.train","lassoCL.scKFO_sc.test",
         "PDS.42_sc.train",  "PDS.42_sc.test" ) ] ,
      ncol = 6 )
    cowplot::plot_grid( plotlist= list( g1,g2 ) , ncol = 1 )
     }
 
### density
    {
      ggplot2::ggplot( toPlot4, aes(   x=value , color=condition) ) + 
        scale_color_colorblind() +  
        geom_density(size=1.5 ) + theme_bw() + 
        ylab("density")+xlab("damage score")+
        ggtitle("\nall cells")+  
        theme( text = element_text(size = 24), legend.position = "none") +
        geom_vline( data=grp.mean, 
                    aes( xintercept=value, color=condition), linetype="dashed") +
        geom_text_repel( data = grp.mean , size=10 ,
                         aes(x=-0.2, y=0, label = t.test.padj  )) + 
        coord_cartesian(xlim = c(-0.4,0.15))+
        facet_grid(rows=vars(variable),scales = "free_y"  ) 
    }
    
  
 
}

#### Robustness test ####
  {
 
  # percent of studies randomized
  pert.percent <- c( 0 , 10, 20, 30, 40 , 50 , 60 , 70 , 80 , 90, 100)

  
    # run test for subsampled  Wt1hetdel KFO dataset
    # suset data 50 times to avoid dataset bias
    PDS_pert.lvl <- Reduce( cbind.data.frame , lapply( 1:50, function(kk)
      {
      
      # print(paste("randomisation round", kk ))
      
      datt <- DEdat
      # choose dataset to randomize
      whichR <- sample( 1:length(datt), length(datt)*( pert.percent[ii]/100) )
      # randomise pvalue and lfc of chosen datasets
      datt[whichR]  <- lapply( datt[whichR]  , function(XX){
        XX$log2FoldChange <- sample( XX$log2FoldChange)
        XX$pvalue <- sample( XX$pvalue)
        return(XX) })
      
      PDS <- damage_signature.func( datt , thresh = 0.75, 
                                    geneFilter = SCSNpodogenes )
      
      # calculate disease score for one dataset, bulk KFO Wt1
      XXX <- DS_calc.func( exprMatrices = exprDat , ceilThrsh = ceilThrsh , 
                    genesUP = PDS$gene_symbol[1:ntop][
                      which( PDS$direction_foldchange[1:ntop]==1)],
                    genesDOWN = PDS$gene_symbol[1:ntop][
                      which( PDS$direction_foldchange[1:ntop]==-1)] )
      return(XXX)
    } ) )
    
    # prepare for ploting with ggplot
    PDS_pert.lvl$gtype <- sub(".*__" , "", rownames(PDS_pert.lvl))
    XX <- reshape::melt( data=PDS_pert.lvl )
    XX$random.percent <- as.factor(pert.percent[ii])
    return(XX)
   
    
  saveRDS( perturb_test, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Pertrubation/RandomTest_bulkWt1.KFO.rda" )
  
  # make boxplot of results
  toPlot <- Reduce( rbind, lapply( perturb_test , function(X) {
    X$value <- scale(X$value)
    return(X)
  }))
  ggp <- ggplot( toPlot , aes( x=random.percent , y=value , color=gtype)) +
    geom_boxplot(outlier.shape = NA) + theme_bw() +  
    coord_cartesian(ylim =  c(-2.5, 2.5)) +
    scale_color_colorblind() + theme( text = element_text(size=18)) +
    labs(x = "percent of studies randomised" , y = "scaled PDS", 
         title = "randomisation effect on PDS\n tested on bulk Wt1.het.del")
  
  # save boxplots  
  pdf(height = 6, width = 6, file =   "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Pertrubation/RandomTest_bulkWt1.KFO_plot.pdf")
    print(ggp)
  dev.off()
  
  png(height = 400, width = 400, file =   "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Pertrubation/RandomTest_bulkWt1.KFO__plot.png" )
    print(ggp)
  dev.off()
}



#### cross-validate damage types #### 
annotations <- readRDS(file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/PDS/explist_44.ma.bulk.sc_annot.rda")
DS_all<- read.table(header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")
explist <- readRDS(  file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/PDS/explist_44.ma.bulk.scPbulk.rda")
XX <- readRDS(  file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/PDS/explist_44.ma.bulk.sc.rda")
names(explist) <- names(XX)

### clreate a list of models    
  {
  CV_model.list <- list(
    
    # diabetic datasets
    diab_models = c( "GSE131266", "GSE134327", "GSE123853", "GSE77717", "GSE79291", 
                     "GSE36209" ,"GSE106841", "GSE20844", "GSE168676",
                     "btbr_sc","GSE146912.BTBR_podo"),
    # slit diaphragm datasets
    sltdph_models = c( "GSE63272", "GSE123179" , "GSE110092" , "KFO.Nphs2",
                       "Nphs2_sc", "cd2ap_sc", "KFO.Nphs2_podo",
                       "GSE146912.Cd2ap_podo"),
    # toxic damage datasets
    txcdmg_models = c("GSE108629","KFO.adria", "GSE154955", 
                      "GSE174013_sc","doxo_sc",
                      "GSE174013_podo" , "GSE146912.doxo_podo") ,
    # TF models
    tf_models = c( "GSE96044","GSE18358","GSE117571","GSE17709","KFO.Wt1",
                   "GSE174102_sc","Wt1_sc","KFO.Wt1_podo","GSE174102_podo"))
  
  CV_model.list_DS <- lapply(seq(CV_model.list), function(ii){
    models <- CV_model.list[[ii]]
    print( names(DE_list)[(names(DE_list)%in% models)])
    damage_signature.func( DE_list[ !names(DE_list)%in% models], 
                           thresh = 0.75, geneFilter = allPodoGenes )
    
  })
  names(CV_model.list_DS) <- names(CV_model.list)
  saveRDS(CV_model.list_DS, file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/models/CV_model.list_DS.rda")
  
}
      
### calculate PDS for left out datasets
  {
  CV_model.list_PDS42 <- lapply( seq(CV_model.list), 
                                 function(ii){
    print(names(CV_model.list)[ii])
    # names of the left out expression tabs
    dataSets <- CV_model.list[[ii]][ CV_model.list[[ii]]%in% 
                                       names(explist)] 
    
    # calculate PDS for left out datasets
    ll<- lapply(seq(dataSets) , function(jj){
      print(jj)
      expr <- explist[[dataSets[jj]]]
      # remove genes expressed in less than 1% of cells from sc.sn
      if( grepl("podo",dataSets[jj])){
        expr <- expr[ rowSums( round(expr)>0)> 0.01*ncol(expr),]
      }
      
      DS_calc.func( exprMatrices = expr ,
                    ceilThrsh = 0.05 ,
                    DSignature= DS_all , 
                    ntop=150 , wghtd=T )
      
    })
    names(ll) <- dataSets
    return(ll)
  })
  names( CV_model.list_PDS42 ) <- names(CV_model.list)
  saveRDS( CV_model.list_PDS42, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/Dev.Valid/CV_model.list_PDS42.rda")
 
  # aggregate by sample 
  CV_model.list_PDS42.agg <- lapply( seq(CV_model.list), 
                 function(ii){
                   print(names(CV_model.list)[ii])
                   dataSets <- CV_model.list_PDS42[[ii]]
                   # calculate PDS for left out datasets
                   ll <- lapply(seq(dataSets) , 
                                              function(jj){
                     print(jj)
                     datt <- dataSets[[jj]]
                     # add annotation
                     if( !grepl("podo",names(dataSets)[jj])){
                       datt <- data.frame( PDS.42= datt,
                                           sample= sub( "__.*", "",names(datt)),
                                           gtypeDE= sub( ".*__", "",names(datt)),
                                           dataset= names(dataSets)[jj],
                                           CV_model = names(CV_model.list)[ii])
                     }  else  {
                       datt <- data.frame( PDS.42= datt,
                                           sample= annotations$annotlist_scSample[[ 
                                             names(dataSets)[jj] ]],
                                           gtypeDE= annotations$annotlist_gtype[[
                                             names(dataSets)[jj] ]] ,
                                           dataset = names(dataSets)[jj] ,
                                           CV_model = names(CV_model.list)[ii]) 
                       # aggregate by samples
                       datt <- aggregate( .~sample+gtypeDE+dataset+CV_model, data=datt, FUN=mean)
                       
                     } 
                     return(datt)
                     })
                   names(ll) <- names(dataSets)
                   return(ll) 
                  
})
  
  # fix GSE168676 annotation!!!
  CV_model.list_PDS42.agg[[1]]$GSE168676$gtypeDE <- ifelse( bulk_annlist$GSE168676$groups=="Control (baseline)", 
                                                            "control","experiment")
  ## fix 117571 annot
  CV_model.list_PDS42.agg[[4]]$GSE117571$gtypeDE <- ifelse( CV_model.list_PDS42.agg[[4]]$GSE117571$gtypeDE =="control",
                                                            "control","experiment")
  # select gloms
  CV_model.list_PDS42.agg[[4]]$GSE117571 <- 
    CV_model.list_PDS42.agg[[4]]$GSE117571[
      CV_model.list_PDS42.agg[[4]]$GSE117571$sample%in%c(
        "GSM3304015", "GSM3304016", "GSM3304017", "GSM3304018"), ]
  names(CV_model.list_PDS42.agg) <- names(CV_model.list_PDS42)
  saveRDS( CV_model.list_PDS42.agg, file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/Dev.Valid/CV_model.list_PDS42.sampAgg.rda")
  
}

### plot CV results
  {
    GSE134327ann <- read.table( header = T, sep = ",","/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/models/GSE134327_SraRunTable.csv")
    GSE168676ann <- read.table( header = T, sep = ",","/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/models/GSE168676_SraRunTable.csv")
    CV_model.list_PDS42.agg <- readRDS( file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/Dev.Valid/CV_model.list_PDS42.sampAgg.rda")
    
    # scale PDS
    toPlot.list <- lapply(seq(CV_model.list_PDS42.agg), function(ii){
      print(ii)
      toPlot <- CV_model.list_PDS42.agg[[ii]]
    Reduce( rbind, lapply( seq(toPlot), function(jj){
      print(jj)
         datt <- toPlot[[jj]]
         if (names(toPlot)[jj]=="GSE134327") {
           datt$sample <- GSE134327ann$Sample.Name[ match(
             datt$sample  , GSE134327ann$Run 
           )]
           datt <- aggregate( PDS.42 ~ . , data=datt, FUN=mean)
           datt$PDS.42 <- scale(datt$PDS.42)
           datt <- datt[ ,c( "PDS.42","sample","gtypeDE","dataset",  "CV_model")]
         } else if  (names(toPlot)[jj]=="GSE168676") {
            datt$sample <- GSE168676ann$Sample.Name[ match(
              datt$sample , GSE168676ann$Run 
            )]
            datt <- aggregate( PDS.42 ~ . , data=datt, FUN=mean)
            datt$PDS.42 <- scale(datt$PDS.42)
            datt <- datt[ ,c( "PDS.42","sample","gtypeDE","dataset",  "CV_model")]
       } else {
         datt$PDS.42 <- scale(datt$PDS.42)
       }
     
        return(datt)
      }))
    }) 
names(toPlot.list) <- names(CV_model.list_PDS42.agg)
    
    ### plot individual models
datt <-  Reduce( rbind, lapply(seq(toPlot.list), function(ii){
  toPlot.list[[ii]]
  
  }))
     gg0 <-  ggplot2::ggplot( datt , 
                              aes( y= PDS.42, x=gtypeDE, color=gtypeDE)) +
        scale_color_colorblind() +  geom_boxplot(lwd=1.5,outlier.shape = NA) + 
        geom_jitter(size=5,width = 0.1, alpha=0.5)+
        theme_bw() + theme( text = element_text(size = 24) ,
                            axis.title.x=element_blank(),
                            axis.text.x=element_blank(),
                            axis.ticks.x=element_blank()) +
        # stat_compare_means()+
        stat_compare_means( size = 6   )+
        facet_wrap(facets = vars(dataset),scales = "free")
    
     gg0
        pdf( height = 20, width = 20, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/models/PDS_cv.models_boxplot_Indv.test.pdf" )
        print( gg0 )
        dev.off()
        
        png( height = 1800, width = 1500, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/models/PDS_cv.models_boxplot_Indv.png" )
        print( ppglist )
        dev.off()
        
        ### make one plot
        toPlot <- Reduce( rbind, lapply(seq(toPlot.list), function(ii){
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
          stat_compare_means( size = 6 ) +
          facet_grid(cols = vars(CV_model),scales = "free")
        
        pdf( height = 6, width = 12, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/models/PDS_cv.models_boxplot_agg.pdf" )
        print( gg )
        dev.off()
        
        png( height = 600, width = 1000, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/models/PDS_cv.models_boxplot_agg.png" )
        print( gg )
        dev.off()
  }
  

#### calculate PDS for all datasets ####
{
  # use 42 genes and all_DS signature
  PDS.42 <-  lapply( seq(explist) , function(ii )
  {
    print(ii)
    DS_calc.func( exprMatrices = explist[[ii]], ceilThrsh = 0.05 , 
                  DSignature= DS_all , ntop=42 )
  })
  # use all  genes and all_DS signature
  PDS.all <-  lapply( seq(explist) , function(ii )
  {
    print(ii)
    DS_calc.func( exprMatrices = explist[[ii]], ceilThrsh = 0.05 , 
                  DSignature= DS_all , ntop=nrow(DS_all) )
  })
  # name and save
  names(PDS.all) <-   names(PDS.42) <-     names(explist)
  saveRDS( PDS.42, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/AUcell.42_PDSall.10.06.2022.rda")
  saveRDS( PDS.all , "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/AUcell.all_PDSall.10.06.2022.rda")
  
  # # pseudotime
  # PDS.ptime <-  lapply( seq(explist) , function(ii )
  # {
  #   print(ii)
  #   DS_calc.func( exprMatrices = explist[[ii]], ceilThrsh = 0.05 , 
  #                 DSignature= DS_cv.list[["ptime_PDS"]] , ntop=42 )
  # })
  

  
  PDS_plot <-lapply( seq(PDS.all), function(ii )
  {
    GSE <-  names(PDS.all)[ii]  # get GSE id of a current experiment
    toPlot <- data.frame( score= PDS.all[[ii]] , gt =annotlist[[ii]])
    
    if ( length(PDS.all[[ii]])<50) {
      gg <- ggplot2::ggplot( toPlot , aes( y=score, x=gt, color=gt)) +
        scale_color_colorblind() +  geom_jitter(width = 0.1, size=3) + 
        theme_bw() + theme( text = element_text(size = 18) , legend.position = "none") +
        ggtitle(GSE) +
        stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean,
                     geom = "crossbar", width = .2, color = "red")  
    } else {
      # make a density plot
      gg <-  ggplot2::ggplot( toPlot , aes( x=score, color=gt)) + 
        scale_color_colorblind() +  geom_density(size=1.5) + 
        theme_bw() +  ggtitle( GSE ) +
        theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
               text = element_text(size = 18)) +
        geom_vline(data=ddply( toPlot, "gt",
                               summarise, grp.mean=mean(score)), 
                   aes(xintercept=grp.mean, color=gt), linetype="dashed")
    }
    return(gg)
    
  })
  
  ppg1 <- cowplot::plot_grid( plotlist=PDS_plot[1:14], nrow = 4)
  ppg2 <- cowplot::plot_grid( plotlist=PDS_plot[15:32],nrow = 5)
  ppg3  <- cowplot::plot_grid( plotlist=PDS_plot[33:44], nrow = 4)
  
  pdf(height = 12, width = 12, file =  "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/platforms/PDS.allPlatforms_all449_plot.pdf")
  print( ppg1 )
  print( ppg2 )
  print( ppg3 )
  
  dev.off()
  
  png(height = 1000, width = 1000, file = paste( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/platforms/",
                                                 names(PDS.platforms)[jj],"_plot.png", sep = ""))
  print( ppg1 )
  print( ppg2 )
  print( ppg3 )
  dev.off()
  
}
#### visualise datasets with different stages ####
  {
    
  ### prepare annotations
    {
      ### GSE106841__obob_MA
      library( "annotationTools")
      GSE106841 <- GEOquery::getGEO( "GSE106841" , GSEMatrix=T )
      GSE106841_des <- GSE106841$GSE106841_series_matrix.txt.gz@phenoData@data["source_name_ch1"]
      GSE106841_des$groups <- sub( "Glomerulus, ", "", GSE106841_des$source_name_ch1)
      GSE106841_des$groups <-factor(GSE106841_des$groups, levels = c("WT 4 weeks","ob 4 weeks",  "WT 8 weeks","ob 8 weeks",
                                                                     "WT 16 weeks", "ob 16 weeks", "WT 24 weeks", "ob 24 weeks" )) 
      ###  GSE108629__NEP25_MA
      GSE108629 <- GEOquery::getGEO("GSE108629")
      GSE108629_des <-  GSE108629$GSE108629_series_matrix.txt.gz@phenoData@data[1]
      GSE108629_des$ctype <- sub( "_.*", "",GSE108629_des$title)
      GSE108629_des$groups <- as.factor( sub(".*_(.+)_.*", "\\1",GSE108629_des$title))
      GSE108629_des$groups <- relevel(  GSE108629_des$groups , ref="normal")
      
      ###  GSE43061 aging
      GSE43061 <- GEOquery::getGEO("GSE43061")
      GSE43061_des <- GSE43061$GSE43061_series_matrix.txt.gz@phenoData@data[c("age:ch1", "genotype/variation:ch1")]
      colnames(GSE43061_des) <- c("age","gtype")
      GSE43061_des$groups <- paste( GSE43061_des$gtype,GSE43061_des$age ,sep = ", " )
      GSE43061_des$groups <- factor(GSE43061_des$groups, 
                                    levels = c("wild type, 4 weeks","Ercc1[-/Δ], 4 weeks",
                                               "wild type, 14 weeks","Ercc1[-/Δ], 14 weeks",
                                               "wild type, 96 weeks")) 
      ### GSE154955 
      GSE154955_des <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/PDS/bulkRNAseq/expr_tabs/GSE154955_counts.rda")[[2]]
      GSE154955_des$groups <- factor(GSE154955_des$groups, levels = c("PBS", 
                      "D9 after ADR injection", "D14 after ADR injection" ))
      rownames(GSE154955_des) <- GSE154955_des$Run
      
      ### KFO.Nphs2 
      annot <- read.table( sep = "\t", header = T, "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_NPHS2/NPHS2_Sample_NamesAll.csv")
      rownames(annot) <- annot$CCG.Sample.ID
      annot$gtype_vec <- ifelse( annot$gtype=="wt" , "control","experiment")
      annot$groups <- paste(annot$gtype , annot$age , sep = "_")
      annot$groups <-  factor( annot$groups, 
                              levels = c("wt_4w","mut_4w","wt_12w","mut_12w")) 
      KFO.Nphs2_des <- annot
      
      ### KFO.Wt1
      annot <- read.table( sep = "\t", header = T, "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/bulk_Wt1/annot_tab.csv")
      rownames(annot) <- annot$sample
      # make sure we compare experiment VS control
      annot$groups <- paste(annot$gtype, annot$age , sep = "_")
      annot$gtype_vec <- ifelse( annot$gtype=="wt" , "control","experiment")
      annot$groups <-  factor( annot$groups, 
                               levels = c("wt_4w","ko_4w","wt_12w","ko_12w")) 
      KFO.Wt1_des <- annot
        
        # combine
      annot_bulkStage <- list( GSE106841_des, GSE108629_des , GSE43061_des ,
                               GSE154955_des , KFO.Nphs2_des, KFO.Wt1_des )
      names(annot_bulkStage) <- c("GSE106841","GSE108629","GSE43061", "GSE154955","KFO.Nphs2","KFO.Wt1")
    }
 
  ### plot MA and bulk datasets with stages
    {
      # use PDS calulated with the "all inclusive" signature
      PDS.50_list <- PDS.platforms[[5]]

      plotList <- lapply( seq(annot_bulkStage) , function( ii ){
        id <- names(annot_bulkStage)[ii]
        score <- PDS.50_list[[which(names(explist)==id)]]
        annot <- annot_bulkStage[[id]]$groups[ 
          match( sub( "__.*", "", names(score)), 
          rownames( annot_bulkStage[[id]] ))]
        toPlot <- data.frame(score=score, groups= annot)
        gg <- ggplot2::ggplot( toPlot , aes( y=score, x=groups, color=groups)) +
          scale_color_colorblind() +  geom_jitter(width = 0.1, size=5) + 
          theme_bw() + theme( text = element_text(size = 18) , 
                              axis.text.x = element_text(angle = 45, hjust = +1) ,
                              legend.position = "none") +
          ggtitle( id ) +
          stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean,
                       geom = "crossbar", width = .2, color = "red") 
        return(gg)
      })
      
      ppg <- cowplot::plot_grid( plotlist = plotList , nrow = 2 )
      # save the figure
      pdf(height = 8, width = 10, file ="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Disease_stages/bulkStages_jitterPlot.pdf" )
      print(ppg)
      dev.off()
      
      png(height = 600, width = 800, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Disease_stages/bulkStages_jitterPlot.png")
      print( ppg )
      dev.off( )
    }
  
  }

#### correlate PDS with Albumin/creatinine ####
  {
  
  #### bulk ####
    
      stud <- c( "GSE117571", "GSE108629", "GSE17709" ,
                 "GSE117987", "GSE126217" ,"GSE154955",
                 "GSE112116", "GSE131266", "GSE134327", 
                 "GSE110092" ,"GSE77717" , "KFO.Wt1"   )
      
      annot_bulkStage <- readRDS(file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Disease_stages/annot_bulkStage.rda")
      
      # calculate PDS
      blkData <- explist[ names(explist) %in% stud]
      PDSsize.bulk <-  lapply( seq(blkData) , function(ii , ceilThrsh = 0.05  )
      {
        DS_calc.func( exprMatrices = blkData[[ii]], 
                      ceilThrsh = ceilThrsh , 
                      wghtd = F, useOrder = "mean_rank",
                      DSignature= DS_all, ntop = 42)
      })
      names(PDSsize.bulk) <- names(blkData)
      # order like in the stud vector
      PDSsize.bulk <- PDSsize.bulk[stud]
      # GSE117571 use only gloms
      PDSsize.bulk[["GSE117571"]] <- PDSsize.bulk[["GSE117571"]][1:4]
      
      # combine in one df
      ll<- lapply( seq(PDSsize.bulk), function(jj)
      {
        score <- PDSsize.bulk[[jj]]
        id <- names(PDSsize.bulk)[jj]
        # combine score and annotation, then aggregate the score by annotation
        # if a multistaged study - use the prepared annotation list,
        # otherwise extract annotation from sample names
        if( id %in% names(annot_bulkStage)) {
          annot <- annot_bulkStage[[id]]$groups[ 
            match( sub( "__.*", "", names(score)), 
                   rownames( annot_bulkStage[[id]] ))]
        } else {
          annot <- sub( ".*__", "", names(score))
        }
        
        XX <- data.frame(score =score, study  = id,
                         groups= as.factor(annot) )
        XX <- aggregate( .~groups+study, data=XX , FUN=median)
        
        return(XX)
      })
      names(ll) <- stud 
      datMean  <- Reduce( rbind, ll)
      
       # ad proteinuria
        mgmg.Bulk <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsProtoneuria/bulk/mgmg.Bulk.rda")
        datMean$mgmg <-  unlist(mgmg.Bulk) # group  Wt1h.d. KFO
        
        # add platform type annotation
        datMean$platform <- c( rep( "MA" , 2) , rep( "MA" ,3 ) ,  rep( "MA" , 2) , 
                               rep(  "bulk" , 2) , rep( "bulk", 2) ,  rep( "bulk" , 3),
                               rep(  "MA" , 2) , rep( "MA", 2) , rep( "bulk", 2) ,  
                               rep( "bulk", 2) , rep(  "bulk", 2), rep(  "bulk", 4))
        ## add stage
        datMean$stage <- ifelse( datMean$groups %in% c("experiment","mutant", "ko_4w","D9 after ADR injection","LMB2day4"),
                                 "stage.1", ifelse(datMean$groups %in% c( "ko_12w", "D14 after ADR injection","LMB2day7"), 
                                                   "stage.2", "control"))
          
        # make a plot 
        gg <- ggplot2::ggplot( data = datMean, aes( x=score , y=log(mgmg) )) +
          geom_point( aes( color=study , shape=stage), size=6 )+ theme_bw() +  
          ggtitle(paste("42 genes damage signature",sep = ""))+
          theme( text = element_text(size = 22) )  + 
          # geom_text(hjust=0, vjust=0)+
          geom_smooth(method='lm', se = FALSE) + stat_cor(size=7, method = "spearman") 
        gg
        
        # save the plot
        pdf(height = 6, width = 9, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsProtoneuria/bulk/PDS.42vsAlbCr_bul.pdf")
        print(gg)
        dev.off()
        # save the plot
        png(height = 500, width = 800, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsProtoneuria/bulk/PDS.42vsAlbCr_bul.png")
        print(gg)
        dev.off()
      
       
      
  #### sn KFO ####
    
     #plot correlation between PDS and albCr
      source("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/cell-damage-score/AUCell_script.r")
      ##  load Damage signature
      DS_all <- read.table( header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")
      # load data
        listPodo <- c("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_Nphs2/Seurat/decontX.allcells.scDblFind_Seur_podo.rda",
                      "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/snRNAseq_WT1hetdel/Seurat/decontX.allcells.scDblFind_Seur_podo.rda",
                      "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/Cem_data/Seurat/decontX.allcells_Seur_podo.pdss2.rda" )
                      
        listKFO <- lapply(seq(listPodo), function(ii)
          {
          datt <- readRDS(listPodo[ii])
        })
        # names(listKFO )<- c("Nphs2","Wt1","Pdss2")
       
        ## calculate PDS
        listKFO.podo_seur.PDS <- lapply( seq(listKFO), 
                                         function(ii){
          print(ii)
          newSeu <- listKFO[[ii]]
          # subset 
          Idents(newSeu) <- newSeu$sample
          newSeu <- subset(newSeu , downsample=1000)
          Idents(newSeu) <- newSeu$group
          newSeu <- subset(newSeu , downsample=1000)
          
          # get counts
          expr <- newSeu@assays$RNA@counts
          expr <- expr[ rowSums(expr> 0) > 0.01*ncol(expr) , ]
          
          ## equal number of pos and neg
          newSeu$PDS.42 <-   DS_calc.func( exprMatrices = expr ,
                                                 ceilThrsh = 0.05 ,
                                                 DSignature= DS_all , 
                                                 ntop=42 , wghtd=F )
          newSeu$PDS.42.wghtd <-   DS_calc.func( exprMatrices = expr ,
                                                 ceilThrsh = 0.05 ,
                                                DSignature= DS_all , 
                                                 ntop=42 , wghtd=T )
          return(newSeu)
          
        })
        names(listKFO.podo_seur.PDS )<-   c("Nphs2","Wt1","Pdss2")
        saveRDS(listKFO.podo_seur.PDS,file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listKFO.podo_seur.PDS" )
        # saveRDS(listKFO.podo_seur.PDS,file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listKFO.podo_seur.19.09.PDS" )

        # laod annotation
        library(ggpubr)
        annot_tab <- read.table(sep = "\t",header = T, "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/Sample_Names_KFO.csv")
        annot_tab$group <- paste(annot_tab$Genotype,annot_tab$Age_weeks,sep = "_")
        
        # combine metadata from 3 experiments 
        datt <- Reduce( rbind , lapply( listKFO.podo_seur.PDS, function(X) {
          XX <- X@meta.data 
          return(XX[,c( "group","sample","gtype", 
                        grep("PDS",names(XX), value = T)
                        )])
        }))
        
        # treat carefully 21 week Pdss2 samples since they have only per group measurements
        datt1 <- datt[ datt$sample %in% c( "146985", "146986", "143485" , "143486") ,]
        # aggregate
        aggPDS1 <- aggregate( .~group , FUN = mean , 
                              data= datt1[ , c("group", 
                  grep("PDS",names(datt1), value = T))] )
        aggPDS1 <- aggPDS1[rep(seq_len(nrow(aggPDS1)), each = 2), ]
        aggPDS1$group <-  c( "146985", "146986", "143485" , "143486")
        colnames(aggPDS1)[1] <- "sample"
        # the rest of samples
        datt2 <- datt[ !(datt$sample %in%  c("146985", "146986", "143485" , "143486")), ]
        aggPDS2 <- aggregate( .~sample , FUN=mean,
                              data= datt2[ , c( "sample",
                              grep("PDS",names(datt2), value = T))] ) 
                              
        aggPDS <- rbind(aggPDS2, aggPDS1)
        
        aggPDS$AlbCrRatio <- annot_tab$AlbCrRatio[ match( sub("SID","" ,aggPDS$sample) , 
                                                          annot_tab$CCG_Sample_ID)]
        aggPDS$group <- annot_tab$group[ match( sub("SID","" ,aggPDS$sample ), 
                                                annot_tab$CCG_Sample_ID)]
        aggPDS$gtype <- as.factor(annot_tab$Genotype[ match( sub("SID","" ,aggPDS$sample ), 
                                                             annot_tab$CCG_Sample_ID)])
        
        PDSvec <- grep("PDS",names(datt1), value = T)
        gglist <- Reduce(c,  lapply( seq( PDSvec ), 
                                     function(ii){
              # plot lm and correlation
              gg1 <- ggplot2::ggplot( data = aggPDS, aes( x=aggPDS[,PDSvec[ii]], y=log(AlbCrRatio))) +
                geom_point(  size=6, aes(col=gtype)) +
                theme_bw() +  theme( text = element_text(size = 22)) + 
                geom_smooth(method='lm', se = FALSE) + ggtitle( PDSvec[ii])+ 
                stat_cor( size=7, method = "spearman" ) 
              # geom_text(aes(label = sample  ), size=6, position = "dodge")
              
              # density plots 
              gg2 <- ggplot2::ggplot( data = datt, 
                                      aes( x=datt[,PDSvec[ii]], 
                                           color=gtype)) +
                geom_density(size=1.5) + ggtitle( PDSvec[ii])+ 
                theme_bw() +  theme( text = element_text(size = 22)) 
              
              # dotplot for samples of Nphs2mut
              aggPDS.nphs2 <- aggPDS[aggPDS$sample%in% listKFO.podo_seur.PDS$Nphs2$sample,]
              gg3 <- ggplot2::ggplot( data = aggPDS.nphs2, 
                                      aes( y=aggPDS.nphs2[,PDSvec[ii]], 
                                           x=group, 
                                           color=group)) +
                geom_jitter(size=1.5) + ggtitle( PDSvec[ii])+ 
                theme_bw() +  theme( text = element_text(size = 22)) +
                geom_label(aes(label=sample))
              
              # density plots for nphs2
              datt.nphs2 <- datt[datt$sample%in% listKFO.podo_seur.PDS$Nphs2$sample,]
              gg4 <- ggplot2::ggplot( data = datt.nphs2, 
                                      aes( x=datt.nphs2[,PDSvec[ii]], 
                                           color=group)) +
                geom_density(size=1.5) + ggtitle( PDSvec[ii])+ 
                theme_bw() +  theme( text = element_text(size = 22)) 
          
          return( list(gg1, gg2, gg3 ,gg4))
        }))
       cowplot::plot_grid(plotlist = gglist , ncol=2)
        
        
        
}
    
#### correlate PDS with SD length ####
### sic! its not possible for the same mice, but can be done for matched age/disease stages
# nphs2 
# SD.nphs2.meta <- readxl::read_xlsx(sheet = 1,"/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsSDlength/SD_length_per_image_Nphs2mut.2metadata.xlsx")

  SD.nphs2 <- readxl::read_xlsx(sheet = 1,"/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsSDlength/SD_length_per_image__Nphs2mut.xlsx")
  SD.nphs2$age <- sub( "_mut.*" ,"" , SD.nphs2$animal)
  SD.nphs2$strain <- gsub( ".*_mut_| .*" ,"" , SD.nphs2$animal)
  SD.nphs2$mouse <- gsub( ".*_mut_|_img.*| " ,"" , SD.nphs2$animal)
  SD.nphs2$average_SD <- SD.nphs2$`sd length`
  SD.nphs2$GT <- "Nphs2mut"
  SD.nphs2$age_simple <- ifelse(SD.nphs2$age=="04w", 4, "8_9")
  
  # pdss2 
  SD.pdss2 <- read.csv(skip = 2,"/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsSDlength/STED_Makro_Werte__Pdss2mut.csv")
  SD.pdss2 <- SD.pdss2[,1:6]
  SD.pdss2 <- SD.pdss2[ !is.na(SD.pdss2$average_SD) , ]
  SD.pdss2 <- SD.pdss2[SD.pdss2$GT %in% c("pdss2 mid","pdss2 hom","wt", "wt mid"),]
  colnames(SD.pdss2)[1] <- "animal"
  SD.both <- rbind(SD.pdss2[,c("animal", "GT","age_simple","average_SD")] ,
                   SD.nphs2[, c("animal","GT","age_simple","average_SD")])
  SD.both$gtype = ifelse(SD.both$GT%in% c("wt", "wt mid") , "wt", 
                         ifelse(SD.both$GT %in% c("pdss2 hom", "pdss2 mid"), "Pdss2","Nphs2" ))
  SD.both$dataset <- ifelse( SD.both$gtype== "Nphs2", "Nphs2", "Pdss2")
  # gg1 <- ggplot(SD.nphs2, aes(y=average, x=age))+ geom_boxplot(outlier.size = 0)+
  #   geom_point( position = position_dodge(0.1), alpha = 0.5, 
  #                # aes( color=strain) ,
  #               size=3)+
  #   stat_compare_means()
  # gg2 <- ggplot(SD.pdss2, aes(y=average, x=GT))+ geom_boxplot(outlier.size = 0)+
  #   geom_point( position = position_dodge(0.1), alpha = 0.5, 
  #               aes( color=age_simple) ,
  #               size=3)+
  #   stat_compare_means()
  gg0 <- ggplot(SD.both, aes(y=average_SD, x=gtype))+ geom_boxplot(outlier.size = 0)+
      geom_point( position = position_dodge(0.1), alpha = 0.5,
                   aes( color=age_simple) ,
                  size=3)+
      # stat_compare_means()+
    theme_bw() + theme( 
                         text = element_text( size = 14))
  
  
  # calculate average PDS per mice 

  PDS_tab <- Reduce( rbind, lapply( c(1,3), function(ii){
    datt<- listSCSN.1K.sampl[[ii]]@meta.data[ , c("PDS", "gtype", "age", "sample", "group")]
    datt$dataset <- names(listSCSN.1K.sampl)[ii]
    return(datt)
    }))
  
  PDS_tab <- PDS_tab[ PDS_tab$sample %in% 
                        names(table(PDS_tab$sample)[
                          table(PDS_tab$sample)>10]), ]
  aggPDS <- aggregate( .~sample , FUN = mean , na.rm=T,
                        data= PDS_tab[, c(1,4)] )
  aggPDS <- cbind(aggPDS , PDS_tab[ match(aggPDS$sample , PDS_tab$sample),
                                    c("gtype", "age", "group","dataset")])
  aggPDS$age_simple <- aggPDS$age
  
  aggPDS$age_simple[aggPDS$age_simple==8] <- "8_9"
  aggPDS$age_simple[aggPDS$age_simple==6] <- "6_7"
  aggPDS$age_simple[aggPDS$age_simple%in%c(12,14)] <- "12_14"

  aggPDS$age_simple <- as.factor(aggPDS$age_simple)

  gg3 <- ggplot(aggPDS, aes(y=PDS, x=gtype))+ geom_boxplot(outlier.size = 0)+
    geom_point( position = position_dodge(0.1), alpha = 0.5, 
                aes( color=age_simple, shape=dataset) ,
                size=3)+
    # stat_compare_means()+
    theme_bw() + theme( 
                          text = element_text( size = 14))
  
  ### combine SD and PDS
  aggPDS.SD <- merge( aggregate( .~age_simple+gtype+dataset , FUN = mean , na.rm=T,
                    data= aggPDS[,c("age_simple","gtype","PDS","dataset")] ) , 
         aggregate( .~age_simple+gtype+dataset , FUN = mean , na.rm=T,
                    data= SD.both[,c("age_simple","gtype","average_SD","dataset")] ),  )
  
  write.table(aggPDS.SD , sep = "\t", row.names = F, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/WRITING/PDS_manuscript/Revision1/KFO.snRNAseq_aggPDS.SD.tsv")
  
  
  gg4 <- ggplot(aggPDS.SD, aes(y=average_SD, x=PDS ))+ 
    geom_point(  alpha = 0.5, 
                aes( color=gtype , shape=dataset) ,
                size=3)+stat_cor(method="spearman")+
    geom_smooth(method = "lm")+ 
    theme_bw() + theme(  text = element_text( size = 14))
  
  gll <- cowplot::plot_grid( plotlist = list( gg0, gg3, gg4), nrow = 1)
  
  pdf( height = 6, width = 16, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSvsSDlength/PDSvsSD_Nphs2Pdss2.pdf")
  print(gll)
  dev.off()
  
# 
  write.table(aggPDS.SD , sep = "\t", row.names = F, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/WRITING/PDS_manuscript/Revision1/KFO.snRNAseq_aggPDS.vs.SD.tsv")
  write.table(SD.both , sep = "\t", row.names = F, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/WRITING/PDS_manuscript/Revision1/KFO.snRNAseq_SD.tsv")
  write.table(aggPDS , sep = "\t", row.names = F, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/WRITING/PDS_manuscript/Revision1/KFO.snRNAseq_aggPDS.tsv")
  



#### explore RBO instead of AUcell ####
#### Result: gives good results but too slow
  {
    ### prepare PDS ranked marker list
    # read PD signatures
    # DS_cv.list <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/platforms/DS_platform.list.tsv")
    DS_all <- read.table( sep = "\t", header = T, 
                          file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/CrossValidation/platforms/DamageSignatures/DS_all.02.06.2022.tsv")
      
    s <- DS_all$mean_rank
    names(s) <- DS_all$gene_symbol
    spos <- s[which(DS_all$direction_foldchange==1)]
    sneg <- s[which(DS_all$direction_foldchange==-1)]
   
    
    # define rbo testing fuction
    rbo_DS.func <- function( tt , spos , sneg , depth=NULL )
      {
      
      require(gespeR)
      
      # tt is a ranked gene expression list
      rbo_pos <- apply( tt , 2,  function ( list2 ){
        names(list2) <- rownames(tt)
        rbo(spos, list2, p=1, side = "top",  k =depth) 
      })
      
      rbo_neg <- apply(tt , 2, function ( list2) {
        names(list2) <- rownames(tt)
        rbo(sneg, list2, p=1, side = "top", k =depth )
      })
      
      rbo_both <- rbo_pos-rbo_neg
      return(rbo_both)
    } 
    
    #### for all datasets$
    # load bulk and pseudo-bulk
    explist_scPbulk <- readRDS( file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/PDS/explist_44.ma.bulk.scPbulk.rda")
    # calculate RBO
    rbo_PDS <- lapply( seq(explist_scPbulk) , function(ii, ddepth=1000) 
      {
      # t0<- Sys.time()
      print(ii)
        expr <- explist_scPbulk[[ii]]
        
        X <- rbo_DS.func( tt=expr , spos=spos, sneg=sneg , depth=ddepth ) 
     # t1 <- Sys.time()-t0   
      # print(t1)  
        return(X)
      })
    names(rbo_PDS) <- names(explist)
    
    saveRDS( rbo_PDS , file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/rbo_PDSall.20.04.2022.rda")
    
    # compute rbo for podocin
    rbo.nphs2 <- readRDS( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Disease_stages/rbo.PDS_Nphs2.snRNAseq_podo.rda" )
    toPlot <- data.frame( score= rbo.nphs2 , gt = annotlist[[34]])
    ggplot2::ggplot( toPlot , aes( x=score, color=gt)) + 
      scale_color_colorblind() +  geom_density(size=1.5) + 
      theme_bw() +  ggtitle( GSE ) +
      theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
             text = element_text(size = 18)) +
      geom_vline(data=ddply( toPlot, "gt",
                             summarise, grp.mean=mean(score)), 
                 aes(xintercept=grp.mean, color=gt), linetype="dashed")
    
    ## plots for pbulk
    rbo_scPlot <- lapply( 33:44 , function(ii, bulk=F) {  
      
      print(ii)
      toPlot <- data.frame( rboScore= rbo_PDS[[ii]] , 
                            gtype =annotlist[[ii]] )
      
      if(isTRUE(bulk)){
   
      gg <- ggplot( toPlot , aes( x=gtype , y=rboScore , color=gtype)) +
        geom_jitter(width = 0.1, size=3) + theme_bw() +
        scale_color_colorblind() +
        theme( text = element_text(size=18)) +
        theme( text = element_text(size = 18) , legend.position = "none") +
        labs(x = "" , y = "rbo Score",
             title = names(rbo_PDS)[ii] ) +
        stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean,
                     geom = "crossbar", width = .2, color = "red")
    } else { 

      gg <- ggplot2::ggplot( toPlot , aes( x=rboScore, color=gtype)) + 
      scale_color_colorblind() +  geom_density(size=1.5) + 
      theme_bw() +  ggtitle(names(rbo_PDS)[ii]) +
      theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
             text = element_text(size = 18)) +
      geom_vline(data=ddply( toPlot, "gtype",
                             summarise, grp.mean=mean(rboScore)), 
                 aes(xintercept=grp.mean, color=gtype), linetype="dashed")
    }
    return(gg)
  })
  
    ppg <- cowplot::plot_grid( rbo_scPlot[[1]], rbo_scPlot[[2]], rbo_scPlot[[3]], rbo_scPlot[[4]],
                               rbo_scPlot[[5]],rbo_scPlot[[6]],rbo_scPlot[[7]],rbo_scPlot[[8]],
                               rbo_scPlot[[9]],rbo_scPlot[[10]],rbo_scPlot[[11]],rbo_scPlot[[12]],
                               nrow = 4 )
    # save the figure
   pdf(height = 12, width = 15, file ="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/ScoreMethod/rbo_DSall_scPlot.pdf" )
    print(ppg)
    dev.off()
    
    png(height = 1000, width = 1200, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/ScoreMethod/rbo_DSall_scPlot.png")
    print(ppg)
    dev.off( )
    
      # pdf(height = 6, width = 8, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDSvsProtoneuria/PDSvsMgMg_bulk_rboScore.pdf")
      # ggplot2::ggplot( data = datMean, aes( x=rboScore , y=log(AlbCrRatio))) +
      #   geom_point( aes( color=study, shape=as.factor(gt) ), size=4) +
      #   theme_bw() +  geom_smooth(method='lm', se = FALSE) +
      #   stat_cor() 
      # 
      # dev.off()
      
    }

#### compare PDS ptime, trajectory ptime and rbo ####
{
  library( scater )
  library( TSCAN )
  # # merge and process data
  # # this is done with Seurat v5
  # listSCSN.1K_integrated <- IntegrateLayers(object = listSCSN.1K_merged, 
  #                         method = CCAIntegration,
  #                         orig.reduction = "pca", 
  #                         new.reduction = "integrated.cca",
  #                         verbose = FALSE)
  # # process integrated data
  # listSCSN.1K_integrated <- FindNeighbors( 
  #   listSCSN.1K_integrated , reduction = "integrated.cca", dims = 1:10)
  # listSCSN.1K_integrated <- FindClusters(
  #   listSCSN.1K_integrated, resolution = 2, cluster.name = "cca_clusters")
  # listSCSN.1K_integrated <- RunUMAP(
  #   listSCSN.1K_integrated, reduction = "integrated.cca", dims = 1:10, 
  #   reduction.name = "umap.cca")
  # ElbowPlot( listSCSN.1K_integrated , ndims = 50 ) 
  # DimPlot( listSCSN.1K_integrated, group.by = "gtypeDE" ,label = TRUE)
  
  # saveRDS( listSCSN.1K_integrated , file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN.1K_integratedLayers.rda")
  
  
  sce_clust_SCE <- listSCSN.1K_integrated
  sce_clust_SCE[["RNA"]] <- as(sce_clust_SCE[["RNA"]], Class="Assay")
  
  # convert to sce
  sce_clust_SCE <- as.SingleCellExperiment( sce_clust_SCE )
  
  plotReducedDim(  sce_clust_SCE, colour_by=("group"),
                   dimred= "INTEGRATED.CCA",
                   text_by="group", text_colour="red")
  

  library(slingshot)
  # fit a single principal curve 
  sce_clust_SCE$groupSling <- sce_clust_SCE$group
  sce_clust_SCE$groupSling <- ifelse(sce_clust_SCE$gtypeDE !="control", sce_clust_SCE$groupSling , "control")
  sce.sling <- slingshot::slingshot( sce_clust_SCE, reducedDim='UMAP.CCA',
                                     clusterLabels =sce_clust_SCE$groupSling )
  saveRDS(sce.sling ,  "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listSCSN.1K_sling.ptime_sce.rda" )
  
  embedded <- embedCurves(sce.sling, "UMAP.CCA")
  embedded <- slingCurves(embedded)[[1]] # only 1 path.
  embedded <- data.frame(embedded$s[embedded$ord,])
  
  # plot with pseudotime
  gg1 <- plotReducedDim( sce.sling, colour_by=("slingPseudotime_1"),
                         text_by="group", text_colour="black",
                         dimred= "UMAP.CCA") + ggtitle("trajectory ptime") +
    geom_path(data=embedded, aes(x=umapcca_1, y=umapcca_2), size=1.2) 
  # xlim(-13,-8)
  # plot with PDS
  gg2 <-plotReducedDim( sce_clust_SCE , colour_by=("PDS"),
                        text_by="group", text_colour="red",
                        dimred= "UMAP.CCA")+ ggtitle("aucell.42 PDS")
  # plot conditions
  gg0 <-  plotReducedDim( sce_clust_SCE , colour_by=("gtypeDE"),
                          text_by="group", text_colour="red",
                          dimred= "UMAP.CCA")
  # xlim(-13,-8)+ggtitle("Wt1het.del. podocytes")
  ggl <- cowplot::plot_grid( gg0, gg1, gg2, nrow = 1)
  
  # save plots
  pdf(height = 6, width = 15, file = "Supl.Fig1/listSCSN.1K_sling.ptimeDE_dimRed.pdf")
  ggl
  dev.off()
  png(height = 600, width = 1500, file = "Supl.Fig1/listSCSN.1K_sling.ptimeDE_dimRed.png")
  ggl
  dev.off() 
  
  ### DE along the trajectory
  pseudo <- TSCAN::testPseudotime( sce.sling, 
                                   pseudotime=sce.sling$slingPseudotime_1)
  pseudo$Gene.Symbol <- rownames( pseudo)
  pseudo <- pseudo[order(pseudo$p.value),]
  
  colnames(pseudo)[1:2] <-  c("log2FoldChange", "pvalue" )
  # saveRDS( pseudo, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/WRITING/PDS_manuscript/Figures/Figure1/Supl.Fig1/listSCSN.1K_sling.ptimeDE.rda")
  
    ### plot pseudotime score
    toPlot <- data.frame( score= sce.sling$slingPseudotime_1 , gt = sce.sling$groupSling)
    toPlot <- toPlot[toPlot$gt %in% c( "control", "Wt1het.del._12",
                                       "Wt1het.del._25"),]
    toPlot <- toPlot[toPlot$gt %in% c( "control", "Nphs2_4",
                                       "Nphs2_8","Nphs2_12"),]
     ggplot2::ggplot( toPlot , aes( x=score, color=gt)) + 
      # scale_color_colorblind() +  
      geom_density(size=1.5) + theme_bw() +  
      theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
             text = element_text(size = 18)) +
      geom_vline(data=ddply( toPlot, "gt",
                             summarise, grp.mean=mean(score)), 
                 aes(xintercept=grp.mean, color=gt), linetype="dashed")
  
  
  
  ### compare with PDS and other measures
  {
    
    ## enrichment plot
    { 
      X <- pseudo[ !(is.na(pseudo$pvalue)|pseudo$pvalue==1) ,]
    sstats <- 1/seq(X$Gene.Symbol)
    names(sstats) <- X$Gene.Symbol
    ppathways <- list(DS_all.all=DS_all$gene_symbol, DS_all.42=DS_all$gene_symbol[1:42])
    
    library(fgsea)
    fgsea::fgseaSimple(pathways= ppathways, stats=sstats, nperm=10000 )
    gg1 <- plotEnrichment(ppathways[[1]], sstats)+ggtitle("all 449 PD markers")
    gg2 <- plotEnrichment(ppathways[[2]], sstats)+ggtitle("42 PD markers")
    cowplot::plot_grid( gg1, gg2, nrow = 1 )
    
    }
   
    ## correlation with proteinuria for Nphs2 study
    {

      ## correlate with proteinuria
      annot_mgmg <- read.table(header = T , sep = "\t", "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/Sample_Names_KFO.csv" ,fill = T)
      
      # load rbo for sc
      X <- readRDS( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Disease_stages/rbo.PDS_Nphs2.snRNAseq_podo.rda" )
             
      
      ### collect differnt DS measures
      # calculate kendal
      kenCor <- apply(  explist[["KFO.Nphs2_sc"]] , 2,  function(X){
        X1 <- X[ match( names(spos), names(X))]
        k1 <- cor( X1 , spos , method = "spearman" , use = "pairwise.complete.obs")
        X2 <- X[ match( names(sneg), names(X))]
        k2 <- cor( X2 , sneg , method = "spearman", use = "pairwise.complete.obs" )
        kk <- k2 - k1
        # print(kk)
        return(kk)
      })
      # add sc score
      mean_tab <- data.frame( PtimeSling = sce.sling$slingPseudotime_1 ,
                              # KendalCor = kenCor , 
                              aucell.42=PDS.42[["KFO.Nphs2_sc"]] ,
                              aucell.all= PDS.all[["KFO.Nphs2_sc"]],
                              rbo.all = X,
                              sample = sce.sling$orig.ident )
      # aggregate
      mean_tab <- aggregate( .~ sample , data=mean_tab , FUN = mean )
      

      ### add score for pseudobulk
      # 42 markers
      PDS.42.pbulk <- DS_calc.func( exprMatrices = explist_scPbulk[["KFO.Nphs2_sc"]], ceilThrsh = 0.05 , 
                    DSignature= DS_all , ntop=42 )
      mean_tab$aucell.42.pbulk <- PDS.42.pbulk[mean_tab$sample]
      # all markers
      PDS.all.pbulk <- DS_calc.func( exprMatrices = explist_scPbulk[["KFO.Nphs2_sc"]], ceilThrsh = 0.05 , 
                                    DSignature= DS_all , ntop=449 )
      mean_tab$aucell.all.pbulk <- PDS.all.pbulk[mean_tab$sample]
      ### add rbo 
      # 42 markers
      # rbo.42_PDS <- readRDS(  file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/rbo.42_PDSall.20.04.2022.rda")
      mean_tab$rbo.42.pbulk <- rbo.42_PDS[["KFO.Nphs2_sc"]][mean_tab$sample]
      # all markers
      # rbo_PDS <- readRDS(  file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/rbo_PDSall.20.04.2022.rda")
      mean_tab$rbo.all.pbulk <- rbo_PDS[["KFO.Nphs2_sc"]][mean_tab$sample]
     
        
        
      # add metadata
      mean_tab$study = annot_mgmg$Genotype[ match( mean_tab$sample, annot_mgmg$CCG_Sample_ID)]
      mean_tab$batch = as.factor( annot_mgmg$run[ match( mean_tab$sample, annot_mgmg$CCG_Sample_ID)])
      mean_tab$AlbCrRatio  <- annot_mgmg$AlbCrRatio[ match( mean_tab$sample, annot_mgmg$CCG_Sample_ID)]
      
      # make a plot
      gg_list <- lapply( c( 1,4 , 2:3,5:8), function(ii){
        gg <- ggplot2::ggplot( data = mean_tab, aes( x= mean_tab[,1+ii], y=log(AlbCrRatio))) +
          geom_point( aes( color=study , shape= batch), size=6) +
          theme_bw() +  theme( text = element_text(size = 18)) + 
          geom_smooth(method='lm', se = FALSE) + ggtitle(colnames(mean_tab)[1+ii]) +
          stat_cor( size=5 )  +  geom_text(aes(label = sample ), size=4, position = "dodge")
        return(gg)
      })
      # 
      cowplot::plot_grid( plotlist = gg_list, ncol = 2)
    }
    
    ## correlation with proteinuria for all studies
    {
      ## correlate with proteinuria
      annot_mgmg <- read.table(header = T, sep = "\t","/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/Sample_Names_KFO.csv" ,fill = T)
      
        ### collect scores for sc
      rr <- readRDS("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/Disease_stages/rbo.PDS_KFO.snRNAseq_podo.rda")
      # 42 and markers
      mean_tab <- cbind.data.frame(
        PDS.42=unlist(PDS.42[c("KFO.Wt1_sc","KFO.Nphs2_sc","KFO.Pdss2_sc")]) ,
        PDS.all= unlist( PDS.all[c("KFO.Wt1_sc","KFO.Nphs2_sc","KFO.Pdss2_sc")]),
        rbo.all= unlist(rr),                   
        sample=unlist( annotlist_scSample[c("KFO.Wt1_sc","KFO.Nphs2_sc","KFO.Pdss2_sc")])) 
      mean_tab <- aggregate( .~ sample , data=mean_tab , FUN = mean )

      ### add PDS for pseudobulk
      # 42 markers
      PDS.42.pbulk <- c( DS_calc.func( exprMatrices = explist_scPbulk[["KFO.Wt1_sc"]], ceilThrsh = 0.05 ,
                                    DSignature= DS_all , ntop=42 ) , 
                         DS_calc.func( exprMatrices = explist_scPbulk[["KFO.Nphs2_sc"]], ceilThrsh = 0.05 ,
                                       DSignature= DS_all , ntop=42 ),
                         DS_calc.func( exprMatrices = explist_scPbulk[["KFO.Pdss2_sc"]], ceilThrsh = 0.05 ,
                                       DSignature= DS_all , ntop=42 ))
      mean_tab$PDS.42.pbulk <- PDS.42.pbulk[mean_tab$sample]
      # all markers
      PDS.all.pbulk <- c( DS_calc.func( exprMatrices = explist_scPbulk[["KFO.Wt1_sc"]], ceilThrsh = 0.05 ,
                                        DSignature= DS_all , ntop=nrow(DS_all) ) , 
                          DS_calc.func( exprMatrices = explist_scPbulk[["KFO.Nphs2_sc"]], ceilThrsh = 0.05 ,
                                        DSignature= DS_all , ntop=nrow(DS_all)  ),
                          DS_calc.func( exprMatrices = explist_scPbulk[["KFO.Pdss2_sc"]], ceilThrsh = 0.05 ,
                                        DSignature= DS_all , ntop=nrow(DS_all)  ))
      mean_tab$PDS.all.pbulk <- PDS.all.pbulk[mean_tab$sample]
     
      ### add rbo 
      # 42 markers
      rbo.42_PDS <- readRDS(  file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/rbo.42_PDSall.20.04.2022.rda")
      XX <- unlist( rbo.42_PDS[c("KFO.Wt1_sc","KFO.Nphs2_sc","KFO.Pdss2_sc")] )
      names(XX) <- sub(".*\\.","", names(XX))
      mean_tab$rbo.42.pbulk <- XX[mean_tab$sample]
      # all markers
      rbo_PDS <- readRDS(  file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/rbo_PDSall.20.04.2022.rda")
        XX <- unlist( rbo_PDS[c("KFO.Wt1_sc","KFO.Nphs2_sc","KFO.Pdss2_sc")] )
      names(XX) <- sub(".*\\.","", names(XX))
      mean_tab$rbo.all.pbulk <- XX[mean_tab$sample]
      
      # add metadata
      mean_tab$study = annot_mgmg$Genotype[ match( mean_tab$sample, annot_mgmg$CCG_Sample_ID)]
      mean_tab$run = as.factor( annot_mgmg$run[ match( mean_tab$sample, annot_mgmg$CCG_Sample_ID)] )
      

      mean_tab$AlbCrRatio  <- annot_mgmg$AlbCrRatio[ match( mean_tab$sample, annot_mgmg$CCG_Sample_ID)]
      mean_tab$group <- paste( mean_tab$study , mean_tab$run , sep = "_")
        
      # make a plot
      gg_list <- lapply( c(1,2,4,5,7,3), function(ii){
        gg <- ggplot2::ggplot( data = mean_tab, aes( x= mean_tab[,1+ii], y=log(AlbCrRatio))) +
          geom_point( aes( color=run  ), size=6) +
          theme_bw() +  theme( text = element_text(size = 18)) + 
          geom_smooth(method='lm', se = FALSE) + ggtitle(colnames(mean_tab)[1+ii]) +
          stat_cor( size=5 )  +  geom_text(aes(label = sample ), size=6, position = "dodge")
        return(gg)
      })
      # 
      cowplot::plot_grid( plotlist = gg_list, ncol = 2)
    }
  }
  
  ### normalise by wtype samples
  {
    summary(SCSNdata_list_sub[[2]]$orig.ident)
    SCSNdata_list_sub[[2]]$PDS.42 <- PDS.42[[34]]
    
    # normalise the score by wt samples
    ll <- c("Nphs2mut_12",  "Nphs2mut_4-5" , "Nphs2mut_6" , "Nphs2mut_8")
    Nphs2_PDS.42_wt.norm <- Reduce(c , lapply(seq(ll), function(ii) {
      # PDS.42[[34]][SCSNdata_list_sub[[2]]$groups==ll[ii]] - 
        # mean(PDS.42[[34]][SCSNdata_list_sub[[2]]$groups==sub("Nphs2mut","control",ll[ii])])
    print(mean(PDS.42[[34]][SCSNdata_list_sub[[2]]$groups==sub("Nphs2mut","control",ll[ii])]))
      } ))
     Nphs2_PDS.42_wt.norm <- Nphs2_PDS.42_wt.norm[order(Nphs2_PDS.42_wt.norm)]
    
     Nphs2_PDS.all_wt.norm <- Reduce(c , lapply(seq(ll), function(ii) {
       PDS.all[[34]][SCSNdata_list_sub[[2]]$groups==ll[ii]] - 
        mean(PDS.all[[34]][SCSNdata_list_sub[[2]]$groups==sub("Nphs2mut","control",ll[ii])])
     } ))
     Nphs2_PDS.all_wt.norm <- Nphs2_PDS.all_wt.norm[order(Nphs2_PDS.all_wt.norm)]
     
    mean_tab <- data.frame( PtimeSling = sce.sling$slingPseudotime_1[ match(  names(Nphs2_PDS.42_wt.norm) , colnames(sce.sling))] ,
                            # KendalCor = kenCor , 
                            aucell.42 = PDS.42[["KFO.Nphs2_sc"]][ match(  names(Nphs2_PDS.42_wt.norm) , names(PDS.42[["KFO.Nphs2_sc"]]))] ,
                            aucell.42.wtypeNorm = Nphs2_PDS.42_wt.norm ,
                            aucell.all = PDS.all[["KFO.Nphs2_sc"]][ match(  names(Nphs2_PDS.42_wt.norm) , names(PDS.all[["KFO.Nphs2_sc"]]))] ,
                            sample = sce.sling$orig.ident[ match(  names(Nphs2_PDS.42_wt.norm) , colnames(sce.sling))] )
    mean_tab <- aggregate( .~sample, data=mean_tab, FUN=mean)
    # add metadata
    mean_tab$study = annot_mgmg$Genotype[ match( mean_tab$sample, annot_mgmg$CCG_Sample_ID)]
    mean_tab$batch = as.factor( annot_mgmg$run[ match( mean_tab$sample, annot_mgmg$CCG_Sample_ID)])
    mean_tab$AlbCrRatio  <- annot_mgmg$AlbCrRatio[ match( mean_tab$sample, annot_mgmg$CCG_Sample_ID)]
    mean_tab$gtype  <- annot_mgmg$Genotype[ match( mean_tab$sample, annot_mgmg$CCG_Sample_ID)]
    mean_tab$age  <- annot_mgmg$Age_weeks[ match( mean_tab$sample, annot_mgmg$CCG_Sample_ID)]
    
    
    gg_list <- lapply( c( 1:4), function(ii){
      gg <- ggplot2::ggplot( data = mean_tab, aes( x= mean_tab[,1+ii], y=log(AlbCrRatio))) +
        geom_point( aes( color=study , shape= batch), size=6) +
        theme_bw() +  theme( text = element_text(size = 18)) + 
        geom_smooth(method='lm', se = FALSE) + ggtitle(colnames(mean_tab)[1+ii]) +
        stat_cor( size=5 )  +  geom_text(aes(label = sample ), size=4, position = "dodge")
      return(gg)
    })
    cowplot::plot_grid( plotlist = gg_list, ncol=2)
    
    
  }
 
  
}

#### charachterisation of PDS signature genes ####
  {
    DS_all <- read.table( header = T, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/DamageSignatures/DS_all.20.09.2023.tsv")
    allPodoGenes <- readRDS( file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_pseudotime/SCSN_allPodoGenes.rda")
    
    ### top 42 markers
    DS_all.42 <- DS_all[ 1:42, ]
    # all genes UP
    upp <- DS_all.42$gene_symbol[DS_all.42$direction_foldchange==1 ]
    upp  <-  c( toupper( upp ), 
                tx2gene$uniprotsptrembl[ ( tx2gene$external_gene_name ) %in%  upp ] , 
                tx2gene$uniprotswissprot[ ( tx2gene$external_gene_name ) %in%  upp ] )
    # all genes DOWN
    downn <- DS_all.42$gene_symbol[DS_all.42$direction_foldchange==-1 ]
    downn <- c( toupper( downn ) ,
                tx2gene$uniprotsptrembl[ ( tx2gene$external_gene_name ) %in%  downn ] , 
                tx2gene$uniprotswissprot[ ( tx2gene$external_gene_name ) %in%  downn ] )
    
    ### DE heatmap of markers
    {
      DEall_LFC <- Reduce( cbind , lapply( seq(DE_list), function(ii){
        # select 50 genes and column with LFC values
        X <- DE_list[[ii]]
        X$log2FoldChange[!is.finite(X$log2FoldChange)] <- max(X$log2FoldChange[is.finite(X$log2FoldChange)])
        X$log2FoldChange <- scale( X$log2FoldChange, center = F )
         X <- X[ match( DS_all.42$gene_symbol , 
                        X$Gene.Symbol ) , "log2FoldChange" ]
         return(X)
         }))
      colnames (DEall_LFC ) <- names(DE_list)
      rownames( DEall_LFC ) <- DS_all.42$gene_symbol
      
      # tame LFC outliers
      toPlot <- DEall_LFC
      toPlot[toPlot< -3] <- -3
      toPlot[toPlot> 3] <- 3
      # chose color scheme
      paletteLength <- 20
      myColor <-  colorspace::divergingx_hcl(palette = 'RdBu', 
                                             rev = T, n=20)
      #  save plots 
      library(viridis)
      setwd("/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_charachterise")
      heatmaply::heatmaply(toPlot, colors =myColor , na.value	=  "grey25", 
                                       file = c("heatmaply_plot.pdf", "heatmaply_plot.png"),
                  width=1500, height= 1000)
      
    }
    
    ### protein expression of markers
    {
      library(DEP)
      ll <- list.files( pattern = "proteinGroups", full.names = T , 
                        path = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/ProteinValidation/data/")
  
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
                          path = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/ProteinValidation/data/")
        
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
        # DE analysis
        podoGlom_proteome_DE <- lapply( seq(podoGlom_proteome_filt), function( ii ){
          proteome <- podoGlom_proteome_filt[[ii]]
          data_diff <- test_diff( proteome, type = "all")
          data_diff <- add_rejections(data_diff, alpha = 0.1, lfc=0)
          return(data_diff)
        } )
        
        # make plots
        lapply( podoGlom_proteome_DE , function(X) {
          plot_pca( X ) 
          plot_cor( X ,  pal = "Reds", lower = 0)
          plot_volcano( X , contrast = "all", label_size = 2, add_names = TRUE)
          
        })
        
        ### boxplots of proteins of interest
        lapply( seq(podoGlom_proteome_DE) , function(ii){
          print(ii)
         
          # check if any proteins (up) of interest are detected and if yes - plot
          if( length( intersect( upp, podoGlom_proteome_DE[[ii]]@NAMES))!=0){
            pp1 <- plot_single( podoGlom_proteome_DE[[ii]], type= "contrast",
                                proteins = upp )
          } 
          # check if any proteins (down) of interest are detected and if yes - plot
          if( length( intersect( downn, podoGlom_proteome_DE[[ii]]@NAMES))!=0) {
            pp2 <- plot_single( podoGlom_proteome_DE[[ii]],  type= "contrast",
                              proteins =  downn )
          }
          # plot 
          if( exists("pp1") && exists("pp2")) {
            gg <- cowplot::plot_grid( pp1, pp2 , labels = c("UP","DOWN") )
            } else if( exists("pp1")){
              gg <- cowplot::plot_grid( pp1, labels = "UP" )
            } else if( exists("pp2")){
              gg <- cowplot::plot_grid( pp2, labels = "DOWN" )
            } else gg <- NULL
          # save plots
          pdf(width = 10, height = 10, file = paste( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/ProteinValidation/figures/",
                                                    sub("_.*","",basename(ll)[ii]),"_PDS.50_contrasts.pdf"))
            print( gg )
          dev.off()
          png(width = 800, height = 800, file = paste( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/ProteinValidation/figures/",
                                                     sub("_.*","",basename(ll)[ii]),"_PDS.50_contrasts.png"))
            print( gg )
          dev.off()
        } )
    
        ### custom jitterplot
        lapply( seq(podoGlom_proteome_se) , function( ii )
          {
          
          ### prepare the data
          XX <- podoGlom_proteome_se[[ii]]@assays@data[[1]]
          XX <- XX[  rownames(XX) %in% c(upp, downn), ]
          # remove genes with less than 3 measurments 
          XX <- as.data.frame(  XX[rowSums( !is.na(XX) )>=3, ] )
          
          XX$name <- rownames(XX)
          
          XX.m <- reshape2::melt( XX )
          ## add metadata
          XX.m$condition <- sub( "_.*", "", XX.m$variable)
          XX.m$direction <- ifelse( XX.m$name %in%  downn ,
                                    "DOWN", "UP")
          
          # plot 
          gg <- ggplot2::ggplot( data = XX.m ,  
                           aes( x = value , 
                                y = reorder(name, value, FUN = mean, na.rm=T) ) ) + 
            geom_jitter(shape=16, position=position_jitter(width =0.1, height = 0.1),
                        size = 5, aes( color=condition)) + 
            facet_grid(rows = vars(direction) ,scales = "free_y", space = "free_y") + 
            scale_colour_colorblind(name = "proteomics\ndataset") +
            theme_bw() + labs(title = sub("_.*","",basename(ll)[ii]) ) +
            labs(y = "Gene symbol", x = "LFQ itensity") +
            theme(text=element_text(size=20 )) +
            guides(color = guide_legend(override.aes = list(size=5)))
          
          # save plot
          pdf(width = 10, height =  nrow( XX)*0.4, file = paste( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/ProteinValidation/figures/",
                                                     sub("_.*","",basename(ll)[ii]),"_PDS.50_jitterplotNOTscaled.pdf", sep = ""))
            print(gg)
          dev.off()
          png(width = 800, height =  nrow( XX)*40, file = paste( "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/ProteinValidation/figures/",
                                                                                                          sub("_.*","",basename(ll)[ii]),"_PDS.50_jitterplotNOTscaled.png", sep = ""))
            print(gg)
          dev.off()
        })
        
        ### average expression in podocytes
        controlIDs <- c("control", "c_", "Con_","Podo","control")
        podoGlom_exprLVL <- Reduce( function( x , y ) merge( x, y, by="Row.names", all=T) , 
                                    lapply( seq(podoGlom_proteome_se), function(ii){
          XX <- podoGlom_proteome_se[[ii]]@assays@data[[1]]
          XX <- scale(XX ,center = F )
          XX <- XX[  rownames(XX) %in% c(upp, downn) , ]
          XX <- XX[ rowSums( is.na(XX)) < ncol(XX)-1,]
          print( dim(XX))
          XX <- data.frame( value = rowMedians( XX, na.rm = T) ,
                               row.names = rownames(XX))
          XX$Row.names <- rownames(XX) 
          return(XX)
        }) )
        rownames(podoGlom_exprLVL) <- podoGlom_exprLVL$Row.names
        podoGlom_exprLVL <- podoGlom_exprLVL[,-1]
       
        X <- rowMeans( podoGlom_exprLVL, na.rm = T)
        X <- X[order(X)]
        par(mar=c(4,6,2,2))
        barplot( X[order(X)], las=2, horiz = T, xpd = FALSE ,xlim = c(0.8, max(X)),
                 main="average protein expression (scaled)\nof PDS markers")

        
        
    }
    
    ### gene expression of markers in bulk KFO
    {
      # load expression
      Wt1bulk <- readRDS( "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/bulk_Wt1/expr_IntExon_glvl.rda")
      Wt1bulk_annot <- read.table(header=T , "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/bulk_Wt1/annot_tab.csv")
      Nphs2bulk <- readRDS( "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/bulk_NPHS2/Nphs2expr_IntExon_glvl.rda")
      Nphs2bulk_annot <- read.table(header=T , "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/bulk_NPHS2/NPHS2_Sample_NamesAll.csv")
      colnames(Nphs2bulk) <- sub( "_.*", "",  colnames(Nphs2bulk))
      
      Wt1bulk.podo <- readRDS( "/home/tim_nevelsk/PROJECTS/PODOCYTE/RNAseq/wtPodo_bulkRNAseq/bulkRNAseq_podo_ExInt.rda")
      
    
      ### 
      top50_ID <- tx2gene$ensembl_gene_id[ match( top50PDS.mod, toupper( tx2gene$external_gene_name) )]
      top50PDS_bulk <- cbind( 
        Wt1bulk[ match( top50_ID , rownames(Wt1bulk) ) , 
              colnames(Wt1bulk) %in% Wt1bulk_annot$sample[Wt1bulk_annot$gtype=="wt"] ] ,
      Nphs2bulk[ match( top50_ID , rownames(Nphs2bulk) ) ,
                colnames(Nphs2bulk) %in% Nphs2bulk_annot$CCG.Sample.ID[Nphs2bulk_annot$gtype=="wt"] ] 
      # Wt1bulk.podo[ which( rownames(Wt1bulk.podo) %in% top50_ID) ,  1:3]
      )
      
      top50PDS_bulk <- as.data.frame( apply( top50PDS_bulk, 2, function(x) (x/sum(x))*1000000) )
        
      top50PDS_bulk$Gene.symbol <- top50PDS.mod
        
      top50PDS_bulk.melt <-  reshape2::melt(top50PDS_bulk, value.name = "cpm")
      top50PDS_bulk.melt$Experiment <- c( rep( "Wt1 wt gloms", 300),
                        rep( "Nphs2 wt podo",300))
      
      ggplot2::ggplot( data = top50PDS_bulk.melt , 
                       aes( x= log(cpm) , 
                            y = reorder(Gene.symbol,cpm, FUN = median) ) ) + 
        geom_boxplot()+
        geom_jitter(shape=16, position=position_jitter(0.02),size = 1,
                    aes( color=Experiment)) + 
        theme_bw() + labs(y = "Gene symbol") + 
        theme(axis.text=element_text(face="bold"), legend.title=element_text(size=16),
              axis.title=element_text(size=16), legend.text=element_text(size=16)) +
        scale_colour_colorblind() + guides(color = guide_legend(override.aes = list(size=5)))
      
    }
    
    ### gene expression of markers in snRNAseq 
    {
      
      # podocytes
      {
        sn_subsPodo_wtype <- lapply(list( sce_subsPodo_nphs2 , sce_subsPodo_wt1 ,sce_subsPodo_cem ), function(datt){
        XX <- subset( datt , subset = gtype %in% c("wildtype", "Podocin.wt.wt", "Wt1.wt.wt."))
        
        XX <- subset( XX ,cells = sample( colnames(XX), 1000) )
        return(XX)
      })
      
      # merge all KFO wt podo
      sn_subsPodo_wtype <- merge( sn_subsPodo_wtype[[1]] , 
                                  y = c( sn_subsPodo_wtype[[2]],
                                         sn_subsPodo_wtype[[3]]), merge.data = T,
                                  add.cell.ids = c("Nphs2___", "Wt1___", "Cem___"))
      # sn_subsPodo_wtype <- ScaleData( sn_subsPodo_wtype )
      # downsample to 1000 cells
      sn_subsPodo_wt.PDS50 <- subset(sn_subsPodo_wtype , cells = sample( colnames(sn_subsPodo_wtype), 1000) )
      # select PDS genes availale in the dataset
      ggenes <- stringr::str_to_title(union(top50PDS.mod, top50PDS))
      ggenes <- ggenes[ ggenes %in% rownames(sn_subsPodo_wt.PDS50@assays$RNA@data)]
      # subset expr data to 50 PDS markers
      sn_subsPodo_wt.PDS50 <- as.data.frame( sn_subsPodo_wt.PDS50@assays$RNA@data[ ggenes ,] )
      # sn_subsPodo_wt.PDS50 <- as.data.frame( sn_subsPodo_wt.PDS50@assays$RNA@counts[ which( rownames(sn_subsPodo_wt.PDS50@assays$RNA@counts) %in% stringr::str_to_title(union(top50PDS.mod, top50PDS))),] )
      sn_subsPodo_wt.PDS50$Gene.symbol <- rownames(sn_subsPodo_wt.PDS50)
      sn_subsPodo_wt.PDS50.melt <- reshape2::melt( sn_subsPodo_wt.PDS50 )
      sn_subsPodo_wt.PDS50.melt$experiment <- sub( "___.*","",sn_subsPodo_wt.PDS50.melt$variable)
      # sn_subsPodo_wt.PDS50$log10value <- log10(sn_subsPodo_wt.PDS50$value+1)
      # sn_subsPodo_wt.PDS50$log10value[sn_subsPodo_wt.PDS50$log10value==0] <- NA
      
      XX <- sn_subsPodo_wt.PDS50.melt[sample(1:nrow(sn_subsPodo_wt.PDS50.melt)),]
      ggplot2::ggplot( data = XX ,  
                       aes( x = value , 
                            y = reorder(Gene.symbol,value, FUN = mean) ) ) + 
        geom_jitter(shape=16, position=position_jitter(0.08),size = 0.8,
                    aes( color=experiment)) + 
        scale_colour_colorblind(name = "snRNAseq\nexperiment") +
        theme_bw() + geom_boxplot( outlier.size = 0, alpha = 0.0)+
        labs(y = "Gene symbol", x = "log( counts_libNorm )") +
        theme(axis.text=element_text(size=12 , face="bold"), legend.title=element_text(size=16),
              axis.title=element_text(size=16), legend.text=element_text(size=16)) +
        guides(color = guide_legend(override.aes = list(size=5)))
      }
      
      # all clusts
      {
        sce_nphs2 <- readRDS("/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_Nphs2/Seurat/snRNAseq_Nphs2_premRNA.Seurat.SubsTOviz.rda")
        sce_subset_nphs2 <- subset(sce_nphs2 , idents=c(0,1,2,3,5))
        sce_subset_nphs2$celltype <- droplevels(revalue( sce_subset_nphs2$clust ,  c("0"="endothelial", "1"="mesangial","2"="podocytes",
                                                                                 "3"="proxim.tubule","5"="podocytes")))
        
        sce_wt1 <- readRDS( "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/snRNAseq_WT1hetdel/Seurat/snRNAseq_Wt1hd_premRNA.Seurat.SubsTOviz.rda" )
        sce_subset_wt1 <- subset(sce_wt1 , idents=0:3)
        sce_subset_wt1$celltype <- droplevels(revalue( sce_subset_wt1$clust ,  c("0"="endothelial", "1"="mesangial","2"="podocytes",
                                                                                  "3"="proxim.tubule")))
        
        sce_cem <- readRDS(  "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Cem_data/Seurat/snRNAseq_Cem_premRNA.Seurat.rda" )
        sce_subs_cem <- subset(sce_cem , downsample=500, idents=0:3)
        sce_subs_cem$celltype <- droplevels(revalue( sce_subs_cem$seurat_clusters ,  c("0"="endothelial", "1"="mesangial","2"="podocytes",
                                                                             "3"="proxim.tubule")))
          
        sn_subs_wtype <- merge( sce_subset_nphs2 , 
                                y = c( sce_subset_wt1 ,
                                       sce_subs_cem ), merge.data = TRUE )
        # select PDS genes availale in the dataset
        ggenes <- stringr::str_to_title(union(top50PDS.mod, top50PDS))
        ggenes <- ggenes[ ggenes %in% rownames(sn_subs_wtype@assays$RNA@data)]
        # subset exprData to only PDS markers
        sn_subs_wt.PDS50 <- as.data.frame( sn_subs_wtype@assays$RNA@data[ ggenes,] )
        colnames(sn_subs_wt.PDS50) <- c( paste("nphs2___", sce_subset_nphs2$celltype,"__",colnames(sce_subset_nphs2), sep = ""), 
                                         paste("wt1___",sce_subset_wt1$celltype,"__", colnames(sce_subset_wt1) ,sep = ""), 
                                         paste("cem___",sce_subs_cem$celltype,"__", colnames(sce_subs_cem) ,sep = "")
        )
        sn_subs_wt.PDS50 <- sn_subs_wt.PDS50[,sample( colnames(sn_subs_wt.PDS50),1000)]
        
        # sn_subsPodo_wt.PDS50 <- as.data.frame( sn_subsPodo_wt.PDS50@assays$RNA@counts[ which( rownames(sn_subsPodo_wt.PDS50@assays$RNA@counts) %in% stringr::str_to_title(union(top50PDS.mod, top50PDS))),] )
        sn_subs_wt.PDS50$Gene.symbol <- rownames(sn_subs_wt.PDS50)
        sn_subs_wt.PDS50.melt <- reshape2::melt( sn_subs_wt.PDS50 )
        sn_subs_wt.PDS50.melt$experiment <- sub( "___.*","",sn_subs_wt.PDS50.melt$variable)
        sn_subs_wt.PDS50.melt$cluster <- sub( "__.*", "", sub( ".*___","",sn_subs_wt.PDS50.melt$variable) )
        # sn_subsPodo_wt.PDS50$log10value <- log10(sn_subsPodo_wt.PDS50$value+1)
        # sn_subsPodo_wt.PDS50$log10value[sn_subsPodo_wt.PDS50$log10value==0] <- NA
        
        XX <- sn_subs_wt.PDS50.melt[sample(1:nrow(sn_subs_wt.PDS50.melt)),]
        ggplot2::ggplot( data = XX ,  
                         aes( x = value , 
                              y = reorder(Gene.symbol,value, FUN = mean) ) ) + 
          scale_colour_colorblind(name = "snRNAseq\nexperiment") +
          theme_bw() + 
          geom_jitter(shape=16, position=position_jitter(0.08),size = 0.8,
                      aes( color=cluster)) +
          geom_boxplot( outlier.size = 0, alpha = 0.0,  aes( color=cluster))+
          theme_bw() + labs(y = "Gene symbol", x = "log( counts_libNorm )")+
          theme(axis.text=element_text(face="bold"), legend.title=element_text(size=16),
                axis.title=element_text(size=16), legend.text=element_text(size=16)) +
          guides(color = guide_legend(override.aes = list(size=5)))
        
      }
      
      
     ### make vln plots 
      listKFO.podo_seur.PDS <- readRDS(file="/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/disease.score/listKFOseur.PDS.rda")
      source('~/PROJECTS/myCode/stacked_vlnPlot_seurat.R')
      
      # vln plots for individual genes
      # StackedVlnPlot( 
      #   obj = sn_subsPodo_wtype, 
      #   features = stringr::str_to_title( DS_all$gene_symbol[1:5]))

      ### all up VS all dow
      allPodoGenes.inter <- Reduce( intersect,  
                                lapply( seq(listSCSN.PDS) ,
                                        function(ii){
                                          # all genes expressed in 1%of cells in at least one group
                                          rownames(listSCSN.PDS[[ii]])[ rowSums( 
                                            sapply( unique(listSCSN.PDS[[ii]]$group), 
                                                    function(XX){
                                                      rowSums(round(listSCSN.PDS[[ii]]@assays$RNA@counts[,listSCSN.PDS[[ii]]$group==XX])>0) > 
                                                        0.01*ncol(listSCSN.PDS[[ii]]@assays$RNA@counts[,listSCSN.PDS[[ii]]$group==XX])
                                                    })) >0 ]
                                          }))
      # calculate means
      scsnRNAseq_means <- lapply( seq(listSCSN.PDS) , function(ii){
        DD<- data.frame( mean_expr= rowMeans( 
          listSCSN.PDS[[ii]]@assays$RNA@data[allPodoGenes.inter,] ),
          dataset=names(listSCSN.PDS)[ii],
          FC_drctn=DS_all$direction_foldchange[ match( 
            allPodoGenes.inter, DS_all$gene_symbol)] )
        DD$scsn.Rank <- rank(-DD$mean_expr)
        return(DD)
      })
      scsnRNAseq_mean.rank <- rowMeans( Reduce(cbind , 
                                     lapply(scsnRNAseq_means, 
                                            function(X) X$scsn.Rank )))
      names(scsnRNAseq_mean.rank)<-allPodoGenes.inter
      scsnRNAseq_mean.rank <- scsnRNAseq_mean.rank[order(scsnRNAseq_mean.rank)]
      saveRDS(scsnRNAseq_mean.rank, "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/scsnRNAseq.podo_gene.inter_mean.rank.rda")
      
      # prepare for plot
      scsnRNAseq_means.melt <- Reduce( rbind, scsnRNAseq_means)
      scsnRNAseq_means.melt$FC_drctn <- ifelse( is.na(scsnRNAseq_means.melt$FC_drctn),
      "not.DS",ifelse(scsnRNAseq_means.melt$FC_drctn<0, "DS.DOWN", "DS.UP"))
      scsnRNAseq_means.melt$dataset <- factor( as.factor(scsnRNAseq_means.melt$dataset),
                                          levels=names(listSCSN.PDS))
      
      # plot
      ggplot(data = scsnRNAseq_means.melt, aes(color=FC_drctn,
                                          x = log10(mean_expr)))+
        geom_density( lwd=2) + theme_bw()+ 
        theme(text = element_text(size=20))+
        facet_grid(rows=vars(dataset), scales = "free_y" )
    }
    
    ### gene expression of markers in scRNAseq
    {
      GSE146912_doxo <- readRDS( file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/GSE146912_RAW/Seurat/GSE146912_doxo_seurat.rda")
      GSE146912_doxo_podo <- subset( GSE146912_doxo , cells = colnames(GSE146912_doxo)[ colnames(GSE146912_doxo)  %in% rownames(GSE146912_doxo@meta.data)[ GSE146912_doxo@meta.data$celltype.stim =="1_contol" ]] )
      Idents(GSE146912_doxo) <- GSE146912_doxo$seurat_clusters
      GSE146912_doxo_subset <- subset(GSE146912_doxo, downsample=500,
                                      idents = 0:3 )
      GSE146912_doxo_subset$celltype <- droplevels(revalue( GSE146912_doxo_subset$celltype ,  c("0"="endothelial", "1"="podocytes","2"="mesangial", "3"="SMC")))
      GSE146912_cd2ap <- readRDS( file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/RNAseq/public/PDS/sc/GSE146912_scRNAseq/Seurat/GSE146912_cd2ap_seurat.rda")
      GSE146912_cd2ap_podo <- subset( GSE146912_cd2ap , cells = colnames(GSE146912_cd2ap)[ colnames(GSE146912_cd2ap)  %in% rownames(GSE146912_cd2ap@meta.data)[ GSE146912_cd2ap@meta.data$celltype.stim =="2_control" ]] )
      Idents(GSE146912_cd2ap) <- GSE146912_cd2ap$seurat_clusters
      GSE146912_cd2ap_subset <- subset(GSE146912_cd2ap,  downsample=500,
                                       idents = 0:2 )
      GSE146912_cd2ap_subset$celltype <- droplevels(revalue( GSE146912_cd2ap_subset$celltype ,  c("0"="endothelial", "1"="mesangial","2"="podocytes")))
      GSE146912_btbr <- readRDS( file = "/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/GSE146912_RAW/Seurat/GSE146912_btbr_seurat.rda")
      GSE146912_btbr_podo <- subset( GSE146912_btbr , cells = colnames(GSE146912_btbr)[ colnames(GSE146912_btbr)  %in% rownames(GSE146912_btbr@meta.data)[ GSE146912_btbr@meta.data$celltype.stim =="2_control" ]] )
      Idents(GSE146912_btbr) <- GSE146912_btbr$seurat_clusters
      GSE146912_btbr_subset <- subset(GSE146912_btbr, downsample=500,
                                      idents = 0:2 )
      GSE146912_btbr_subset$celltype <- droplevels(revalue( GSE146912_btbr_subset$celltype ,  c("0"="endothelial", "1"="mesangial","2"="podocytes")))
      
      ## podocytes
      {
        GSE146912_list <- lapply(list(GSE146912_doxo_podo , GSE146912_cd2ap_podo, GSE146912_btbr_podo),
                                 function(XX) subset(XX, cells = sample( colnames(XX), 500) ))
        sc_subsPodo_wtype <- merge( GSE146912_list[[1]] , 
                                    y = c( GSE146912_list[[2]],
                                           GSE146912_list[[3]] ), merge.data = TRUE,
                                    add.cell.ids = c("doxo___", "cd2ap___", "btbr___"))
        # downsample to 100 cells
        sc_subsPodo_wt.PDS50 <- subset(sc_subsPodo_wtype , cells = sample( colnames(sc_subsPodo_wtype), 1000))
        # select PDS genes availale in the dataset
        ggenes <- stringr::str_to_title(union(top50PDS.mod, top50PDS))
        ggenes <- ggenes[ ggenes %in% rownames(sc_subsPodo_wt.PDS50@assays$data@data)]
        # subset exprData to only PDS markers
        sc_subsPodo_wt.PDS50 <- as.data.frame( sc_subsPodo_wt.PDS50@assays$data@data[ ggenes,] )
        # sn_subsPodo_wt.PDS50 <- as.data.frame( sn_subsPodo_wt.PDS50@assays$RNA@counts[ which( rownames(sn_subsPodo_wt.PDS50@assays$RNA@counts) %in% stringr::str_to_title(union(top50PDS.mod, top50PDS))),] )
        sc_subsPodo_wt.PDS50$Gene.symbol <- rownames(sc_subsPodo_wt.PDS50)
        sc_subsPodo_wt.PDS50.melt <- reshape2::melt( sc_subsPodo_wt.PDS50 )
        sc_subsPodo_wt.PDS50.melt$experiment <- sub( "___.*","",sc_subsPodo_wt.PDS50.melt$variable)
        # sn_subsPodo_wt.PDS50$log10value <- log10(sn_subsPodo_wt.PDS50$value+1)
        # sn_subsPodo_wt.PDS50$log10value[sn_subsPodo_wt.PDS50$log10value==0] <- NA
        
        XX <- sc_subsPodo_wt.PDS50.melt[sample(1:nrow(sc_subsPodo_wt.PDS50.melt)),]
        ggplot2::ggplot( data = XX ,  
                         aes( x = value , 
                              y = reorder(Gene.symbol,value, FUN = mean) ) ) + 
          geom_jitter(shape=16, position=position_jitter(0.08),size = 0.8,
                      aes( color=experiment)) +
          scale_colour_colorblind(name = "scRNAseq\nexperiment") +
          theme_bw() + geom_boxplot( outlier.size = 0, alpha = 0.0)+
          theme_bw() + labs(y = "Gene symbol", x = "log( counts_libNorm )")+
          theme(axis.text=element_text(face="bold"), legend.title=element_text(size=16),
                axis.title=element_text(size=16), legend.text=element_text(size=16)) +
          guides(color = guide_legend(override.aes = list(size=5)))
        
      }
     
      ## all cells
      {

        sc_subs_wtype <- merge( GSE146912_doxo_subset , 
                                    y = c( GSE146912_cd2ap_subset ,
                                           GSE146912_btbr_subset ), merge.data = TRUE )
        # select PDS genes availale in the dataset
        ggenes <- stringr::str_to_title(union(top50PDS.mod, top50PDS))
        ggenes <- ggenes[ ggenes %in% rownames(sc_subs_wtype@assays$data@data)]
        # subset exprData to only PDS markers
        sc_subs_wt.PDS50 <- as.data.frame( sc_subs_wtype@assays$data@data[ ggenes,] )
        colnames(sc_subs_wt.PDS50) <- c( paste("doxo___", GSE146912_doxo_subset$celltype,"__",colnames(GSE146912_doxo_subset), sep = ""), 
           paste("cd2ap___",GSE146912_cd2ap_subset$celltype,"__", colnames(GSE146912_cd2ap_subset) ,sep = ""), 
           paste("btbr___",GSE146912_btbr_subset$celltype,"__", colnames(GSE146912_btbr_subset) ,sep = "")
        )
        sc_subs_wt.PDS50 <- sc_subs_wt.PDS50[,sample(colnames(sc_subs_wt.PDS50),1000)]
        
        # sn_subsPodo_wt.PDS50 <- as.data.frame( sn_subsPodo_wt.PDS50@assays$RNA@counts[ which( rownames(sn_subsPodo_wt.PDS50@assays$RNA@counts) %in% stringr::str_to_title(union(top50PDS.mod, top50PDS))),] )
        sc_subs_wt.PDS50$Gene.symbol <- rownames(sc_subs_wt.PDS50)
        sc_subs_wt.PDS50.melt <- reshape2::melt( sc_subs_wt.PDS50 )
        sc_subs_wt.PDS50.melt$experiment <- sub( "___.*","",sc_subs_wt.PDS50.melt$variable)
        sc_subs_wt.PDS50.melt$cluster <- sub( "__.*", "", sub( ".*___","",sc_subs_wt.PDS50.melt$variable) )
        # sn_subsPodo_wt.PDS50$log10value <- log10(sn_subsPodo_wt.PDS50$value+1)
        # sn_subsPodo_wt.PDS50$log10value[sn_subsPodo_wt.PDS50$log10value==0] <- NA
        
        XX <- sc_subs_wt.PDS50.melt[sample(1:nrow(sc_subs_wt.PDS50.melt)),]
        ggplot2::ggplot( data = XX ,  
                         aes( x = value , 
                              y = reorder(Gene.symbol,value, FUN = mean) ) ) + 
          geom_jitter(shape=16, position=position_jitter(0.08),size = 0.8,
                       aes( color=cluster)) +
          scale_colour_colorblind(name = "scRNAseq\nexperiment") +
          theme_bw() + geom_boxplot( outlier.size = 0, alpha = 0.0,  aes( color=cluster))+
          theme_bw() + labs(y = "Gene symbol", x = "log( counts_libNorm )")+
          theme(axis.text=element_text(face="bold"), legend.title=element_text(size=16),
                axis.title=element_text(size=16), legend.text=element_text(size=16)) +
          guides(color = guide_legend(override.aes = list(size=5)))
        
      }
      }

    ### functional annotation
    {
      gene.bckgrnd_eID <- tx2gene$entrezgene_id[ match( allPodoGenes, tx2gene$external_gene_name)]
      gene.bckgrnd_eID <- as.character( gene.bckgrnd_eID[!is.na(gene.bckgrnd_eID)] )
      tx2gene <- getBM(attributes=c( "entrezgene_id", "external_gene_name"),  mart = mart_mouse)
      DS.42.eID <- as.character( tx2gene$entrezgene_id[ match( DS_all$gene_symbol[1:42],
                                                               tx2gene$external_gene_name)] )
      
      # test with cluster profiler
      DS.42_cProfiler <- cProfiler.GKR( ggenes.eID=DS.42.eID , 
                                        gene.bckgrnd_eID=gene.bckgrnd_eID )
      
      
      # make a df for visualising
      toPlot <- Reduce( rbind , lapply(DS.42_cProfiler, function(X) X@result))
        
      # toPlotGO <- toPlot[ toPlot$ID %in% toPlot$ID[toPlot$qvalue<0.01] &
      #                       toPlot$ONTOLOGY%in%c("BP","MF","CC"), ]
      # toPlotGO$Description <- sapply( toPlotGO$Description , wrap_text, 40)
      toPlotPTH <- unique( toPlot[(toPlot$pvalue< 0.1 & toPlot$ONTOLOGY%in%c("KEGG","REACT")) | 
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
                                 scales = "free_y", space = "free" ) + theme( text = element_text(size=20)) + 
        ylab("pathway names")
      
      ### cluster terms
      gg<- enrichplot::emapplot_cluster(enrichplot::pairwise_termsim(
        DS.42_cProfiler$GO.enrich), color = "pvalue", cex_label_group=1.2,
        showCategory = 15)
      # p3 <-  enrichplot::emapplot_cluster(enrichplot::pairwise_termsim(
      #   datt$KEGG.enrich), color = "pvalue")
      # p2  <- enrichplot::emapplot_cluster(enrichplot::pairwise_termsim(
      #   DS.42_cProfiler$REACT.enrich), color = "pvalue", cex_label_group=1.2)
      # cowplot::plot_grid( plotlist =  list(p1,p2), nrow = 1, labels = c("GO","REACT"))
     pdf(height = 4, width = 6, file="/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/PDS_charachterise/DS.42_cProfiler_GOclustered")
      print(gg)
     dev.off()
      
      }
   
  }

