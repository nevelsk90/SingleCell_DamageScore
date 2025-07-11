# ### show and test alignment of distribution of PDS in ctrl VS xprmnt ### #

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

#### plot densitites ####
## density Plot of distributions
PDS.list.tab <-Reduce( rbind, lapply( seq(listSCSN.1K.sampl), function(ii){
  datt <-  listSCSN.1K.sampl[[ii]]@meta.data[ , c( "gtypeDE", "PDS","sample")]
  datt$dataSet <- names(listSCSN.1K.sampl)[ii]
  
  return(datt)
}) )
PDS.list.tab$seqType <- ifelse(PDS.list.tab$dataSet%in% c(
  "Nphs2","Wt1","Lmx1b","Pdss2") , "sn","sc" )

toPlot <- PDS.list.tab[PDS.list.tab$dataSet!="Lmx1b",]
toPlot$group <- paste0(toPlot$gtypeDE, toPlot$dataSet)

gg00 <- ggplot(toPlot ,
               aes(x=PDS, color=gtypeDE, group= sample))+
  geom_line( stat="density", adjust=1.2,
             linewidth=1.2, alpha=0.5)+
  theme_minimal() + theme( legend.position = "none") + 
  scale_color_colorblind()

gg0 <- ggplot(toPlot , 
              aes(x=PDS, color=gtypeDE, group= sample))+
  geom_line( stat="density", adjust=1.2,
             linewidth=1.2, alpha=0.5)+
  facet_wrap( vars(seqType), ncol = 1)+
  
  theme_minimal()+ theme( legend.position = "bottom") + scale_color_colorblind()


pdf( height = 5, width = 5, file = "/home/tim_nevelsk/PROJECTS/PODOCYTE/DiseaseScore/Development_Validation/PDSctrlVSxprmnt_density.Smpls.v2.pdf")
cowplot::plot_grid(plotlist =  list(gg00 , gg0), ncol = 1,
                   rel_heights =c(0.5,1, 4) )
dev.off()

#### Kolmogorov-smirnov distances between distributions ####
PDS.list <- unlist( lapply( seq(listSCSN.1K.sampl), function(ii){
  datt <-  listSCSN.1K.sampl[[ii]]@meta.data[ , c( "gtypeDE", "PDS")]
  datt$dataSet <- names(listSCSN.1K.sampl)[ii]
  datt.ct <- datt[datt$gtypeDE =="control" , ]
  datt.xp <- datt[datt$gtypeDE !="control", ]
  dattSS <- list(datt.ct$PDS , datt.xp$PDS)
  names(dattSS) <- c( paste0( names(listSCSN.1K.sampl)[ii], "_ctrl") ,
                      paste0(names(listSCSN.1K.sampl)[ii], "_xprmt"))
  return(dattSS)
}) ,recursive = F)

PDS.list <- PDS.list[lengths(PDS.list) > 0]

combinations <- combn(names(PDS.list), 2, simplify = FALSE)
PDS.kstest <- Reduce( rbind ,lapply(combinations, function(pair)
{
  x <- PDS.list[[pair[1]]]
  y <- PDS.list[[pair[2]]]
  test <- tryCatch(  ks.test(x, y),  error = function(e) NA)
  if( is.na(test)) {
    sstat <- NA 
    ppval <- NA
  } else {
    sstat <- test$statistic 
    ppval <- test$p.value
  }
  data.frame(
    var1 = pair[1],
    var2 = pair[2],
    gtypeVar1 = sub(".*_","",pair[1]),
    gtypeVar2 = sub(".*_","",pair[2]),
    
    statistic = sstat,
    p_value = ppval,
    stringsAsFactors = FALSE
  ) } ) )

PDS.kstest$seqType1 <- ifelse( sub("_.*","",PDS.kstest$var1) %in% c(
  "Nphs2","Wt1","Lmx1b","Pdss2") , "sn","sc" )
PDS.kstest$seqType2 <- ifelse(sub("_.*","",PDS.kstest$var2) %in% c(
  "Nphs2","Wt1","Lmx1b","Pdss2") , "sn","sc" )

wilcox.test( PDS.kstest$statistic[ PDS.kstest$gtypeVar1=="ctrl" &
                                     PDS.kstest$gtypeVar2=="ctrl"   ] ,
             PDS.kstest$statistic[PDS.kstest$gtypeVar1!="ctrl" &
                                    PDS.kstest$gtypeVar2!="ctrl"],
             alternative =  "less")


# ## plot a ditance heatmap
# PDS.kstest.matrixP <- triangle_vec_to_matrix( PDS.kstest$p_value, )
# PDS.kstest.matrixS <- triangle_vec_to_matrix( PDS.kstest$statistic )
# 
# colnames(PDS.kstest.matrixS) <- rownames(PDS.kstest.matrixS)<-
#   colnames(PDS.kstest.matrixP) <- rownames(PDS.kstest.matrixP)<- names(PDS.list)
# 
# PDS.kstest.matrixS[PDS.kstest.matrixS==0] <- NA
# pheatmap::pheatmap( -log10(PDS.kstest.matrixS) ,
#          cluster_rows = T,
#          cluster_cols = T,
#           # display_numbers = PDS.kstest.matrixS  ,
#          color = colorRampPalette(c("white", "red"))(100),
#          main = "Pairwise KS D-statistic Heatmap")

PDS.list.sampl <- unlist( lapply( seq(listSCSN.1K.sampl), function(ii){
  datt <-  listSCSN.1K.sampl[[ii]]@meta.data[ , c( "gtypeDE", "PDS","sample")]
  datt$dataSet <- names(listSCSN.1K.sampl)[ii]
  
  lst <- split(datt, datt$sample)
  
  names(lst) <- paste(  names(listSCSN.1K.sampl)[ii], names(lst),
                        datt$gtypeDE[ match( names(lst), datt$sample)],
                        sep="__")
  return(lst)
}) ,recursive = F)


combinations.S <- combn(names(PDS.list.sampl), 2, simplify = FALSE)
PDS.kstest.S <- Reduce( rbind ,lapply(combinations.S, function(pairS)
{
  x <- PDS.list.sampl[[pairS[1]]]$PDS
  y <- PDS.list.sampl[[pairS[2]]]$PDS
  test <- tryCatch(  ks.test(x, y),  error = function(e) NA)
  if( is.na(test)) {
    sstat <- NA 
    ppval <- NA
  } else {
    sstat <- test$statistic 
    ppval <- test$p.value
  }
  data.frame(
    var1 = pairS[1],
    var2 = pairS[2],
    gtypeVar1 = sub(".*__","",pairS[1]),
    gtypeVar2 = sub(".*__","",pairS[2]),
    
    statistic = sstat,
    p_value = ppval,
    stringsAsFactors = FALSE
  ) }
) )

PDS.kstest.S$seqType1 <- ifelse( sub("__.*","",PDS.kstest.S$var1) %in% c(
  "Nphs2","Wt1","Lmx1b","Pdss2") , "sn","sc" )
PDS.kstest.S$seqType2 <- ifelse(sub("__.*","",PDS.kstest.S$var2) %in% c(
  "Nphs2","Wt1","Lmx1b","Pdss2") , "sn","sc" )

wilcox.test( PDS.kstest.S$statistic[ PDS.kstest.S$gtypeVar1=="control" &
                                       PDS.kstest.S$gtypeVar2=="control"   ] ,
             PDS.kstest.S$statistic[PDS.kstest.S$gtypeVar1!="control" &
                                      PDS.kstest.S$gtypeVar2!="control"],
             alternative =  "less")

wilcox.test( PDS.kstest.S$statistic[ PDS.kstest.S$gtypeVar1=="control" &
                                       PDS.kstest.S$gtypeVar2=="control"&
                                       PDS.kstest.S$seqType1=="sn" &
                                       PDS.kstest.S$seqType2=="sn" ] ,
             PDS.kstest.S$statistic[PDS.kstest.S$gtypeVar1!="control" &
                                      PDS.kstest.S$gtypeVar2!="control" &
                                      PDS.kstest.S$seqType1=="sn" &
                                      PDS.kstest.S$seqType2=="sn" ],
             alternative =  "less")
wilcox.test(PDS.kstest.S$statistic[ PDS.kstest.S$gtypeVar1=="control" &
                                      PDS.kstest.S$gtypeVar2=="control"&
                                      PDS.kstest.S$seqType1!="sn" &
                                      PDS.kstest.S$seqType2!="sn" ] ,
            PDS.kstest.S$statistic[PDS.kstest.S$gtypeVar1!="control" &
                                     PDS.kstest.S$gtypeVar2!="control" &
                                     PDS.kstest.S$seqType1!="sn" &
                                     PDS.kstest.S$seqType2!="sn" ],
            alternative =  "less")

