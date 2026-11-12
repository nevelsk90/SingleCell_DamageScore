# ###################################################### #
#  prepare SCSN podocyte data for PDS analysis, select TFs, calculate PDS #
# ###################################################### #
mallinfo::mallinfo()
mallinfo::malloc.trim()
gc()

options(connectionObserver = NULL)
library( Matrix)
library( AUCell )
library( org.Mm.eg.db)
library( ggplot2 )
library( reshape2)
library( plyr )
library( Seurat )
library( viridis )
library( ggthemes )
library( ggpubr )
# library( biomaRt )
# mart_homo <- useMart( "ensembl",dataset="hsapiens_gene_ensembl" , host="www.ensembl.org")
# mart_mouse <- useMart( "ensembl",dataset="mmusculus_gene_ensembl")
# tx2gene <- getBM( attributes=c('ensembl_gene_id', 'external_gene_name',"entrezgene_id"),  mart = mart_mouse)
# tx2prot <-  getBM( attributes=c("uniprotswissprot","uniprot_gn_id", 'external_gene_name',"mgi_symbol"),  mart = mart_mouse)

setwd("/data/user/tpadvits/PROJECTS/PDS")

#### load functions
# calculate damage score
source("AUCell_script.r")
# source("myCode/usefulRfunc.r")



# load the podo damage signature
DS_all <- read.table( header = T, "DS_all.20.09.2023.tsv")
# load KFO and GSE146912 podocyte data
listSCSN <- readRDS( "listSCSN.22.12.23.rda" )


### calculate PDS for all cells
  listSCSN.PDS <- lapply( seq(listSCSN), function(ii)
    {
    print(ii)
    newSeu <- listSCSN[[ii]]


    ## exclude genes non.expressed in more than certain percentage of cells
    expr <- newSeu@assays$RNA@counts
    expr <- expr[ rowSums( round(expr) > 0 ) > 0.005*ncol(expr), ]

    set.seed(42)
    # # calculate damage signatures
    newSeu@meta.data$PDS <- DS_calc.func( exprMatrices = expr ,
                           DSignature = DS_all ,
                           ntop = 42 , wghtd = T,
                           ceilThrsh =  0.05 )
    

    # adjust ceilThrsh based on a total number of genes in the matrix!
    # the top should include ~ 1K genes

  return(newSeu)
})
 names(listSCSN.PDS) <- names(listSCSN )
  ## save
  saveRDS( listSCSN.PDS , file="listSCSN.PDS_22.12.23.rda")
  
 toPlot <-listSCSN.PDS$Nphs2@meta.data[ listSCSN.PDS$Nphs2@meta.data$sample %in% 
                                          c("140739", "140738", "140740", "140741", "139919",
                                            "139917", "139921", "139913", "139915", "139911"),] 
ggplot( toPlot , aes(x=sample, y=PDS , fill=group))+
  geom_boxplot()+theme_bw() # + ggtitle("filt.0.03, thrshld=0.05")

