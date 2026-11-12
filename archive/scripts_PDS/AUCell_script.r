##############################################################
### calculating damage score, using a gene set and AUCell ####
########################################################## ####

DS_calc.func <- function(exprMatrices , DSignature , 
                         geneIDname = "gene_symbol", 
                         useOrder=c( "mean_rank", "median_rank"),
                         ntop= 42 , ceilThrsh=0.05 , progStat =F ,
                         wghtd=T )
  {
  require(AUCell)
  require(GSEABase)
  
  ## exprMatrices - expression matrix with genes in rows and samples/cells
  # in columns, row names ID type should be the same as in Signature[geneIDname,]
  ## DSignature - an output of damage_signature.func - a dataframe containing columns gene_symbol, mean_rank and direction_foldchange
  ## geneIDname - a charachter or a numeric value, specifying which 
  # DSignature column to use for gene names 
  # useOrder - whether to use mean_rank or median_rank columns of DSignature for ordering genes
  # nTop - how many genes from the top of a gene signature should be used to calculate the score
  # ceilThrsh - a threshold to calculate the AUC
  # progStat - to print progress messages or not
  # genesUP and genesDOWN are charachter vectors containing the names of genes UP and DOWN regulated upon damage 
  ## wghtd - logical indicating whether score for up and down-regulated genes 
  # should be weghted by the number of genes
  
  ## evaluate choices
  useOrder <- match.arg(useOrder)
  
  # order genes according to mean or median
  DSignature <- DSignature[ order(DSignature[[useOrder]] ),]
  
  ### make a collection of gene sets  to use with AUCell
  genesUP = DSignature[ 1:ntop, geneIDname][ 
    which( DSignature$direction_foldchange[1:ntop]==1)]
  genesDOWN =  DSignature[ 1:ntop, geneIDname][ 
    which( DSignature$direction_foldchange[1:ntop]== -1)]
  
geneSets <-  GSEABase::GeneSetCollection( list( GeneSet( genesUP , 
                                               setName= "UP" ) , 
                                      GeneSet( genesDOWN , 
                                               setName= "DOWN" )) )

### uppercase row. names
exprMatrices <- as.matrix(exprMatrices)

###  1. Build gene-expression rankings for each cell  
cellRanks <- AUCell_buildRankings(  exprMatrices, nCores=1, plotStats=F, 
                                    verbose = progStat )

###  2. Calculate enrichment for the gene signatures (AUC)
# gene sets should include UP and DOWN regulated genes seperately
cells_AUC <- AUCell_calcAUC( geneSets, cellRanks , verbose = progStat,
                             # aucMaxRank parameter is a tricky one,
                             # for single cell 0.05 is a reasonable 
                             # but for bulk a bigger, e.g. 0.2 can be used
                             aucMaxRank = ceiling( ceilThrsh * nrow(cellRanks))
                             )
cells_AUC <- getAUC( cells_AUC) #this will extract results 

### calculate the damage score
# indeces specify elements of geneSets that contain UP and DOWN regulated genes, corespondingly
index1 <- 1
index2 <- 2
# account for possible empty UP or DOWN gene sets
if( isTRUE(wghtd)) {
  coef1 <- sum( genesUP%in%rownames(exprMatrices))/ntop
  coef2 <- sum( genesDOWN%in%rownames(exprMatrices))/ntop
} else coef1 <- coef2 <- 1
  
if( nrow(cells_AUC)==2 ) { 
  Dscore = cells_AUC[index1,] * coef1 - cells_AUC[index2,] * coef2
  } else if(rownames(cells_AUC)=="UP") {
    Dscore= cells_AUC[index1,] } else { Dscore= - cells_AUC[index1,] }

return(Dscore)
}
