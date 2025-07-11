### check PDS enrichment in clinically relevant variants ###
# setwd("/data/user/tpadvits/PROJECTS/PDS/HDS_genomic.variants/")
setwd("/data/user/tpadvits/PROJECTS/PDS/PDS_genomic.variants/")

library(GenomicRanges)
library(ggthemes)
library(BSgenome)
library(biomaRt)
mart  <- useEnsembl("ensembl",
                    dataset="hsapiens_gene_ensembl")

source("~/PROJECTS/scripts/extendGenes.coord_script.r")

# PDS
DS_all <- read.table(header = T, "/data/user/tpadvits/PROJECTS/PDS/DSpodo_musTOhomo_oneTOone.tsv")

#### prepare the data ####
### PDS
DS42.coords <- getBM(
  attributes = c("hgnc_symbol","chromosome_name",
                 "start_position","end_position","strand"),
  filters    = "hgnc_symbol",
  values     = DS_all$Human_Symbol[1:42], mart = mart)

DS42.GR <- GRanges(
  seqnames =  DS42.coords$chromosome_name ,
  ranges = IRanges( start = DS42.coords$start_position, 
                    end = DS42.coords$end_position, names =DS42.coords$hgnc_symbol ),
  strand = DS42.coords$strand,
  gName = (DS42.coords$hgnc_symbol)) 

### extend to include romoters and nearby regions
si <- Seqinfo(genome="hg38")
seqlevelsStyle(si) <- "NCBI"
DS42.GR.plus <- extendGenes( DS42.GR,
                             prm_up        = 2000,
                             upstream_max  = 250000,
                             downstream_max  = 5000 )

### HDS
HDS_all <- read.csv(
  file = "/cellfile/datapublic/pungerav/cell-damage-score/hepatocyte-damage-score/Data/Output/HDAGTop42mappedToHumanGenesManuallyModified",
  sep = '\t')
HDS_all$Human_Symbol <- HDS_all$HumanGeneID

HDS42.coords <- getBM(
  attributes = c("hgnc_symbol","chromosome_name",
                 "start_position","end_position","strand"),
  filters    = "hgnc_symbol",
  values     = HDS_all$Human_Symbol[1:42], mart = mart)

HDS42.GR <- GRanges(
  seqnames =  HDS42.coords$chromosome_name ,
  ranges = IRanges( start = HDS42.coords$start_position, 
                    end = HDS42.coords$end_position, names =HDS42.coords$hgnc_symbol ),
  strand = HDS42.coords$strand,
  gName = (HDS42.coords$hgnc_symbol)) 

### extend to include romoters and nearby regions
si <- Seqinfo(genome="hg38")
seqlevelsStyle(si) <- "NCBI"
HDS42.GR.plus <- extendGenes( HDS42.GR,
                             prm_up        = 2000,
                             upstream_max  = 250000,
                             downstream_max  = 5000 )
# write.table(coords, "genes.bed",
#             sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

#### extract snps from ClinVar db, intersecting with gene ranges #### 
# ClinVar.vars <- read.table(header = F, sep = "\t",
#                            fill = T,
#                            gzfile("/cellfile/datapublic/tpadvits/global_data/ClinVar.2025_variant_summary.txt.gz"))
# hhd <- readxl::read_xlsx("/cellfile/datapublic/tpadvits/global_data/ClinVar.2025_variant_summary.firstRow.xlsx")
# colnames(ClinVar.vars) <- colnames(hhd)
# saveRDS(ClinVar.vars , file="/cellfile/datapublic/tpadvits/global_data/ClinVar.2025_variant_summary.Rda")


### Granges intersect with variant_summary
{
  ClinVar.vars <- readRDS( file="/cellfile/datapublic/tpadvits/global_data/ClinVar.2025_variant_summary.Rda")
  
  ClinVar.vars.inGenes <- ClinVar.vars[ (ClinVar.vars$GeneSymbol %in%  DS_all$Human_Symbol[1:42]) &
                                          ClinVar.vars$ClinicalSignificance %in% c(
                                            "Uncertain significance","Pathogenic/Likely pathogenic",
                                            "Pathogenic"
                                          ),]
  
  ClinVar.vars.GR <- ClinVar.vars[ClinVar.vars$ClinSigSimple==1,]
  ClinVar.vars.GR$Start <-as.numeric(ClinVar.vars.GR$Start)
  ClinVar.vars.GR$Stop <-as.numeric(ClinVar.vars.GR$Stop)
  
  ClinVar.vars.GR <- ClinVar.vars.GR[!is.na( ClinVar.vars.GR$Start) & 
                                       !is.na( ClinVar.vars.GR$Stop),]
  # ClinVar.vars <- read.table(header = T, gzfile("variant_summary.txt.gz"))
  
  ClinVar.vars.GR <- GRanges(
    seqnames =  ClinVar.vars.GR$Chromosome ,
    ranges = IRanges( start =  ClinVar.vars.GR$Start, 
                      end =  ClinVar.vars.GR$Stop, names =ClinVar.vars.GR$VariationID ),
    genomeV=ClinVar.vars.GR$Assembly,
    phenotype=ClinVar.vars.GR$PhenotypeList,
    rsid=ClinVar.vars.GR$`RS# (dbSNP)`,
    Type = ClinVar.vars.GR$Type,
    RCVaccession=ClinVar.vars.GR$RCVaccession
  ) 
  
  # intersect with PDS genes
  fo <- findOverlaps(PDS42.GR.plus, ClinVar.vars.GR)
  II <- pintersect(PDS42.GR.plus[queryHits(fo)], ClinVar.vars.GR[subjectHits(fo)])
  
  II2 <- pintersect(ClinVar.vars.GR[subjectHits(fo)], PDS42.GR.plus[queryHits(fo)])
  II2$gName <- II$gName
  names(II2) <- make.unique(names(II2))
  II2$withinGeneBody <- ifelse( II2$RCVaccession %in% ClinVar.vars.inGenes$RCVaccession, "TRUE", "FALSE")
  # write.table( II2 , sep = "\t", file = "DS.42gene.Extnd_ClinVar.tsv", row.names = F)
  write.table( II2 , sep = "\t", file = "HDS.DS.42gene.Extnd_ClinVar.tsv", row.names = F)
  
}

### download variants for 42 PDS genes (querry terms) from https://www.ncbi.nlm.nih.gov/clinvar
{
  ClinVar.vars_PDS42 <- read.table( header = T, fill = T,  sep = "\t", "/data/user/tpadvits/PROJECTS/PDS/PDS_genomic.variants/DS42_clinvar_result.txt")
  
  ClinVar.vars_PDS42_Kidney <- ClinVar.vars_PDS42[  grepl( 
    "Nephr|Glom|Kidney|Podo",ClinVar.vars_PDS42$Condition.s.,ignore.case = T),]

  # chose (likely) pathogenic variants
  ClinVar.vars_PDS42_Kidney.sel <- ClinVar.vars_PDS42_Kidney[
    ClinVar.vars_PDS42_Kidney$Germline.classification %in% c(
      "Likely pathogenic", "Pathogenic","Pathogenic/Likely pathogenic"
    ),
  ]
  ClinVar.vars_PDS42_Kidney.likely <- intersect( DS_all$Human_Symbol[1:42] , 
                                                 unique(unlist(strsplit(ClinVar.vars_PDS42_Kidney.sel$Gene.s., split="\\|"))))
  # chose uncertain significance variants
  ClinVar.vars_PDS42_Kidney.sel2 <- ClinVar.vars_PDS42_Kidney[
    ClinVar.vars_PDS42_Kidney$Germline.classification %in% c("Uncertain significance") , ]
  
  ClinVar.vars_PDS42_Kidney.uncertain <- intersect( DS_all$Human_Symbol[1:42] , 
             unique(unlist(strsplit(ClinVar.vars_PDS42_Kidney.sel2$Gene.s., split="\\|"))))
  
  }

### download variants for 42 HDS genes (querry terms) from https://www.ncbi.nlm.nih.gov/clinvar
{
  ClinVar.vars_HDS42 <- read.table( header = T, fill = T,  sep = "\t", "/data/user/tpadvits/PROJECTS/PDS/HDS_genomic.variants/HDS.DS42_clinvar_result.txt")
  
  
  liver_key <- c("liver", "hepatic", "cirrhosis", "steatosis",
                 "cholestasis", "bilirubin", "alanine aminotransferase",
                 "aspartate aminotransferase", "gamma glutamyl","bilirubin")
  # no significant relevant hits
  ClinVar.vars_HDS42_Liver <- ClinVar.vars_HDS42[  grepl(paste(liver_key, collapse="|"),
                                                         ClinVar.vars_HDS42$Condition.s., ignore.case = T), ]
  
  
  
  
  }


#### GWAS Catalog via gwasrapidd ####
library(gwasrapidd)
    library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
     library(BSgenome.Hsapiens.NCBI.GRCh38)

### PDS
{
  ### search relevant terms
  ckd.nephr.onc_hits <- get_associations(efo_trait = c( "chronic kidney disease",
                                                        "glomerulonephritis","nephritis",
                                                        "kidney disease","kidney cancer",
                                                        "kidney failure", 
                                                        "acute kidney failure",
                                                        "renal fibrosis", 
                                                        "diabetic nephropathy", 
                                                        "nephrotic syndrome", 
                                                        "minimal change disease", 
                                                        "focal segmental glomerulosclerosis", 
                                                        "membranous glomerulonephritis", 
                                                        "IgA glomerulonephritis", 
                                                        "lupus nephritis", 
                                                        "glomerular filtration rate measurement", 
                                                        "serum creatinine measurement", 
                                                        "urinary albumin-to-creatinine ratio"
  ))
  ### extract genes and/or coordinates of variants
  withinGeneBody <- intersect(ckd.nephr.onc_hits@genes$gene_name ,
                              DS_all$Human_Symbol[1:42])
  
  ckd.nephr.onc_var <- unique( ckd.nephr.onc_hits@risk_alleles$variant_id )
  
  # search stractural variants
  # ckd.nephr.onc_varSV <- ckd.nephr.onc_var[ -grep("rs", ckd.nephr.onc_var) ]
## analyse rs only
  ckd.nephr.onc_varSNP<- grep("rs", ckd.nephr.onc_var, value = T)
  ckd.nephr.onc_varSNP.GR <- snpsById( SNPlocs.Hsapiens.dbSNP155.GRCh38,
                                       ckd.nephr.onc_varSNP ,
                                       ifnotfound="drop", genome=BSgenome.Hsapiens.NCBI.GRCh38)
  
  # vars <- get_variants(variant_id = ckd.nephr.onc_var)
  # 
  # ## The table we need lives in the “variants” slot:
  # vars_tbl <- vars@variants %>%
  #   select(variant_id,
  #          chr  = chromosome_name,
  #          pos  = chromosome_position)      # GRCh38 coordinates
  # 
  # ##  Merge coordinates back onto the association table
  # assoc_with_coords <- ckd.nephr.onc_var@risk_alleles %>%
  #   left_join(vars_tbl, by = "variant_id")
  # # create Granges
  # ckd.nephr.onc_var.GR <-  GRanges(
  #   seqnames =  assoc_with_coords$ ,
  #   ranges = IRanges( start = HDS42.coords$start_position, 
  #                     end = HDS42.coords$end_position, names =HDS42.coords$hgnc_symbol ),
  #   strand = HDS42.coords$strand,
  #   gName = (HDS42.coords$hgnc_symbol)) 
  # 
  # overlaps with PDS 42
  fo <- findOverlaps(DS42.GR.plus, ckd.nephr.onc_varSNP.GR)
  X <- pintersect( DS42.GR.plus[queryHits(fo)], ckd.nephr.onc_varSNP.GR[subjectHits(fo)])
  PDS42.plus.GWAS <- pintersect( GRanges(ckd.nephr.onc_varSNP.GR)[subjectHits(fo)], 
                                 DS42.GR.plus[queryHits(fo)])
  PDS42.plus.GWAS$gName <- X$gName
  PDS42.plus.GWAS$withinGeneBody <- ifelse( PDS42.plus.GWAS$gName %in% withinGeneBody , "TRUE", "FALSE")
  PDS42.plus.GWAS$alt_alleles <- sapply( PDS42.plus.GWAS$alt_alleles, paste, collapse = ",")
  write.table(  as.data.frame( PDS42.plus.GWAS) , sep = "\t", file = "DS.42gene.Extnd_GWAS.tsv", row.names = F)
  
  
  ### test significance of enrichment
  {

    ### test only hits within gene bodies
    {
      genes_of_interest <- DS_all$Human_Symbol[1:42]
      gwas_genes <- unique(ckd.nephr.onc_hits@genes$gene_name)   # catalog gene symbols
      
      
      # overlap counts
      a <- sum(genes_of_interest %in% gwas_genes)              # in G and in GWAS
      b <- sum(!(genes_of_interest %in% gwas_genes))           # in G, not in GWAS
      c <- sum(setdiff(gwas_genes, genes_of_interest) %in% PDSuniverse.h)  # GWAS but not in G
      d <- length(PDSuniverse.h) - a - b - c                        # neither
      
      fisher_mat <- matrix(c(a, b, c, d), nrow = 2,
                           dimnames = list(GWAS = c("yes","no"),
                                           GeneSet = c("yes","no")))
      fisher.test(fisher_mat, alternative = "greater")  # one-sided enrichment
      
      
    }
    
    ### test hits within extended Granges
    {
      #### create extended ranges for all podo genes ####
      {
        ## all sc/snRNA-seq podocyte expressed genes
        allPodoGenes <- readRDS("/data/user/tpadvits/PROJECTS/PDS/SCSN_allPodoGenes.rda")
        library(gprofiler2)
        PDSuniverse.h <- gorth(
          query   = allPodoGenes,
          source_organism = "mmusculus", 
          target_organism = "hsapiens",
          filter_na       = TRUE,
          mthreshold      = Inf,        # return all matches
          numeric_ns      = ""
        )
        PDSuniverse.h <- unique(PDSuniverse.h$ortholog_name)     
        
        # create GR for all podo genes
        allPodoGenes.h.coords <- getBM(
          attributes = c("hgnc_symbol","chromosome_name",
                         "start_position","end_position","strand"),
          filters    = "hgnc_symbol",
          values     = PDSuniverse.h, mart = mart)
        
        allPodoGenes.h.GR <- GRanges(
          seqnames =  allPodoGenes.h.coords$chromosome_name ,
          ranges = IRanges( start = allPodoGenes.h.coords$start_position, 
                            end = allPodoGenes.h.coords$end_position, names =allPodoGenes.h.coords$hgnc_symbol ),
          strand = allPodoGenes.h.coords$strand,
          gName = (allPodoGenes.h.coords$hgnc_symbol)) 
        allPodoGenes.h.GR <- keepStandardChromosomes(allPodoGenes.h.GR,pruning.mode = "coarse")
        allPodoGenes.h.GR <- unique(allPodoGenes.h.GR)
        allPodoGenes.h.GR <- allPodoGenes.h.GR[ !allPodoGenes.h.GR$gName %in% 
                                                  allPodoGenes.h.GR$gName[
                                                    duplicated(allPodoGenes.h.GR$gName)] ]
        ### extend to include romoters and nearby regions
        si <- Seqinfo(genome="hg38")
        seqlevelsStyle(si) <- "NCBI"
        allPodoGenes.h.GR.plus <- extendGenes( allPodoGenes.h.GR,
                                               prm_up        = 2000,
                                               upstream_max  = 250000,
                                               downstream_max  = 5000 )
      }
   
      ### 1. overlaps with podoGenes
      fo <- findOverlaps( allPodoGenes.h.GR.plus, ckd.nephr.onc_varSNP.GR)
      X <- pintersect( allPodoGenes.h.GR.plus[queryHits(fo)], ckd.nephr.onc_varSNP.GR[subjectHits(fo)])
      allPodoGenes.h.GR.plus.GWAS <- pintersect( GRanges(ckd.nephr.onc_varSNP.GR)[subjectHits(fo)], 
                                     allPodoGenes.h.GR.plus[queryHits(fo)])
      allPodoGenes.h.GR.plus.GWAS$gName <- X$gName
      allPodoGenes.h.GR.plus.GWAS$withinGeneBody <- ifelse( allPodoGenes.h.GR.plus.GWAS$gName %in% withinGeneBody , "TRUE", "FALSE")
      allPodoGenes.h.GR.plus.GWAS$alt_alleles <- sapply( allPodoGenes.h.GR.plus.GWAS$alt_alleles, paste, collapse = ",")

      ### 2 Collapse to gene-level counts
      gwas_gene_counts <- as.data.frame( allPodoGenes.h.GR.plus.GWAS) %>%
        count(gName, name = "nHits")  
      
      ### 3 Define your gene set and the background (universe)
      df <- tibble(gene = PDSuniverse.h) %>%
        left_join( gwas_gene_counts, by = c("gene" = "gName")) %>%
        mutate(nHits = replace_na(nHits, 0),
               inSet = gene %in% genes_of_interest)
      
      ### 4 Binary enrichment (≥1 hit per gene)
      # Interpretation: “Do my genes contain any GWAS hit more often than random genes do?”
      tab <- table(inSet = df$inSet, hasHit = df$nHits > 0)
      fisher_res <- fisher.test(tab, alternative = "greater")
      print(tab)
      print(fisher_res$p.value)
     

    }
    
    }
 
}

### HDS
{
  ### search relevant terms
  liver_hits <- get_associations(efo_trait = c( 
    "chronic liver disease",
                                                "liver disease",
                                                "liver enzyme measurement",
                                                "liver steatosis",
                                                "fatty liver disease",
                                                "hepatitis",
                                                "liver fibrosis",
                                                "non-alcoholic steatohepatitis",
                                                "non-alcoholic fatty liver disease", 
                                               "alanine aminotransferase measurement",
                                                 "aspartate aminotransferase measurement",
                                                "serum gamma-glutamyl transferase measurement",
                                                "liver fibrosis measurement", 
                                                "serum albumin change measurement",
                                                "alcoholic liver disease", 
                                                "alcoholic liver cirrhosis",
                                                "hepatocellular carcinoma"
                                                
  ))
  
  
  ### extract genes and/or coordinates of variants
  withinGeneBody <- intersect( liver_hits@genes$gene_name ,
                             HDS_all$Human_Symbol[1:42] )
  
  liver_hits_var <- unique( liver_hits@risk_alleles$variant_id )
  
  # search stractural variants
  liver_varSV <- liver_hits[ -grep("rs", liver_hits_var) ]
  
  liver_varSNP<- grep("rs", liver_hits_var, value = T)
  liver_varSNP.GR <- snpsById( SNPlocs.Hsapiens.dbSNP155.GRCh38, 
                               liver_varSNP ,
                               ifnotfound="drop", 
                               genome=BSgenome.Hsapiens.NCBI.GRCh38 )

    
  # overlaps with PDS 42
  HDS42.GR.plus <- HDS42.GR.plus[ unique(HDS42.GR.plus$gName),]
  fo <- findOverlaps( HDS42.GR.plus, liver_varSNP.GR )
  X <- pintersect(HDS42.GR.plus[queryHits(fo)], liver_varSNP.GR[subjectHits(fo)])
  HDS42.plus.GWAS <- pintersect( GRanges(liver_varSNP.GR)[subjectHits(fo)], HDS42.GR.plus[queryHits(fo)])
  HDS42.plus.GWAS$gName <- X$gName
  withinGeneBody <- intersect( liver_hits@genes$gene_name ,
                               HDS_all$Human_Symbol[1:42] )
  HDS42.plus.GWAS$withinGeneBody <- ifelse( HDS42.plus.GWAS$gName%in% withinGeneBody , "TRUE", "FALSE")
  HDS42.plus.GWAS$alt_alleles <- sapply(HDS42.plus.GWAS$alt_alleles, paste, collapse = ",")
  write.table(  as.data.frame(HDS42.plus.GWAS) , sep = "\t", file = "HDS.DS.42gene.Extnd_GWAS.tsv", row.names = F)
  
  ### test significance of enrichment
  {
    ###  enrichment within gene bodies
    {
      HDSgenes_of_interest <- HDS_all$Human_Symbol[1:42]
      HDSgwas_genes <- unique(liver_hits@genes$gene_name)   # catalog gene symbols
      
      # overlap counts
      a <- sum(HDSgenes_of_interest %in% HDSgwas_genes)              # in G and in GWAS
      b <- sum(!(HDSgenes_of_interest %in% HDSgwas_genes))           # in G, not in GWAS
      c <- sum(setdiff(HDSgwas_genes, HDSgenes_of_interest) %in% HDSuniverse.h)  # GWAS but not in G
      d <- length(HDSuniverse.h) - a - b - c                        # neither
      
      fisher_mat <- matrix(c(a, b, c, d), nrow = 2,
                           dimnames = list(GWAS = c("yes","no"),
                                           GeneSet = c("yes","no")))
      fisher.test(fisher_mat, alternative = "greater")  # one-sided enrichment
    }  
    
    ### enrichment within extended Granges
    {
      #### create extended ranges for all podo genes ####
      {
        
        ## all sc/snRNA-seq podocyte expressed genes
        allHepatoGenes <- read.table(sep = ",",header = T,row.names = 1,
                                     "/data/user/tpadvits/PROJECTS/PDS/GenesExpressedInHepatocytes.csv" )
        allHepatoGenes <- allHepatoGenes$Genes
        HDSuniverse.h <- gorth(
          query   = allHepatoGenes,
          source_organism = "mmusculus", 
          target_organism = "hsapiens",
          filter_na       = TRUE,
          mthreshold      = Inf,        # return all matches
          numeric_ns      = ""
        )
        HDSuniverse.h <- unique(HDSuniverse.h$ortholog_name)     
        
        
        # create GR for all podo genes
        allHepatoGenes.h.coords <- getBM(
          attributes = c("hgnc_symbol","chromosome_name",
                         "start_position","end_position","strand"),
          filters    = "hgnc_symbol",
          values     = HDSuniverse.h, mart = mart)
        
        allHepatoGenes.h.GR <- GRanges(
          seqnames =  allHepatoGenes.h.coords$chromosome_name ,
          ranges = IRanges( start = allHepatoGenes.h.coords$start_position, 
                            end = allHepatoGenes.h.coords$end_position, names =allHepatoGenes.h.coords$hgnc_symbol ),
          strand = allHepatoGenes.h.coords$strand,
          gName = allHepatoGenes.h.coords$hgnc_symbol ) 
        allHepatoGenes.h.GR <- keepStandardChromosomes(allHepatoGenes.h.GR,pruning.mode = "coarse")
        # allHepatoGenes.h.GR <- unique(allHepatoGenes.h.GR)
        # allHepatoGenes.h.GR <- allHepatoGenes.h.GR[ !allHepatoGenes.h.GR$gName %in%
        #                                               allHepatoGenes.h.GR$gName[
        #                                             duplicated(allHepatoGenes.h.GR$gName)] ]
        ### extend to include romoters and nearby regions
        si <- Seqinfo(genome="hg38")
        seqlevelsStyle(si) <- "NCBI"
        allHepatoGenes.h.GR.plus <- extendGenes( allHepatoGenes.h.GR,
                                                 prm_up        = 2000,
                                                 upstream_max  = 250000,
                                                 downstream_max  = 5000 )
      }
      
      
      ### 1. overlaps with liver Genes
      fo <- findOverlaps( allHepatoGenes.h.GR.plus, liver_varSNP.GR )
      X <- pintersect( allHepatoGenes.h.GR.plus[queryHits(fo)], liver_varSNP.GR[subjectHits(fo)])
      allHepatoGenes.h.GR.plus.GWAS <- pintersect( GRanges(liver_varSNP.GR)[subjectHits(fo)], 
                                                   allHepatoGenes.h.GR.plus[queryHits(fo)])
      allHepatoGenes.h.GR.plus.GWAS$gName <- X$gName
      withinGeneBody <- intersect(liver_hits@genes$gene_name ,
                                  HDSgenes_of_interest )
      allHepatoGenes.h.GR.plus.GWAS$withinGeneBody <- ifelse( allHepatoGenes.h.GR.plus.GWAS$gName %in% withinGeneBody , "TRUE", "FALSE")
      allHepatoGenes.h.GR.plus.GWAS$alt_alleles <- sapply( allHepatoGenes.h.GR.plus.GWAS$alt_alleles, paste, collapse = ",")
      
      
      # overlap counts
      a <- sum(HDSgenes_of_interest %in%  HDS42.plus.GWAS$gName )              # in G and in GWAS
      b <- sum(!(HDSgenes_of_interest %in% HDS42.plus.GWAS$gName))           # in G, not in GWAS
      c <- sum(setdiff( allHepatoGenes.h.GR.plus.GWAS$gName, HDSgenes_of_interest) %in% HDSuniverse.h)  # GWAS but not in G
      d <- length(HDSuniverse.h) - a - b - c                        # neither
      
      fisher_mat <- matrix(c(a, b, c, d), nrow = 2,
                           dimnames = list(GWAS = c("yes","no"),
                                           GeneSet = c("yes","no")))
      fisher.test(fisher_mat, alternative = "greater")  # one-sided enrichment
    }
  }
  

    
  }



#### Szustak eGFRcrea GWAS (2.2 million individuals) and Kidney Disease http://www.susztaklab.com/  ####
## load Genome-wide Kidney Disease Genetic Scorecard
eGFR.GWAS <- readxl::read_xlsx("Genome_wide.Kidney.Disease.Genetic.Scorecard.xlsx")
eGFR.GWAS.601Gene <- readxl::read_xlsx(col_names = T,skip = 2, sheet = "S25 601.Genes.Reg&CDS", "science.adp4753_tables_s1_to_s31.xlsx")

# intersect with experimentally val Nephropathy genes
nephr.Genes <- readxl::read_xlsx(col_names = T,skip = 1, sheet = "S22 Neph.Genes", "science.adp4753_tables_s1_to_s31.xlsx")
nephr.Genes$DS42 <- ifelse(nephr.Genes$Gene %in% DS_all$Human_Symbol[1:42],TRUE , FALSE )

# intersect with 601 gene
eGFR.GWAS.601Gene <- readxl::read_xlsx(col_names = T,skip = 2, sheet = "S25 601.Genes.Reg&CDS", "science.adp4753_tables_s1_to_s31.xlsx")

# focus on Glomerular data
eGFR.GWAS.glom <- eGFR.GWAS[ eGFR.GWAS$`ASE of Glomerlus`==1 ,]
eGFR.GWAS.glom$DS42 <- eGFR.GWAS.glom$`Target gene symbol` %in% DS_all$Human_Symbol[1:42]

X <- eGFR.GWAS.glom[eGFR.GWAS.glom$DS42==TRUE,]

# toPlot <- eGFR.GWAS.glom[ (eGFR.GWAS.glom$DS42==TRUE),]
# ggplot2::ggplot(toPlot, aes(x=`Variant Score`, color=DS42))+
#   geom_histogram()+theme_minimal()+ scale_color_colorblind()
# toPlot2 <- eGFR.GWAS.glom[ (eGFR.GWAS.glom$DS42!=TRUE),]
# 
# ggplot2::ggplot(toPlot2, aes(x=`Variant Score`, color=DS42))+
#   geom_histogram()+theme_minimal()+ scale_color_colorblind()

# #### retrieve all variants from dbSNP using gene coordinates ####
# 
#     library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
#      library(BSgenome.Hsapiens.NCBI.GRCh38)
#      
#      genome <- BSgenome.Hsapiens.NCBI.GRCh38
#      all_snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38
#      
#      ## try to retrieve rsID for the variants using biomart
#      my_snps <- snpsByOverlaps(all_snps, PDS42.GR , genome = genome)
#      my_snps.xtnd <- snpsByOverlaps(all_snps, PDS42.GR.plus , genome = genome)
#      # saveRDS(my_snps.xtnd, file = "my_snps.xtnd.rda")
#      
#      ### get variant annotation from myvariant.info
#      my_snps.clinvar <- queryVariants( my_snps$RefSNP_id  , species = "human")
#      
