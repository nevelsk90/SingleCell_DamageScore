
#### generate pre-mRNA indeces to use with Salmon ####
# https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/ 
library(eisaR)




grl <- eisaR::getFeatureRanges(
  gtf = "/cellfile/datapublic/tpadvits/global_data/Genome_annot/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz",
  featureType = c("unspliced"), 
  verbose = TRUE
)


genome <- Biostrings::readDNAStringSet(
  "/cellfile/datapublic/ypaul1/genome_assemblies/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
)


### clean chromosome names
names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)
### make sure you sequences and transcrip coordinates are in harmony

# 1. Get sequence lengths from DNAStringSet
seqlens <- setNames(width(genome), names(genome))  # named vector: seqnames → length
# 2. Keep ranges that are on chromosomes present in the DNAStringSet
gr_valid <- grl[seqnames(grl) %in% names(seqlens)]
# 3. Remove ranges that exceed the sequence length
within_bounds <- start(gr_valid) >= 1 & end(gr_valid) <= seqlens[as.character(seqnames(gr_valid))]
# Filter GRanges
gr_valid <- gr_valid[within_bounds]

# extract pre-mRNA sequences , and write to a fasta file for later indexing with Salmon
seqs <- GenomicFeatures::extractTranscriptSeqs(
  x = genome, 
  transcripts = gr_valid
)

Biostrings::writeXStringSet(
  seqs, filepath = "/cellfile/datapublic/tpadvits/global_data/Genome_annot/mus_musculus/Mus_musculus.GRCm39.113.expanded.fa"
)

eisaR::exportToGtf(
  gr_valid, 
  filepath = "/cellfile/datapublic/tpadvits/global_data/Genome_annot/mus_musculus/Mus_musculus.GRCm39.113.annotation.expanded.gtf"
)


write.table(
  metadata(gr_valid)$corrgene, 
  file = "/cellfile/datapublic/tpadvits/global_data/Genome_annot/mus_musculus/Mus_musculus.GRCm39.113.annotation.expanded.features.tsv",
  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
)


df <- eisaR::getTx2Gene(
  gr_valid, filepath = "/cellfile/datapublic/tpadvits/global_data/Genome_annot/mus_musculus/Mus_musculus.GRCm39.113.annotation.expanded.tx2gene.tsv"
)


#### Isoform Switch analyser ####
library("IsoformSwitchAnalyzeR")

# read annotation
allDesign <-read.table(sep = "\t", header = T,
                       "/cellfile/datapublic/tpadvits/PROJECTS/PodocytePJ/RNA_seq/bulk_NPHS2mut./Sample_NamesAll.csv")


myDesign <- data.frame( sampleID= allDesign$CCG.Sample.ID ,
                        condition =  allDesign$gtype ) 

### Import Salmon example data in R package
salmonQuant <- importIsoformExpression(
  parentDir = "/cellfile/datapublic/tpadvits/PROJECTS/PodocytePJ/RNA_seq/bulk_NPHS2mut./Salmon/",
  addIsofomIdAsColumn = TRUE
)
# clean sample names
colnames(salmonQuant$counts)[ 2:15] <- colnames(salmonQuant$abundance)[ 2:15] <- colnames(salmonQuant$length)[ 2:15] <-
  gsub("_S.*|A006......_","",colnames(salmonQuant$counts)[ 2:15])

### Create switchAnalyzeRlist
GTFannot <- rtracklayer::import("/cellfile/datapublic/ypaul1/genome_assemblies/Mus_musculus.GRCm39.113.gtf.gz")


# ## create sequence file for transcripts
# fasta_file <- "/cellfile/datapublic/ypaul1/genome_assemblies/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
# genome_sequences <- Biostrings::readDNAStringSet(fasta_file)
# txdb <- GenomicFeatures::makeTxDbFromGFF("/cellfile/datapublic/ypaul1/genome_assemblies/Mus_musculus.GRCm39.113.gtf.gz", format = "gtf")
# tx_sequences <- GenomicFeatures::extractTranscriptSeqs(genome_sequences, txdb, use.names=TRUE)
# Biostrings::writeXStringSet( tx_sequences, filepath = "/cellfile/datapublic/tpadvits/global_data/Genome_annot/mus_musculus/Mus_musculus.GRCm39.113.transcriptsDNA.fa.gz", compress=T)

#
SwitchList <- importRdata(
  isoformCountMatrix = salmonQuant$counts[ salmonQuant$counts$isoform_id %in% 
                                            paste( GTFannot$transcript_id , GTFannot$transcript_version, sep = ".") , ] ,
  isoformRepExpression = salmonQuant$abundance[ salmonQuant$abundance$isoform_id %in% 
                                                  paste( GTFannot$transcript_id , GTFannot$transcript_version, sep = ".")  , ] ,
  designMatrix         = myDesign  ,
  ignoreAfterPeriod=T, ignoreAfterSpace=T, 
  isoformExonAnnoation = "/cellfile/datapublic/ypaul1/genome_assemblies/Mus_musculus.GRCm39.113.gtf.gz" ,
  isoformNtFasta       = "/cellfile/datapublic/tpadvits/global_data/Genome_annot/mus_musculus/Mus_musculus.GRCm39.cdna.all.fa.gz" ,
  showProgress = T
)


p01 <-switchPlotTranscript( SwitchList , gene = 'Pkm',
                            condition1="mut",condition2 = "wt",
                            
                            optimizeForCombinedPlot = T) # Visualizes the transcripts and their annotation
p01

#### visualise isoform expr
datt <- SwitchList 

tpm_matrix <- datt$isoformRepExpression

# Also, your design matrix should include sample identifiers and their conditions:
# For example, a design matrix with a column "sample_id" and a column "condition"
design_matrix <- datt$designMatrix
design_matrix$time.hr <-as.factor( allDesign$Time_.hr.[ match(design_matrix$sampleID , allDesign$Run)])

# Choose a gene of interest to visualize (replace 'GeneX' with your gene of interest)
gene_id <- "Pklr"

# Identify transcripts corresponding to the gene
transcripts_of_gene <- datt$isoformFeatures$isoform_id[
  datt$isoformFeatures$gene_name == gene_id  ]

# Subset the TPM matrix for the transcripts of the gene
tpm_gene <- tpm_matrix[ tpm_matrix$isoform_id  %in% transcripts_of_gene, ,drop = FALSE]
rownames(tpm_gene)<- tpm_gene$isoform_id 
tpm_gene[ is.na(tpm_gene)] <- 0

# Transpose the matrix so that rows are samples and columns are transcripts
tpm_df <- as.data.frame(t(tpm_gene[,-1]))
tpm_df$sampleID <- rownames(tpm_df)

# Convert the data frame to long format for ggplot2
tpm_melt <- reshape2::melt(tpm_df, id.vars = "sampleID", 
                           variable.name = "transcript", value.name = "TPM")

# Merge the melted expression data with the design matrix to incorporate condition info
# (adjust the column names in design_matrix if necessary)
tpm_melt <- merge(tpm_melt, design_matrix, by.x = "sampleID", by.y = "sampleID")

# Plot transcript expression across conditions
p02 <-  ggplot(tpm_melt, aes(x = transcript , y = TPM, fill = condition)) +
  geom_boxplot(outlier.size = 0) +   
  geom_point(position=position_jitterdodge(jitter.width = 0.1), 
             aes(fill=condition), size=2)+
  labs(title = paste("Transcript Expression for", gene_id),
       x = "Transcript",
       y = "TPM") +
  theme_minimal()  +   scale_shape_manual(values = c(21, 24)) +
  # facet_grid(rows   = vars(study), scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.1)) 





## Testing for Isoform Switches via isoformSwitchTestSatuRn
iso_filt <- preFilter(
  switchAnalyzeRlist = SwitchList,
  removeSingleIsoformGenes = TRUE,
  # e.g. keep isoforms with >= 10 counts in >= 3 samples:
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 3
)

# Perform test
exampleSwitchListAnalyzed <- isoformSwitchTestSatuRn(
  switchAnalyzeRlist = iso_filt,
  reduceToSwitchingGenes=F
)
p03 <- switchPlot( exampleSwitchListAnalyzed , gene = 'Pkm', 
            condition1= "mut", condition2 = "wt" )

pdf(width = 8 , height = 8, file="/data/user/tpadvits/PROJECTS/PDS/Manuscript_submission/Nphs2mut.bulkRNAseq_Pklr.expr.TPM.boxplots.pdf" )
p02
dev.off()

pdf(width = 8 , height = 8, file="/data/user/tpadvits/PROJECTS/PDS/Manuscript_submission/Nphs2mut.bulkRNAseq_Pklr.switchPlot.pdf" )
switchPlot( exampleSwitchListAnalyzed , gene = 'Pklr', 
            condition1= "mut", condition2 = "wt" )
dev.off()


# # Summarize switching features
# extractSwitchSummary(exampleSwitchListAnalyzed)
# 
# 
# ### Predicting Alternative Splicing
# exampleSwitchListAnalyzed <- analyzeAlternativeSplicing(
#   switchAnalyzeRlist = exampleSwitchListAnalyzed,onlySwitchingGenes = F,
#   quiet=TRUE
# )
# 
# # overview of number of intron retentions (IR)
# table( exampleSwitchListAnalyzed$AlternativeSplicingAnalysis$IR )
# 
# 
#### visualise
# switchPlotTranscript(exampleSwitchListFiltered, gene = 'Parp11') # Visualizes the transcripts and their annotation
# switchPlotGeneExp(exampleSwitchListFiltered, gene = 'Parp11')    # Visualizes the gene expression
# switchPlotIsoExp(exampleSwitchListFiltered, gene = 'Parp11')     # Visualizes the isoform expression
# switchPlotIsoUsage(exampleSwitchListFiltered, gene = 'Parp11')   # Visualizes the isoform usage


