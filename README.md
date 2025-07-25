# SingleCell_DamageScore

Repository that contains files for the Single-cell Damage Score manuscript, PDS analysis.

*  podocyteDamageSignature.tsv - podocyte damage signature consisting of 381 markers, columns **mean_rank** and **direction_foldchange** shows average rank of differential expression (based on p-value) and consistent direction of change in experimental compared to control samples, across the tested studies.

*  SCSN_allPodoGenes.rda - an R file containing a list of genes expressed in healthy and damaged podocytes in single cell RNAseq datasets

### Figure_scripts

folder that contains scripts for making figures from preprocessed data.

### Analysis_scripts

folder that contains scripts for analysing data using the damage score. 

* code_KNNfilt - a folder with KNN filtration code (python) to filter sporadic cells and identify circular clusters of glomerular cells in spatial transcriptomics data. Modified from <https://github.com/marshalljamie/Kidney-Slide-seq> 
* snRNAseq_analysis_Nphs2mut.R - analysis of Nphs2 mut. single-nuclei RNAseq data. Same code was used to analyse Wt1het.del and Pdss2 mut snRNAseq datasets generated for the study. 
* xenium.human_PDSanalysis.r - analysis of human kidney Xenium 5K spatial transcriptomics data 
* Supervised.vs.Unsupervised_DamageSignature.r - comparison of the supervised and unsupervised (trajectory inference) damage signatures 
* nphs2mut_isoformAnalysis.r - script to analyse isoform lvl expression in murine glomeruli 
* DiseaseScore_TrlvlExpr.r - visualise expression of a gene of interest in podocytes across 3 murine snRNAseq datasets, generated for the manuscript, and 4 public scRNAseq datasets used in the manuscript. 
* DiseaseScore_tempo.res.prep.r - processing results of running tempo algotythm that estimates circadian time of individual cells 
* DiseaseScore_spatial.r - PDS analysis in mouse and human spatial transcriptomics datasets 
* DiseaseScore_Proteomics.r - PDS analysis of proteomics data 
* DiseaseScore_PDS.DataPrep.r - combine 3 snRNAseq datasets, generated for the manuscript, and 4 public scRNAseq datasets of murine podocyte damage in one object, that can be used for PDS analysis and visualisation. 
* DiseaseScore_PDS.application.r - use PDS to explore common and model-specific changes in activity of biological molecules and processes associated with the cellular damage. Interpret results of PDS analysis with podocyte-specifc gene regulatory network, constructed using podocyte ATACseq and TF motifs, to charachterise transcription factors involved in the damage. 
* DiseaseScore_PAanalysis.r - pathway activity calculation and visualisation
* DiseaseScore_Human.r - PDS analysis of human single-cell RNAseq data 
* DiseaseScore_GWAS.analysis.r - annotation of the podocyte and the liver damage markers with results of GWAS studies 
* DiseaseScore_Dev_SigSize.test.r - damage signature size test 
* DiseaseScore_Dev_PDS.aligment.test.r - testing if PDS distributions are statistically more simmilar in control compared to experimental samples. 
* DiseaseScore_Dev_CrossValidation.r - cross-validation of the damage score across varios podocyte damage models 
* DiseaseScore_Dev.r - analysing public datasets, generating damage signature, testing varios approaches to damage score calculation. 
* DiseaseScore_circad.analyis.r - analysis of bulk RNAseq circadian dataset. Analysis of circadian (dys)regulation in podocyte damage with CRD and tempo.
