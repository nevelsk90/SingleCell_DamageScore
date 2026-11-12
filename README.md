# SingleCell_DamageScore

Code and gene lists for the Single-cell Damage Score manuscript (podocyte damage score, PDS).

## Repository layout

| Folder / file | Content |
|---|---|
| `gene_sets_and_tables/` | Damage signature, podocyte gene list and derived tables (e.g. the signature translated to human orthologs) |
| `preprocessing/` | Preprocessing of public datasets, one script per dataset |
| `Analysis_scripts/` | Analyses using the damage score |
| `Figure_scripts/` | Figures of the manuscript, made from preprocessed data |
| `archive/` | Older scripts kept for reference; some are still sourced by the analysis scripts (e.g. `archive/scripts_PDS/AUCell_script.r`) |

## gene_sets_and_tables

* `gene_sets_and_tables/podocyteDamageSignature.tsv` - podocyte damage signature consisting of 381 markers. Columns **mean_rank** and **direction_foldchange** show the average rank of differential expression (based on p-value) and the consistent direction of change in experimental compared to control samples across the tested studies.
* `gene_sets_and_tables/SCSN_allPodoGenes.rda` - R file with a list of genes expressed in healthy and damaged podocytes in single-cell RNAseq datasets.
* `gene_sets_and_tables/DSpodo_musTOhomo_oneTOone.tsv` - the podocyte damage signature translated from mouse to human genes, using one-to-one orthologs only.

## preprocessing

Scripts that turn raw public data into the objects used by the analysis scripts. Each script is named after the dataset accession.

* `GSE146912.R` - preprocessing of the public mouse glomerular scRNAseq dataset GEO GSE146912, used for the damage score validation.

## Analysis_scripts

* `code_KNNfilt` - KNN filtration code (Python) that filters sporadic cells and identifies circular clusters of glomerular cells in spatial transcriptomics data. Modified from <https://github.com/marshalljamie/Kidney-Slide-seq>.
* `snRNAseq_analysis_Nphs2mut.R` - analysis of the Nphs2 mut. single-nuclei RNAseq data. The same code was used for the Wt1het.del and Pdss2 mut. snRNAseq datasets generated for the study.
* `xenium.human_PDSanalysis.r` - analysis of human kidney Xenium 5K spatial transcriptomics data.
* `Supervised.vs.Unsupervised_DamageSignature.r` - comparison of the supervised and unsupervised (trajectory inference) damage signatures.
* `nphs2mut_isoformAnalysis.r` - isoform-level expression in murine glomeruli.
* `DiseaseScore_TrlvlExpr.r` - expression of a gene of interest in podocytes across the 3 murine snRNAseq datasets generated for the manuscript and the 4 public scRNAseq datasets used in it.
* `DiseaseScore_tempo.res.prep.r` - processing of the results of the Tempo algorithm, which estimates the circadian time of individual cells.
* `DiseaseScore_spatial.r` - PDS analysis of mouse and human spatial transcriptomics datasets.
* `DiseaseScore_Proteomics.r` - PDS analysis of proteomics data.
* `DiseaseScore_PDS.DataPrep.r` - combines the 3 snRNAseq datasets generated for the manuscript and the 4 public scRNAseq datasets of murine podocyte damage into one object for PDS analysis and visualisation.
* `DiseaseScore_PDS.application.r` - uses PDS to explore common and model-specific changes in the activity of biological molecules and processes associated with cellular damage. Results are interpreted with a podocyte-specific gene regulatory network, built from podocyte ATACseq and TF motifs, to characterise transcription factors involved in the damage.
* `DiseaseScore_PAanalysis.r` - pathway activity calculation and visualisation.
* `DiseaseScore_Human.r` - PDS analysis of human single-cell RNAseq data.
* `DiseaseScore_GWAS.analysis.r` - annotation of the podocyte and liver damage markers with GWAS results.
* `DiseaseScore_Dev_SigSize.test.r` - damage signature size test.
* `DiseaseScore_Dev_PDS.aligment.test.r` - tests whether PDS distributions are statistically more similar in control than in experimental samples.
* `DiseaseScore_Dev_CrossValidation.r` - cross-validation of the damage score across podocyte damage models.
* `DiseaseScore_Dev.r` - analysis of public datasets, generation of the damage signature, and comparison of approaches to damage score calculation.
* `DiseaseScore_circad.analyis.r` - analysis of the bulk RNAseq circadian dataset, and of circadian (dys)regulation in podocyte damage with CRD and Tempo.

## Figure_scripts

Scripts that make the main figures (`figures_scripts_fig1-4.r`) and supplementary figures (`Supplimentary_figures/`) from preprocessed data.

## Data and paths

The scripts read input data from, and write results to, absolute paths on the Beyer lab file server (`/cellfile/DATASETS/...` for data, `/cellfile/PROJECTS/...` for results). Results are not stored in this repository. To run the code elsewhere, download the public datasets by their accessions and adjust the paths at the top of each script.
