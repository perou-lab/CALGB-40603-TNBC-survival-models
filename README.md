# CALGB 40603 TNBC survival models
Code and data used to train and evaluate the TNBC survival models in "Prognostic and molecular multi-platform analysis of CALGB 40603 (Alliance) and public triple-negative breast cancer datasets" ([Felsheim et al. NPJ Breast Cancer 2025, PMID: 40057511](https://www.nature.com/articles/s41523-025-00740-z))

## Overview

**Method used for model training:** Cox proportional hazards regression with elastic net regularization. Regularization parameters were selected via sequential bootstrapping over an alpha/lambda grid; the alpha/lambda combination with the lowest average out-of-bag deviance was chosen to fit the final model.

**Types of models:** A total of 7 models were trained, one with each combination of clinical (stage), RNA, and DNA input feature types. Each model was trained to predict overall survival.

**Training set:** CALGB 40603 (n = 238)

**Testing sets:** FUSCC (n = 157), METABRIC (n = 90), and TCGA (n = 133)

## Repository contents

1. **`train_elastic_net.R`:** The base R script used to train and evaluate the elastic net models.     
2. **`model_scripts/`:** Shell scripts specifying the parameters used to run `train_elastic_net.R` for each type of model.
3. **`input_data/`:** The input data needed to train and evaluate the elastic net models.     
    a. **`combined_molecular_feature_data.rds`:** An RDS file containing the processed RNA and DNA input feature values for all samples considered in the CALGB 40603, FUSCC, METABRIC, and TCGA datasets. This contains the exact sample IDs used for training/testing from each dataset in the 'Sample' and 'Set' columns. NOTE: overall survival and stage data are needed to train the models, but cannot be provided here due to data sharing restrictions. These data must be separately obtained from their respective sources (see clinical data availability section below) and added to the input RDS file as 'OS_event' (0 - no event, 1 - event), 'OS_years', and 'Stage' (0 - stage II, 1 - stage III) columns to run `train_elastic_net.R`.   
    b. **`model_predictor_types.tsv`:** A TSV file that maps the names of all features considered by the model to their category (Stage, RNA, or DNA).     
4. **`final_model_coefficients/`:** Contains the scaled and unscaled coefficient values for each trained model.     
5. **`gene_lists/`:** The gene lists used to calculate molecular features considered as features during model training.     
   a. **`RNA_expression_signature_gene_annotations.xlsx`:** Contains the gene lists used to calculate the RNA expression signature features.     
   b. **`copy_number_segment_gene_annotations.rds`:** Contains the gene lists used to calculate the DNA copy number segment features.     

## Clinical data availability

- **CALGB 40603:** clinical data can be requested via the NCBI database of Genotypes and Phenotypes (dbGaP): phs003801.v1.p1
- **FUSCC:** clinical data are available in [Jiang et al. Cancer Cell 2019 (PMID: 30853353)](https://www.sciencedirect.com/science/article/pii/S1535610819300960). Overall survival data was provided by request to the first author.
- **METABRIC:** clinical data are available via cBioPortal: https://www.cbioportal.org/study/summary?id=brca_metabric.     
- **TCGA:** clinical data are available via the NCI Genomic Data Commons (GDC) Data Portal: https://portal.gdc.cancer.gov/projects/TCGA-BRCA

## RNA expression signature features

RNA expression signature features were calculated from log2 transformed, upper-quartile normalized RNA-seq counts for all samples from a given dataset, aligned to hg38 with [GENCODE v36](https://www.gencodegenes.org/human/release_36.html) annotations mapped to Entrez IDs. RNA-seq counts were aligned with STAR v2.7.6a and quantified with Salmon v1.4.0, as detailed [here](https://sc.unc.edu/lbg/workflows/nextflow/rna-seq-quant/-/tree/rna-seq-quant-2.5?ref_type=tags). The normalized count matrix for each set was filtered to include only those with non-zero expression in at least 70% of samples, and was subsequently median-centered and column standardized. From these processed expression matrices, each RNA expression signature was calculated with its corresponding method and gene list as noted in `gene_lists/RNA_expression_signature_gene_annotations.xlsx`. Genes present in the signature's gene list but not the processed expression matrix were ignored. This follows the methods from [Fan et al. BMC Med Genomics 2011, PMID: 21214954](https://bmcmedgenomics.biomedcentral.com/articles/10.1186/1755-8794-4-3). Calculated RNA expression signatures for each dataset were subset to include only the stage II-III TNBC samples listed in `input_data/combined_molecular_feature_data.rds` and were median-centered after subsetting. These RNA expression signature scores were considered as features in the models.

## DNA copy number segment features

DNA copy number segment features were calculated from GISTIC2 gene-level output calculated from ASCAT log2 ratio copy number scores. Gene-level copy number output was collapsed to the level of 534 recurrent copy number segments, listed in `gene_lists/copy_number_segment_gene_annotations.rds`. Segment-level copy number scores were calculated for each segment as the mean GISTIC2 copy number score of the genes listed for that segment in this file. This follows the methods from [Xia et al. Nat Commun 2019, PMID:  31827079](https://www.nature.com/articles/s41467-019-13588-2). Segment-level copy number scores > 0.3 were annotated as a copy number gain, and segment-level copy number scores < −0.3 were annotated as a copy number loss. Both copy number gains and losses were considered as features in the models.

