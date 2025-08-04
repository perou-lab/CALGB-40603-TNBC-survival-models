#!/bin/bash
Rscript train_elastic_net.R \
	--featureDataFile input_data/combined_molecular_clinical_feature_data.rds \
	--predictorTypeFile input_data/model_predictor_types.tsv \
	--includeTypes DNA,RNA \
	--nBootstraps 100 \
	--alphas 0.1,0.9,0.1 \
	--lambdas 2,-2,20 \
	--seed 822 \
	--outputDir output/RNA_DNA
