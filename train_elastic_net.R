#' train_elastic_net.R
#'  
#' Bootstrapping workflow to train elastic net model on clinical and/or 
#' molecular features using CALGB 40603 training set, with evaluation on FUSCC,
#' METABRIC, and TCGA test sets
#' 
#' Prognostic overall survival model of stage II+III TNBC

library(glmnet)
library(tidyverse)
library(sampling)
library(survival)
library(survminer)
library(maftools)
library(data.table)
library(gridExtra)
library(matrixStats)
library(Hmisc)
library(caret)
library(optparse)
library(foreach)
library(doParallel)
library(parallelly)

# Input parameters -------------------------------------------------------------

option_list <- list(
	make_option("--featureDataFile", type="character", default=NULL,
							help="Feature data RDS file path, the feature data is a data frame
							with a 'Sample' column, overall survival columns ('OS_years', 
							'OS_events'), and feature columns present in predictor types file",
							metavar="character"),
	make_option("--predictorTypeFile", type="character", default=NULL,
							help="Predictor types tsv file path, containing 'name' column
							(required - feature name) and 'type' column (optional - data type)",
							metavar="character"),
	make_option("--includeTypes", type="character", default=NULL,
							help="Data types to include (must be a category in the 'type'
							column in the predictor type file. Default is NULL (all features
							in predictor type file are kept in model). To include multiple
							data types, separate with a comma (e.g. 'Stage,RNA')",
							metavar="character"),
	make_option("--nBootstraps", type="numeric", default=100,
							help="Number of bootstraps to use",
							metavar="character"),
	make_option("--alphas", type="character", default="0.1,0.9,0.1",
							help="Alpha values to use, will generate based on 
							seq(from = from, to = to, by = by). Specify 'from,to,by' with 
							comma separation. Default is 0.1,0.9,0.1",
							metavar="character"),
	make_option("--lambdas", type="character", default="2,-2,20",
							help="Alpha values to use, will generate based on 
							10^(seq(from = from, to = to, length = length)). Specify 
							'to,from,length' with comma separation. Default is 2,-2,20",
							metavar="character"),
	make_option("--seed", type="numeric", default=123,
							help="Seed to use for reproducibility",
							metavar="character"),
	make_option("--outputDir", type="character", default=NULL,
							help="Path to output directory",
							metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

feature_data_file <- opt$featureDataFile
predictor_type_file <- opt$predictorTypeFile
include_types <- unlist(str_split(opt$includeTypes, ","))
n_bootstraps <- opt$nBootstraps
alpha_params <- as.numeric(unlist(str_split(opt$alphas, ",")))
alphas <- seq(alpha_params[1], alpha_params[2], by = alpha_params[3])
lambda_params <- as.numeric(unlist(str_split(opt$lambdas, ",")))
lambdas <- 10^(seq(lambda_params[1], lambda_params[2], length = lambda_params[3]))
seed <- opt$seed
out_dir <- opt$outputDir

med_colors <- c("Low" = "#000000", 
								"High" = "#b34e3f")

tert_colors <- c("Low" = "#0172b1",
								 "Med" = "#000000", 
								 "High" = "#b34e3f")

set_colors <- c("CALGB40603" = "#2baad3",
								"FUSCC" = "#ca6551",
								"METABRIC" = "#81a24c",
								"TCGA" = "#a466ba",
								"Combined Test" = "#555555")

test_colors <- c("FUSCC" = "#ca6551",
								 "METABRIC" = "#81a24c",
								 "TCGA" = "#a466ba",
								 "Combined" = "#555555")

dir_colors <- c("Negative" = "#0172b1", 
								"Positive" = "#b34e3f")

# Prepare data -----------------------------------------------------------------

set.seed(seed)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

predictor_types <- read_tsv(predictor_type_file) 

if (!is.null(include_types)) {
	predictor_types <- predictor_types %>% 
		filter(type %in% include_types)
}

metadata <- readRDS(feature_data_file) 

predictors <- metadata[,predictor_types$name]

# Initial filtering of predictors ----------------------------------------------

x <- as.matrix(predictors) # predictors
y <- Surv(metadata$OS_years, metadata$OS_event) # predict OS

c03_idx <- as.numeric(row.names(metadata)[metadata$Set == "CALGB40603"])
fuscc_idx <- as.numeric(row.names(metadata)[metadata$Set == "FUSCC"])
mb_idx <- as.numeric(row.names(metadata)[metadata$Set == "METABRIC"])
tcga_idx <- as.numeric(row.names(metadata)[metadata$Set == "TCGA"])
test_idx <- c(fuscc_idx, mb_idx, tcga_idx)

c03_metadata <- metadata[c03_idx,]
fuscc_metadata <- metadata[fuscc_idx,]
mb_metadata <- metadata[mb_idx,]
tcga_metadata <- metadata[tcga_idx,]
test_metadata <- metadata[test_idx,]

# Training set

c03_x <- x[c03_idx,]
c03_y <- y[c03_idx]

# Remove sparse features

keep_features <- which(colSums(c03_x == 0)/nrow(c03_x) < 0.95)
x <- x[,keep_features]
c03_x <- c03_x[,keep_features]

# Keep only one feature among highly correlated segments

M <- cor(c03_x)
diag(M) <- 0

M_high <- M == 1

features <- row.names(M_high)[rowSums(M_high) == 0]
correlated <- row.names(M_high)[rowSums(M_high) > 0]

sink(paste0(out_dir, "/correlation_feature_selection.txt"))
while(length(correlated) > 0) {
	columns <- c(correlated[1], colnames(M_high)[M_high[correlated[1],]])
	print("Correlated:")
	print(columns)
	values <- c03_x[,columns]
	if (sum(grepl("wholearm", columns)) > 0) {
		keep <- columns[grepl("wholearm", columns)]
	} else {
		keep <- colnames(values)[which.max(colSds(values))]
	}
	print("Keep:")
	print(keep)
	features <- c(features, keep)
	correlated <- correlated[!(correlated %in% columns)]
}
sink()

# Features to keep after filtering

keep_features <- unique(features)

x <- x[,keep_features]
c03_x <- c03_x[,keep_features]

fuscc_x <- x[fuscc_idx,]
fuscc_y <- y[fuscc_idx]
mb_x <- x[mb_idx,]
mb_y <- y[mb_idx]
tcga_x <- x[tcga_idx,]
tcga_y <- y[tcga_idx]
test_x <- x[test_idx,]
test_y <- y[test_idx]

input_features <- predictor_types %>% 
	filter(name %in% colnames(x))
write_tsv(input_features, 
					paste0(out_dir, "/input_features.txt"))

sink(paste0(out_dir, "/input_feature_counts.txt"))
print("Total input features:")
print(nrow(input_features))
if (!is.null(include_types)) {
	print("Total features by type:")
	print(table(input_features$type))
}
sink()

# Train glmnet model using bootstrapping approach ------------------------------

bootstrapped_samples <- lapply(1:n_bootstraps, function(i) sample(1:length(c03_y), replace = T))
oob_samples <- sapply(bootstrapped_samples, function(x) setdiff(1:length(c03_y), x))

bootstrap_table <- expand.grid(boot_iter = 1:n_bootstraps,
															 alpha = alphas, 
															 lambda = lambdas,
															 train_dev = NA,
															 train_C = NA,
															 oob_dev = NA,
															 oob_C = NA,
															 n_coef = NA,
															 model_coefs = NA)

print(paste0("# bootstraps: ", n_bootstraps, 
						 ", # alphas: ", length(alphas), 
						 ", # lambdas: ", length(lambdas)))

cores <- availableCores(omit = 2)
cl <- makeCluster(cores) 
registerDoParallel(cl)

print(paste0("Running in parallel with ", cores, " cores"))

bootstrap_table <- foreach(i = 1:nrow(bootstrap_table),
													 .packages=c('glmnet', 'stats'),
													 .combine = rbind) %dopar% {
													 	
													 	boot_iter <- bootstrap_table[i, "boot_iter"]
													 	alpha <- bootstrap_table[i, "alpha"]
													 	lambda <- bootstrap_table[i, "lambda"]
													 	
													 	boot_c03_x <- c03_x[bootstrapped_samples[[boot_iter]],]
													 	boot_c03_y <- c03_y[bootstrapped_samples[[boot_iter]]]
													 	
													 	oob_c03_x <- c03_x[oob_samples[[boot_iter]],]
													 	oob_c03_y <- c03_y[oob_samples[[boot_iter]]]
													 	
													 	# Train models on bootstrapped data
													 	model_obj <- glmnet(boot_c03_x, boot_c03_y, family = "cox",
													 											alpha = alpha, lambda = lambda)
													 	
													 	# Apply bootstrapped models to training/oob data
													 	pred_train <- predict(model_obj, 
													 												newx = boot_c03_x)
													 	pred_oob <- predict(model_obj, 
													 											newx = oob_c03_x)
													 	
													 	# Calculate C index on training/oob data
													 	bootstrap_table[i, "train_C"] <- Cindex(pred = pred_train, y = boot_c03_y)
													 	bootstrap_table[i, "oob_C"] <- Cindex(pred = pred_oob, y = oob_c03_y)
													 	
													 	# Calculate deviance on training/oob data
													 	bootstrap_table[i, "train_dev"] <- coxnet.deviance(pred = pred_train, y = boot_c03_y)
													 	bootstrap_table[i, "oob_dev"] <- coxnet.deviance(pred = pred_oob, y = oob_c03_y)
													 	
													 	# Save coefficient info
													 	model_coefs <- coef(model_obj)[,1]
													 	model_coefs <- model_coefs[model_coefs != 0]
													 	bootstrap_table[i, "n_coef"] <- length(model_coefs)
													 	bootstrap_table[i, "model_coefs"] <- paste(names(model_coefs), collapse = ",")
													 	
													 	bootstrap_table[i,]
													 	
													 }

write_tsv(bootstrap_table, paste0(out_dir, "/bootstrap_table.txt"))

# Train a final model based on the alpha/lambda combination that had best
# average performance across bootstrapping parameters

tuning_results <- bootstrap_table %>% 
	group_by(alpha, lambda) %>% 
	dplyr::summarize(avg_train_dev = mean(train_dev),
									 train_dev_CI_lo = quantile(train_dev, 0.025),
									 train_dev_CI_hi = quantile(train_dev, 0.975),
									 avg_oob_dev = mean(oob_dev),
									 oob_dev_CI_lo = quantile(oob_dev, 0.025),
									 oob_dev_CI_hi = quantile(oob_dev, 0.975),
									 avg_optimism_dev = mean(train_dev - oob_dev),
									 avg_train_C = mean(train_C),
									 train_C_CI_lo = quantile(train_C, 0.025),
									 train_C_CI_hi = quantile(train_C, 0.975),
									 avg_oob_C = mean(oob_C),
									 oob_C_CI_lo = quantile(oob_C, 0.025),
									 oob_C_CI_hi = quantile(oob_C, 0.975),
									 avg_optimism_C = mean(train_C - oob_C),
									 avg_n_coef = mean(n_coef),
									 min_n_coef = min(n_coef),
									 max_n_coef = max(n_coef))

write_tsv(tuning_results, 
					paste0(out_dir, "/tuning_grid_bootstrap_summary_metrics.txt"))

# Optimal alpha, lambda as combination with lowest OOB deviance

best_params <- tuning_results %>%
	dplyr::arrange(avg_oob_dev) %>% 
	head(n = 1)

best_alpha <- best_params$alpha
best_lambda <- best_params$lambda
best_table <- bootstrap_table %>% 
	filter(alpha == best_alpha) %>% 
	filter(lambda == best_lambda)
best_avg_oob_dev <- best_params$avg_oob_dev
best_avg_oob_C <- best_params$avg_oob_C

write_tsv(best_table, 
					paste0(out_dir, "/training_bootstrap_table_for_optimal_alpha_lambda.txt"))

final_model <- glmnet(c03_x, c03_y, family = "cox",
											alpha = best_alpha, lambda = best_lambda)

saveRDS(final_model, paste0(out_dir, "/final_model.RDS"))

# Model characterization / visualization ---------------------------------------

# Model coefficients

unscaled_coef <- coef(final_model)[,1]
unscaled_coef <- unscaled_coef[unscaled_coef != 0]
scaled_coef <- coef(final_model)[,1] *
	apply(c03_x, 2, function(x) sqrt(sum((x - mean(x)) ^ 2) / length(c03_y)))
scaled_coef <- scaled_coef[scaled_coef != 0 & !is.na(scaled_coef)]
model_coef <- data.frame(name = names(unscaled_coef), unscaled_coef = unscaled_coef,
												 scaled_coef = scaled_coef)

write_tsv(model_coef, paste0(out_dir, "/final_model_coefficients.txt"))

# Plot model coefficients

coefficients <- model_coef %>% 
	dplyr::arrange(scaled_coef) %>% 
	mutate(name = factor(name, levels = name))

coefficients$direction <- factor(ifelse(coefficients$scaled_coef > 0, "Positive", 
																				"Negative"), levels = c("Negative", "Positive"))

plot <- ggplot(coefficients, aes(x = name, y = scaled_coef, fill = direction)) +
	geom_bar(stat = "identity") +
	coord_flip() +
	guides(fill = "none") +
	theme_minimal(base_size = 18) +
	ylab("Scaled coefficient value") +
	xlab(NULL) +
	scale_fill_manual(values = dir_colors)

pdf(paste0(out_dir, "/scaled_coefficient_value_plot_dir_colored.pdf"),
		width = 12, height = 1 + nrow(coefficients)*0.5)
print(plot)
dev.off()

# Final model evaluation on test data ------------------------------------------

pred_c03 <- predict(final_model, newx = c03_x)
pred_fuscc <- predict(final_model, newx = fuscc_x)
pred_mb <- predict(final_model, newx = mb_x)
pred_tcga <- predict(final_model, newx = tcga_x)
pred_test <- predict(final_model, newx = test_x)

# Add prediction values to metadata

metadata$prediction <- NA
metadata$prediction[c03_idx] <- unlist(pred_c03)
metadata$prediction[fuscc_idx] <- unlist(pred_fuscc)
metadata$prediction[mb_idx] <- unlist(pred_mb)
metadata$prediction[tcga_idx] <- unlist(pred_tcga)

# Final model evaluation metrics (C-index, deviance)

evaluation_metrics <- data.frame(Set = c("CALGB 40603 (Training)", "FUSCC (Test)", 
																				 "METABRIC (Test)", "TCGA (Test)",
																				 "Combined Test"),
																 C_index = c(Cindex(pred = pred_c03, y = c03_y), 
																 						Cindex(pred = pred_fuscc, y = fuscc_y),
																 						Cindex(pred = pred_mb, y = mb_y),
																 						Cindex(pred = pred_tcga, y = tcga_y),
																 						Cindex(pred = pred_test, y = test_y)),
																 Deviance = c(coxnet.deviance(pred = pred_c03, y = c03_y), 
																 						 coxnet.deviance(pred = pred_fuscc, y = fuscc_y),
																 						 coxnet.deviance(pred = pred_mb, y = mb_y),
																 						 coxnet.deviance(pred = pred_tcga, y = tcga_y),
																 						 coxnet.deviance(pred = pred_test, y = test_y)))

write_tsv(evaluation_metrics, 
					paste0(out_dir, "/final_model_evaluation_metrics.txt"))

eval_metrics <- evaluation_metrics %>% 
	mutate(C_index = round(C_index, digits = 2)) %>% 
	rename(`C-index` = C_index) %>% 
	mutate(Set = sub(" \\(.*", "", Set)) %>% 
	mutate(Set = sub("B ", "B", Set)) %>% 
	mutate(Set = factor(Set, levels = c("CALGB40603", "FUSCC", "METABRIC", 
																			"TCGA", "Combined Test"))) %>% 
	select(Set, `C-index`)

# Final model C-index barplots

pdf(paste0(out_dir, "/final_model_C_index_barplots.pdf"),
		width = 5, height = 5)
print(eval_metrics %>% 
				ggplot(aes(x = Set, y = `C-index`, fill = Set)) +
				geom_bar(stat = 'identity') +
				ylim(c(0,1)) +
				scale_fill_manual(values = set_colors) + 
				theme_minimal(base_size = 14) + 
				geom_text(aes(label = `C-index`), vjust = -0.2, size = 6) +
				theme(legend.position = "none") 
) 
dev.off()

pdf(paste0(out_dir, "/final_model_C_index_barplots_test_only.pdf"),
		width = 5, height = 5)
print(eval_metrics %>% 
				filter(Set != "CALGB40603") %>% 
				mutate(Set = sub(" Test", "", Set)) %>% 
				mutate(Set = factor(Set, levels = c("FUSCC", "METABRIC", 
																						"TCGA", "Combined"))) %>% 
				ggplot(aes(x = Set, y = `C-index`, fill = Set)) +
				geom_bar(stat = 'identity') +
				ylim(c(0,1)) +
				scale_fill_manual(values = test_colors) + 
				theme_minimal(base_size = 18) + 
				geom_text(aes(label = `C-index`), vjust = -0.2, size = 6) +
				theme(legend.position = "none") 
) 
dev.off()

# Likelihood ratio test - stage + model risk

sink(paste0(out_dir, "/final_model_LRT_output.txt"))

print("CALGB 40603 -------------------")
print(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ Stage"),
						data = metadata[c03_idx,]))
print(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ prediction"),
						data = metadata[c03_idx,]))
print(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ Stage + prediction"),
						data = metadata[c03_idx,]))

print("FUSCC -------------------------")
print(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ Stage"),
						data = metadata[fuscc_idx,]))
print(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ prediction"),
						data = metadata[fuscc_idx,]))
print(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ Stage + prediction"),
						data = metadata[fuscc_idx,]))

print("METABRIC ----------------------")
print(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ Stage"),
						data = metadata[mb_idx,]))
print(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ prediction"),
						data = metadata[mb_idx,]))
print(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ Stage + prediction"),
						data = metadata[mb_idx,]))

print("TCGA --------------------------")
print(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ Stage"),
						data = metadata[tcga_idx,]))
print(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ prediction"),
						data = metadata[tcga_idx,]))
print(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ Stage + prediction"),
						data = metadata[tcga_idx,]))

print("Combined test -----------------")
print(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ Stage"),
						data = metadata[test_idx,]))
print(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ prediction"),
						data = metadata[test_idx,]))
print(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ Stage + prediction"),
						data = metadata[test_idx,]))
sink()

test_stage <- summary(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ Stage"),
														data = metadata[test_idx,]))
stage_lr <- test_stage$logtest[1]
stage_p <- test_stage$coefficients[5]
stage_pc <- ifelse(stage_p < 0.001, "p < 0.001",
									 paste0("p = ", signif(stage_p, 2)))

test_risk <- summary(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ prediction"),
													 data = metadata[test_idx,]))
risk_lr <- test_risk$logtest[1]
risk_p <- test_risk$coefficients[5]
risk_pc <- ifelse(risk_p < 0.001, "p < 0.001",
									paste0("p = ", signif(risk_p, 2)))

test_both <- summary(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ Stage + prediction"),
													 data = metadata[test_idx,]))
both_lr <- test_both$logtest[1]
both_stage_p <- test_both$coefficients[1,5]
both_risk_p <- test_both$coefficients[2,5]
both_stage_pc <- ifelse(both_stage_p < 0.001, "p < 0.001",
												paste0("p = ", signif(both_stage_p, 2)))
both_risk_pc <- ifelse(both_risk_p < 0.001, "p < 0.001",
											 paste0("p = ", signif(both_risk_p, 2)))

lrt_data <- data.frame(LR = c(stage_lr, both_lr-stage_lr, risk_lr, both_lr-risk_lr),
											 p = c(stage_pc, both_risk_pc, risk_pc, both_stage_pc),
											 type = c("Stage", "Model risk", "Model risk", "Stage"),
											 group = c("Order 1", "Order 1", "Order 2", "Order 2"),
											 ordering = c(1,2,1,2)) %>%
	group_by(group) %>%
	mutate(cumulative_LR = cumsum(LR),
				 mid_LR = cumulative_LR - 0.5 * LR)

pdf(paste0(out_dir, "/LRT_plot_stage_model_risk_test.pdf"),
		width = 5, height = 5)
lrt_data %>% 
	ggplot(aes(x = group, y = LR, fill = reorder(type,-ordering))) + 
	geom_bar(stat = "identity", position = "stack", data = lrt_data %>% filter(group == "Order 1")) +
	geom_bar(stat = "identity", position = "stack", data = lrt_data %>% filter(group == "Order 2")) +
	geom_text(aes(y = mid_LR, label = p), 
						color = "black", 
						size = 5) +
	theme_minimal(base_size = 18) +
	ylab("Order of addition to model\n\nChange in LR statistic") +
	xlab("") +
	theme(legend.title=element_blank()) +
	scale_fill_manual(values = c("Stage" = "#edb2d5", "Model risk" = "#34b399")) +
	theme(legend.position="top")
dev.off()

# Likeihood ratio test - stage + model risk (stratified by set)

sink(paste0(out_dir, "/final_model_LRT_output_set_strata.txt"))

print("Combined test -----------------")
print(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ Stage + strata(Set)"),
						data = metadata[test_idx,]))
print(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ prediction + strata(Set)"),
						data = metadata[test_idx,]))
print(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ Stage + prediction + strata(Set)"),
						data = metadata[test_idx,]))
sink()

test_stage <- summary(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ Stage + strata(Set)"),
														data = metadata[test_idx,]))
stage_lr <- test_stage$logtest[1]
stage_p <- test_stage$coefficients[5]
stage_pc <- ifelse(stage_p < 0.001, "p < 0.001",
									 paste0("p = ", signif(stage_p, 2)))

test_risk <- summary(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ prediction + strata(Set)"),
													 data = metadata[test_idx,]))
risk_lr <- test_risk$logtest[1]
risk_p <- test_risk$coefficients[5]
risk_pc <- ifelse(risk_p < 0.001, "p < 0.001",
									paste0("p = ", signif(risk_p, 2)))

test_both <- summary(coxph(formula = as.formula("Surv(OS_years, OS_event) ~ Stage + prediction + strata(Set)"),
													 data = metadata[test_idx,]))
both_lr <- test_both$logtest[1]
both_stage_p <- test_both$coefficients[1,5]
both_risk_p <- test_both$coefficients[2,5]
both_stage_pc <- ifelse(both_stage_p < 0.001, "p < 0.001",
												paste0("p = ", signif(both_stage_p, 2)))
both_risk_pc <- ifelse(both_risk_p < 0.001, "p < 0.001",
											 paste0("p = ", signif(both_risk_p, 2)))

lrt_data <- data.frame(LR = c(stage_lr, both_lr-stage_lr, risk_lr, both_lr-risk_lr),
											 p = c(stage_pc, both_risk_pc, risk_pc, both_stage_pc),
											 type = c("Stage", "Model risk", "Model risk", "Stage"),
											 group = c("Order 1", "Order 1", "Order 2", "Order 2"),
											 ordering = c(1,2,1,2)) %>%
	group_by(group) %>%
	mutate(cumulative_LR = cumsum(LR),
				 mid_LR = cumulative_LR - 0.5 * LR)

pdf(paste0(out_dir, "/LRT_plot_stage_model_risk_test_set_strata.pdf"),
		width = 5, height = 5)
lrt_data %>% 
	ggplot(aes(x = group, y = LR, fill = reorder(type,-ordering))) + 
	geom_bar(stat = "identity", position = "stack", data = lrt_data %>% filter(group == "Order 1")) +
	geom_bar(stat = "identity", position = "stack", data = lrt_data %>% filter(group == "Order 2")) +
	geom_text(aes(y = mid_LR, label = p), 
						color = "black", 
						size = 5) +
	theme_minimal(base_size = 18) +
	ylab("Order of addition to model\n\nChange in LR statistic") +
	xlab("") +
	theme(legend.title=element_blank()) +
	scale_fill_manual(values = c("Stage" = "#edb2d5", "Model risk" = "#34b399")) +
	theme(legend.position="top")
dev.off()

# Final model KM plots (median)

metadata$risk_groups_med <- NA
metadata$risk_groups_med[c03_idx] <- ifelse(metadata$prediction[c03_idx] >= median(metadata$prediction[c03_idx]), "High", "Low")
metadata$risk_groups_med[fuscc_idx] <- ifelse(metadata$prediction[fuscc_idx] >= median(metadata$prediction[fuscc_idx]), "High", "Low")
metadata$risk_groups_med[mb_idx] <- ifelse(metadata$prediction[mb_idx] >= median(metadata$prediction[mb_idx]), "High", "Low")
metadata$risk_groups_med[tcga_idx] <- ifelse(metadata$prediction[tcga_idx] >= median(metadata$prediction[tcga_idx]), "High", "Low")
metadata$risk_groups_med <- factor(metadata$risk_groups_med,  c("Low", "High"))

splots <- list()
formula <- as.formula("Surv(OS_years, OS_event) ~ risk_groups_med")
fit <- surv_fit(formula, data = metadata[c03_idx,])
splots[[1]] <- ggsurvplot(fit, data = metadata[c03_idx,], 
													risk.table = TRUE, pval = TRUE,
													legend.labs = names(med_colors),
													palette = med_colors,
													legend.title = "Predicted risk",
													xlab = "Time (years)", ylab = "Overall survival",
													title = "CALGB 40603 (Training)",
													ggtheme = theme_classic2(base_size=16),
													xlim = c(0,7), break.x.by = 1, legend = "none")
summary_fit_7 <- summary(fit, times = 7)
sink(paste0(out_dir, "/final_model_median_surv_7yr_C03_training.txt"))
print(summary_fit_7)
sink()

fit <- surv_fit(formula, data = metadata[fuscc_idx,])
splots[[2]] <- ggsurvplot(fit, data = metadata[fuscc_idx,], risk.table = TRUE, pval = TRUE,
													legend.labs = names(med_colors),
													palette = med_colors,
													legend.title = "Predicted risk",
													xlab = "Time (years)", ylab = "Overall survival",
													title = "FUSCC (Test)",
													ggtheme = theme_classic2(base_size=16),
													xlim = c(0,7), break.x.by = 1, legend = "none")
summary_fit_7 <- summary(fit, times = 7)
sink(paste0(out_dir, "/final_model_median_surv_7yr_FUSCC.txt"))
print(summary_fit_7)
sink()

fit <- surv_fit(formula, data = metadata[mb_idx,])
splots[[3]] <- ggsurvplot(fit, data = metadata[mb_idx,], risk.table = TRUE, pval = TRUE,
													legend.labs = names(med_colors),
													palette = med_colors,
													legend.title = "Predicted risk",
													xlab = "Time (years)", ylab = "Overall survival",
													title = "METABRIC (Test)",
													ggtheme = theme_classic2(base_size=16),
													xlim = c(0,7), break.x.by = 1, legend = "none")
summary_fit_7 <- summary(fit, times = 7)
sink(paste0(out_dir, "/final_model_median_surv_7yr_METABRIC.txt"))
print(summary_fit_7)
sink()

fit <- surv_fit(formula, data = metadata[tcga_idx,])
splots[[4]] <- ggsurvplot(fit, data = metadata[tcga_idx,], risk.table = TRUE, pval = TRUE,
													legend.labs = names(med_colors),
													palette = med_colors,
													legend.title = "Predicted risk",
													xlab = "Time (years)", ylab = "Overall survival",
													title = "TCGA (Test)",
													ggtheme = theme_classic2(base_size=16),
													xlim = c(0,7), break.x.by = 1, legend = "none")
summary_fit_7 <- summary(fit, times = 7)
sink(paste0(out_dir, "/final_model_median_surv_7yr_TCGA.txt"))
print(summary_fit_7)
sink()

plot <- arrange_ggsurvplots(splots, print = TRUE,
														ncol = 4, nrow = 1, risk.table.height = 0.33)
pdf(paste0(out_dir, "/final_model_median_KMs_7yr.pdf"), width = 20, height = 5)
print(plot)
dev.off()

splots <- list()
formula <- as.formula("Surv(OS_years, OS_event) ~ risk_groups_med")
fit <- surv_fit(formula, data = metadata[test_idx,])
splots[[1]] <- ggsurvplot(fit, data = metadata[test_idx,], risk.table = TRUE, pval = TRUE,
													legend.labs = names(med_colors),
													palette = med_colors,
													legend.title = "Predicted risk",
													xlab = "Time (years)", ylab = "Overall survival",
													title = "Combined Test",
													ggtheme = theme_classic2(base_size=16),
													xlim = c(0,7), break.x.by = 1, legend = "none")

plot <- arrange_ggsurvplots(splots, print = TRUE,
														ncol = 1, nrow = 1, risk.table.height = 0.33)
pdf(paste0(out_dir, "/final_model_median_KM_7yr_combined_test.pdf"),
		width = 5, height = 5)
print(plot)
dev.off()

summary_fit_7 <- summary(fit, times = 7)
sink(paste0(out_dir, "/final_model_median_surv_7yr_combined_test.txt"))
print(summary_fit_7)
sink()

# KM plots (tertile)

metadata$risk_groups_tert <- NA
metadata$risk_groups_tert[c03_idx] <- ifelse(metadata$prediction[c03_idx] >= quantile(metadata$prediction[c03_idx], 2/3), "High", 
																						 ifelse(metadata$prediction[c03_idx] <= quantile(metadata$prediction[c03_idx], 1/3), "Low", "Med"))
metadata$risk_groups_tert[fuscc_idx] <- ifelse(metadata$prediction[fuscc_idx] >= quantile(metadata$prediction[fuscc_idx], 2/3), "High", 
																							 ifelse(metadata$prediction[fuscc_idx] <= quantile(metadata$prediction[fuscc_idx], 1/3), "Low", "Med"))
metadata$risk_groups_tert[mb_idx] <- ifelse(metadata$prediction[mb_idx] >= quantile(metadata$prediction[mb_idx], 2/3), "High", 
																						ifelse(metadata$prediction[mb_idx] <= quantile(metadata$prediction[mb_idx], 1/3), "Low", "Med"))
metadata$risk_groups_tert[tcga_idx] <- ifelse(metadata$prediction[tcga_idx] >= quantile(metadata$prediction[tcga_idx], 2/3), "High", 
																							ifelse(metadata$prediction[tcga_idx] <= quantile(metadata$prediction[tcga_idx], 1/3), "Low", "Med"))
metadata$risk_groups_tert <- factor(metadata$risk_groups_tert,  c("Low", "Med", "High"))

splots <- list()
formula <- as.formula("Surv(OS_years, OS_event) ~ risk_groups_tert")
fit <- surv_fit(formula, data = metadata[c03_idx,])
splots[[1]] <- ggsurvplot(fit, data = metadata[c03_idx,], 
													risk.table = TRUE, pval = TRUE,
													legend.labs = names(tert_colors),
													palette = tert_colors,
													legend.title = "Predicted risk",
													xlab = "Time (years)", ylab = "Overall survival",
													title = "CALGB 40603 (Training)",
													ggtheme = theme_classic2(base_size=16),
													xlim = c(0,7), break.x.by = 1, legend = "none")
summary_fit_7 <- summary(fit, times = 7)
sink(paste0(out_dir, "/final_model_tertile_surv_7yr_C03_training.txt"))
print(summary_fit_7)
sink()
fit <- surv_fit(formula, data = metadata[fuscc_idx,])
splots[[2]] <- ggsurvplot(fit, data = metadata[fuscc_idx,], risk.table = TRUE, pval = TRUE,
													legend.labs = names(tert_colors),
													palette = tert_colors,
													legend.title = "Predicted risk",
													xlab = "Time (years)", ylab = "Overall survival",
													title = "FUSCC (Test)",
													ggtheme = theme_classic2(base_size=16),
													xlim = c(0,7), break.x.by = 1, legend = "none")
summary_fit_7 <- summary(fit, times = 7)
sink(paste0(out_dir, "/final_model_tertile_surv_7yr_FUSCC.txt"))
print(summary_fit_7)
sink()
fit <- surv_fit(formula, data = metadata[mb_idx,])
splots[[3]] <- ggsurvplot(fit, data = metadata[mb_idx,], risk.table = TRUE, pval = TRUE,
													legend.labs = names(tert_colors),
													palette = tert_colors,
													legend.title = "Predicted risk",
													xlab = "Time (years)", ylab = "Overall survival",
													title = "METABRIC (Test)",
													ggtheme = theme_classic2(base_size=16),
													xlim = c(0,7), break.x.by = 1, legend = "none")
summary_fit_7 <- summary(fit, times = 7)
sink(paste0(out_dir, "/final_model_tertile_surv_7yr_METABRIC.txt"))
print(summary_fit_7)
sink()
fit <- surv_fit(formula, data = metadata[tcga_idx,])
splots[[4]] <- ggsurvplot(fit, data = metadata[tcga_idx,], risk.table = TRUE, pval = TRUE,
													legend.labs = names(tert_colors),
													palette = tert_colors,
													legend.title = "Predicted risk",
													xlab = "Time (years)", ylab = "Overall survival",
													title = "TCGA (Test)",
													ggtheme = theme_classic2(base_size=16),
													xlim = c(0,7), break.x.by = 1, legend = "none")
summary_fit_7 <- summary(fit, times = 7)
sink(paste0(out_dir, "/final_model_tertile_surv_7yr_TCGA.txt"))
print(summary_fit_7)
sink()

plot <- arrange_ggsurvplots(splots, print = TRUE,
														ncol = 4, nrow = 1, risk.table.height = 0.33)
pdf(paste0(out_dir, "/final_model_tertile_KMs_7yr.pdf"),
		width = 20, height = 5)
print(plot)
dev.off()

splots <- list()
formula <- as.formula("Surv(OS_years, OS_event) ~ risk_groups_tert")
fit <- surv_fit(formula, data = metadata[test_idx,])
splots[[1]] <- ggsurvplot(fit, data = metadata[test_idx,], risk.table = TRUE, pval = TRUE,
													legend.labs = names(tert_colors),
													palette = tert_colors,
													legend.title = "Predicted risk",
													xlab = "Time (years)", ylab = "Overall survival",
													title = "Combined Test",
													ggtheme = theme_classic2(base_size=16),
													xlim = c(0,7), break.x.by = 1, legend = "none")

plot <- arrange_ggsurvplots(splots, print = TRUE,
														ncol = 1, nrow = 1, risk.table.height = 0.33)
pdf(paste0(out_dir, "/final_model_tertile_KM_7yr_combined_test.pdf"),
		width = 5, height = 5)
print(plot)
dev.off()

summary_fit_1 <- summary(fit, times = 1)
sink(paste0(out_dir, "/final_model_tertile_surv_1yr_combined_test.txt"))
print(summary_fit_1)
sink()

summary_fit_2 <- summary(fit, times = 2)
sink(paste0(out_dir, "/final_model_tertile_surv_2yr_combined_test.txt"))
print(summary_fit_2)
sink()

summary_fit_3 <- summary(fit, times = 3)
sink(paste0(out_dir, "/final_model_tertile_surv_3yr_combined_test.txt"))
print(summary_fit_3)
sink()

summary_fit_4 <- summary(fit, times = 4)
sink(paste0(out_dir, "/final_model_tertile_surv_4yr_combined_test.txt"))
print(summary_fit_4)
sink()

summary_fit_5 <- summary(fit, times = 5)
sink(paste0(out_dir, "/final_model_tertile_surv_5yr_combined_test.txt"))
print(summary_fit_5)
sink()

summary_fit_6 <- summary(fit, times = 6)
sink(paste0(out_dir, "/final_model_tertile_surv_6yr_combined_test.txt"))
print(summary_fit_6)
sink()

summary_fit_7 <- summary(fit, times = 7)
sink(paste0(out_dir, "/final_model_tertile_surv_7yr_combined_test.txt"))
print(summary_fit_7)
sink()

saveRDS(metadata, paste0(out_dir, "/combined_metadata_with_final_model_scores.RDS"))
