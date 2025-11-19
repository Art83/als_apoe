suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(caret))


setwd("D:/AZ/als_apoe")
df_raw <- read.csv("raw/transcriptomics/4_matrix/AnswerALS-290-T-v1-release5_rna-counts-matrix.csv", check.names = FALSE)

mapper <- read.csv("raw/transcriptomics/4_matrix/Sample Mapping Information/Sample Mapping File Feb 2024.csv", check.names = FALSE)
meta <- read.csv("proc/meta.csv")

rownames(df_raw) <- df_raw[,1]
df_raw <- df_raw[,-1]

mapper <- mapper %>% filter(Sample_ID %in% colnames(df_raw))

meta <- meta %>%
  inner_join(mapper[,c("Participant_ID", "Sample_ID", "Batch")], by="Participant_ID") %>%
  filter(SUBJECT_GROUP %in% c("ALS", "Healthy Control"))

idx <- which(meta$Batch %in% c(1,8,9,35,37,38))
meta <- meta[-idx,]

df_raw <- df_raw[, meta$Sample_ID]

stopifnot(all(colnames(df_raw) == meta$Sample_ID))

# basic gene filters (unsupervised)
df_raw <- df_raw[rowSums(df_raw) > 0, ]
df_raw <- df_raw[apply(df_raw, 1, var) > 0, ]
keep <- rowMeans(df_raw > 1) >= 0.2
df_raw <- df_raw[keep, ]



dds <- DESeq2::DESeqDataSetFromMatrix(df_raw, meta, design = ~ 1)
dds <- DESeq2::estimateSizeFactors(dds)
vst_mat <- DESeq2::vst(dds, blind = TRUE)
expr <- SummarizedExperiment::assay(vst_mat)


# Batch correction on training
design_bc <- model.matrix(~ SUBJECT_GROUP + sex + age, data = meta)
expr_bc <- limma::removeBatchEffect(expr, batch = meta$Batch, design = design_bc)

# transpose to samples x genes
dat <- as.data.frame(t(expr_bc))
dat$outcome <- meta$apoe_status

# ML
set.seed(42)
library(caret)

# Ensure outcome is factor with consistent labels
dat$outcome <- factor(ifelse(dat$outcome == "APOE4+", "pos", "neg"), levels = c("pos" , "neg"))

train_idx <- createDataPartition(dat$outcome, p = 0.7, list = FALSE)

train_df <- dat[train_idx, ]
test_df  <- dat[-train_idx, ]


# FS only on training

# remove outcome column from feature names
feat_names <- colnames(train_df)[colnames(train_df) != "outcome"]

mi_scores_train <- furrr::future_map_dfr(
  feat_names,
  ~ FSelectorRcpp::information_gain(
    as.formula(paste("outcome ~", .x)),
    train_df[, c(.x, "outcome")]
  )
)

mi_scores_train <- mi_scores_train[order(mi_scores_train$importance, decreasing = TRUE), ]

# choose features; you can keep your threshold or adjust to top N
mi_scores_not_zero <- mi_scores_train[mi_scores_train$importance > 0, ]
best <- mi_scores_not_zero$attributes[mi_scores_not_zero$importance > 0.06]


# Construct ML data frames with selected features
training <- train_df[, best, drop = FALSE]
training$target <- train_df$outcome

test <- test_df[, best, drop = FALSE]
test$target <- test_df$outcome 


# SMOTE only on training data
smote_res <- smotefamily::SMOTE(
  X = training[, -ncol(training)],
  target = training$target,
  K = 5
)

new_training <- smote_res$data
# smotefamily::SMOTE usually names the target column "class"
new_training$class <- factor(new_training$class, levels = c("pos", "neg"))



# ------------------ STEP 8: Model training ------------------
new_training <- smotefamily::SMOTE(training[,-ncol(training)], training$target, K = 5)

new_training <- new_training$data


control <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 5,
  summaryFunction = twoClassSummary,
  classProbs = TRUE,
  savePredictions = "final"
)

set.seed(42)
model <- train(
  class ~ .,
  data = new_training,
  method = "rf",
  metric = "ROC",
  tuneLength = 100,
  trControl = control,
  preProcess = c("center", "scale")  # caret will estimate params on training folds only
)

# Predict on test set (no resampling here)
pred_probs <- predict(model, newdata = test, type = "prob")[, "pos"]
pred_class <- predict(model, newdata = test)

cm <- confusionMatrix(pred_class, test$target, positive = "pos")

roc_obj <- pROC::roc(response = test$target, predictor = pred_probs, levels = c("neg", "pos"))



library(limma)

# expr_bc: genes x samples, we want genes x samples for limma
design <- model.matrix(~ apoe_status + sex + age + Batch, data = meta)

fit <- lmFit(expr_bc, design)
fit <- eBayes(fit)
res_limma <- topTable(fit, coef = "apoe_statusAPOE4+", n = Inf)



# Save metrics
metrics_all <- list(
  Type = "APOE+",
  Sensitivity = cm$byClass["Sensitivity"],
  Specificity = cm$byClass["Specificity"],
  PPV = cm$byClass["Pos Pred Value"],
  NPV = cm$byClass["Neg Pred Value"],
  AUC_test = as.numeric(pROC::auc(roc_obj)),
  AUC = mean(model$resample$ROC),
  SD = sd(model$resample$ROC)
)

jsonlite::write_json(metrics_all, "output/transcriptomics/final_ML_test_metrics.json", pretty = TRUE)
jsonlite::write_json(model$bestTune, "output/transcriptomics/ML_best_hyperparameters_somascan.json", pretty = TRUE)
write.csv(res_limma, "output/transcriptomics/limma_results.csv", row.names = T)