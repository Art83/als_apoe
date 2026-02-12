# MS

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(here))


df <- read.csv(paste0(here(), "/raw/proteomics/4_matrix/AnswerALS-436-P_proteomics-protein-matrix_correctedImputed.txt"), 
               sep='\t', check.names = F)
mapper <- read.csv(paste0(here(), "/raw/transcriptomics/4_matrix/Sample Mapping Information/Sample Mapping File Feb 2024.csv"), check.names = F)
meta <- read.csv(paste0(here(), "/proc/meta.csv"))


df <- df[,!colnames(df) %in% c("nFragment", "nPeptide")]
colnames(df)[colnames(df) != "Protein"] <- gsub("(-P).+$","\\1",  colnames(df)[colnames(df) != "Protein"])
colnames(df)[colnames(df) != "Protein"] <- gsub("_","-",  colnames(df)[colnames(df) != "Protein"])
df$Protein <-  sapply(strsplit(df$Protein, "[|]"),"[", 2)
df <- df[!is.na(df$Protein), ]
row.names(df) <- df$Protein
df <- df[,-1]


mapper <- mapper %>% filter(Sample_ID %in% colnames(df))

meta_P <- meta %>%
  inner_join(mapper[,c("Participant_ID", "Sample_ID", "Batch")], by="Participant_ID") %>%
  filter(SUBJECT_GROUP %in% c("ALS", "Healthy Control")) %>%
  filter(Sample_ID %in% colnames(df))

df <- df[,colnames(df) %in% meta_P$Sample_ID]

stopifnot(all(colnames(df) == meta_P$Sample_ID))



##------------------Step 1. Preprocessing---------------------
offset <- min(df[df > 0]) / 2
df <- df[rowSums(df) > 0, ]
df <- df[apply(df, 1, var) > 0, ]
keep <- rowMeans(df > 1) >= 0.2
df <- df[keep, ]
df_t <- log2(df+offset)


design_bc <- model.matrix(~ SUBJECT_GROUP + sex + age, data = meta_P)
expr_bc <- limma::removeBatchEffect(df_t, batch = meta_P$Batch, design = design_bc)

dat <- as.data.frame(t(expr_bc))
dat$outcome <- meta_P$apoe_status

##---------------------Step 2.ML prep------------------------
set.seed(42)
library(caret)


dat$outcome <- factor(ifelse(dat$outcome == "APOE4+", "pos", "neg"), levels = c("pos" , "neg"))
train_idx <- createDataPartition(dat$outcome, p = 0.7, list = FALSE)
train_df <- dat[train_idx, ]
test_df  <- dat[-train_idx, ]

feat_names <- colnames(train_df)[colnames(train_df) != "outcome"]

mi_scores_train <- furrr::future_map_dfr(
  feat_names,
  ~ FSelectorRcpp::information_gain(
    as.formula(paste("outcome ~", .x)),
    train_df[, c(.x, "outcome")]
  )
)

mi_scores_train <- mi_scores_train[order(mi_scores_train$importance, decreasing = TRUE), ]
mi_scores_not_zero <- mi_scores_train[mi_scores_train$importance > 0, ]
best <- mi_scores_not_zero$attributes

training <- train_df[, best, drop = FALSE]
training$target <- train_df$outcome

test <- test_df[, best, drop = FALSE]
test$target <- test_df$outcome 



# ------------------ STEP 3.Model training ------------------
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
  preProcess = c("center", "scale")
)

# Predict on test set (no resampling here)
pred_probs <- predict(model, newdata = test, type = "prob")[, "pos"]
pred_class <- predict(model, newdata = test)
cm <- confusionMatrix(pred_class, test$target, positive = "pos")
roc_obj <- pROC::roc(response = test$target, predictor = pred_probs, levels = c("neg", "pos"))




## --------------Step 4. Limma ------------------------------------
design <- model.matrix(~ apoe_status + sex + age, data = meta_P)
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

volc <- res_limma |>
  mutate(
    negLog10FDR = -log10(adj.P.Val),
    sig = adj.P.Val < 0.05
  )


dir <- paste0(here(), "/output/ms")
if(!dir.exists(dir)){
  dir.create(dir, recursive = T)
}

tiff(paste0(here(), "/output/ms/volcano.tiff"), width = 4, height = 4, units = "in", res = 300)
ggplot(volc[volc$logFC < 0.6, ], aes(x = logFC, y = negLog10FDR)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +  # FDR = 0.05
  theme_bw() +
  xlab("log2 FC (APOE4- vs APOE4+)") +
  ylab("-log10 FDR")
dev.off()


jsonlite::write_json(metrics_all, paste0(here(), "/output/ms/final_ML_test_metrics.json"), pretty = TRUE)
jsonlite::write_json(model$bestTune, paste0(here(), "/output/ms/ML_best_hyperparameters_somascan.json"), pretty = TRUE)
write.csv(res_limma, paste0(here(), "/output/ms/limma_results.csv"), row.names = T)
