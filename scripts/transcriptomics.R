library(dplyr)

# Transcriptomics
df_raw <- read.csv("D:/AZ/als/raw/transcriptomics/4_matrix/AnswerALS-290-T-v1-release5_rna-counts-matrix.csv", check.names = F)

mapper <- read.csv("D:/AZ/als/raw/transcriptomics/4_matrix/Sample Mapping Information/Sample Mapping File Feb 2024.csv", check.names = F)
meta <- read.csv("D:/AZ/als/proc/meta.csv")


row.names(df_raw) <- df_raw[,1]
df_raw <- df_raw[,-1]

mapper <- mapper %>% filter(Sample_ID %in% colnames(df_raw))

meta <- meta %>%
  inner_join(mapper[,c("Participant_ID", "Sample_ID", "Batch")], by="Participant_ID") %>%
  filter(SUBJECT_GROUP %in% c("ALS", "Healthy Control")) 


ggplot(meta, aes(x = Batch, fill = apoe_status)) + geom_bar(position = "fill") +
  ylab("Proportion") + ggtitle("APOE distribution across batches")

idx <- which(meta$Batch %in% c(1,8,9,35,37,38))
meta <- meta[-idx,]


mt_apoe <- as.integer(meta$apoe_status == "APOE4+")
cor.test(mt_apoe, as.integer(as.factor(meta$Batch)))



df_raw <- df_raw[,meta$Sample_ID]

all(colnames(df_raw) == meta$Sample_ID)

df_raw <- df_raw[rowSums(df_raw) > 0, ]
df_raw <- df_raw[apply(df_raw, 1, var) > 0, ]
keep <- rowMeans(df_raw > 1) >= 0.2  # expressed in >20% of samples
df_raw<- df_raw[keep, ]






dds_train <- DESeq2::DESeqDataSetFromMatrix(df_raw, meta, design=~1)
dds_train <-  DESeq2::estimateSizeFactors(dds_train)
vst_train <-  DESeq2::vst(dds_train, blind=TRUE)
expr_train <-  SummarizedExperiment::assay(vst_train)




expr_train <- expr_train[rowSums(expr_train) > 0, ] # drop 0 count
expr_train <- expr_train[apply(expr_train, 1, var) > 0, ] # drop low variance 
keep <- rowMeans(expr_train > 1) >= 0.2  # expressed in >20% of samples
expr_train <- expr_train[keep, ]



# Batch correction on training
combat_train <- sva::ComBat(expr_train, batch=meta$Batch,
                       mod=model.matrix(~ SUBJECT_GROUP+sex+age+apoe_status, data=meta))

combat_train <- limma::removeBatchEffect(expr_train, batch=meta$Batch,
                            mod=model.matrix(~ SUBJECT_GROUP+sex+age+apoe_status, data=meta))


combat_train <- limma::normalizeBetweenArrays(combat_train)

combat_train <- as.data.frame(t(combat_train))

combat_train$outcome <- meta$apoe_status[match(row.names(combat_train), meta$Sample_ID)]




# Surrogates
null_model <- model.matrix(~ Batch + age + sex, data = meta)
num_sv <- sva::num.sv(expr_train, model.matrix(~ Batch + SUBJECT_GROUP+sex+age+apoe_status, data = meta), method = "leek")
sva_obj <- sva::sva(expr_train, 
                    mod = model.matrix(~ Batch + age + sex + SUBJECT_GROUP + apoe_status, data = meta), 
                    mod0 = model.matrix(~ Batch + age + sex + SUBJECT_GROUP, data = meta))




expression_corrected <- limma::removeBatchEffect(expr_train, 
                                                 covariates = sva_obj$sv, 
                                                 design = model.matrix(~ Batch + age + sex + SUBJECT_GROUP + apoe_status, data = meta))

expression_corrected <- limma::normalizeBetweenArrays(expression_corrected)
combat_train <- as.data.frame(t(expression_corrected))
combat_train$outcome <- meta$apoe_status[match(row.names(combat_train), meta$Sample_ID)]

colnames(combat_train) <- make.names(colnames(combat_train))
feat_names <- colnames(combat_train)[colnames(combat_train) != "outcome"]

mi_scores <- furrr::future_map_dfr(
  feat_names,
  ~ FSelectorRcpp::information_gain(as.formula(paste("outcome ~", .x)), combat_train)
)


mi_scores <- mi_scores[order(mi_scores$importance,decreasing = T),]
mi_scores_not_zero <- mi_scores[mi_scores$importance > 0, ]      
best <- mi_scores_not_zero$attributes[mi_scores_not_zero$importance > 0.05]




set.seed(42)
split <- sample(2, nrow(combat_train), prob = c(0.7, 0.3), replace = TRUE)
train_df <- combat_train[split == 1, ]
test_df  <- combat_train[split == 2, ]


training <- train_df[, best]
training$target <- factor(ifelse(train_df$outcome == "APOE4+", "pos", "neg"), levels = c("pos", "neg"))

test <- test_df[, best]
test$target <- factor(ifelse(test_df$outcome == "APOE4+" , "pos", "neg"), levels = c("pos", "neg"))

# ------------------ STEP 8: Model training ------------------
new_training <- smotefamily::SMOTE(training[,-ncol(training)], training$target, K = 5)

new_training <- new_training$data

control <- trainControl(
  method = "repeatedcv", number = 5, repeats = 5,
  summaryFunction = twoClassSummary,
  classProbs = TRUE,
  sampling = "up",
  savePredictions = "final"
)

model <- train(class ~ ., data = new_training, method = "xgbTree",
               metric = "ROC", tuneLength = 100, trControl = control,
               preProcess = c("center", "scale"))

preds_class <- predict(model, test)
confusionMatrix(preds_class, test$target, positive = "pos")


varImp(model)


# Get predicted probabilities on test data
probs <- predict(model, newdata = test, type = "prob")

# Lower threshold: e.g. 0.3 instead of 0.5
preds <- ifelse(probs$pos > 0.17, "pos", "neg") |> factor(levels = c("pos","neg"))

confusionMatrix(preds, test$target, positive = "pos")


probs <- predict(model, newdata = test, type = "prob")[,1]
pROC::roc(response = test$target, predictor = probs, plot=T)



















find_best_thresholds <- function(fit, positive_class = "pos") {
  
  # caret stores CV predictions in fit$pred
  probs <- fit$pred[[positive_class]]
  true  <- fit$pred$obs
  
  # --- ROC-based optimal cutoff (Youden's J) ---
  roc_obj <- pROC::roc(response = true, predictor = probs, levels = rev(levels(true)))
  roc_thr <- pROC::coords(roc_obj, "best", ret = "threshold", transpose = FALSE)
  
  # --- Custom search for F1 & Balanced Accuracy ---
  f1_score <- function(thr) {
    preds <- ifelse(probs > thr, positive_class, setdiff(levels(true), positive_class))
    cm <- confusionMatrix(factor(preds, levels = levels(true)), true, positive = positive_class)
    as.numeric(cm$byClass["F1"])
  }
  
  bal_acc <- function(thr) {
    preds <- ifelse(probs > thr, positive_class, setdiff(levels(true), positive_class))
    cm <- confusionMatrix(factor(preds, levels = levels(true)), true, positive = positive_class)
    as.numeric(cm$byClass["Balanced Accuracy"])
  }
  
  thr_seq <- seq(0.05, 0.95, 0.01)
  
  best_f1_thr  <- thr_seq[which.max(sapply(thr_seq, f1_score))]
  best_bal_thr <- thr_seq[which.max(sapply(thr_seq, bal_acc))]
  
  list(
    roc_best = roc_thr,
    f1_best  = best_f1_thr,
    bal_best = best_bal_thr
  )
}

find_best_thresholds(model)






find_best_thresholds <- function(fit, positive_class = "pos") {
  # caret stores CV predictions in fit$pred
  probs <- fit$pred[[positive_class]]
  true  <- fit$pred$obs
  
  # sanity checks
  if (nlevels(true) != 2) stop("Only binary classification supported")
  negative_class <- setdiff(levels(true), positive_class)
  
  # --- ROC-based optimal cutoff (Youden's J) ---
  roc_obj <- pROC::roc(response = true, predictor = probs, levels = rev(levels(true)))
  roc_thr <- pROC::coords(roc_obj, x = "best", best.method = "youden",
                          ret = "threshold", transpose = FALSE)
  
  # --- helper: compute confusion matrix from threshold ---
  get_cm <- function(thr) {
    preds <- factor(ifelse(probs > thr, positive_class, negative_class),
                    levels = levels(true))
    table(preds, true)
  }
  
  # --- F1 score ---
  f1_score <- function(cm) {
    TP <- cm[positive_class, positive_class]
    FP <- cm[positive_class, negative_class]
    FN <- cm[negative_class, positive_class]
    precision <- ifelse((TP + FP) == 0, 0, TP / (TP + FP))
    recall    <- ifelse((TP + FN) == 0, 0, TP / (TP + FN))
    if ((precision + recall) == 0) return(0)
    2 * precision * recall / (precision + recall)
  }
  
  # --- Balanced Accuracy ---
  bal_acc <- function(cm) {
    TP <- cm[positive_class, positive_class]
    TN <- cm[negative_class, negative_class]
    FP <- cm[positive_class, negative_class]
    FN <- cm[negative_class, positive_class]
    sens <- ifelse((TP + FN) == 0, 0, TP / (TP + FN))  # recall
    spec <- ifelse((TN + FP) == 0, 0, TN / (TN + FP))
    (sens + spec) / 2
  }
  
  # --- Threshold search ---
  thr_seq <- sort(unique(probs))  # search only over unique probability values
  
  f1_vals  <- sapply(thr_seq, function(thr) f1_score(get_cm(thr)))
  bal_vals <- sapply(thr_seq, function(thr) bal_acc(get_cm(thr)))
  
  # pick median in case of ties
  best_f1_thr  <- median(thr_seq[f1_vals == max(f1_vals)])
  best_bal_thr <- median(thr_seq[bal_vals == max(bal_vals)])
  
  # --- return both thresholds and metric values ---
  list(
    roc_best  = list(threshold = roc_thr),
    f1_best   = list(threshold = best_f1_thr, value = max(f1_vals)),
    bal_best  = list(threshold = best_bal_thr, value = max(bal_vals))
  )
}



# Diagnostics

form <- ~ (1|Batch) + (1|apoe_status)
vp <- variancePartition::fitExtractVarPartModel(expr_train, form, meta)
variancePartition::plotVarPart(vp)



combat_train %>%
  select(ENSG00000171530) %>%
  cbind(apoe = meta$apoe_status) %>%
  ggplot(aes(x=apoe,y=ENSG00000171530)) +
  geom_boxplot()

