#Somalogic data

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(here))


df <- read.csv(paste0(here(), "/raw/somascan/normalized/csv/SS-2342615_v4.1_anmlSMP__AALS-v1.0-participant-samples.csv"), check.names = F)
mapper <- read.csv(paste0(here(), "/raw/als-metadata/answerals/aals_dataportal_datatable.csv"), check.names = F)
meta <- read.csv(paste0(here(), "/proc/meta.csv"))
annotation <- read.csv(paste0(here(),"/raw/somascan/documents/annotation.csv"))


misc_columns <- c("Barcode2d", "SubjectUID", "TimePoint")
colnames(df)[!colnames(df) %in% misc_columns] <- annotation$UniProt.ID[match(colnames(df)[!colnames(df) %in% misc_columns], 
                                                       annotation$SeqId)]

df <- df[,colnames(df) != "" & !is.na(colnames(df))]

meta_cols <- c("Participant_ID", "GUID", "Number_of_Visits", "SEX", "SUBJECT_GROUP", "Site_of_Onset", "AGE_AT_SYMPTOM_ONSET",
               "AGE_AT_DEATH","ALSFRS_R_Baseline_Value","ALSFRS_R_Latest_Value","ALSFRS_R_PROGRESSION_SLOPE",   
               "CBS_Baseline_Value", "CBS_Latest_Value", "COGNITIVE_PROGRESSION_SLOPE")

df1 <- df %>%
  dplyr::filter(TimePoint == 1) %>%
  dplyr::rename(GUID = SubjectUID) %>%
  dplyr::inner_join(mapper[,meta_cols], by="GUID") %>%
  dplyr::inner_join(meta[, c("Participant_ID", "apoe_status")], by="Participant_ID") %>%
  dplyr::select(-Barcode2d, -TimePoint)

tb <- table(df1$SUBJECT_GROUP)
cat("Somalogic sttudy\n")
cat(" Number of ALS cases:", tb["ALS"],"\n",
    "Number of non-ALS MND cases:", tb["Non-ALS MND"], "\n",
    "Number of Healthy Control cases:", tb["Healthy Control"])


meta_P <- df1[,meta_cols]
X <- df1[,!colnames(df1) %in% meta_cols]


##--------- Step 1. Feature Selection ---------------
idx <- which(meta_P$SUBJECT_GROUP == "Healthy Control")
df2 <- X[idx, ]

df2[,!colnames(df2) %in% c("apoe_status")] <- log2(df2[,!colnames(df2) %in% c("apoe_status")])
colnames(df2) <- make.names(colnames(df2))
feat_names <- colnames(df2)[colnames(df2) != "apoe_status"]



mi_scores <- furrr::future_map_dfr(
  feat_names,
  ~ FSelectorRcpp::information_gain(as.formula(paste("apoe_status ~", .x)), df2),
  .options = furrr::furrr_options(seed = 1234)
)
mi_scores <- mi_scores[order(mi_scores$importance,decreasing = T),]
mi_scores_not_zero <- mi_scores[mi_scores$importance > 0, ] 
mi_scores <- mi_scores[order(mi_scores$importance, decreasing = TRUE), ]
best <- mi_scores_not_zero$attributes


means_dif <- sapply(df2, function(x) {
  m <- tapply(x, df2$apoe_status, mean, na.rm=T)
} )

res_dir <- as.data.frame(t(means_dif)) %>%
  dplyr::mutate(UniProt = row.names(.)) %>%
  dplyr::mutate(direction = ifelse(`APOE4+` - `APOE4-` < 0, "down", "up") ) %>%
  dplyr::select(UniProt, direction)

best_direction <- data.frame(UniProt = mi_scores_not_zero$attributes,
                             importance = mi_scores_not_zero$importance)%>%
  dplyr::mutate(name = annotation$Target.Name[match(UniProt, annotation$UniProt.ID)] ) %>%
  dplyr::inner_join(res_dir, by="UniProt") %>%
  dplyr::select(UniProt, name, direction, importance)

dir <- paste0(here(), "/output/somascan")
if(!dir.exists(dir)){
  dir.create(dir, recursive = T)
}

write.csv(best_direction, paste0(here(), "/output/somascan/mi_score_somascan.csv"), row.names = F)


# ------------------ STEP 2: Model training ------------------
colnames(df1) <- make.names(colnames(df1))

train_df <- df1[df1$SUBJECT_GROUP == "Healthy Control",]
test_df  <- df1[df1$SUBJECT_GROUP == "ALS",]
test_df_2  <- df1[df1$SUBJECT_GROUP == "Non-ALS MND",]


training <- log2(train_df[, best])
training$target <- factor(ifelse(train_df$apoe_status == "APOE4+", "pos", "neg"), levels = c("pos", "neg"))

new_training <- smotefamily::SMOTE(training[,-ncol(training)], training$target, K = 5)

new_training <- new_training$data


test <- log2(test_df[, best])
test$target <- factor(ifelse(test_df$apoe_status == "APOE4+" , "pos", "neg"), levels = c("pos", "neg"))


test2 <- log2(test_df_2[, best])
test2$target <- factor(ifelse(test_df_2$apoe_status == "APOE4+" , "pos", "neg"), levels = c("pos", "neg"))

control <- trainControl(
  method = "repeatedcv", number = 5, repeats = 5,
  summaryFunction = twoClassSummary,
  classProbs = TRUE,
  sampling = "up"
)


model <- train(class ~ ., data = new_training, method = "rf",
               metric = "ROC", tuneLength = 100, trControl = control,
               preProcess = c("center", "scale"))

# Test on ALS
preds_class <- predict(model, test)
cm <- confusionMatrix(preds_class, test$target, positive = "pos")


pred_prob <- predict(model, test, type = "prob")[, "pos"]
roc_obj <- pROC::roc(test$target, pred_prob, levels = c("neg", "pos"))


# Test on non-ALS MND
preds_class_2 <- predict(model, test2)
cm_2 <- confusionMatrix(preds_class_2, test2$target, positive = "pos")


pred_prob_2 <- predict(model, test2, type = "prob")[, "pos"]
roc_obj_2 <- pROC::roc(test2$target, pred_prob_2, levels = c("neg", "pos"))

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

metrics_all_2 <- list(
  Type = "APOE+",
  Sensitivity = cm_2$byClass["Sensitivity"],
  Specificity = cm_2$byClass["Specificity"],
  PPV = cm_2$byClass["Pos Pred Value"],
  NPV = cm_2$byClass["Neg Pred Value"],
  AUC_test = as.numeric(pROC::auc(roc_obj_2)),
  AUC = mean(model$resample$ROC),
  SD = sd(model$resample$ROC)
)

jsonlite::write_json(metrics_all, paste0(here(), "/output/somascan/final_ML_test_metrics_somascan_ALS.json"), pretty = TRUE)
jsonlite::write_json(metrics_all_2, paste0(here(), "/output/somascan/final_ML_test_metrics_somascan_nonALSmnd.json"), pretty = TRUE)
jsonlite::write_json(model$bestTune, paste0(here(), "/output/somascan/ML_best_hyperparameters_somascan.json"), pretty = TRUE)

# Save feature importance
importance <- varImp(model)$importance
importance$Feature <- rownames(importance)
write.csv(importance, paste0(here(), "/output/somascan/feature_importance_somascan.csv"), row.names = FALSE)
