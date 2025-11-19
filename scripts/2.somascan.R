suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(caret))

# ------------------ STEP 1: Feature Selection ------------------

setwd("D:/AZ/als_apoe")
df <- read.csv("raw/somascan/normalized/csv/SS-2342615_v4.1_anmlSMP__AALS-v1.0-participant-samples.csv", check.names = F)
mapper <- read.csv("raw/transcriptomics/4_matrix/Sample Mapping Information/Sample Mapping File Feb 2024.csv", check.names = F)
meta <- read.csv("proc/meta.csv")
annotation <- read.csv("raw/somascan/documents/annotation.csv")

meta <- meta %>%
  dplyr::filter(SUBJECT_GROUP %in% c("ALS", "Healthy Control")) 


misc_columns <- c("Barcode2d", "SubjectUID", "TimePoint")
colnames(df)[!colnames(df) %in% misc_columns] <- annotation$UniProt.ID[match(colnames(df)[!colnames(df) %in% misc_columns], 
                                                       annotation$SeqId)]

df <- df[,colnames(df) != "" & !is.na(colnames(df))]


df1 <- df %>%
  dplyr::filter(TimePoint == 1) %>%
  dplyr::rename(GUID = SubjectUID) %>%
  dplyr::inner_join(mapper[,c("GUID", "Participant_ID")], by="GUID") %>%
  dplyr::inner_join(meta[, c("Participant_ID", "SUBJECT_GROUP", "apoe_status")], by="Participant_ID") %>%
  dplyr::select(-Barcode2d, -TimePoint, -Participant_ID, - GUID)


meta_P <- df1[,colnames(df1) %in% c("apoe_status", "SUBJECT_GROUP")]

idx <- which(meta_P$SUBJECT_GROUP == "Healthy Control")

df2 <- df1[idx, ]
df2 <- df2[,colnames(df2) != 'SUBJECT_GROUP']


df2[,!colnames(df2) %in% c("apoe_status")] <- log2(df2[,!colnames(df2) %in% c("apoe_status")])


colnames(df2) <- make.names(colnames(df2))
feat_names <- colnames(df2)[colnames(df2) != "apoe_status"]


mi_scores <- furrr::future_map_dfr(
  feat_names,
  ~ FSelectorRcpp::information_gain(as.formula(paste("apoe_status ~", .x)), df2),
  .options = furrr::furrr_options(seed = TRUE)
)


mi_scores <- mi_scores[order(mi_scores$importance,decreasing = T),]
mi_scores_not_zero <- mi_scores[mi_scores$importance > 0, ] 


# Sorted MI scores (descending)
mi_scores <- mi_scores[order(mi_scores$importance, decreasing = TRUE), ]
best <- mi_scores_not_zero$attributes[mi_scores_not_zero$importance > 0.3]


means_dif <- sapply(df2, function(x) {
  m <- tapply(x, df2$apoe_status, mean, na.rm=T)
} )

res_dir <- as.data.frame(t(means_dif)) %>%
  dplyr::mutate(UniProt = row.names(.)) %>%
  dplyr::mutate(direction = ifelse(`APOE4+` - `APOE4-` < 0, "down", "up") ) %>%
  dplyr::select(UniProt, direction)

best_direction <- data.frame(UniProt = mi_scores_not_zero$attributes[mi_scores_not_zero$importance > 0.1],
                             importance = mi_scores_not_zero$importance[mi_scores_not_zero$importance > 0.1])%>%
  dplyr::mutate(name = annotation$Target.Name[match(UniProt, annotation$UniProt.ID)] ) %>%
  dplyr::inner_join(res_dir, by="UniProt") %>%
  dplyr::select(UniProt, name, direction, importance)


write.csv(best_direction, "output/mi_score_somascan.csv", row.names = F)

# ------------------ STEP 2: Model training ------------------
colnames(df1) <- make.names(colnames(df1))

train_df <- df1[df1$SUBJECT_GROUP == "Healthy Control",]
test_df  <- df1[df1$SUBJECT_GROUP == "ALS",]



training <- train_df[, best]
training$target <- factor(ifelse(train_df$apoe_status == "APOE4+", "pos", "neg"), levels = c("pos", "neg"))

test <- test_df[, best]
test$target <- factor(ifelse(test_df$apoe_status == "APOE4+" , "pos", "neg"), levels = c("pos", "neg"))


control <- trainControl(
  method = "repeatedcv", number = 5, repeats = 5,
  summaryFunction = twoClassSummary,
  classProbs = TRUE,
  sampling = "up"
)

model <- train(target ~ ., data = training, method = "rpart",
               metric = "ROC", tuneLength = 100, trControl = control,
               preProcess = c("center", "scale"))

preds_class <- predict(model, test)
cm <- confusionMatrix(preds_class, test$target, positive = "pos")


pred_prob <- predict(model, test, type = "prob")[, "pos"]
roc_obj <- pROC::roc(test$target, pred_prob, levels = c("neg", "pos"))

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

jsonlite::write_json(metrics_all, "output/somascan/final_ML_test_metrics_somascan.json", pretty = TRUE)
jsonlite::write_json(model$bestTune, "output/somascan/ML_best_hyperparameters_somascan.json", pretty = TRUE)

# Save feature importances
importance <- varImp(model)$importance
importance$Feature <- rownames(importance)
write.csv(importance, "output/somascan/feature_importance_somascan.csv", row.names = FALSE)
