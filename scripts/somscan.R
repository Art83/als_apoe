df <- read.csv("D:/AZ/als/raw/somascan/normalized/csv/SS-2342615_v4.1_anmlSMP__AALS-v1.0-participant-samples.csv")

df1 <- df %>%
  filter(TimePoint == 1) %>%
  dplyr::rename(GUID = SubjectUID) %>%
  inner_join(mapper[,c("GUID", "Participant_ID")], by="GUID") %>%
  inner_join(meta[, c("Participant_ID", "apoe_status")]) %>%
  select(-Barcode2d, -TimePoint, -Participant_ID, - GUID)



df1[colnames(df1)!="apoe_status"] <- log2(df1[colnames(df1)!="apoe_status"])


colnames(df1) <- make.names(colnames(df1))
feat_names <- colnames(df1)[colnames(df1) != "apoe_status"]

mi_scores <- furrr::future_map_dfr(
  feat_names,
  ~ FSelectorRcpp::information_gain(as.formula(paste("apoe_status ~", .x)), df1)
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
control <- trainControl(
  method = "repeatedcv", number = 5, repeats = 5,
  summaryFunction = twoClassSummary,
  classProbs = TRUE,
  sampling = "down"
)

model <- train(target ~ ., data = training, method = "rpart",
               metric = "ROC", tuneLength = 100, trControl = control,
               preProcess = c("center", "scale"))

preds_class <- predict(model, test)
confusionMatrix(preds_class, test$target, positive = "pos")


varImp(model)


