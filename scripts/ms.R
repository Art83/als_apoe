setwd("D:/AZ/als_apoe/")
library(dplyr)



df <- read.csv("raw/proteomics/4_matrix/AnswerALS-436-P_proteomics-protein-matrix_correctedImputed.txt", 
               sep='\t', check.names = F)
mapper <- read.csv("raw/transcriptomics/4_matrix/Sample Mapping Information/Sample Mapping File Feb 2024.csv", check.names = F)
meta <- read.csv("proc/meta.csv")


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

df_t <- log2(df+1e-6)



# Design matrix
meta_P$SUBJECT_GROUP <- factor(meta_P$SUBJECT_GROUP, levels = c("Healthy Control", "ALS"))
meta_P$Batch <- factor(meta_P$Batch)
design <- model.matrix(~ 0 + SUBJECT_GROUP, data = meta_P)
colnames(design) <- gsub("SUBJECT_GROUP", "", colnames(design))
colnames(design) <- make.names(colnames(design))

fit <- limma::lmFit(df_t, design)

# Contrast: ALS vs Healthy Control (adjusted for Batch)
cont.matrix <- limma::makeContrasts(ALS - `Healthy.Control`, levels = design)
fit2 <- limma::contrasts.fit(fit, cont.matrix)
fit2 <- limma::eBayes(fit2,, trend = TRUE, robust = TRUE)

# ---- Results ----
results <- limma::topTable(fit2, number = Inf, adjust.method = "fdr")

ncol(df_t) == nrow(design)

# Fit linear model
fit <- limma::lmFit(df_t, design) |> limma::eBayes()

# For 4w comparison
results_4w <- limma::topTable(fit2, coef = "SUBJECT_GROUPALS", number = Inf, adjust = "fdr")








