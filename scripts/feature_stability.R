library(infotheo)  # for mutual information
library(pROC)      # for AUC
library(dplyr)



suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(infotheo))
suppressPackageStartupMessages(library(caret))


df <- read.csv(paste0(here(), "/raw/somascan/normalized/csv/SS-2342615_v4.1_anmlSMP__AALS-v1.0-participant-samples.csv"), check.names = F)
mapper <- read.csv(paste0(here(), "/raw/als-metadata/answerals/aals_dataportal_datatable.csv"), check.names = F)
meta <- read.csv(paste0(here(),"/proc/meta.csv"))
annotation <- read.csv(paste0(here(), "/raw/somascan/documents/annotation.csv"))


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


meta_P <- df1[,meta_cols]
X <- df1[,!colnames(df1) %in% meta_cols]



idx_control <- which(meta_P$SUBJECT_GROUP == "Healthy Control")
idx_als <- which(meta_P$SUBJECT_GROUP == "ALS")
idx_nonals <- which(meta_P$SUBJECT_GROUP == "Non-ALS MND")


y <- df1$apoe_status
X <- log2(X[,!colnames(X) %in% c("apoe_status", "SUBJECT_GROUP")])


colnames(X) <- make.names(colnames(X))


# Helper: compute MI for each column of X vs a binary outcome
compute_mi_vector <- function(X, y, nbins = 5) {
  # X: n x p matrix/data.frame
  # y: factor or vector, length n
  X <- as.matrix(X)
  y <- as.factor(y)
  
  # discretize all proteins at once (equal-frequency bins)
  X_disc <- infotheo::discretize(X, disc = "equalfreq", nbins = nbins)
  
  p <- ncol(X_disc)
  mi <- numeric(p)
  
  for (j in seq_len(p)) {
    mi[j] <- infotheo::mutinformation(X_disc[, j], y)
  }
  
  names(mi) <- colnames(X)
  mi
}


mi_stability <- function(X, y, train_idx,
                         B = 100,       # number of resamples
                         top_k = 100,  # how many features per resample
                         nbins = 5) {
  
  X_train <- X[train_idx, , drop = FALSE]
  y_train <- y[train_idx]
  
  p <- ncol(X_train)
  prot_names <- colnames(X_train)
  
  # storage
  sel_counts <- numeric(p)
  mi_mat <- matrix(NA_real_, nrow = B, ncol = p,
                   dimnames = list(NULL, prot_names))
  
  n_train <- length(train_idx)
  
  set.seed(1)  # or pass seed as arg if you want
  
  for (b in seq_len(B)) {
    # bootstrap indices within training set
    idx_b <- sample.int(n_train, size = n_train, replace = TRUE)
    
    X_b <- X_train[idx_b, , drop = FALSE]
    y_b <- y_train[idx_b]
    
    mi_b <- compute_mi_vector(X_b, y_b, nbins = nbins)
    mi_mat[b, ] <- mi_b
    
    # select top_k features in this resample
    top_idx <- order(mi_b, decreasing = TRUE)[seq_len(top_k)]
    sel_counts[top_idx] <- sel_counts[top_idx] + 1L
  }
  
  # selection frequency across resamples
  sel_freq <- sel_counts / B
  
  # summary MI across resamples (median, mean, SD etc.)
  mi_median <- apply(mi_mat, 2, median, na.rm = TRUE)
  mi_mean   <- apply(mi_mat, 2, mean,   na.rm = TRUE)
  mi_sd     <- apply(mi_mat, 2, sd,     na.rm = TRUE)
  
  tibble(
    protein   = prot_names,
    mi_median = mi_median,
    mi_mean   = mi_mean,
    mi_sd     = mi_sd,
    sel_freq  = sel_freq
  ) %>%
    arrange(desc(mi_median))
}


stab_res_control <- mi_stability(
  X          = X,
  y          = y,
  train_idx  = idx_control,
  B          = 100,
  top_k      = 100
)

stab_res_als <- mi_stability(
  X          = X,
  y          = y,
  train_idx  = idx_als,
  B          = 100,
  top_k      = 100
)

stab_res_nonals <- mi_stability(
  X          = X,
  y          = y,
  train_idx  = idx_nonals,
  B          = 100,
  top_k      = 100
)

mi_ctrl <- compute_mi_vector(X[idx_control,], y[idx_control])
mi_als <- compute_mi_vector(X[idx_als,], y[idx_als])
mi_nonals <- compute_mi_vector(X[idx_nonals,], y[idx_nonals])

# Combine stability results
combined <- stab_res_control %>%
  select(protein, 
         mi_median_ctrl = mi_median,
         sel_freq_ctrl  = sel_freq) %>%
  full_join(
    stab_res_als %>%
      select(protein,
             mi_median_als = mi_median,
             sel_freq_als  = sel_freq),
    by = "protein"
  ) %>%
  full_join(
    stab_res_nonals %>%
      select(protein,
             mi_median_nonals = mi_median,
             sel_freq_nonals  = sel_freq),
    by = "protein"
  ) %>%
  mutate(
    mi_full_ctrl = mi_ctrl[protein],
    mi_full_als  = mi_als[protein],
    mi_full_nonals  = mi_nonals[protein]
  )

# Optional: deltas (can be handy later)
combined <- combined %>%
  mutate(
    delta_mi_median = mi_median_als - mi_median_ctrl - mi_median_nonals,
    delta_sel_freq  = sel_freq_als  - sel_freq_ctrl - sel_freq_ctrl
  )


thr_mi_ctrl <- quantile(combined$mi_median_ctrl, probs = 0.8, na.rm = TRUE)
thr_mi_als  <- quantile(combined$mi_median_als,  probs = 0.8, na.rm = TRUE)
thr_mi_nonals  <- quantile(combined$mi_median_nonals,  probs = 0.8, na.rm = TRUE)
thr_sel <- 0.8

combined_tagged <- combined %>%
  mutate(
    high_ctrl = !is.na(mi_median_ctrl) & 
      mi_median_ctrl >= thr_mi_ctrl &
      sel_freq_ctrl  >= thr_sel,
    
    high_als  = !is.na(mi_median_als) & 
      mi_median_als >= thr_mi_als &
      sel_freq_als  >= thr_sel,
    
    high_nonals  = !is.na(mi_median_nonals) & 
      mi_median_nonals >= thr_mi_nonals &
      sel_freq_nonals  >= thr_sel
  ) %>%
  mutate(
    category = case_when(
      high_ctrl & high_als & high_nonals ~ "shared",
      high_ctrl & !high_als & !high_nonals ~ "control_specific",
      !high_ctrl & high_nonals & !high_als ~ "nonALS_specific",
      !high_ctrl & !high_nonals & high_als ~ "ALS_specific",
      TRUE ~ "low_or_nonspecific"
    )
  )


combined_tagged$protein[combined_tagged$category == 'shared']

write.csv(combined_tagged, paste0(here(),"output/somascan/proteins_category.csv"), row.names = F)
