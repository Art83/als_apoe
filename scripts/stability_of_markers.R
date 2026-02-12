suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(caret))



df <- read.csv(paste0(here(),"/raw/somascan/normalized/csv/SS-2342615_v4.1_anmlSMP__AALS-v1.0-participant-samples.csv"), check.names = F)
mapper <- read.csv(paste0(here(), "/raw/als-metadata/answerals/aals_dataportal_datatable.csv"), check.names = F)
alsfrs <- read.csv(paste0(here(), "/raw/als-metadata/answerals/clinical/ALSFRS_R.csv"))
clin <- read.csv(paste0(here(), "/proc/clin_data.csv"))
proteins <- read.csv(paste0(here(), "/output/somascan/mi_score_somascan.csv"))
meta <- read.csv(paste0(here(),"/proc/meta.csv"))
annotation <- read.csv(paste0(here(), "/raw/somascan/documents/annotation.csv"))

prot_cols <- proteins$UniProt



misc_columns <- c("Barcode2d", "SubjectUID", "TimePoint")
colnames(df)[!colnames(df) %in% misc_columns] <- annotation$UniProt.ID[match(colnames(df)[!colnames(df) %in% misc_columns], 
                                                                             annotation$SeqId)]

df <- df[,colnames(df) != "" & !is.na(colnames(df))]

meta_cols <- c("Participant_ID", "GUID", "Number_of_Visits", "SEX", "SUBJECT_GROUP", "Site_of_Onset", "AGE_AT_SYMPTOM_ONSET",
               "AGE_AT_DEATH","ALSFRS_R_Baseline_Value","ALSFRS_R_Latest_Value","ALSFRS_R_PROGRESSION_SLOPE",   
               "CBS_Baseline_Value", "CBS_Latest_Value", "COGNITIVE_PROGRESSION_SLOPE")


df1 <- df %>%
  dplyr::select(SubjectUID, TimePoint, prot_cols) %>%
  dplyr::rename(GUID = SubjectUID) %>%
  dplyr::inner_join(mapper[,meta_cols], by="GUID") %>%
  dplyr::inner_join(meta[, c("Participant_ID", "apoe_status")], by="Participant_ID") %>%
  dplyr::select(-Participant_ID, -Number_of_Visits)



alsfrs_proc <- alsfrs %>%
    select(SubjectUID, Visit_Date, alsfrst) %>%
    rename(GUID = SubjectUID) %>%
    group_by(GUID) %>%
    filter(n() >= 2 ) %>%
    arrange(Visit_Date, .by_group = T) %>%
    mutate(time_years = Visit_Date / 365.25) %>%
    mutate(TimePoint = seq_len(n()))

  
  
df1 <- df1 %>%
    inner_join(alsfrs_proc[,c("GUID", "TimePoint", "time_years")])
  
# function to get APOE4 logFC at one timepoint
get_apoe_effect_tp <- function(df_tp) {
  df_tp %>%
    tidyr::pivot_longer(all_of(prot_cols),
                 names_to = "protein", values_to = "expr") %>%
    group_by(protein) %>%
    group_modify(~ {
      fit <- lm(expr ~ apoe_status, data = .x)
      broom::tidy(fit)
    }) %>%
    ungroup() %>%
    filter(term == "apoe_statusAPOE4+") %>%
    select(protein, estimate, std.error, p.value)
}


effects_by_tp <- df1 %>%
  group_by(TimePoint) %>%
  group_modify(~ get_apoe_effect_tp(.x)) %>%
  ungroup()

# Wide matrix of estimates: rows = protein, cols = TimePoint
effects_mat <- effects_by_tp %>%
  select(TimePoint, protein, estimate) %>%
  tidyr::pivot_wider(names_from = TimePoint, values_from = estimate)

# Pairwise correlations of APOE4 effects across timepoints
cor_mat <- cor(as.matrix(effects_mat[,-1]), use = "pairwise.complete.obs")

cat("# Pairwise correlations of APOE4 effects across 5 timepoints\n")
print(cor_mat)


breaks <- seq(0.90, 1.00, length.out = 101)   # explicitly 0.90–1.00
cols   <- colorRampPalette(c("white", "orange", "red"))(length(breaks) - 1)

tiff(paste0(here(), 'output/somascan/pics/heatmap.tiff'), width = 5, height = 5, units = "in", res = 300)
pheatmap::pheatmap(
  cor_mat,
  color           = cols,
  breaks          = breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  main = "Correlation of APOE4 effects across visits"
)
dev.off()

# 1. Fit PCA on baseline (TimePoint == 1) only
baseline <- df1 %>%
  filter(TimePoint == 1)

pc_sig <- prcomp(baseline[, prot_cols], center = TRUE, scale. = TRUE)

# 2. Project ALL timepoints into that PCA space
X_all <- scale(df1[, prot_cols],
               center = pc_sig$center,
               scale  = pc_sig$scale)
scores_all <- X_all %*% pc_sig$rotation
df1$apoe4_PC1 <- scores_all[, 1]

# 3. Mixed model: PC1 ~ time * APOE + random intercept
library(lme4)



m_pc <- lmer(apoe4_PC1 ~ time_years * apoe_status + (1 | GUID),
             data = df1)

cat("Mixed model: PC1 ~ time * APOE + random intercept\n")
print(summary(m_pc))





tp1 <- effects_by_tp %>% filter(TimePoint == 1) %>% select(protein, est1 = estimate)
tp2 <- effects_by_tp %>% filter(TimePoint == 2) %>% select(protein, est2 = estimate)

plot_df <- inner_join(tp1, tp2, by = "protein")

tiff(paste0(here(), '/output/somascan/pics/cor_visit_1_vs_2.tiff'), width = 5, height = 5, units = "in", res = 300)
ggplot(plot_df, aes(x = est1, y = est2)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlab("APOE4 effect at Visit 1") +
  ylab("APOE4 effect at Visit 2") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = paste0("r = ", round(cor(plot_df$est1, plot_df$est2), 2))) +
  theme_bw()
dev.off()


tiff(paste0(here(), '/output/somascan/pics/trajectories_long.tiff'), width = 5, height = 5, units = "in", res = 300)
ggplot(df1, aes(x = time_years, y = apoe4_PC1, group = GUID, colour = apoe_status)) +
  geom_line(alpha = 0.2) +
  stat_smooth(aes(group = apoe_status), method = "lm", se = FALSE, size = 1.2) +
  theme_bw() +
  xlab("Years from baseline") +
  ylab("APOE4 PC1 score")
dev.off()


# Numeric Stability
tp <- as.numeric(colnames(cor_mat))

cor_long <- as.data.frame(cor_mat) %>%
  tibble::rownames_to_column("tp1") %>%
  tidyr::pivot_longer(-tp1, names_to = "tp2", values_to = "cor") %>%
  mutate(tp1 = as.numeric(tp1),
         tp2 = as.numeric(tp2)) %>%
  filter(tp2 > tp1) %>%                  # off-diagonal only
  mutate(lag = tp2 - tp1)

# Overall stability
cor_long %>%
  summarise(
    mean_cor   = mean(cor, na.rm = TRUE),
    median_cor = median(cor, na.rm = TRUE),
    min_cor    = min(cor, na.rm = TRUE),
    max_cor    = max(cor, na.rm = TRUE)
  )

# Correlation vs lag
cor_by_lag <- cor_long %>%
  group_by(lag) %>%
  summarise(
    mean_cor   = mean(cor, na.rm = TRUE),
    median_cor = median(cor, na.rm = TRUE),
    n_pairs    = n()
  )
cor_by_lag




# Site consistency
sign_consistency <- effects_by_tp %>%
  group_by(protein) %>%
  summarise(
    n_tp      = n(),
    frac_pos  = mean(estimate > 0, na.rm = TRUE),
    frac_neg  = mean(estimate < 0, na.rm = TRUE),
    all_same  = (frac_pos == 1 | frac_neg == 1),
    mostly_same = (frac_pos >= 0.8 | frac_neg >= 0.8)
  )

# How many proteins are direction-consistent?
sign_consistency %>%
  summarise(
    n_proteins      = n(),
    n_all_same      = sum(all_same),
    n_mostly_same   = sum(mostly_same)
  )



# Effect of slopes (ES across visits)
effect_slopes <- effects_by_tp %>%
  group_by(protein) %>%
  filter(n_distinct(TimePoint) >= 2) %>%
  summarise(
    slope      = coef(lm(estimate ~ TimePoint))[2],
    intercept  = coef(lm(estimate ~ TimePoint))[1],
    n_tp       = n(),
    mean_est   = mean(estimate, na.rm = TRUE)
  )

# Quick look at distribution of slopes
tiff(paste0(here(),'output/somascan/pics/hist_slopes_long.tiff'), width = 5, height = 5, units = "in", res = 300)
ggplot(effect_slopes, aes(x = slope)) +
  geom_histogram(bins = 50) +
  theme_bw() +
  xlab("Change in APOE4 effect per TimePoint") +
  ylab("Number of proteins")

dev.off()
