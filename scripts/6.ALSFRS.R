suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(lmerTest))
suppressPackageStartupMessages(library(here))


df <- read.csv(paste0(here(),"/raw/somascan/normalized/csv/SS-2342615_v4.1_anmlSMP__AALS-v1.0-participant-samples.csv"), check.names = F)
mapper <- read.csv(paste0(here(), "/raw/als-metadata/answerals/aals_dataportal_datatable.csv"), check.names = F)
alsfrs <- read.csv(paste0(here(), "/raw/als-metadata/answerals/clinical/ALSFRS_R.csv"))
proteins <- read.csv(paste0(here(), "/output/somascan/mi_score_somascan.csv"))
meta <- read.csv(paste0(here(),"/proc/meta.csv"))
annotation <- read.csv(paste0(here(), "/raw/somascan/documents/annotation.csv"))

prot_cols <- proteins$UniProt
misc_columns <- c("Barcode2d", "SubjectUID", "TimePoint")
colnames(df)[!colnames(df) %in% misc_columns] <- annotation$UniProt.ID[match(colnames(df)[!colnames(df) %in% misc_columns], 
                                                                             annotation$SeqId)]
df <- df[,colnames(df) != "" & !is.na(colnames(df))]
colnames(df) <- make.names(colnames(df))


meta_cols <- c("Participant_ID", "GUID", "Number_of_Visits", "SEX", "SUBJECT_GROUP", "Site_of_Onset", "AGE_AT_SYMPTOM_ONSET",
               "AGE_AT_DEATH")


df1 <- df %>%
  dplyr::select(SubjectUID, TimePoint, all_of(prot_cols)) %>%
  dplyr::rename(GUID = SubjectUID) %>%
  dplyr::inner_join(mapper[,meta_cols], by="GUID") %>%
  dplyr::inner_join(meta[, c("Participant_ID", "age", "apoe_status")], by="Participant_ID") %>%
  dplyr::select(-Participant_ID, -Number_of_Visits)


alsfrs_proc <- alsfrs %>%
  select(SubjectUID, Visit_Date, alsfrst) %>%
  rename(GUID = SubjectUID) %>%
  group_by(GUID) %>%
  filter(n() >= 2 ) %>%
  arrange(Visit_Date, .by_group = T) %>%
  mutate(time_years = Visit_Date / 365.25) %>%
  mutate(TimePoint = seq_len(n()))


########### Longitudinal effect of APOE4 ############
apoe_eff <- 
  alsfrs_proc %>% 
  dplyr::inner_join(df1[, c("GUID", "TimePoint", "age", "SEX", "apoe_status")])

m <- lmer(
  alsfrst ~ time_years * apoe_status + age + SEX + (time_years | GUID),
  data = apoe_eff,
  REML = TRUE
)


m_alsfrs <- lmer(
  alsfrst ~ time_years + (time_years | GUID),
  data = alsfrs_proc 
)

als_slopes <- ranef(m_alsfrs)$GUID %>%
  tibble::rownames_to_column("GUID") %>%
  rename(
    ALSFRS_int  = `(Intercept)`,
    ALSFRS_slope = time_years
  )

dir <- paste0(here(), "/output/somascan/pics")
if(!dir.exists(dir)){
  dir.create(dir, recursive = T)
}

tiff( paste0(here(), '/output/somascan/pics/alsfrs_apoe.tiff'), width = 5, height = 5, units = "in", res = 300)
ggplot(apoe_eff, aes(time_years, alsfrst, colour = apoe_status)) +
  geom_point(position = position_jitter(width = 0.03, height = 0),
             alpha = 0.25, size = 1) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, linewidth = 1) +
  labs(x = "Years from baseline", y = "ALSFRS-R", colour = NULL) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA)
  )
dev.off()


als_df <- alsfrs_proc %>%
  inner_join(df1[df1$SUBJECT_GROUP == "ALS", c("GUID", "TimePoint", prot_cols)])


als_df[,!colnames(als_df) %in% c("GUID", "TimePoint", "Visit_Date", "alsfrst", "time_years")] <- 
  log2(als_df[,!colnames(als_df) %in% c("GUID", "TimePoint", "Visit_Date", "alsfrst", "time_years")])


prot_summ <- als_df %>%
  group_by(GUID) %>%
  group_modify(~ {
    dat <- .x

    base_vals <- dat %>%
      filter(time_years == min(time_years)) %>%
      select(all_of(prot_cols)) %>%
      slice(1)
    
    slope_vals <- purrr::map_dfc(
      prot_cols,
      ~ {
        coef(lm(dat[[.x]] ~ dat$time_years))[2]
      }
    )
    names(slope_vals) <- paste0(prot_cols, "_slope")
    
    bind_cols(base_vals, slope_vals)
  }) %>%
  ungroup()


als_meta <- df1 %>%
  filter(SUBJECT_GROUP == "ALS") %>%
  group_by(GUID) %>%
  summarise(
    apoe_status = first(apoe_status),
    sex = first(SEX),
    onset = first(Site_of_Onset),
    age = first(age),
    age_onset = first(AGE_AT_SYMPTOM_ONSET),
    age_death = first(AGE_AT_DEATH),
    .groups = "drop"
  )


features <- als_slopes %>%
  inner_join(prot_summ, by = "GUID") %>%
  inner_join(als_meta, by = "GUID")


# Slopes
res_slope <- features %>%
  rename(ALSFRS_slope_mm = ALSFRS_slope) %>%
  select(GUID, ALSFRS_slope_mm, apoe_status, ends_with("_slope")) %>%
  tidyr::pivot_longer(
    cols = ends_with("_slope"),
    names_to = "protein",
    values_to = "prot_slope"
  ) %>%
  group_by(protein) %>%
  group_modify(~ {
    fit <- lm(prot_slope ~ ALSFRS_slope_mm + apoe_status, data = .x)
    broom::tidy(fit)
  }) %>%
  ungroup()


res_slope_clean <- res_slope %>%
  filter(term %in% c("ALSFRS_slope_mm", "apoe_statusAPOE4+")) %>%
  group_by(term) %>%
  mutate(p_fdr = p.adjust(p.value, method = "fdr")) %>%
  ungroup() %>%
  filter(p_fdr <= 0.05)


# VIS
res_slope_alsfrs <- res_slope %>%
  filter(term == "ALSFRS_slope_mm") %>%
  mutate(
    protein = sapply(strsplit(protein, "_"), "[[", 1),
    p_fdr = p.adjust(p.value, method = "fdr"),
    signif = p_fdr <= 0.05
  ) %>%
  arrange(estimate)

### 1A. Forest-style coefficient plot
tiff(paste0(here(), '/output/somascan/pics/forest_alsfrs.tiff'), width = 5, height = 5, units = "in", res = 300)
ggplot(res_slope_alsfrs,
       aes(x = estimate,
           y = reorder(protein, estimate),
           colour = signif)) +
  geom_point() +
  geom_errorbarh(aes(xmin = estimate - 1.96 * std.error,
                     xmax = estimate + 1.96 * std.error),
                 height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = c("FALSE" = "grey50", "TRUE" = "red")) +
  labs(
    x = "Effect estimate for ALSFRS-R slope",
    y = "Protein slope",
    colour = "FDR < 0.05"
  ) +
  theme_bw()
dev.off()


## Long per-patient table
long_slope <- features %>%
  rename(ALSFRS_slope_mm = ALSFRS_slope) %>%
  select(GUID, ALSFRS_slope_mm, apoe_status, ends_with("_slope")) %>%
  tidyr::pivot_longer(
    cols = ends_with("_slope"),
    names_to = "protein",
    values_to = "prot_slope"
  )

hits <- res_slope_clean$protein


res_for_hits <- long_slope %>%
  filter(protein %in% hits)

res_for_hits$protein <- sapply(strsplit(res_for_hits$protein,"_"),"[[",1)

## Scatter + regression lines, faceted by protein
tiff(paste0(here(), '/output/somascan/pics/3prot_slopes.tiff'), width = 6, height = 5, units = "in", res = 300)
ggplot(res_for_hits,
       aes(x = ALSFRS_slope_mm,
           y = prot_slope)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ protein, scales = "free_y") +
  labs(
    x = "ALSFRS-R slope",
    y = "Protein slope (log2)",
    colour = "APOE status"
  ) +
  theme_bw()
dev.off()


# 3-protein molecular progression score
prog_mat <- features %>%
  select(all_of(hits)) %>%
  scale()

idx_remove <- which(rowSums(apply(prog_mat, 2, is.na)) > 0)

prog_mat_n <- prog_mat[-idx_remove, ]
features_n <- features[-idx_remove, ]

pc <- prcomp(prog_mat_n, center = FALSE, scale. = FALSE)

features_n$molecular_prog_PC1 <- pc$x[,1]

cr <- cor.test(features_n$molecular_prog_PC1, features_n$ALSFRS_slope)

cat("Correlation of PC1 (3 proteins) and ALSFRS slope:\n")
print(cr)

if (cor(features_n$molecular_prog_PC1, features_n$ALSFRS_slope, use="complete.obs") > 0) {
  features_n$molecular_prog_PC1 <- -features_n$molecular_prog_PC1
}


m <- lm(ALSFRS_slope ~ molecular_prog_PC1 + apoe_status + age + sex, data=features_n)
cat("Linear model ALSFRS slope ~ PC1 + apoe + age + sex\n")
summary(m)




tiff(paste0(here(), '/output/somascan/pics/pc1_vs_alsfrs.tiff'), width = 6, height = 5, units = "in", res = 300)
ggplot(features_n,
       aes(x = molecular_prog_PC1,
           y = ALSFRS_slope)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Molecular progression PC1 (protein-slope axis)",
    y = "ALSFRS-R slope (points per year)",
    colour = "APOE status"
  ) +
  theme_bw()
dev.off()


tiff(paste0(here(), '/output/somascan/pics/pc1_vs_alsfrs_apoe.tiff'), width = 6, height = 5, units = "in", res = 300)
ggplot(features_n,
       aes(x = molecular_prog_PC1,
           y = ALSFRS_slope,
           color=apoe_status)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Molecular progression PC1 (protein-slope axis)",
    y = "ALSFRS-R slope (points per year)",
    colour = "APOE status"
  ) +
  theme_bw()
dev.off()


# features_n <- features_n %>%
# mutate(
#   PC1_tertile = cut(
#     molecular_prog_PC1,
#     breaks = quantile(molecular_prog_PC1, probs = c(0, 1/2, 1), na.rm = TRUE),
#     include.lowest = TRUE,
#     labels = c("Low", "High")
#   )
# ) 
# 
# tiff(paste0(here(), '/output/somascan/pics/pc1_vs_alsfrs_boxplots.tiff'), width = 6, height = 5, units = "in", res = 300)
# ggplot(features_n, aes(x = PC1_tertile, y = ALSFRS_slope)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(width = 0.15, alpha = 0.5) +
#   labs(
#     x = "3-protein molecular progression score",
#     y = "ALSFRS-R slope"
#   ) +
#   theme_classic()
# dev.off()
# 
# 
# clin <- clin %>%
#   rename(GUID = SubjectUID)
# 
# 
# features_n %>%
#   inner_join(clin, by="GUID") %>%
#   group_by(PC1_tertile) %>%
#   summarise(
#     n          = n(),
#     ALSFRS_slope_mean = mean(ALSFRS_slope),
#     APOE4_plus = mean(apoe_status == "APOE4+"),
#     bulbar_onset = mean(onset == "Bulbar"),
#     female       = mean(sex == "Female"),
#     age_onset_mean     = mean(age_onset),
#     bmi = mean(bmi,na.rm=T),
#     head_t = mean(headrb == 1, na.rm=T),
#     concus = mean(concus, na.rm=T),
#     mutations = mean(general_mut, na.rm=T),
#     history = mean(history_tot,na.rm=T),
#     milit = mean(milirb == 1,na.rm=T),
#     city_status = mean(PrimaryRUCA,na.rm=T),
#     drink = mean(driavgtb,na.rm=T),
#     smoking = mean(smokerb==1,na.rm=T),
#     .groups = "drop"
#   )
# 
# for_stats <- features_n %>%
#   filter(PC1_tertile != "Mid") %>% 
#   mutate(PC1_tertile = factor(PC1_tertile, levels = c("Low", "High"))) %>% 
#   inner_join(clin, by="GUID") 
# 
# kruskal.test(PC1_tertile ~ , data=for_stats)
# 
# 
# chisq.test(table(for_stats$PC1_tertile, for_stats$general_mut), simulate.p.value = T)



# features_glm <- features_n %>%
#   inner_join(clin, by="GUID") %>%
#   filter(PC1_tertile != "Mid") %>% 
#   mutate(outcome = ifelse(PC1_tertile == "Low", 1,0) )
# 
# 
# 
# 
# m <- glm(outcome ~ onset, family = binomial, data=features_glm)
# summary(m)


# Checks for stats models

library(dplyr)

qs <- quantile(features_n$ALSFRS_slope, probs = c(1/3, 2/3), na.rm = TRUE)

features_n <- features_n %>%
  mutate(
    prog_group = case_when(
      ALSFRS_slope <= qs[1] ~ "fast",
      ALSFRS_slope >= qs[2] ~ "slow",
      TRUE                  ~ "mid"
    )
  )

# Keep only fast vs slow
dat_fs <- features_n %>%
  filter(prog_group != "mid") %>%
  mutate(
    fast = if_else(prog_group == "fast", 1L, 0L)  # 1 = fast, 0 = slow
  )


m_log <- glm(
  fast ~ molecular_prog_PC1 + apoe_status + sex,
  data   = dat_fs,
  family = binomial
)

cat("Checks: logistic regression on  binomial outcome\n")
summary(m_log)

cat("OR\n")
print(exp(cbind(
  OR  = coef(m_log),
  confint.default(m_log)
)))


#m_lin <- lm(ALSFRS_slope ~ molecular_prog_PC1 + apoe_status + sex + age,
#            data = features_n)
#summary(m_lin)

tiff(paste0(here(), 'output/somascan/pics/pc1_vs_alsfrs_violin.tiff'), width = 6, height = 5, units = "in", res = 300)
dat_fs %>% 
  mutate(prog_group = factor(prog_group, levels = c("slow", "fast"))) %>% 
  ggplot(
       aes(x = prog_group, y = molecular_prog_PC1, fill = prog_group)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("slow" = "grey70", "fast" = "black")) +
  labs(
    x = "Clinical progression group",
    y = "Molecular progression PC1"
  ) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size=16))

dev.off()

set.seed(1)
B <- 2000
beta_pc1 <- numeric(B)

for (b in seq_len(B)) {
  idx <- sample(seq_len(nrow(features_n)), replace = TRUE)
  d   <- features_n[idx, ]
  fit <- lm(ALSFRS_slope ~ molecular_prog_PC1 + apoe_status + sex + onset,
            data = d)
  beta_pc1[b] <- coef(fit)["molecular_prog_PC1"]
}

quantile(beta_pc1, c(0.025, 0.5, 0.975))
cat("Checks: mean beta of PC1 in 2000 permutations\n")
m <- mean(beta_pc1 > 0)
cat(m, "\n")



set.seed(2)
B <- 2000
beta_pc1_logit <- numeric(B)

for (b in seq_len(B)) {
  idx <- sample(seq_len(nrow(dat_fs)), replace = TRUE)
  d   <- dat_fs[idx, ]
  fit <- glm(fast ~ molecular_prog_PC1 + apoe_status + sex + onset,
             data = d, family = binomial)
  beta_pc1_logit[b] <- coef(fit)["molecular_prog_PC1"]
}

OR_boot <- exp(beta_pc1_logit)
quantile(OR_boot, c(0.025, 0.5, 0.975))
m2 <- mean(OR_boot > 1)
cat("Checks: mean OR > 1 in 2000 permutations\n")
cat(m2, "\n")



set.seed(3)
B <- 2000
beta_perm <- numeric(B)

beta_obs <- coef(m_lin)["molecular_prog_PC1"]

for (b in seq_len(B)) {
  d <- features_n
  d$ALSFRS_slope <- sample(d$ALSFRS_slope)  # permute outcome
  fit <- lm(ALSFRS_slope ~ molecular_prog_PC1 + apoe_status + sex + onset,
            data = d)
  beta_perm[b] <- coef(fit)["molecular_prog_PC1"]
}

p_perm <- mean(abs(beta_perm) >= abs(beta_obs))
cat("Checks: mean(abs(beta_perm) >= abs(beta_obs)) in 2000 permutations\n")
cat(p_perm, "\n")



set.seed(4)
B <- 2000
beta_perm_log <- numeric(B)

beta_obs_log <- coef(m_log)["molecular_prog_PC1"]

for (b in seq_len(B)) {
  d <- dat_fs
  d$fast <- sample(d$fast)
  fit <- glm(fast ~ molecular_prog_PC1 + apoe_status + sex + onset,
             data = d, family = binomial)
  beta_perm_log[b] <- coef(fit)["molecular_prog_PC1"]
}

p_perm_log <- mean(abs(beta_perm_log) >= abs(beta_obs_log))
p_perm_log
cat("Checks: mean(abs(beta_perm_log) >= abs(beta_obs_log)) in 2000 permutations\n")
cat(p_perm_log, "\n")



set.seed(5)
B <- 1000
hits <- matrix(0, nrow = B, ncol = length(prot_cols))
colnames(hits) <- paste0(prot_cols, "_slope")

for (b in seq_len(B)) {
  idx <- sample(seq_len(nrow(features_n)), replace = TRUE)
  d   <- features[idx, ]
  
  # re-run 5A-style per-protein regressions on this bootstrap sample
  res_b <- d %>%
    rename(ALSFRS_slope_mm = ALSFRS_slope) %>%
    select(GUID, ALSFRS_slope_mm, ends_with("_slope")) %>%
    tidyr::pivot_longer(cols = ends_with("_slope"),
                 names_to = "protein",
                 values_to = "prot_slope") %>%
    group_by(protein) %>%
    group_modify(~ broom::tidy(lm(prot_slope ~ ALSFRS_slope_mm, data = .x))) %>%
    ungroup() %>%
    filter(term=="ALSFRS_slope_mm")
  
  # mark top k or FDR<0.05 as "hit"
  top <- res_b %>%
    arrange(p.value) %>%
    slice_head(n = 3) %>%
    pull(protein)
  
  hits[b, top] <- 1
}

cat("Checks: frequency of main hits popping up in 1000 permutations\n")
print(colMeans(hits)[order(colMeans(hits), decreasing = T)][1:3])






# med <- median(features_n$ALSFRS_slope, na.rm = TRUE)
# 
# features_n <- features_n %>%
#   mutate(
#     fast_med = if_else(ALSFRS_slope >= med, 1L, 0L)
#   )
# 
# m_log_med <- glm(
#   fast_med ~ molecular_prog_PC1 + apoe_status + sex,
#   data = features_n,
#   family = binomial
# )
# summary(m_log_med)





# Comparison of PC1 in apoe4+ vs apoe4- carriers  
df2 <- df %>%
  filter(TimePoint == 1) %>%
  dplyr::select(SubjectUID, TimePoint, prot_cols) %>%
  dplyr::rename(GUID = SubjectUID) %>%
  dplyr::inner_join(mapper[,meta_cols], by="GUID") %>%
  dplyr::inner_join(meta[, c("Participant_ID", "apoe_status", 'age')], by="Participant_ID") %>%
  dplyr::select(-Participant_ID, -Number_of_Visits) %>%
  filter(SUBJECT_GROUP == "Healthy Control") %>%
  select(P07196, O43175, O75347, apoe_status, SEX, age)



pc_ctrl <- prcomp(df2[, c("P07196", "O43175", "O75347")],
                  center = T, scale. = T)
df2$prot3_PC1 <- pc_ctrl$x[,1]

cat("Comparison of PC1 in apoe4+ vs apoe4- carriers in Control group\n")
print(summary(lm(prot3_PC1 ~ apoe_status + SEX + age, data = df2)))



df3 <- df %>%
  filter(TimePoint == 1) %>%
  dplyr::select(SubjectUID, TimePoint, prot_cols) %>%
  dplyr::rename(GUID = SubjectUID) %>%
  dplyr::inner_join(mapper[,meta_cols], by="GUID") %>%
  dplyr::inner_join(meta[, c("Participant_ID", "apoe_status", 'age')], by="Participant_ID") %>%
  dplyr::select(-Participant_ID, -Number_of_Visits) %>%
  filter(SUBJECT_GROUP == "ALS") %>%
  select(P07196, O43175, O75347, apoe_status, SEX, age)




pc_als <- prcomp(df3[, c("P07196", "O43175", "O75347")],
                  center = T, scale. = T)
df3$prot3_PC1 <- pc_als$x[,1]

cat("Comparison of PC1 in apoe4+ vs apoe4- carriers in ALS group\n")
print(summary(lm(prot3_PC1 ~ apoe_status + SEX + age, data = df3)))


if (mean(df2$prot3_PC1[df2$apoe_status=="APOE4+"]) <
    mean(df2$prot3_PC1[df2$apoe_status=="APOE4-"])) {
  df2$prot3_PC1 <- -df2$prot3_PC1
}

if (mean(df3$prot3_PC1[df3$apoe_status=="APOE4+"]) <
    mean(df3$prot3_PC1[df3$apoe_status=="APOE4-"])) {
  df3$prot3_PC1 <- -df3$prot3_PC1
}


als_df_apoe <- als_df %>%
  inner_join(als_meta[,c("GUID", "apoe_status")])

pc_all <- prcomp(als_df[, c("P07196", "O43175", "O75347")],
                 center = T, scale. = T)
als_df_apoe$prot3_PC1 <- pc_all$x[,1]

m_PC1_mm <- lmer(
  prot3_PC1 ~ time_years * apoe_status + (1 | GUID),
  data = als_df_apoe
)
#summary(m_PC1_mm)




ctrl_plot <- df2 %>%
  mutate(group = "Control") %>%
  select(group, prot3_PC1, apoe_status)

als_plot <- df3 %>%
  mutate(group = "ALS") %>%
  select(group, prot3_PC1, apoe_status)

plot_dat <- bind_rows(ctrl_plot, als_plot)



tiff(paste0(here(), 'output/somascan/pics/pc1_apoe.tiff'), width = 6, height = 5, units = "in", res = 300)
ggplot(plot_dat,
       aes(x = apoe_status, y = prot3_PC1, fill = apoe_status)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  facet_wrap(~ group, scales = "free_y") +
  labs(
    x = "APOE status",
    y = "3-protein PC1 score"
  ) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size=16),
        strip.text.x = element_text(size = 14))

dev.off()



# df_all <- df %>%
#   filter(TimePoint == 1) %>%
#   dplyr::select(SubjectUID, TimePoint, prot_cols) %>%
#   dplyr::rename(GUID = SubjectUID) %>%
#   dplyr::inner_join(mapper[,meta_cols], by="GUID") %>%
#   dplyr::inner_join(meta[, c("Participant_ID", "apoe_status", 'age')], by="Participant_ID") %>%
#   dplyr::select(-Participant_ID, -Number_of_Visits) %>%
#   select(P07196, O43175, O75347, apoe_status, SEX, age, SUBJECT_GROUP)
# 
# 
# pc_als <- prcomp(df_all[, c("P07196", "O43175", "O75347")],
#                  center = T, scale. = T)
# df_all$prot3_PC1 <- pc_als$x[,1]
# 
# summary(lm(prot3_PC1 ~ SUBJECT_GROUP + apoe_status + SUBJECT_GROUP:apoe_status + age + SEX, data = df_all))
# 
# 
# if (mean(df2$prot3_PC1[df2$apoe_status=="APOE4+"]) <
#     mean(df2$prot3_PC1[df2$apoe_status=="APOE4-"])) {
#   df2$prot3_PC1 <- -df2$prot3_PC1
# }