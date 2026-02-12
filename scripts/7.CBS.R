suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(lmerTest))
suppressPackageStartupMessages(library(here))




df <- read.csv(paste0(here(), "/raw/somascan/normalized/csv/SS-2342615_v4.1_anmlSMP__AALS-v1.0-participant-samples.csv"), check.names = F)
mapper <- read.csv(paste0(here(), "/raw/als-metadata/answerals/aals_dataportal_datatable.csv"), check.names = F)
cbs <- read.csv(paste0(here(),"/raw/als-metadata/answerals/clinical/ALS_CBS.csv"))
meta <- read.csv(paste0(here(),"/proc/meta.csv"))
annotation <- read.csv(paste0(here(), "/raw/somascan/documents/annotation.csv"))
proteins <- read.csv(paste0(here(), "/output/somascan/mi_score_somascan.csv"))


prot_cols <- proteins$UniProt
misc_columns <- c("Barcode2d", "SubjectUID", "TimePoint")
colnames(df)[!colnames(df) %in% misc_columns] <- annotation$UniProt.ID[match(colnames(df)[!colnames(df) %in% misc_columns], 
                                                                             annotation$SeqId)]

df <- df[,colnames(df) != "" & !is.na(colnames(df))]

meta_cols <- c("Participant_ID", "GUID", "Number_of_Visits", "SEX", "SUBJECT_GROUP", "Site_of_Onset", "AGE_AT_SYMPTOM_ONSET",
               "AGE_AT_DEATH","ALSFRS_R_Baseline_Value","ALSFRS_R_Latest_Value","ALSFRS_R_PROGRESSION_SLOPE",   
               "CBS_Baseline_Value", "CBS_Latest_Value", "COGNITIVE_PROGRESSION_SLOPE")
df1 <- df %>%
  dplyr::select(SubjectUID, TimePoint, all_of(prot_cols)) %>%
  dplyr::rename(GUID = SubjectUID) %>%
  dplyr::inner_join(mapper[,meta_cols], by="GUID") %>%
  dplyr::inner_join(meta[, c("Participant_ID", "age", "apoe_status")], by="Participant_ID") %>%
  dplyr::select(-Participant_ID, -Number_of_Visits)


cbs_proc <- cbs %>%
  filter(Child_Name == "ALS CBS") %>%
  select(SubjectUID, Visit_Date, conscr, iniscr, trkscr, cbstot, attscr) %>%
  rename(GUID = SubjectUID) %>%
  group_by(GUID) %>%
  filter(n() >= 2 ) %>%
  arrange(Visit_Date, .by_group = T) %>%
  mutate(time_years = Visit_Date / 365.25) %>%
  mutate(TimePoint = seq_len(n()))



########### Longitudinal effect of APOE4 ############
apoe_eff <- 
  cbs_proc %>% 
  filter(cbstot != 0) %>% 
  dplyr::inner_join(df1[, c("GUID", "TimePoint", "age", "SEX", "apoe_status")])

m <- lmer(
  cbstot ~ time_years * apoe_status + age + SEX + (1 | GUID),
  data = apoe_eff,
  REML=T
)



tiff(paste0(here(), 'output/somascan/pics/cbs_apoe.tiff'), width = 5, height = 5, units = "in", res = 300)
ggplot(apoe_eff, aes(time_years, cbstot, colour = apoe_status)) +
  geom_point(position = position_jitter(width = 0.03, height = 0),
             alpha = 0.25, size = 1) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, linewidth = 1) +
  labs(x = "Years from baseline", y = "CBS", colour = NULL) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA)
  )
dev.off()



m_alsfrs <- lmer(
  cbstot ~ time_years + (time_years | GUID),
  data = cbs_proc 
)

als_slopes <- ranef(m_alsfrs)$GUID %>%
  tibble::rownames_to_column("GUID") %>%
  rename(
    ALSFRS_int  = `(Intercept)`,
    ALSFRS_slope = time_years
  )


als_df <- cbs_proc %>%
  inner_join(df1[df1$SUBJECT_GROUP == "ALS", c("GUID", "TimePoint", prot_cols)])



als_df[,!colnames(als_df) %in% c("GUID", "TimePoint", "Visit_Date", "cbstot", "conscr", "iniscr", "trkscr", "attscr", "time_years")] <- 
  log2(als_df[,!colnames(als_df) %in% c("GUID", "TimePoint", "Visit_Date", "cbstot", "conscr", "iniscr", "trkscr","attscr", "time_years")])



prot_summ <- als_df %>%
  group_by(GUID) %>%
  group_modify(~ {
    dat <- .x
    # baseline = first visit (time 0)
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
    age = first(age),
    onset = first(Site_of_Onset),
    age_onset = first(AGE_AT_SYMPTOM_ONSET),
    age_death = first(AGE_AT_DEATH),
    .groups = "drop"
  )


features <- als_slopes %>%
  inner_join(prot_summ, by = "GUID") %>%
  inner_join(als_meta, by = "GUID")



res_slope <- features %>%
  slice(-176) %>% 
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
tiff(paste0(here(),'/output/somascan/pics/forest_cbs.tiff'), width = 5, height = 5, units = "in", res = 300)
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
    x = "Effect estimate for CBS slope",
    y = "Protein slope",
    colour = "FDR < 0.05"
  ) +
  theme_bw()
dev.off()

### 1B. Volcano plot version (optional)

# ggplot(res_slope_alsfrs,
#        aes(x = estimate, y = -log10(p_fdr), colour = signif)) +
#   geom_point() +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
#   scale_colour_manual(values = c("FALSE" = "grey50", "TRUE" = "black")) +
#   labs(
#     x = "Effect of ALSFRS-R slope on protein slope",
#     y = expression(-log[10]("FDR-adjusted p-value")),
#     colour = "FDR < 0.05"
#   ) +
#   theme_bw()


## Long per-patient table
long_slope <- features %>%
  slice(-176) %>% 
  rename(ALSFRS_slope_mm = ALSFRS_slope) %>%
  select(GUID, ALSFRS_slope_mm, apoe_status, ends_with("_slope")) %>%
  tidyr::pivot_longer(
    cols = ends_with("_slope"),
    names_to = "protein",
    values_to = "prot_slope"
  )

## Proteins that passed FDR
hits <- c("P56704_slope", "P78325_slope")

res_for_hits <- long_slope %>%
  filter(protein %in% hits)

## Scatter + regression lines, faceted by protein
ggplot(res_for_hits,
       aes(x = ALSFRS_slope_mm,
           y = prot_slope)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ protein, scales = "free_y") +
  labs(
    x = "ALSFRS-R slope (points per year)",
    y = "Protein slope (log2 intensity per year)",
    colour = "APOE status"
  ) +
  theme_bw()

ggplot(res_for_hits,
       aes(x = ALSFRS_slope_mm,
           y = prot_slope)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ apoe_status, scales = "free_y") +
  labs(
    x = "ALSFRS-R slope (points per year)",
    y = "Protein slope (log2 intensity per year)",
    colour = "APOE status"
  ) +
  theme_bw()



res_int <- features %>%
  rename(ALSFRS_slope_mm = ALSFRS_slope) %>%
  select(GUID, ALSFRS_slope_mm, apoe_status, ends_with("_slope")) %>%
  tidyr::pivot_longer(
    cols = ends_with("_slope"),
    names_to = "protein",
    values_to = "prot_slope"
  ) %>%
  group_by(protein) %>%
  group_modify(~ {
    fit <- lm(prot_slope ~ ALSFRS_slope_mm * apoe_status, data = .x)
    broom::tidy(fit)
  }) %>%
  ungroup() %>%
  filter(term == "ALSFRS_slope_mm:apoe_statusAPOE4+") %>%
  mutate(p_fdr = p.adjust(p.value, method = "fdr"))


# 3-protein molecular progression score
prog_mat <- features %>%
  select(hits) %>%
  scale()

idx_remove <- which(rowSums(apply(prog_mat, 2, is.na)) > 0)

prog_mat_n <- prog_mat[-idx_remove, ]
features_n <- features[-idx_remove, ]

pc <- prcomp(prog_mat_n, center = FALSE, scale. = FALSE)

features_n$molecular_prog_PC1 <- pc$x[,1]



cor.test(features_n$molecular_prog_PC1, features_n$ALSFRS_slope)



m <- lm(ALSFRS_slope ~ molecular_prog_PC1 + apoe_status + age + sex, data=features_n)
summary(m)



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






# checks
# qs <- quantile(features_n$ALSFRS_slope, probs = c(1/3, 2/3), na.rm = TRUE)
# 
# features_n <- features_n %>%
#   mutate(
#     prog_group = case_when(
#       ALSFRS_slope <= qs[1] ~ "slow",
#       ALSFRS_slope >= qs[2] ~ "fast",
#       TRUE                  ~ "mid"
#     )
#   )

# Keep only fast vs slow
# dat_fs <- features_n %>%
#   filter(prog_group != "mid") %>%
#   mutate(
#     fast = if_else(prog_group == "fast", 1L, 0L)  # 1 = fast, 0 = slow
#   )
# table(dat_fs$prog_group)




# m_log <- glm(
#   fast ~ molecular_prog_PC1 + apoe_status + sex,
#   data   = dat_fs,
#   family = binomial
# )
# summary(m_log)



# exp(cbind(
#   OR  = coef(m_log),
#   confint.default(m_log)   # use confint(m_log) for profile CI, but slower
# ))




# m_lin <- lm(ALSFRS_slope ~ molecular_prog_PC1 + apoe_status + sex + onset,
#             data = features)
# summary(m_lin)
# 
# 
# 
# ggplot(dat_fs,
#        aes(x = prog_group, y = molecular_prog_PC1, fill = prog_group)) +
#   geom_violin(trim = FALSE, alpha = 0.5) +
#   geom_boxplot(width = 0.15, outlier.shape = NA) +
#   geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
#   scale_fill_manual(values = c("slow" = "grey70", "fast" = "black")) +
#   labs(
#     x = "Clinical progression group",
#     y = "Molecular progression PC1"
#   ) +
#   theme_bw() +
#   theme(legend.position = "none")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# set.seed(1)
# B <- 2000
# beta_pc1 <- numeric(B)
# 
# for (b in seq_len(B)) {
#   idx <- sample(seq_len(nrow(features_n)), replace = TRUE)
#   d   <- features_n[idx, ]
#   fit <- lm(ALSFRS_slope ~ molecular_prog_PC1 + apoe_status + sex + onset,
#             data = d)
#   beta_pc1[b] <- coef(fit)["molecular_prog_PC1"]
# }
# 
# quantile(beta_pc1, c(0.025, 0.5, 0.975))
# mean(beta_pc1 > 0)   # proportion > 0
# 
# 
# 
# set.seed(2)
# B <- 2000
# beta_pc1_logit <- numeric(B)
# 
# for (b in seq_len(B)) {
#   idx <- sample(seq_len(nrow(dat_fs)), replace = TRUE)
#   d   <- dat_fs[idx, ]
#   fit <- glm(fast ~ molecular_prog_PC1 + apoe_status + sex + onset,
#              data = d, family = binomial)
#   beta_pc1_logit[b] <- coef(fit)["molecular_prog_PC1"]
# }
# 
# OR_boot <- exp(beta_pc1_logit)
# quantile(OR_boot, c(0.025, 0.5, 0.975))
# mean(OR_boot > 1)
# 
# 
# 
# set.seed(3)
# B <- 2000
# beta_perm <- numeric(B)
# 
# beta_obs <- coef(m_lin)["molecular_prog_PC1"]
# 
# for (b in seq_len(B)) {
#   d <- features
#   d$ALSFRS_slope <- sample(d$ALSFRS_slope)  # permute outcome
#   fit <- lm(ALSFRS_slope ~ molecular_prog_PC1 + apoe_status + sex + onset,
#             data = d)
#   beta_perm[b] <- coef(fit)["molecular_prog_PC1"]
# }
# 
# p_perm <- mean(abs(beta_perm) >= abs(beta_obs))
# p_perm
# 
# 
# 
# set.seed(4)
# B <- 2000
# beta_perm_log <- numeric(B)
# 
# beta_obs_log <- coef(m_log)["molecular_prog_PC1"]
# 
# for (b in seq_len(B)) {
#   d <- dat_fs
#   d$fast <- sample(d$fast)
#   fit <- glm(fast ~ molecular_prog_PC1 + apoe_status + sex + onset,
#              data = d, family = binomial)
#   beta_perm_log[b] <- coef(fit)["molecular_prog_PC1"]
# }
# 
# p_perm_log <- mean(abs(beta_perm_log) >= abs(beta_obs_log))
# p_perm_log
# 
# 
# 
# set.seed(5)
# B <- 1000
# hits <- matrix(0, nrow = B, ncol = length(prot_cols))
# colnames(hits) <- paste0(prot_cols, "_slope")
# 
# for (b in seq_len(B)) {
#   idx <- sample(seq_len(nrow(features)), replace = TRUE)
#   d   <- features[idx, ]
#   
#   # re-run 5A-style per-protein regressions on this bootstrap sample
#   res_b <- d %>%
#     rename(ALSFRS_slope_mm = ALSFRS_slope) %>%
#     select(GUID, ALSFRS_slope_mm, ends_with("_slope")) %>%
#     tidyr::pivot_longer(cols = ends_with("_slope"),
#                         names_to = "protein",
#                         values_to = "prot_slope") %>%
#     group_by(protein) %>%
#     group_modify(~ broom::tidy(lm(prot_slope ~ ALSFRS_slope_mm, data = .x))) %>%
#     ungroup() %>%
#     filter(term=="ALSFRS_slope_mm")
#   
#   # mark top k or FDR<0.05 as "hit"
#   top <- res_b %>%
#     arrange(p.value) %>%
#     slice_head(n = 3) %>%
#     pull(protein)
#   
#   hits[b, top] <- 1
# }
# 
# colMeans(hits)[order(colMeans(hits), decreasing = T)][1:3]
# 
# 
# 
# 
# 
# 
# med <- median(features$ALSFRS_slope, na.rm = TRUE)
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



