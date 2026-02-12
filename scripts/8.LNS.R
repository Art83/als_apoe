suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(lmerTest))
suppressPackageStartupMessages(library(here))


df <- read.csv(paste0(here(), "/raw/somascan/normalized/csv/SS-2342615_v4.1_anmlSMP__AALS-v1.0-participant-samples.csv"), check.names = F)
mapper <- read.csv(paste0(here(), "/raw/als-metadata/answerals/aals_dataportal_datatable.csv"), check.names = F)
lns_df <- read.csv(paste0(here(), "/raw/als-metadata/answerals/clinical/CNS_Lability_Scale.csv"))
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
  dplyr::select(SubjectUID, TimePoint, all_of(prot_cols)) %>%
  dplyr::rename(GUID = SubjectUID) %>%
  dplyr::inner_join(mapper[,meta_cols], by="GUID") %>%
  dplyr::inner_join(meta[, c("Participant_ID", "age", "apoe_status")], by="Participant_ID") %>%
  dplyr::select(-Participant_ID, -Number_of_Visits)



lns_proc <- lns_df %>%
  select(SubjectUID, Visit_Date, cnslstot) %>%
  rename(GUID = SubjectUID) %>%
  group_by(GUID) %>%
  filter(n() >= 2 ) %>%
  arrange(Visit_Date, .by_group = T) %>%
  mutate(time_years = Visit_Date / 365.25) %>%
  mutate(TimePoint = seq_len(n()))


########### Longitudinal effect of APOE4 ############
apoe_eff <- 
  lns_proc %>% 
  dplyr::inner_join(df1[, c("GUID", "TimePoint", "age", "SEX", "apoe_status")])

m <- lmer(
  cnslstot ~ time_years * apoe_status + age + SEX + (1 | GUID),
  data = apoe_eff,
  REML=T
)

m_alsfrs <- lmer(
  cnslstot ~ time_years + (time_years | GUID),
  data = lns_proc 
)

als_slopes <- ranef(m_alsfrs)$GUID %>%
  tibble::rownames_to_column("GUID") %>%
  rename(
    ALSFRS_int  = `(Intercept)`,
    ALSFRS_slope = time_years
  )

tiff(paste0(here(),'/output/somascan/pics/lns_apoe.tiff'), width = 5, height = 5, units = "in", res = 300)
ggplot(apoe_eff, aes(time_years, cnslstot, colour = apoe_status)) +
  geom_point(position = position_jitter(width = 0.03, height = 0),
             alpha = 0.25, size = 1) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, linewidth = 1) +
  labs(x = "Years from baseline", y = "LNS-LS", colour = NULL) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA)
  )
dev.off()


als_df <- lns_proc %>%
  inner_join(df1[df1$SUBJECT_GROUP == "ALS", c("GUID", "TimePoint", prot_cols)])



als_df[,!colnames(als_df) %in% c("GUID", "TimePoint", "Visit_Date", "cnslstot", "time_years")] <- 
  log2(als_df[,!colnames(als_df) %in% c("GUID", "TimePoint", "Visit_Date", "cnslstot", "time_years")])



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
tiff(paste0(here(), 'output/somascan/pics/forest_lns.tiff'), width = 5, height = 5, units = "in", res = 300)
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
    x = "Effect estimate for LNS-LS slope",
    y = "Protein slope",
    colour = "FDR < 0.05"
  ) +
  theme_bw()
dev.off()



