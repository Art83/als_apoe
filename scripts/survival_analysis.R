meta <- read.csv("D:/AZ/als/raw/als-metadata/answerals/aals_dataportal_datatable.csv")

apoe_vars <- read.csv("D:/AZ/als/proc/apoe_status.csv")



meta <- meta %>%
  select(Participant_ID, NYGC_CGND_ID, SEX, SUBJECT_GROUP,REVISED_EL_ESCORIAL_CRITERIA) %>%
  inner_join(apoe_vars[,c("NYGC_CGND_ID", "variant")], by="NYGC_CGND_ID")



surv_curves <- meta %>%
  select(Participant_ID, NYGC_CGND_ID, AGE_AT_SYMPTOM_ONSET, AGE_AT_DEATH) %>%
  inner_join(apoe_vars[,c("NYGC_CGND_ID", "variant")], by="NYGC_CGND_ID") %>%
  mutate(apoe_status = ifelse(variant %in% c("44", "34"), 1, 0)) %>%
  mutate(time = AGE_AT_DEATH - AGE_AT_SYMPTOM_ONSET,
         event = 1)


surv_curves <- surv_curves[!is.na(surv_curves$AGE_AT_DEATH) & !is.na(surv_curves$AGE_AT_SYMPTOM_ONSET), ]

surv_curves$apoe_status <- ifelse(surv_curves$variant == "44", 2, 
                                  ifelse(surv_curves$variant == "34", 1, 0))

surv_curves1 <- surv_curves[surv_curves$apoe_status %in% c(0, 2), ]


surv_obj <- survival::Surv(time = surv_curves1$time, event = surv_curves1$event)

km_fit <- survival::survfit(surv_obj ~ apoe_status, data = surv_curves1)

library(survminer)

ggsurvplot(
  km_fit,
  data = surv_curves1,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  palette = c("#E7B800","#2E9FDF"),
  legend.labs = c("APOE 0", "APOE 1")
)


cox_fit <- survival::coxph(surv_obj ~ apoe_status, data = surv_curves1)
summary(cox_fit)



survival::survdiff(Surv(time, event) ~ apoe_status, data = surv_curves1, rho = 1) 