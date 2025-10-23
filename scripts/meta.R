


meta <- read.csv("D:/AZ/als/raw/als-metadata/answerals/aals_dataportal_datatable.csv")
demo <- read.csv("D:/AZ/als/raw/als-metadata/answerals/clinical/Demographics.csv")
apoe_vars <- read.csv("D:/AZ/als/proc/apoe_status.csv")



meta <- meta %>%
  select(Participant_ID, NYGC_CGND_ID, SEX, SUBJECT_GROUP,REVISED_EL_ESCORIAL_CRITERIA) %>%
  inner_join(apoe_vars[,c("NYGC_CGND_ID", "variant")], by="NYGC_CGND_ID") %>%
  inner_join(demo[,c("Participant_ID", "age")], by="Participant_ID") %>%
  mutate(apoe_status = ifelse(variant %in% c("24", "34", "44"), "APOE4+", "APOE4-") ) %>%
  rename(sex= SEX)

write.csv(meta, "D:/AZ/als/proc/meta.csv", row.names = F)
