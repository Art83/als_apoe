library(dplyr)

df_1 <- read.csv("D:/AZ/als/raw/als-metadata/answerals/clinical/AALSDXFX.csv")
meta <- read.csv("D:/AZ/als/proc/meta.csv")


df_1_proc <- df_1 %>%
  select(-Form_Name, -Visit_Date, -Visit_Name, -alsdx2, -alsdx3, -alsdxdt, -blbelmn, -lleelmn, -lueelmn, -rleelmn,
         -rueelmn, -trnkclmn, -trnkcumn, -trnkelmn, -elescrlr) %>%
  mutate(across(alsdx1:ruecumn, ~case_when(
    .x == 90 ~ NA_integer_,
    .x == 2 ~ 0,
    TRUE ~ .x
  )) )


apply(apply(df_1_proc, 2, is.na),2, sum)


df_stat <- df_1_proc %>%
  inner_join(meta[,c("Participant_ID", "apoe_status")], by="Participant_ID") %>%
  select(-Participant_ID, -SubjectUID) %>%
  na.omit()



for(i in 1:(ncol(df_stat)-1)){
  chsk <- chisq.test(table(df_stat$apoe_status, df_stat[,i]), simulate.p.value = T)
  cat(colnames(df_stat)[i], chsk$p.value, "\n")
}


df_2 <- read.csv("D:/AZ/als/raw/als-metadata/answerals/clinical/AALSHXFX.csv")

df_2_proc <- df_2 %>%
  select(-Form_Name, -Visit_Date, -Visit_Name, -alsdxloc, -diagdt, -hxotsp, -onsetdt)


df_stat <- df_2_proc %>%
  inner_join(meta[,c("Participant_ID", "apoe_status")], by="Participant_ID") %>%
  select(-Participant_ID, -SubjectUID) %>%
  na.omit()



df_3 <- read.csv("D:/AZ/als/raw/als-metadata/answerals/clinical/ALS_CBS.csv")


df_3_self <- df_3 %>%
  filter(Child_Name == "ALS CBS" & Visit_Name == "ANSWER-ALS Screening Visit" ) %>%
  select(-Child_Name, -Visit_Name, -Visit_Date, -carbeh:-carbeh15, -cgcuranx:-cgqrel, 
         -cbsdn:-cbswrite, -ini01:-ini20, -sourcesp, -source)



df_stat <- df_3_self %>%
  inner_join(meta[,c("Participant_ID", "apoe_status")], by="Participant_ID") %>%
  select(-Participant_ID, -SubjectUID) %>%
  na.omit() %>%
  mutate(attb1 = ifelse(attb1 == 5, 1, 0),
         attb2 = ifelse(attb2 == 7, 1, 0) )



cont_vars <- c("attscr", "conscr", "iniscr", "trkscr")
cat_vars <- colnames(df_stat)[!colnames(df_stat) %in% c("apoe_status", cont_vars)]

for(i in cat_vars){
  chsk <- chisq.test(table(df_stat$apoe_status, df_stat[,i]), simulate.p.value = T)
  cat(i, chsk$p.value, "\n")
}


t.test(df_stat$trkscr[df_stat$apoe_status == "APOE4+"],
            df_stat$trkscr[df_stat$apoe_status == "APOE4-"])


chisq.test(table(df_stat$apoe_status, df_stat$trkscr), simulate.p.value = T)



