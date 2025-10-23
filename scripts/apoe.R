

library(dplyr)


df <- read.csv("D:/AZ/als/APOE_2025-04-08_1506.csv")



tt <- df %>%
  select(rsID, CASE.NEUJY426MBU:CASE.NEUWT774BN5) %>%
  filter(rsID %in% c("rs429358", "rs7412")) %>%
  tidyr::pivot_longer(
    cols = -rsID,
    names_to = "Sample",
    values_to = "Genotype"
  ) %>%
  tidyr::pivot_wider(
    names_from = rsID,
    values_from = Genotype
  ) %>%
  mutate(gene = "APOE") %>%
  mutate(als = ifelse(grepl("CASE", Sample), 1, 0)) %>%
  mutate(variant = case_when(rs429358 == "'0/0" & rs7412 == "'1/1" ~ "22",
                             rs429358 == "'0/0" & rs7412 == "'0/1" ~ "23",
                             rs429358 == "'0/1" & rs7412 == "'0/1" ~ "24",
                             rs429358 == "'0/0" & rs7412 == "'0/0" ~ "33",
                             rs429358 == "'0/1" & rs7412 == "'0/0" ~ "34",
                             rs429358 == "'1/1" & rs7412 == "'0/0" ~ "44")) %>%
  rename(sample_id_Answer_ALS = Sample) %>%
  select(sample_id_Answer_ALS, als, gene, variant)



table(tt$als)


write.csv(tt, "D:/AZ/als/als_apor_proc.csv", row.names = F)
