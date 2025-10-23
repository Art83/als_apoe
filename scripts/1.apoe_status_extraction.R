# Step 1. Extracting the APOE status from genomics data


# Genotypes for further analysis below were extracted in bash via:
# bcftools query -r chr19:44908684 -f '%CHROM\t%POS[\t%GT]\n' -H ALS.chr19.vcf.gz > rs429358.genotypes.tsv
# bcftools query -r chr19:44908822 -f '%CHROM\t%POS[\t%GT]\n' -H ALS.chr19.vcf.gz > rs7412.genotypes.tsv



setwd("D:/AZ/als_apoe/")

# Read the files
rs429358 <- read.csv("proc/rs429358.genotypes.tsv", sep='\t', header = TRUE, stringsAsFactors = FALSE, check.names = F)
rs7412   <- read.csv("proc/rs7412.genotypes.tsv",sep='\t', header = TRUE, stringsAsFactors = FALSE, check.names = F)

meta <- read.csv("raw/als-metadata/answerals/aals_dataportal_datatable.csv")
meta_sbs <- meta[,c("Participant_ID", "NYGC_CGND_ID")]


# Remove CHROM and POS (keep only genotypes)
remove <- c("#[1]CHROM", "[2]POS")
gts_429358 <- rs429358[ , !colnames(rs429358) %in% remove]
gts_7412   <- rs7412[ , !colnames(rs7412) %in% remove]

# cleaning the column names
colnames(gts_429358) <- gsub("^\\[[0-9]+\\]|-b38(?=:GT)|:GT$", "", colnames(gts_429358), perl = TRUE)
colnames(gts_7412) <- gsub("^\\[[0-9]+\\]|-b38(?=:GT)|:GT$", "", colnames(gts_7412), perl = TRUE)


# Transpose so samples are rows
gts_429358 <- t(gts_429358)
gts_7412   <- t(gts_7412)

if(!all(rownames(gts_429358) == rownames(gts_7412))){
  stop("Samples in datasets don't match")
}

# Put into a data frame with sample names
geno_df <- data.frame(
  NYGC_CGND_ID = rownames(gts_429358),
  rs429358 = gts_429358,
  rs7412   = gts_7412,
  stringsAsFactors = FALSE
)

suppressPackageStartupMessages(library(dplyr))

geno_df_apoe <- geno_df %>%
  dplyr::inner_join(meta_sbs, by="NYGC_CGND_ID") %>%
  dplyr::mutate(variant = case_when(rs429358 == "0/0" & rs7412 == "1/1" ~ "22",
                             rs429358 == "0/0" & rs7412 == "0/1" ~ "23",
                             rs429358 == "0/1" & rs7412 == "0/1" ~ "24",
                             rs429358 == "0/0" & rs7412 == "0/0" ~ "33",
                             rs429358 == "0/1" & rs7412 == "0/0" ~ "34",
                             rs429358 == "1/1" & rs7412 == "0/0" ~ "44"))


write.csv(geno_df_apoe, "proc/apoe_status1.csv", row.names = F)


demo <- read.csv("raw/als-metadata/answerals/clinical/Demographics.csv")
apoe_vars <- read.csv("proc/apoe_status.csv")



meta <- meta %>%
  select(Participant_ID, NYGC_CGND_ID, SEX, SUBJECT_GROUP,REVISED_EL_ESCORIAL_CRITERIA) %>%
  inner_join(apoe_vars[,c("NYGC_CGND_ID", "variant")], by="NYGC_CGND_ID") %>%
  inner_join(demo[,c("Participant_ID", "age")], by="Participant_ID") %>%
  mutate(apoe_status = ifelse(variant %in% c("24", "34", "44"), "APOE4+", "APOE4-") ) %>%
  rename(sex= SEX)

write.csv(meta, "proc/meta1.csv", row.names = F)
