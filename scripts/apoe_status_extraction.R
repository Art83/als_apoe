
# Genotypes for analysis below were extracted via:
# bcftools query -r chr19:44908684 -f '%CHROM\t%POS[\t%GT]\n' -H ALS.chr19.vcf.gz > rs429358.genotypes.tsv
# bcftools query -r chr19:44908822 -f '%CHROM\t%POS[\t%GT]\n' -H ALS.chr19.vcf.gz > rs7412.genotypes.tsv



# Read the files
rs429358 <- read.csv("F:/container/apoe_als/rs429358.genotypes.tsv", sep='\t', header = TRUE, stringsAsFactors = FALSE, check.names = F)
rs7412   <- read.csv("F:/container/apoe_als/rs7412.genotypes.tsv",sep='\t', header = TRUE, stringsAsFactors = FALSE, check.names = F)

meta <- read.csv("D:/AZ/als/raw/als-metadata/answerals/aals_dataportal_datatable.csv")
meta <- meta[,c("Participant_ID", "NYGC_CGND_ID")]



# Remove CHROM and POS (keep only genotypes)
gts_429358 <- rs429358[ , -(1:2)]
gts_7412   <- rs7412[ , -(1:2)]


colnames(gts_429358) <- gsub("^\\[[0-9]+\\]|-b38(?=:GT)|:GT$", "", colnames(gts_429358), perl = TRUE)
colnames(gts_7412) <- gsub("^\\[[0-9]+\\]|-b38(?=:GT)|:GT$", "", colnames(gts_7412), perl = TRUE)

# Transpose so samples are rows
gts_429358 <- t(gts_429358)
gts_7412   <- t(gts_7412)

all(rownames(gts_429358) == rownames(gts_7412))

# Put into a data frame with sample names
geno_df <- data.frame(
  NYGC_CGND_ID = rownames(gts_429358),
  rs429358 = gts_429358[,1],
  rs7412   = gts_7412[,1],
  stringsAsFactors = FALSE
)

library(dplyr)

geno_df_apoe <- geno_df %>%
  inner_join(meta, by="NYGC_CGND_ID") %>%
  mutate(variant = case_when(rs429358 == "0/0" & rs7412 == "1/1" ~ "22",
                             rs429358 == "0/0" & rs7412 == "0/1" ~ "23",
                             rs429358 == "0/1" & rs7412 == "0/1" ~ "24",
                             rs429358 == "0/0" & rs7412 == "0/0" ~ "33",
                             rs429358 == "0/1" & rs7412 == "0/0" ~ "34",
                             rs429358 == "1/1" & rs7412 == "0/0" ~ "44"))


gt <- readClipboard()  
vars <- readClipboard()


for_test <- data.frame(Participant_ID = gt,
                       variant1 = vars)


geno_df_apoe$Participant_ID <- make.names(geno_df_apoe$Participant_ID)


test <- geno_df_apoe %>%
  inner_join(for_test, by="Participant_ID")


all(test$variant == test$variant1)



write.csv(geno_df_apoe, "D:/AZ/als/proc/apoe_status.csv", row.names = F)
