#Retrieving filtered peptide ids from DDA search

#libraries
library("tidyverse")
library("readxl")
library("data.table")
library("readr")

#loading tsv files and creating a column with their file name
# 1. store file names in a vector
data = list.files(path = "data/filtered", full.names = TRUE)
# 2. import multiple files stored in data vector
files = lapply(data, read_xlsx)
# 3. truncate the file names
material <- c("Lymph node", "NAT", "PBMC", "Serum", "Tumor")
# 4. combine (c) lists (mapply) each argument (df) of list of df and its file name
all_lists <- mapply(c, files, material, SIMPLIFY = FALSE)
# 5. bind rows of all dfs within the list
all <- rbindlist(all_lists, fill = TRUE)
# 6. rename the column
names(all)[9] <- "Material"

# overlap of all identified peptides
ids_a <- all %>% filter(Material == "NAT") %>% pull(Peptide) %>% unique()
ids_b <- all %>% filter(Material == "Tumor") %>% pull(Peptide) %>% unique()
ids_c <- all %>% filter(Material == "PBMC") %>% pull(Peptide) %>% unique()
ids_d <- all %>% filter(Material == "Lymph node") %>% pull(Peptide) %>% unique()
ids_e <- all %>% filter(Material == "Serum") %>% pull(Peptide) %>% unique()

library(eulerr)
all_overlap <- euler(list("Lung" = ids_a,
           "Tumor" = ids_b,
           "PBMC" = ids_c,
           "LN" = ids_d,
           "Serum" = ids_e)) %>% 
  plot(quantities = TRUE, 
    #   main = "All identified peptides overlap", 
       labels=FALSE, 
       legend = NULL)

ids <- all %>% group_by(Material, Sample) %>% summarise(Count = n())

ids <- ids %>% ungroup() %>% mutate(Patient = case_when(
  grepl(Sample, pattern = "p1[abcdE]|C1$") ~ "1",
  grepl(Sample, pattern = "p2[abcdE]|C2$") ~ "2",
  grepl(Sample, pattern = "p3[abcdE]|C3$") ~ "3",
  grepl(Sample, pattern = "p4[abcdE]|C4$") ~ "4",
  grepl(Sample, pattern = "p5[abcdE]|C5$") ~ "5",
  grepl(Sample, pattern = "p6[abcdE]|C6$") ~ "6",
  grepl(Sample, pattern = "p7[abcdE]|C7$") ~ "7",
  grepl(Sample, pattern = "p8[abcdE]|C8$") ~ "8",
  grepl(Sample, pattern = "p9[abcdE]|C9$") ~ "9",
  grepl(Sample, pattern = "p10|C10") ~ "10",
  grepl(Sample, pattern = "p11|C11") ~ "11",
  grepl(Sample, pattern = "p12|C12") ~ "12",
  grepl(Sample, pattern = "p13|C13") ~ "13",
  grepl(Sample, pattern = "p14|C14") ~ "14",
  grepl(Sample, pattern = "p15|C15") ~ "15",
  grepl(Sample, pattern = "p16|C16") ~ "16",
  grepl(Sample, pattern = "p17|C17") ~ "17",
  grepl(Sample, pattern = "p18|C18") ~ "18"
)) 

# 

# boxplot of peptide identification of all patients per material
library(ggpubr)

ids_p <- ggboxplot(ids, 
          x = "Material", 
          y = "Count",
          color = "Material", 
          palette = "npg",
          add = "jitter",
          order = c("Lymph node", "PBMC", "NAT", "Tumor", "Serum"),
          legend = "none",
          xlab = FALSE,
          ylab = "All identified peptide count")
