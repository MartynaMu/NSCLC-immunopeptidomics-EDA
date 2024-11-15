#Retrieving filtered peptide ids from DDA search

#libraries
library("tidyverse")
library("dplyr")
library("writexl")
library("readxl")
library("data.table")
library("readr")

#loading tsv files and creating a column with their file name
# 1. store file names in a vector
data = list.files(path = "data/filtered", full.names = TRUE)
# 2. import multiple files stored in data vector
files = lapply(data, read_xlsx)
# 3. truncate the file names
material <- c("Lymph node", "NAT", "PBMC", "Tumor")
# 4. combine (c) lists (mapply) each argument (df) of list of df and its file name
all_lists <- mapply(c, files, material, SIMPLIFY = FALSE)
# 5. bind rows of all dfs within the list
all <- rbindlist(all_lists, fill = T)
# 6. rename the column
names(all)[9] <- "Material"

AB <- filter(all, Material == c("NAT", "Tumor"))
CD <- filter(all, Material == c("PBMC", "Lymph node"))

#_______________________________________________________________________________
#altogether
#sample-wise
ggplot(AB, aes(Peptide.Length)) + 
  geom_bar(width = 0.7) + 
  facet_grid(Material ~ Sample, scales = "free") +
  theme_bw(base_size = 15) +
  labs(title = "Length distribution of identified peptides", 
       subtitle = "Acid elution, all patients",
       x = "Peptide length (AA)")

#___________________________________________________________________________
#Retrieving 8-14AA for every condition
short_peps <- filtered_slim %>% filter((Peptide.Length <= 14) & (Peptide.Length >= 8))
write_xlsx(short_peps, "data/short_peps.xlsx")
