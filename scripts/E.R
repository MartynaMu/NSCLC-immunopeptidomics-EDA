#Retrieving filtered peptide ids from DDA search
# serum peptides analysis
#libraries
library("tidyverse")
library("writexl")
library("data.table")
library("readr")
library("ggplot2")


#loading tsv files and creating a column with their file name
# 1. store file names in a vector
data = list.files(path = "data/E", full.names = TRUE,pattern = ".tsv")
# 2. import multiple files stored in data vector
files = lapply(data, read.delim)
# 3. truncate the file names
data <- str_extract(data, "[:digit:]+")
data <- paste0("p",data,"E")
# 4. combine (c) lists (mapply) each argument (df) of list of df and its file name
all_lists <- mapply(c, files, data, SIMPLIFY = FALSE)
# 5. bind rows of all dfs within the list
all <- rbindlist(all_lists, fill = T)
all <- select(all, !c(Protein.Start, Protein.End))
# 6. rename the column
names(all)[18] <- "Sample"

#filtering out 1% FDR, contams, iRTs and revs
ftr <- function(x) {
  x %>% filter(Probability >= 0.99) %>%
    filter(!grepl("Biognosys", Protein)) %>%
    filter(!grepl("contam", Protein)) %>%
    filter(!grepl("rev", Protein))
}

filtered <- ftr(all)
ids <- filtered %>% group_by(Sample) %>% summarise(ID.nr = n())

#saving into summary altogether
write_xlsx(ids, "data/ids_e.xlsx")

#_______________________________________________________________________________
#retrieving peptide length distribution
filtered_slim <- filtered[,c(1,4,11,12,13,14,15,18)]
write_xlsx(filtered_slim, "data/filtered/serum.xlsx")

#altogether
ggplot(filtered_slim, aes(Peptide.Length)) + 
  geom_bar() + 
  theme_bw(base_size = 15) +
  labs(title = "Length distribution of all identified peptides", 
       subtitle = "Acid elution of serum, all patients", 
       x = "Peptide length (AA)")

#sample-wise
ggplot(filtered_slim, aes(Peptide.Length)) + 
  geom_bar(width = 0.7) + 
  facet_wrap(~ Sample, nrow = 2, scales = "free") +
  theme_bw(base_size = 15) +
  labs(title = "Length distribution of identified peptides", 
       subtitle = "Acid elution of serum, all patients",
       x = "Peptide length (AA)")

#___________________________________________________________________________
#Retrieving 8-14AA for every condition
short_peps <- filtered_slim %>% filter((Peptide.Length <= 14) & (Peptide.Length >= 8))
clipr::write_clip(unique(short_peps$Peptide))
write_xlsx(short_peps, "data/short_peps.xlsx")
