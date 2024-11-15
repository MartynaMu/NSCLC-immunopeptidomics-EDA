#libraries
library("tidyverse")
library("dplyr")
library("writexl")
library("readxl")
library("data.table")
library("readr")

#loading tsv files and creating a column with their file name
# 1. store file names in a vector
data = list.files(path = "data/ids", full.names = TRUE)
# 2. import multiple files stored in data vector
files = lapply(data, read_excel)
# 3. truncate the file names
material <- c("NAT", "Tumor", "PBMC", "Lymph node")
# 4. combine (c) lists (mapply) each argument (df) of list of df and its file name
all_lists <- mapply(c, files, material, SIMPLIFY = FALSE)
# 5. bind rows of all dfs within the list
all <- rbindlist(all_lists, fill = T)
# 6. rename the column
names(all)[3] <- "Material"

#Tidying sample nr
all$Sample = str_extract(all$Sample, "\\d+([d]|[a]|[b])|[C]\\d+")
all <- arrange(all, Sample)

#median
medians <- all %>% group_by(Material) %>% summarise(median=median(ID.nr))


#barplot
ggplot(all) +
  geom_bar(aes(x=Sample, y=ID.nr, fill=Material), stat = "identity", 
           width = 0.5) + 
  theme_bw(base_size = 15) +
  facet_wrap(~ Material, scales = "free") +
  geom_hline(aes(yintercept=median), medians, linetype=2)+
  geom_text(aes(16,median,label = median, vjust = -1.5), data = medians, size = 5)+
  labs(title="Unique peptides identified", y="Count", x="Patient ID", 
       caption = "NAT - normal adjacent tissue \n line - median value per group")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#99E555"))
