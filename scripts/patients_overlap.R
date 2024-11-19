# libs
library(tidyverse)
library(ComplexHeatmap)


# overlap of peptides per patient
all <- all %>% ungroup() %>% mutate(Patient = case_when(
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

# lung
ids_a_patient <- list()
for(i in seq.int(1,18,1)){
  ids_a_patient[[i]] <- all %>% filter(Material == "NAT", Patient == deparse(i)) %>% pull(Peptide) %>% unique()
}
names(ids_a_patient) <- seq.int(1,18,1) %>% sapply(deparse) %>% sapply(paste0,"NAT") %>% as.vector()
temp <- table(stack(ids_a_patient))
a_ht <- ComplexHeatmap::Heatmap(temp, show_row_names = FALSE, show_row_dend = FALSE, col = c("white","black"))

ids_b_patient <- list()
for(i in seq.int(1,18,1)){
  ids_b_patient[[i]] <- all %>% filter(Material == "Tumor", Patient == deparse(i)) %>% pull(Peptide) %>% unique()
}
names(ids_b_patient) <- seq.int(1,18,1) %>% sapply(deparse) %>% sapply(paste0,"Tumor") %>% as.vector()
temp <- table(stack(ids_b_patient))
b_ht <- ComplexHeatmap::Heatmap(temp, show_row_names = FALSE, show_row_dend = FALSE, col = c("white","black"))

ids_c_patient <- list()
for(i in seq.int(1,18,1)){
  ids_c_patient[[i]] <- all %>% filter(Material == "PBMC", Patient == deparse(i)) %>% pull(Peptide) %>% unique()
}
names(ids_c_patient) <- seq.int(1,18,1) %>% sapply(deparse) %>% sapply(paste0,"PBMC") %>% as.vector()
temp <- table(stack(ids_c_patient))
c_ht <- ComplexHeatmap::Heatmap(temp, show_row_names = FALSE, show_row_dend = FALSE, col = c("white","black"))

ids_d_patient <- list()
for(i in seq.int(1,18,1)){
  ids_d_patient[[i]] <- all %>% filter(Material == "Lymph node", Patient == deparse(i)) %>% pull(Peptide) %>% unique()
}
names(ids_d_patient) <- seq.int(1,18,1) %>% sapply(deparse) %>% sapply(paste0,"LN") %>% as.vector()
temp <- table(stack(ids_d_patient))
d_ht <- ComplexHeatmap::Heatmap(temp, show_row_names = FALSE, show_row_dend = FALSE, col = c("white","black"))

ids_e_patient <- list()
for(i in seq.int(1,18,1)){
  ids_e_patient[[i]] <- all %>% filter(Material == "Serum", Patient == deparse(i)) %>% pull(Peptide) %>% unique()
}
names(ids_e_patient) <- seq.int(1,18,1) %>% sapply(deparse) %>% sapply(paste0,"Serum") %>% as.vector()
temp <- table(stack(ids_e_patient))
e_ht <- ComplexHeatmap::Heatmap(temp, show_row_names = FALSE, show_row_dend = FALSE, col = c("white","black"))


# combined heatmap with material and patient annotation
# ADD TUMOR SUBTYPE ANNOTATION
library(ComplexHeatmap)
temp <- Reduce(append,x = list(ids_a_patient,ids_b_patient,ids_c_patient,ids_d_patient,ids_e_patient))
temp <- table(stack(temp))


ht_annot <- HeatmapAnnotation(show_legend = TRUE, 
                              show_annotation_name = FALSE,
                              Patient = anno_text(rep(seq.int(1,18),times=5), rot = 0, gp = gpar(fontsize=8)),
                              Material = anno_simple(rep(c("NAT","Tumor", "PBMC", "LN", "Serum"),each=18),
                                                     col = c("NAT" = "#00A087FF", 
                                                             "Tumor" = "#3C5488FF", 
                                                             "PBMC" = "#4DBBD5FF", 
                                                             "LN" = "#E64B35FF", 
                                                             "Serum" = "#F39B7FFF")))

abcd_ht <- ComplexHeatmap::Heatmap(temp, 
                                   show_row_names = FALSE, 
                                   show_row_dend = FALSE, 
                                   show_column_names = FALSE,
                                   col = c("white","black"),
                                   top_annotation = ht_annot,
                                   border = TRUE,
                                   show_heatmap_legend = FALSE)



svg(file="figures/occ-ht.svg",width=8, height=6)
draw(abcd_ht)
dev.off()

# frequency in all samples combined
temp %>% 
