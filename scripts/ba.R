libs <- c("tidyverse", "ggplot2", "ggpubr", "readxl", "smartPep", "eulerr")
lapply(libs, require, character.only = TRUE)

ba_cancer <- read_xlsx("data/cancer_all_peptides.xlsx", col_names = FALSE)
ba_cancer <- cleanup_ranks(ba_cancer)
ba_cancer[2:13] <- lapply(ba_cancer[2:13], as.numeric)
ba_cancer <- binding_category(ba_cancer)
ba_cancer_L <- pivot_longer(ba_cancer, 2:14,
                        values_to = "Binding category", 
                        names_to = "HLA supertype")

ba_lung <- read_xlsx("data/lung_all_peptides.xlsx", col_names = FALSE)
ba_lung <- cleanup_ranks(ba_lung)
ba_lung[2:13] <- lapply(ba_lung[2:13], as.numeric)
ba_lung <- binding_category(ba_lung)
ba_lung_L <- pivot_longer(ba_lung, 2:14,
                        values_to = "Binding category", 
                        names_to = "HLA supertype")

ba_pbmc <- read_xlsx("data/pbmc_all_peptides.xlsx", col_names = FALSE)
ba_pbmc <- cleanup_ranks(ba_pbmc)
ba_pbmc[2:13] <- lapply(ba_pbmc[2:13], as.numeric)
ba_pbmc <- binding_category(ba_pbmc)
ba_pbmc_L <- pivot_longer(ba_pbmc, 2:14,
                          values_to = "Binding category", 
                          names_to = "HLA supertype")

ba_ln <- read_xlsx("data/ln_all_peptides.xlsx", col_names = FALSE)
ba_ln <- cleanup_ranks(ba_ln)
ba_ln[2:13] <- lapply(ba_ln[2:13], as.numeric)
ba_ln <- binding_category(ba_ln)
ba_ln_L <- pivot_longer(ba_ln, 2:14,
                          values_to = "Binding category", 
                          names_to = "HLA supertype")

ba_s <- read_xlsx("data/serum_all_peptides.xlsx", col_names = FALSE)
ba_s <- cleanup_ranks(ba_s)
ba_s[2:13] <- lapply(ba_s[2:13], as.numeric)
ba_s <- binding_category(ba_s)
ba_s_L <- pivot_longer(ba_s, 2:14,
                        values_to = "Binding category", 
                        names_to = "HLA supertype")

ba <- bind_rows(ba_lung_L, ba_cancer_L, ba_pbmc_L, ba_ln_L, ba_s_L, .id = "Material")
ba <- ba %>% mutate(Material = case_when(
  Material == "1" ~ "NAT",
  Material == "2" ~ "Tumor",
  Material == "3" ~ "PBMC", 
  Material == "4" ~ "LN",
  Material == "5" ~ "Serum"
))

filtr <- ba %>% pull(Peptide) %>% unique()
temp <- all %>% filter(!(Peptide %in% filtr)) %>% select(Material, Peptide) %>% distinct()
temp$`HLA supertype` <- "Any HLA allele"
temp$`Binding category` <- "Non-binder"
ba <- bind_rows(ba,temp)

ba_p <- ba %>% filter(`HLA supertype` == "Any HLA allele") %>%
  group_by(Material, `Binding category`) %>%
  summarise(Count=n()) %>%
  ggbarplot(x="Material", 
            y="Count",
            color = "Binding category",
            fill = "Binding category",
            palette = c("black", "gray"),
            label = TRUE,  
            lab.col = "white",
            lab.pos = "out",
            position = position_fill(.02),
            xlab = FALSE,
            ylab = "Proportion")+
  theme_pubr(border = TRUE)+
  theme(axis.text.x = element_text(vjust = -7),
        plot.margin = margin(b=1,unit="cm"),
        axis.title.x = element_blank())+
  annotate("rect", xmin=x-.35,xmax=x+.35, ymin=-.2, ymax=-.1, 
           fill = rep(c("#00A087FF","#3C5488FF","#4DBBD5FF","#E64B35FF","#F39B7FFF"),2), 
           color = "black")+
  coord_cartesian(ylim=c(0,1.02), clip = "off")
ggsave("figures/binders-proportion.svg", device="svg", width=4, height=3, dpi=100, scale=1)

lung_binders <- ba_lung %>% filter(`Any HLA allele` == 'Binder') %>% pull(Peptide)
cancer_binders <- ba_cancer %>% filter(`Any HLA allele` == 'Binder') %>% pull(Peptide)
pbmc_binders <- ba_pbmc %>% filter(`Any HLA allele` == 'Binder') %>% pull(Peptide)
ln_binders <- ba_ln %>% filter(`Any HLA allele` == 'Binder') %>% pull(Peptide)
s_binders <- ba_s %>% filter(`Any HLA allele` == 'Binder') %>% pull(Peptide)

euler_binders <- euler(list(
  "NAT" = lung_binders, 
  "Tumor" = cancer_binders,
  "PBMCs" = pbmc_binders,
  "LN" = ln_binders,
  "Serum" = s_binders))

binder_overlap <- plot(euler_binders, 
                       quantities = TRUE, 
    # main = "Predicted HLA-I-binders overlap",
     labels=FALSE, 
     legend = list(side="right"),
    fills = list(fill = c("#00A087FF","#3C5488FF","#4DBBD5FF","#E64B35FF","#F39B7FFF"), alpha = .7),
    edges = FALSE)
binder_overlap$vp$height <- unit(0.9, units = "npc")

library(patchwork)
layout <- "
AAABBB
AAABBB
CCDDDD
CCDDDD
CCDDDD
"
wrap_plots(list(ba_p, ids_p, all_overlap, binder_overlap), design = layout)
ggsave("figures/multi-panel.svg", device="svg", width=9, height=7, dpi=100, scale=1)



