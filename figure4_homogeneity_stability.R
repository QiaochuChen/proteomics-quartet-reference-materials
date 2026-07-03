## Title: Homogeneity & Stability.
## Author: Qiaochu Chen
## Date: Jun 23rd, 2026

library(readxl)
library(data.table)
library(rstatix)
library(pbapply)
library(parallel)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(lme4)
library(RColorBrewer)


options(digits = 4)


## Figure 4c & Supp Figure 10a: 标称特性值的均匀性检验结果 -----------------------
all_files <- list.files("./results/tables", full.names = TRUE)
homo_files <- all_files[grepl("5_qualiprop.+\\.rds", all_files)]
homo_labels <- c("DDA", "DIA", "Orbitrap", "TOF", "All")

stat_pep_tables <- pblapply(3:5, function(i) {
  
  homo_tables <- readRDS(homo_files[i])
  
  stat_sw <- homo_tables %>%
    rbindlist(., idcol = "group") %>%
    dplyr::rename(p = `fisher.p`) %>%
    mutate_at("p", ~ ifelse(is.na(.), 0, .)) %>%
    group_by(peptide_sequence, group) %>%
    summarise_at("p", max) %>%
    mutate(pass = ifelse((p > .05) & (!is.na(p)), "Yes", "No")) %>%
    reshape2::dcast(., group ~ pass, value.var = "peptide_sequence",
                    fun.aggregate = length) %>%
    mutate(n = Yes + No) %>%
    mutate(label = homo_labels[i])
  
  return(stat_sw)
})
stat_pep_df <- stat_pep_tables %>%
  rbindlist %>%
  filter(!group %in% c("HeLa","HEK293T")) %>%
  reshape2::melt(., id = c(1, 4:5))

pep_thres <- stat_pep_df %>%
  filter(variable %in% "Yes") %>%
  mutate(prop = value/n * 100)

stat_pro_tables <- pblapply(3:5, function(i) {
  
  homo_tables <- readRDS(homo_files[i])
  
  stat_sw <- homo_tables %>%
    rbindlist(., idcol = "group") %>%
    dplyr::rename(p = `fisher.p`) %>%
    mutate_at("p", ~ ifelse(is.na(.), 0, .)) %>%
    group_by(protein_id, group) %>%
    summarise_at("p", max) %>%
    mutate(pass = ifelse((p > .05) & (!is.na(p)), "Yes", "No")) %>%
    reshape2::dcast(., group ~ pass, value.var = "protein_id",
                    fun.aggregate = length) %>%
    mutate(n = Yes + No) %>%
    mutate(label = homo_labels[i])
  
  return(stat_sw)
})
stat_pro_df <- stat_pro_tables %>%
  rbindlist %>%
  filter(!group %in% c("HeLa","HEK293T")) %>%
  reshape2::melt(., id = c(1, 4:5))

pro_thres <- stat_pro_df %>%
  filter(variable %in% "Yes") %>%
  mutate(prop = value/n * 100)

p_sfigure10a <- ggplot(stat_pro_df, aes(x = label, y = value)) +
  geom_bar(aes(fill = variable), stat = "identity",
           width = .7, color = "black") +
  geom_text(aes(x = label, y = value, label = sprintf("%.2f%%", prop)),
            data = pro_thres, size = 4, color = "black", vjust = -1) +
  scale_fill_manual(values = c("#A6CEE3", "#1F78B4"), name = "Pass") +
  facet_grid(cols = vars(group)) +
  scale_y_continuous(name = "Number of Proteins",
                     expand = expansion(mult = c(0, 0.1)),
                     n.breaks = 10) +
  theme_bw() +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16));p_sfigure10a

ggsave("results/figures/supp_figure10a.pdf", p_sfigure10a, width = 9, height = 4)


p_figure4c <- ggplot(stat_pep_df, aes(x = label, y = value)) +
  geom_bar(aes(fill = variable), stat = "identity", width = .8, color = "black") +
  geom_text(aes(x = label, y = value, label = sprintf("%.2f%%", prop)),
            data = pep_thres, size = 3, color = "white", vjust = 1.5) +
  scale_fill_manual(values = c("#A6CEE3", "#1F78B4"), name = "Pass") +
  facet_grid(cols = vars(group)) +
  scale_y_continuous(name = "Number of Peptides",
                     expand = expansion(mult = c(0, 0.1)),
                     n.breaks = 5) +
  theme_bw() +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16));p_figure4c

ggsave("results/figures/figure4c.pdf", p_figure4c, width = 9.5, height = 5)


## Figure 4d & Supp Figure 10b: 特性量值的均匀性检验结果 -----------------------
all_files <- list.files("./results/tables", full.names = TRUE)
homo_files <- all_files[grepl("5_quantprop.+\\.rds", all_files)]
homo_labels <- c("DDA", "DIA", "Orbitrap", "TOF", "All")
stat_pep_tables <- pblapply(3:5, function(i) {
  
  homo_tables <- readRDS(homo_files[i])
  
  stat_sw <- homo_tables %>%
    rbindlist(., idcol = "group") %>%
    dplyr::rename(p = `Chisq'p`) %>%
    mutate_at("p", ~ ifelse(is.na(.), 0, .)) %>%
    group_by(peptide_sequence, group) %>%
    summarise_at("p", max) %>%
    mutate(pass = ifelse((p > .05) & (!is.na(p)), "Yes", "No")) %>%
    reshape2::dcast(., group ~ pass, value.var = "peptide_sequence",
                    fun.aggregate = length) %>%
    mutate(n = Yes + No) %>%
    mutate(label = homo_labels[i])
  
  return(stat_sw)
})
stat_pep_df <- stat_pep_tables %>%
  rbindlist %>%
  filter(!group %in% "HeLa/HEK293T") %>%
  reshape2::melt(., id = c(1, 4:5))

pep_thres <- stat_pep_df %>%
  filter(variable %in% "Yes") %>%
  mutate(prop = value/n * 100)

stat_pro_tables <- pblapply(3:5, function(i) {
  
  homo_tables <- readRDS(homo_files[i])
  
  stat_sw <- homo_tables %>%
    rbindlist(., idcol = "group") %>%
    dplyr::rename(p = `Chisq'p`) %>%
    mutate_at("p", ~ ifelse(is.na(.), 0, .)) %>%
    group_by(protein_id, group) %>%
    summarise_at("p", max) %>%
    mutate(pass = ifelse((p > .05) & (!is.na(p)), "Yes", "No")) %>%
    reshape2::dcast(., group ~ pass, value.var = "protein_id",
                    fun.aggregate = length) %>%
    mutate(n = Yes + No) %>%
    mutate(label = homo_labels[i])
  
  return(stat_sw)
})
stat_pro_df <- stat_pro_tables %>%
  rbindlist %>%
  filter(!group %in% "HeLa/HEK293T") %>%
  reshape2::melt(., id = c(1, 4:5))

pro_thres <- stat_pro_df %>%
  filter(variable %in% "Yes") %>%
  mutate(prop = value/n * 100)

p_sfigure10b <- ggplot(stat_pro_df, aes(x = label, y = value)) +
  geom_bar(aes(fill = variable), stat = "identity",
           width = .7, color = "black") +
  geom_text(aes(x = label, y = value, label = sprintf("%.2f%%", prop)),
            data = pro_thres, size = 3.5, color = "white", vjust = 1.3) +
  scale_fill_manual(values = c("#A6CEE3", "#1F78B4"), name = "Pass") +
  facet_grid(cols = vars(group)) +
  scale_y_continuous(name = "Number of Proteins",
                     expand = expansion(mult = c(0, 0.1)),
                     n.breaks = 10) +
  theme_bw() +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16));p_sfigure10b

ggsave("results/figures/supp_figure10b.pdf", p_sfigure10b, width = 9, height = 4)

p_figure4d <- ggplot(stat_pep_df, aes(x = label, y = value)) +
  geom_bar(aes(fill = variable), stat = "identity", width = .8, color = "black") +
  geom_text(aes(x = label, y = value, label = sprintf("%.2f%%", prop)),
                vjust = 1.3, color = "white", data = pep_thres, size = 3) +
  scale_fill_manual(values = c("#A6CEE3", "#1F78B4"), name = "Pass") +
  facet_grid(cols = vars(group)) +
  scale_y_continuous(name = "Number of Peptides",
                     expand = expansion(mult = c(0, 0.1)),
                     n.breaks = 5) +
  guides(color = "none") +
  theme_bw() +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16));p_figure4d

ggsave("results/figures/figure4d.pdf", p_figure4d, width = 9.5, height = 5)


## Figure 4e : 特性量值的均匀性检验结果举例 -----------------------
homo_tables <- readRDS("./results/tables/5_quantprop_homotest_reml.rds")

ratiobyd6_tables <- readRDS("./results/tables/2_quantdata_list_pep_ratiobyd6_2025_10labs.rds")
ratiobyd6_tables <- ratiobyd6_tables[c(1, 3:4)]
names(ratiobyd6_tables) <- names(homo_tables)[1:3]

df_homo <- homo_tables[1:3] %>% rbindlist(., idcol = "group")

## 挑选明星蛋白质: P00558/P26038/P01903/P10809
vip_peptides <- df_homo %>%
  filter(protein_id %in% c("P00558", "P26038", "P01903", "P10809")) %>%
  filter(protein_id %in% "P01903") %>%
  filter(peptide_sequence %in% "NGKPVTTGVSETVFLPR") %>%
  filter(is.homogeneous)

df_ratio_vip <- ratiobyd6_tables %>%
  rbindlist(., idcol = "group") %>%
  filter(!lab_id %in% "zhejiang_university") %>%
  inner_join(., vip_peptides, by = c("peptide_sequence", "group"), relationship = "many-to-many")

df_ratio_vip_intra <- df_ratio_vip %>%
  group_by(peptide_sequence, protein_id, group, lab_id, tube) %>%
  summarise_at("value", ~ . - mean(.)) %>%
  mutate(class = "Intra-vial")

df_ratio_vip_inter <- df_ratio_vip %>%
  group_by(peptide_sequence, protein_id, group, lab_id, tube) %>%
  summarise_at("value", mean) %>%
  group_by(peptide_sequence, protein_id, group, lab_id) %>%
  mutate(lab_mean = mean(value)) %>%
  mutate_at("value", ~ . - lab_mean) %>%
  mutate(class = "Inter-vial")

df_ratio_vip_final <- df_ratio_vip_inter %>%
  rbind(., df_ratio_vip_intra) %>%
  # filter(peptide_sequence %in% c("ELSLAGNELGDEGAR")) %>%
  mutate_at("group", ~ factor(., levels = c("HeLa/HEK293T", "D5/D6", "F7/D6", "M8/D6")))

colors.class <- c("#54278F", "#9E9AC8")
names(colors.class) <- c("Inter-vial", "Intra-vial")
colors.group <- c("#4CC3D9", "#FFC65D", "#F16745", "#E7298A")
names(colors.group) <- c("D5/D6", "F7/D6", "M8/D6", "HeLa/HEK293T")

set.seed(1)
p <- ggplot() +
  geom_boxplot(aes(fill = group, alpha = class, x = class, y = value),
               data = df_ratio_vip_final, outliers = FALSE) +
  # geom_jitter(aes(fill = group, alpha = class, x = class, y = value),
  #             data = df_ratio_vip_final, width = 0.3, shape = 21, color = "black") +
  geom_text(aes(x = "Inter-vial", y = Inf, label = sprintf("Chisq'p (REML) = %.2f", `Chisq'p`)),
            data = vip_peptides, size = 3.5, hjust = .25, vjust = 1.5) +
  facet_grid(protein_id + peptide_sequence ~ group, scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10),
        legend.position = "none",
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(name = "Residual (log2 transformed)",
                     expand = expansion(mult = c(0.05, 0.15))) +
  scale_alpha_manual(values = c(1, .6)) +
  scale_fill_manual(values = colors.group);p

ggsave("./results/figures/figure4e.pdf", p, height = 3, width = 10)



## Figure 4f & Supp Figure 11a: 标称特性值的稳定性检验结果 -----------------------
stab_tables <- readRDS("./results/tables/6_qualiprop_stabtest_fisher.rds")
rm(list = setdiff(ls(), c("stab_tables")))
gc()

stat_stab1_1 <- stab_tables %>%
  rbindlist(., idcol = "group") %>%
  filter(category %in% "Short-term-Transport") %>%
  tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
  mutate_at("p", ~ ifelse(is.na(.), 0, .)) %>%
  group_by(category, peptide_sequence, group, temperature) %>%
  summarise_at("p", max) %>%
  mutate(pass = ifelse((p > .05) & (!is.na(p)), "Yes", "No")) %>%
  reshape2::dcast(., category + group + temperature ~ pass, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(No = 0) %>%
  mutate(n = Yes)

stat_stab1_2 <- stab_tables %>%
  rbindlist(., idcol = "group") %>%
  filter(category %in% "Short-term-Transport") %>%
  tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
  mutate_at("p", ~ ifelse(is.na(.), 0, .)) %>%
  group_by(category, protein_id, group, temperature) %>%
  summarise_at("p", max) %>%
  mutate(pass = ifelse((p > .05) & (!is.na(p)), "Yes", "No")) %>%
  reshape2::dcast(., category + group + temperature ~ pass, value.var = "protein_id",
                  fun.aggregate = length) %>%
  mutate(No = 0) %>%
  mutate(n = Yes)

stat_stab2_1 <- stab_tables %>%
  rbindlist(., idcol = "group") %>%
  filter(category %in% "Short-term-Use") %>%
  tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
  mutate_at("p", ~ ifelse(is.na(.), 0, .)) %>%
  group_by(category, peptide_sequence, group, temperature) %>%
  summarise_at("p", max) %>%
  mutate(pass = ifelse((p > .05) & (!is.na(p)), "Yes", "No")) %>%
  reshape2::dcast(., category + group + temperature ~ pass, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(No = 0) %>%
  mutate(n = Yes)

stat_stab2_2 <- stab_tables %>%
  rbindlist(., idcol = "group") %>%
  filter(category %in% "Short-term-Use") %>%
  tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
  mutate_at("p", ~ ifelse(is.na(.), 0, .)) %>%
  group_by(category, protein_id, group, temperature) %>%
  summarise_at("p", max) %>%
  mutate(pass = ifelse((p > .05) & (!is.na(p)), "Yes", "No")) %>%
  reshape2::dcast(., category + group + temperature ~ pass, value.var = "protein_id",
                  fun.aggregate = length) %>%
  mutate(No = 0) %>%
  mutate(n = Yes)

stat_stab3_1 <- stab_tables %>%
  rbindlist(., idcol = "group") %>%
  filter(category %in% "Long-term") %>%
  tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
  mutate_at("p", ~ ifelse(is.na(.), 0, .)) %>%
  group_by(category, peptide_sequence, group, temperature) %>%
  summarise_at("p", max) %>%
  mutate(pass = ifelse((p > .05) & (!is.na(p)), "Yes", "No")) %>%
  reshape2::dcast(., category + group + temperature ~ pass, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(n = Yes + No)

stat_stab3_2 <- stab_tables %>%
  rbindlist(., idcol = "group") %>%
  filter(category %in% "Long-term") %>%
  tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
  mutate_at("p", ~ ifelse(is.na(.), 0, .)) %>%
  group_by(category, protein_id, group, temperature) %>%
  summarise_at("p", max) %>%
  mutate(pass = ifelse((p > .05) & (!is.na(p)), "Yes", "No")) %>%
  reshape2::dcast(., category + group + temperature ~ pass, value.var = "protein_id",
                  fun.aggregate = length) %>%
  mutate(n = Yes + No)

stat_stab1 <- stat_stab1_1 %>%
  filter(!temperature %in% "-20°C") %>%
  rbind(., stat_stab2_1) %>%
  filter(!temperature %in% "-20°C") %>%
  rbind(., stat_stab3_1) %>%
  reshape2::melt(., id = c(1:3, 6)) %>%
  mutate_at("group", ~ gsub("Quartet ", "", .)) %>%
  mutate_at("variable", ~ factor(., levels = c("No", "Yes"), ordered = TRUE))

pep_thres <- stat_stab1 %>%
  filter(variable %in% "Yes") %>%
  mutate(prop = value/n * 100)

stat_stab2 <- stat_stab1_2 %>%
  filter(!temperature %in% "-20°C") %>%
  rbind(., stat_stab2_2) %>%
  filter(!temperature %in% "-20°C") %>%
  rbind(., stat_stab3_2) %>%
  reshape2::melt(., id = c(1:3, 6)) %>%
  mutate_at("group", ~ gsub("Quartet ", "", .)) %>%
  mutate_at("variable", ~ factor(., levels = c("No", "Yes"), ordered = TRUE))

pro_thres <- stat_stab2 %>%
  filter(variable %in% "Yes") %>%
  mutate(prop = value/n * 100)

p_sfigure11a <- ggplot(stat_stab2, aes(x = category, y = value)) +
  geom_bar(aes(fill = variable), stat = "identity",
           width = .85, color = "black") +
  geom_text(aes(x = category, y = value, label = sprintf("%.2f%%", prop)),
            data = pro_thres, size = 3, color = "white", vjust = 1.5) +
  scale_fill_manual(values = c("#A1D99B", "#006D2C"), name = "Pass") +
  facet_grid(cols = vars(group)) +
  scale_y_continuous(name = "Number of Proteins",
                     expand = expansion(mult = c(0, 0.1)), n.breaks = 10) +
  theme_bw() +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16));p_sfigure11a

ggsave("results/figures/supp_figure11a.pdf", p_sfigure11a, width = 9, height = 5)

p_figure4f <- ggplot(stat_stab1, aes(x = category, y = value)) +
  geom_bar(aes(fill = variable), stat = "identity", width = .8, color = "black") +
  geom_text(aes(x = category, y = value, label = sprintf("%.2f%%", prop)),
            data = pep_thres, size = 3, color = "white", vjust = 1.5) +
  scale_fill_manual(values = c("#A1D99B", "#006D2C"), name = "Pass") +
  facet_grid(cols = vars(group)) +
  scale_y_continuous(name = "Number of Peptides",
                     expand = expansion(mult = c(0, 0.1)),
                     n.breaks = 5) +
  theme_bw() +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16));p_figure4f

ggsave("results/figures/figure4f.pdf", p_figure4f, width = 9.5, height = 5)


## Figure 4g & Supp Figure 11b: 特性量值的稳定性检验结果 -----------------------
stab_tables <- readRDS("./results/tables/6_quantprop_stabtest_linear.rds")

all_meta <- fread("./results/tables/6_outlier_madist_quartet.csv")
norm_tables <- readRDS("./results/tables/4_quantprop_normtest_shapirowilk.rds")
norm_tables <- norm_tables[c(1:2, 4)]
df_uchar <- norm_tables %>% rbindlist(., idcol = "group")

stat_stab1_1 <- stab_tables %>%
  rbindlist(., idcol = "group", use.names=TRUE) %>%
  # filter(category %in% "Short-term-Transport", grepl("14", period), !grepl("30", period)) %>%
  filter(category %in% "Short-term-Transport", period %in% ("0, 3, 7, 14")) %>%
  mutate_at(7:12, as.numeric) %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  distinct(category, peptide_sequence, group, temperature, pass) %>%
  reshape2::dcast(., category + group + temperature ~ pass, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(n = Yes + No)

stat_stab1_2 <- stab_tables %>%
  rbindlist(., idcol = "group", use.names=TRUE) %>%
  # filter(category %in% "Short-term-Transport", grepl("14", period), !grepl("30", period)) %>%
  filter(category %in% "Short-term-Transport", period %in% ("0, 3, 7, 14")) %>%
  mutate_at(7:12, as.numeric) %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  distinct(category, protein_id, group, temperature, pass) %>%
  reshape2::dcast(., category + group + temperature ~ pass, value.var = "protein_id",
                  fun.aggregate = length) %>%
  mutate(n = Yes + No)

stat_stab2_1 <- stab_tables %>%
  rbindlist(., idcol = "group", use.names=TRUE) %>%
  # filter(category %in% "Short-term-Use", grepl("7", period), !grepl("30", period)) %>%
  filter(category %in% "Short-term-Use", period %in% ("0, 1, 2, 4, 7")) %>%
  mutate_at(7:12, as.numeric) %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  distinct(category, peptide_sequence, group, temperature, pass) %>%
  reshape2::dcast(., category + group + temperature ~ pass, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(n = Yes + No)

stat_stab2_2 <- stab_tables %>%
  rbindlist(., idcol = "group", use.names=TRUE) %>%
  # filter(category %in% "Short-term-Use", grepl("7", period), !grepl("30", period)) %>%
  filter(category %in% "Short-term-Use", period %in% ("0, 1, 2, 4, 7")) %>%
  mutate_at(7:12, as.numeric) %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  distinct(category, protein_id, group, temperature, pass) %>%
  reshape2::dcast(., category + group + temperature ~ pass, value.var = "protein_id",
                  fun.aggregate = length) %>%
  mutate(n = Yes + No)

stat_stab3_1 <- stab_tables %>%
  rbindlist(., idcol = "group", use.names=TRUE) %>%
  filter(category %in% "Long-term", grepl("0", period) & grepl("58.03", period)) %>%
  mutate_at(7:12, as.numeric) %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  distinct(category, peptide_sequence, group, temperature, pass) %>%
  reshape2::dcast(., category + group + temperature ~ pass, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(n = Yes + No)

stat_stab3_2 <- stab_tables %>%
  rbindlist(., idcol = "group", use.names=TRUE) %>%
  filter(category %in% "Long-term", grepl("0", period) & grepl("58.03", period)) %>%
  mutate_at(7:12, as.numeric) %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  distinct(category, protein_id, group, temperature, pass) %>%
  reshape2::dcast(., category + group + temperature ~ pass, value.var = "protein_id",
                  fun.aggregate = length) %>%
  mutate(n = Yes + No)

stat_stab1 <- stat_stab1_1 %>%
  filter(!temperature %in% "-20°C") %>%
  rbind(., stat_stab2_1) %>%
  filter(!temperature %in% "-20°C") %>%
  rbind(., stat_stab3_1) %>%
  reshape2::melt(., id = c(1:3, 6)) %>%
  mutate_at("variable", ~ factor(., levels = c("No", "Yes"), ordered = TRUE))

pep_thres <- stat_stab1 %>%
  filter(variable %in% "Yes") %>%
  mutate(prop = value/n * 100)

stat_stab2 <- stat_stab1_2 %>%
  filter(!temperature %in% "-20°C") %>%
  rbind(., stat_stab2_2) %>%
  filter(!temperature %in% "-20°C") %>%
  rbind(., stat_stab3_2) %>%
  reshape2::melt(., id = c(1:3, 6)) %>%
  mutate_at("variable", ~ factor(., levels = c("No", "Yes"), ordered = TRUE))

pro_thres <- stat_stab2 %>%
  filter(variable %in% "Yes") %>%
  mutate(prop = value/n * 100)

p_sfigure11b <- ggplot(stat_stab2, aes(x = category, y = value)) +
  geom_bar(aes(fill = variable), stat = "identity",
           width = .85, color = "black") +
  geom_text(aes(x = category, y = value, label = sprintf("%.2f%%", prop)),
            data = pro_thres, size = 3, color = "white", vjust = 1.5) +
  scale_fill_manual(values = c("#A1D99B", "#006D2C"), name = "Pass") +
  facet_grid(cols = vars(group)) +
  scale_y_continuous(name = "Number of Proteins",
                     expand = expansion(mult = c(0, 0.1)), n.breaks = 10) +
  theme_bw() +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16));p_sfigure11b

ggsave("results/figures/supp_figure11b.pdf", p_sfigure11b, width = 9, height = 5)

p_figure4g <- ggplot(stat_stab1, aes(x = category, y = value)) +
  geom_bar(aes(fill = variable), stat = "identity", width = .8, color = "black") +
  geom_text(aes(x = category, y = value, label = sprintf("%.2f%%", prop)),
            data = pep_thres, size = 3, color = "white", vjust = 1.5) +
  scale_fill_manual(values = c("#A1D99B", "#006D2C"), name = "Pass") +
  facet_grid(cols = vars(group)) +
  scale_y_continuous(name = "Number of Peptides",
                     expand = expansion(mult = c(0, 0.1)),
                     n.breaks = 5) +
  theme_bw() +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16));p_figure4g

ggsave("results/figures/figure4g.pdf", p_figure4g, width = 9.5, height = 5)


## Figure 4h: 稳定性结果举例 -----------------------
stab_tables <- readRDS("./results/tables/6_quantprop_stabtest_linear.rds")

ratiobyd6_tables <- readRDS("./results/tables/6_quantdata_list_pep_ratiobyd6_all_longterm.rds")
ratiobyd6_tables <- ratiobyd6_tables[c(1, 3:4)]
names(ratiobyd6_tables) <- names(stab_tables)[1:3]

df_stab <- stab_tables %>%
  rbindlist(., idcol = "group", use.names=TRUE) %>%
  filter(period %in% ("0, 3, 7, 14") | period %in% ("0, 1, 2, 4, 7") | (grepl("0", period) & grepl("58.03", period)))

## 挑选明星蛋白质: P00558/P26038/P01903/P10809
df_stab_vip <- df_stab %>%
  filter(group %in% c("D5/D6", "F7/D6", "M8/D6")) %>%
  filter(pass %in% "Yes", s_b1 < .02) %>%
  filter(protein_id %in% "P01903") %>%
  filter(peptide_sequence %in% "FHYLPFLPSTEDVYDCR") %>%
  mutate(label_tmp = paste(category, temperature, sep = "_")) %>%
  filter(grepl("-78.5°C", label_tmp)) %>%
  filter(grepl("Long-term|-78.5°C|4°C", label_tmp))

all_meta <- fread("./results/tables/6_outlier_madist_quartet.csv")

## 运输稳定性
df_ratio_vip <- ratiobyd6_tables %>%
  rbindlist(., idcol = "group", use.names=TRUE) %>%
  filter(lab_id %in% "HUP_DIA", time_day >= 1290) %>%
  inner_join(., df_stab_vip, by = c("peptide_sequence", "group"), relationship = "many-to-many") %>%
  distinct(group, peptide_sequence, protein_id, analysis_id, value) %>%
  left_join(., all_meta, by = "analysis_id") %>%
  mutate_at("time_day", ~1320 - .) %>%
  mutate_at("value", ~ 2 ^ (.)) %>%
  mutate_at("temperature", ~ ifelse(. %in% "-80°C", "-20°C_-78.5°C", .)) %>%
  tidyr::separate_rows(temperature, sep = "_") %>%
  filter(temperature %in% "-78.5°C") %>%
  filter(time_day != 30)

time_breaks <- df_ratio_vip %>%
  pull(time_day) %>%
  unique()

colors.group <- c("#4CC3D9", "#FFC65D", "#F16745", "#E7298A")
names(colors.group) <- c("D5/D6", "F7/D6", "M8/D6", "HeLa/HEK293T")

p <- ggplot(df_ratio_vip, aes(x = time_day, y = value)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", color = "black", width = 1) +
  stat_summary(aes(fill = group), fun = mean, geom = "point", color = "black", shape = 21, size = 3) +
  # stat_summary(aes(color = group), fun = mean, geom = "line") +
  geom_smooth(aes(group = group, color = group), linewidth = .5, lty = 2, alpha = .2,
              method = "lm", formula = y ~ x) +
  ggpubr::stat_cor(size = 5) +
  # ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~`,`~")),
  #                               label.y = 1.7, size = 4.5) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14)) +
  facet_grid(cols = vars(group), rows = vars(protein_id, peptide_sequence), scales = "fixed") +
  scale_color_manual(values = colors.group, name = "Sample Pair") +
  scale_fill_manual(values = colors.group, name = "Sample Pair") +
  scale_x_continuous(name = "Time Point (Day)", breaks = time_breaks) +
  scale_y_continuous(name = "Fold changes", limits = c(0, 1.5), n.breaks = 5);p

ggsave("./results/figures/figure4h.pdf", p, height = 2.8, width = 10)


