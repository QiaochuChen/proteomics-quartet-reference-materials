## Title: 将研制报告内容转化为sci文章（figure2标称特性部分）。
## Author: Qiaochu Chen
## Date: Jun 8th, 2026

library(readxl)
library(data.table)
library(rstatix)
library(pbapply)
library(parallel)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(MsCoreUtils)
library(arrow)
library(Biostrings)
library(RColorBrewer)
library(ggh4x)

# 基础映射与颜色配置
labels.lab <- c("Lab-1", "Lab-2", "Lab-0", "Lab-3", "Lab-4", "Lab-5", "Lab-6", "Lab-7", "Lab-8", "Lab-9")
names(labels.lab) <- c("QLB", "NCP", "ZJU", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")

colors.sample <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745", "#E7298A", "#4D9221")
names(colors.sample) <- c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")


## Supplementary Figure 1a: 作图比较鉴定数目 ----------------------
rm(list = setdiff(ls(), c("labels.lab", "colors.sample")))
gc()

all_pep <- fread("./results/tables/3_count.csv")

all_long <- all_pep %>%
  filter(!lab %in% "ZJU") %>%
  select(!Precursors.Identified) %>%
  mutate(source = ifelse(grepl("Lab", lab), "Public", "Local")) %>%
  reshape2::melt(., id = c(1:3, 6:7))

all_long2 <- all_long %>%
  filter(source %in% "Local") %>%
  mutate_at("sample", ~ ifelse(. %in% c("HeLa", "HEK293T"), ., paste("Quartet", .))) %>%
  mutate_at("sample", ~ factor(., levels = c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T"))) %>%
  mutate_at("lab", ~ labels.lab[.])

p1a <- ggplot(all_long2, aes(x = sample, y = value)) +
  stat_summary(aes(fill = sample), fun = mean, geom = "bar", width = .7,
               color = "black", linewidth = .3) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  stat_summary(aes(label = scales::label_comma(accuracy = 1)(after_stat(y)),
                   vjust = ifelse(variable %in% "Proteins.Identified", -2, -3)),
               fun = "mean", geom = "text", size = 4.5) +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        strip.text = element_text(size = 14, face = "bold", color = "black", margin = margin(0.3, 0.3, 0.3, 0.3, "cm")),
        strip.background = element_blank(),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm")) +
  scale_y_continuous(n.breaks = 6, name = "Number of Identifications",
                     expand = expansion(mult = c(0, 0.2)),
                     labels = scales::label_comma()) +
  scale_fill_manual(values = colors.sample) +
  facet_wrap(~ variable, scales = "free_y")

ggsave("./results/figures/supp_figure1a.pdf", p1a, height = 4, width = 14, dpi = 600)


## Supplementary Figure 1b: 作图比较 CV ----------------------
rm(list = setdiff(ls(), c("labels.lab", "colors.sample", "all_long")))
gc()

dda_cv_tables <- readRDS("./results/tables/3_dda_cv.rds")
dia_cv_tables <- readRDS("./results/tables/3_dia_cv.rds")

df_cv2 <- c(dia_cv_tables[[2]], dda_cv_tables[[2]]) %>%
  rbindlist %>%
  filter(!lab %in% "ZJU") %>%
  mutate(source = ifelse(grepl("Lab", lab), "Public", "Local")) %>%
  filter(source %in% "Local") %>%
  reshape2::dcast(., Peptides.Quantified + source + tube + mode ~ sample,
                  value.var = "cv", fun.aggregate = median) %>%
  na.omit %>%
  reshape2::melt(., id = 2:10)

df_cv3 <- c(dia_cv_tables[[3]], dda_cv_tables[[3]]) %>%
  rbindlist %>%
  filter(!lab %in% "ZJU") %>%
  mutate(source = ifelse(grepl("Lab", lab), "Public", "Local")) %>%
  filter(source %in% "Local") %>%
  reshape2::dcast(., Proteins.Quantified + source + tube + mode ~ sample,
                  value.var = "cv", fun.aggregate = median) %>%
  na.omit %>%
  reshape2::melt(., id = 2:10)

df_cv_test <- df_cv2 %>%
  rbind(., df_cv3) %>%
  select(!value) %>%
  select(variable, source, tube, mode, everything()) %>%
  reshape2::melt(., id = 1:4, variable.name = "sample", value.name = "cv") %>%
  mutate_at("sample", ~ ifelse(. %in% c("HeLa", "HEK293T"), as.character(.), paste("Quartet", .))) %>%
  mutate_at("sample", ~ factor(., levels = c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")))

cv_thres <- df_cv_test %>%
  group_by(sample, variable) %>%
  summarise(cv_median = median(cv),
            cv_mean = mean(cv),
            cv_sd = sd(cv), .groups = "drop")

p1b <- ggplot(df_cv_test, aes(x = sample, y = cv)) +
  geom_boxplot(aes(fill = sample), outlier.size = .1, width = .7) +
  geom_text(aes(x = sample, y = cv_median,
                label = sprintf("%.2f%%", cv_median * 100)),
            data = cv_thres, size = 3, color = "white") +
  scale_fill_manual(values = colors.sample) +
  scale_y_continuous(n.breaks = 8, name = "Coefficient of Variation (CV)",
                     labels = scales::percent) +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        strip.text = element_text(size = 14, face = "bold", color = "black", margin = margin(0.3, 0.3, 0.3, 0.3, "cm")),
        strip.background = element_blank(),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm")) +
  facet_wrap(~ variable, scales = "free_y")

ggsave("./results/figures/supp_figure1b.pdf", p1b, height = 4, width = 14, dpi = 600)


## Supplementary Figure 2a: 作图比较 HeLa 鉴定数目 ----------------------
rm(list = setdiff(ls(), c("labels.lab", "colors.sample", "all_long")))
gc()

all_pep <- fread("./results/tables/3_count.csv")

all_long <- all_pep %>%
  filter(sample %in% "HeLa", tube == 1, !lab %in% "ZJU") %>%
  mutate(source = ifelse(grepl("Lab", lab), "Public", "Local")) %>%
  reshape2::melt(., id = c(1:3, 7:8)) %>%
  filter(!variable %in% "Precursors.Identified")

colors.source <- c("Local" = "#A6CEE3", "Public" = "#1F78B4")

p2a <- ggplot(all_long, aes(x = source, y = value)) +
  stat_summary(aes(fill = source), fun = mean, geom = "bar", width = .7,
               color = "black", linewidth = .3) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  stat_summary(aes(label = scales::label_comma(accuracy = 1)(after_stat(y)),
                   vjust = ifelse(grepl("Protein", variable), -1.5, -4)),
               fun = "mean", geom = "text", size = 4.5) +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        strip.text = element_text(size = 14, face = "bold", color = "black", margin = margin(0.3, 0.3, 0.3, 0.3, "cm")),
        strip.background = element_blank(),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm")) +
  scale_y_continuous(n.breaks = 6, name = "Number of Identifications",
                     expand = expansion(mult = c(0, 0.2)),
                     labels = scales::label_comma()) +
  scale_fill_manual(values = colors.source) +
  ggh4x::facet_grid2(~ variable, scales = "free_y")

ggsave("./results/figures/supp_figure2a.pdf", p2a, height = 4, width = 8, dpi = 600)


## Supplementary Figure 2b: 作图比较 HeLa CV ----------------------
rm(list = setdiff(ls(), c("labels.lab", "colors.sample", "all_long")))
gc()

dda_cv_tables <- readRDS("./results/tables/3_dda_cv.rds")
dia_cv_tables <- readRDS("./results/tables/3_dia_cv.rds")

df_cv2 <- c(dia_cv_tables[[2]], dda_cv_tables[[2]]) %>%
  rbindlist %>%
  filter(sample %in% "HeLa", tube == 1) %>%
  filter(!lab %in% "ZJU") %>%
  mutate(source = ifelse(grepl("Lab", lab), "Public", "Local")) %>%
  reshape2::dcast(., Peptides.Quantified + sample + tube + mode ~ source,
                  value.var = "cv", fun.aggregate = median) %>%
  na.omit %>%
  reshape2::melt(., id = 2:6)

df_cv3 <- c(dia_cv_tables[[3]], dda_cv_tables[[3]]) %>%
  rbindlist %>%
  filter(sample %in% "HeLa", tube == 1) %>%
  filter(!lab %in% "ZJU") %>%
  mutate(source = ifelse(grepl("Lab", lab), "Public", "Local")) %>%
  reshape2::dcast(., Proteins.Quantified + sample + tube + mode ~ source,
                  value.var = "cv", fun.aggregate = median) %>%
  na.omit %>%
  reshape2::melt(., id = 2:6)

df_cv_test <- df_cv2 %>%
  rbind(., df_cv3) %>%
  select(variable, sample, tube, mode, Local, Public) %>%
  reshape2::melt(., id = 1:4, variable.name = "source", value.name = "cv") %>%
  mutate_at("variable", ~ gsub("s.Quantified", "-level", .)) %>%
  mutate_at("variable", ~ factor(., levels = c("Precursor-level", "Peptide-level", "Protein-level")))

cv_thres <- df_cv_test %>%
  group_by(source, variable) %>%
  summarise(cv_median = median(cv), .groups = "drop")

p2b <- ggplot(df_cv_test, aes(x = reorder(source, cv), y = cv)) +
  geom_boxplot(aes(fill = source), outlier.size = .5, width = .7) +
  geom_hline(aes(yintercept = cv_median, colour = source), cv_thres, lty = 2) +
  geom_text(aes(x = source, y = cv_median,
                label = sprintf("Median = %.2f%%", cv_median * 100)),
            data = cv_thres, size = 3, vjust = -4) +
  scale_fill_manual(values = c("Local" = "#A6CEE3", "Public" = "#1F78B4")) +
  scale_color_manual(values = c("Local" = "darkred", "Public" = "#276419")) +
  scale_y_continuous(limits = c(0, 2), n.breaks = 8, name = "Coefficient of Variation (CV)",
                     labels = scales::percent) +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        strip.text = element_text(size = 14, face = "bold", color = "black", margin = margin(0.3, 0.3, 0.3, 0.3, "cm")),
        strip.background = element_blank(),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm")) +
  ggh4x::facet_grid2( ~ variable, scales = "fixed")

ggsave("./results/figures/supp_figure2b.pdf", p2b, height = 4, width = 8, dpi = 600)


## Supplementary Figure 2c: 定性验证-检查每个实验室的定性 ----------------------
rm(list = setdiff(ls(), c("labels.lab", "colors.sample", "all_long")))
gc()

meta_hela <- fread("./results/tables/3_outlier_madist_hela.csv")
all_tables <- readRDS("./results/tables/3_quantdata_list_pep_2025_21labs.rds")
all_tables <- all_tables %>% rbindlist %>% split(., ~ lab_id)
all_tables <- all_tables[1:20]

all_peptides <- pblapply(all_tables, function(tmp_table) unique(tmp_table$peptide_sequence))

ref_peptides <- all_peptides[[5]]
for (i in 6:15) ref_peptides <- intersect(all_peptides[[i]], ref_peptides)

inter_peptides <- pblapply(all_peptides[c(1:4, 16:20)], function(tmp_peptides) {
  union_peptides <- union(tmp_peptides, ref_peptides)
  inter_peptides <- intersect(tmp_peptides, ref_peptides)
  
  stat_peptides <- data.frame(n_intersect = length(inter_peptides),
                              n_union = length(union_peptides),
                              n_ref = length(ref_peptides))
  return(stat_peptides)
})

df_inter <- inter_peptides %>%
  rbindlist(., idcol = "lab_id") %>%
  mutate(proportion = n_intersect / n_ref)

labels.lab_2c <- c("Lab-1", "Lab-2", "Lab-0", "Lab-3", "Lab-4", "Lab-5", "Lab-6", "Lab-7", "Lab-8", "Lab-9")
names(labels.lab_2c) <- c("qinglian_bio", "phoenix", "zhejiang_university", "fudan_university",
                          "national_institute_of_methodology", "biotech_pack", "thermofisher_shanghai",
                          "omicsolution", "cas_tianjin", "academy_of_chinese_medical_sciences")

p2c <- ggplot(df_inter, aes(x = reorder(lab_id, -proportion), y = proportion)) +
  geom_col(fill = "#A6CEE3", width = .7, color = "black", linewidth = .3) +
  geom_text(aes(label = scales::percent(round(proportion, digits = 4))), nudge_y = .04) +
  geom_hline(yintercept = .6, lty = 2, col = "red", linewidth = .6) +
  scale_y_continuous(label = scales::percent, n.breaks = 10, name = "Overlapping Proportion") +
  scale_x_discrete(labels = labels.lab_2c) +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        strip.text = element_text(size = 14, face = "bold", color = "black", margin = margin(0.3, 0.3, 0.3, 0.3, "cm")),
        strip.background = element_blank(),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1, vjust = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm"))

ggsave("./results/figures/supp_figure2c.pdf", p2c, width = 8, height = 4, limitsize = FALSE, dpi = 600)


## Supplementary Figure 2d: 定量验证-检查每个实验室的CV ----------------------
rm(list = ls())
gc()

meta_hela <- fread("./results/tables/3_outlier_madist_hela.csv")
all_tables <- readRDS("./results/tables/3_quantdata_list_pep_2025_21labs.rds")
all_tables_bylab <- all_tables %>% rbindlist %>% split(., ~ lab_id)
all_tables_bylab <- all_tables_bylab[1:20]

all_peptides <- pblapply(all_tables_bylab, function(tmp_table) unique(tmp_table$peptide_sequence))

ref_peptides <- all_peptides[[5]]
for (i in 6:15) ref_peptides <- intersect(all_peptides[[i]], ref_peptides)

outliers <- meta_hela$analysis_id[meta_hela$is.outlier]
cv_tables <- pblapply(all_tables_bylab, function(tmp_table) {
  df_cv_i <- tmp_table %>%
    filter(!analysis_id %in% outliers) %>%
    filter(peptide_sequence %in% ref_peptides) %>%
    group_by(peptide_sequence, lab_id, sample, tube) %>%
    summarise(cv = sd(value) / mean(value), .groups = "drop") %>%
    na.omit
  return(df_cv_i)
})

labels.lab_2d <- c("Lab-1", "Lab-2", "Lab-0", "Lab-3", "Lab-4", "Lab-5", "Lab-6", "Lab-7", "Lab-8", "Lab-9",
                   paste("Lab", 1:11, sep = "_"))
names(labels.lab_2d) <- c("qinglian_bio", "phoenix", "zhejiang_university", "fudan_university",
                          "national_institute_of_methodology", "biotech_pack", "thermofisher_shanghai",
                          "omicsolution", "cas_tianjin", "academy_of_chinese_medical_sciences",
                          paste("Lab", 1:11, sep = "_"))

df_cv <- cv_tables %>%
  rbindlist %>%
  mutate_at("lab_id", ~ factor(., levels = names(labels.lab_2d)))

cv_threshold <- df_cv %>%
  filter(grepl("Lab", lab_id)) %>%
  group_by(lab_id) %>%
  summarise_at("cv", median) %>%
  pull(cv) %>%
  max

df_cv_test <- df_cv %>% filter(!grepl("Lab", lab_id))

cv_thres <- df_cv_test %>%
  group_by(lab_id) %>%
  summarise(cv_median = median(cv), .groups = "drop")

p2d <- ggplot(df_cv_test, aes(x = reorder(lab_id, cv, median), y = cv)) +
  geom_boxplot(fill = "#A6CEE3", outliers = FALSE, width = .7) +
  geom_hline(yintercept = cv_threshold, lty = 2, col = "red", linewidth = .6) +
  geom_text(aes(x = lab_id, y = cv_median,
                label = sprintf("%.2f%%", cv_median * 100)),
            data = cv_thres, size = 3.5, vjust = 1.5, color = "black") +
  annotate("text", x = "national_institute_of_methodology", y = .2,
           label = "CV = 14.85%", color = "red", size = 5) +
  scale_y_continuous(label = scales::percent, name = "Coefficient of Variation (CV)",
                     limits = c(0, 1)) +
  scale_x_discrete(labels = labels.lab_2d) +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        strip.text = element_text(size = 14, face = "bold", color = "black", margin = margin(0.3, 0.3, 0.3, 0.3, "cm")),
        strip.background = element_blank(),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1, vjust = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm"))

ggsave("./results/figures/supp_figure2d.pdf", p2d, width = 8, height = 4, limitsize = FALSE, dpi = 600)


## Supplementary Figure 3a: 原始值验证-检查本地9家实验室间的原始值散点图 ----------------------
rm(list = setdiff(ls(), c("labels.lab", "colors.sample")))
gc()

local_meta <- fread("./data/multilab/metadata_2025_10labs.csv")
all_tables <- readRDS("./data/multilab/quantdata_list_pep_2025_10labs.rds")
qualipr_tables <- readRDS("./data/multilab/qualidata_list_pep_2025_10labs.rds")

meta_quartet <- fread("./results/tables/2_outlier_madist_quartet.csv")
outliers <- meta_quartet$analysis_id[meta_quartet$is.outlier]

filtered_tables <- pblapply(1:6, function(i) {
  tmp_quant_table <- qualipr_tables[[i]] %>%
    distinct(analysis_id, peptide_sequence, protein_id) %>%
    inner_join(., all_tables[[i]], by = c("analysis_id", "peptide_sequence"), 
               relationship = "many-to-many") %>%
    filter(!analysis_id %in% outliers)
  return(tmp_quant_table)
})
names(filtered_tables) <- names(all_tables)

all_cor_tables <- pblapply(filtered_tables, function(tmp_table) {
  all_peps1 <- tmp_table %>%
    mutate_at("lab_id", ~ labels.lab[.]) %>%
    mutate_at("value", ~ ifelse(. == 0, NA, log2(.))) %>%
    reshape2::dcast(., peptide_sequence + protein_id ~ lab_id + tube + injection,
                    value.var = "value", fun.aggregate = sum) %>%
    select(!contains("Lab-0"))
  
  all_peps_cor <- all_peps1[, 3:ncol(all_peps1)] %>%
    mutate_all(~ ifelse(. == 0, NA, .)) %>%
    cor_test %>%
    filter(cor != 1, p < .05) %>%
    mutate(label = ifelse(str_extract(var1, "^Lab-\\d+") == str_extract(var2, "^Lab-\\d+"),
                          "Intra-lab",
                          "Inter-lab"))
  return(all_peps_cor)
})

df_cor <- all_cor_tables %>%
  rbindlist(., idcol = "Sample") %>%
  mutate_at("Sample", ~ ifelse(. %in% c("D5", "D6", "F7", "M8"), paste("Quartet", .), .)) %>%
  mutate_at("Sample", ~ factor(., levels = c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")))

fwrite(df_cor, "./results/tables/3_cor.csv")

pcc_thres <- df_cor %>%
  group_by(label, Sample) %>%
  summarise(pcc_median = median(cor), .groups = "drop")

p3a <- ggplot(df_cor, aes(x = label, y = cor, fill = Sample)) +
  geom_boxplot(outliers = FALSE, width = .7, position = position_dodge(width = 1), 
               color = "black", linewidth = .5) +
  geom_text(aes(x = label, y = pcc_median,
                label = sprintf("%.2f", pcc_median)),
            data = pcc_thres, size = 3.5, vjust = -.5, color = "black",
            position = position_dodge(width = 1)) +
  scale_y_continuous(name = "Pearson Correlation Coefficient",
                     limits = c(0, 1), breaks = seq(0, 1, .1)) +
  scale_fill_manual(values = colors.sample) +
  labs(title = "Raw") +
  theme_bw() +
  theme(legend.position = "right",
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        strip.text = element_text(size = 14, face = "bold", color = "black", margin = margin(0.3, 0.3, 0.3, 0.3, "cm")),
        strip.background = element_blank(),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm"),
        plot.title = element_text(size = 16, color = "black", face = "bold"))

ggsave("./results/figures/supp_figure3a.pdf",
       p3a, width = 8, height = 5, limitsize = FALSE, dpi = 600, device = cairo_pdf)


## Supplementary Figure 3b ----------------------
rm(list = setdiff(ls(), c("labels.lab", "colors.sample")))
gc()

local_tables <- readRDS("./results/tables/2_limma_bylab_bytube.rds")

all_cor_tables <- pblapply(local_tables, function(tmp_table) {
  all_peps1 <- tmp_table %>%
    filter(adj.P.Val < .05, abs(logFC) >= 1) %>%
    mutate_at("lab_id", ~ labels.lab[.]) %>%
    reshape2::dcast(., peptide_sequence + protein_id ~ lab_id + tube, value.var = "logFC") %>%
    select(!contains("Lab-0"))
  
  all_peps_cor <- all_peps1[, 3:ncol(all_peps1)] %>%
    cor_test %>%
    filter(cor != 1, p < .05) %>%
    mutate(label = ifelse(str_extract(var1, "^Lab-\\d+") == str_extract(var2, "^Lab-\\d+"),
                          "Intra-lab",
                          "Inter-lab"))
  return(all_peps_cor)
})

df_cor <- all_cor_tables %>%
  rbindlist(., idcol = "Sample Pair") %>%
  mutate_at("Sample Pair", ~ factor(., levels = c("D5/D6", "F7/D6", "M8/D6", "HeLa/HEK293T")))

pcc_thres <- df_cor %>%
  group_by(label, `Sample Pair`) %>%
  summarise(pcc_median = median(cor), .groups = "drop")

colors.ratio <- c("#4CC3D9", "#FFC65D", "#F16745", "#E7298A")
names(colors.ratio) <- c("D5/D6", "F7/D6", "M8/D6", "HeLa/HEK293T")

p3b <- ggplot(df_cor, aes(x = label, y = cor, fill = `Sample Pair`)) +
  geom_boxplot(outliers = FALSE, width = .7, position = position_dodge(width = 1), 
               color = "black", linewidth = .5) +
  geom_text(aes(x = label, y = pcc_median,
                label = sprintf("%.2f", pcc_median)),
            data = pcc_thres, size = 3.5, vjust = -.5, color = "black",
            position = position_dodge(width = 1)) +
  scale_y_continuous(name = "Pearson Correlation Coefficient",
                     limits = c(0, 1), breaks = seq(0, 1, .1)) +
  scale_fill_manual(values = colors.ratio) +
  labs(title = "SRR-scaled") +
  theme_bw() +
  theme(legend.position = "right",
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        strip.text = element_text(size = 14, face = "bold", color = "black", margin = margin(0.3, 0.3, 0.3, 0.3, "cm")),
        strip.background = element_blank(),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm"),
        plot.title = element_text(size = 16, color = "black", face = "bold"))

ggsave("./results/figures/supp_figure3b.pdf",
       p3b, width = 8, height = 4.5, limitsize = FALSE, dpi = 600, device = cairo_pdf)

## Supplementary Figure 3: D5/F7/M8/HeLa/HEK293T展示投票+PEP阈值 --------------
library(tidyverse)
library(data.table)

df_subPEP2 <- curvePEP_tables %>%
  rbindlist(., idcol = "sample") %>%
  mutate_at("sample", ~ ifelse(. %in% c("D5", "D6", "F7", "M8"), paste("Quartet", .), .)) %>%
  filter(!sample %in% c("Quartet D6")) %>%
  mutate(pep_trans = -log10(pep_min)) %>%
  mutate_at("pep_trans", ~ ifelse(is.infinite(.), 400, .))

df_subPEP3 <- df_subPEP2 %>%
  group_by(sample, group) %>%
  summarise(feature_n = length(unique(peptide_sequence)), .groups = "drop")

ymax_pep <- max(df_subPEP2$pep_trans, na.rm = TRUE)
ymax_n   <- max(df_subPEP3$feature_n, na.rm = TRUE)
scale_factor <- ymax_pep / ymax_n

no_thres <- data.frame(
  sample = c("Quartet D5", "Quartet F7", "Quartet M8", "HeLa", "HEK293T"),
  number_thres = c(20817, 23981, 23948, 20983, 23962)
)

colors.sample <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745", "#E7298A", "#4D9221")
names(colors.sample) <- c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")

p_figure3c <- ggplot(df_subPEP2, aes(x = group, y = pep_trans)) +
  geom_boxplot(aes(fill = group), outliers = FALSE) +
  facet_wrap(~ sample, ncol = 1, strip.position = "top") +
  geom_line(data = df_subPEP3, aes(x = group, y = feature_n * scale_factor, group = 1),
            inherit.aes = FALSE, linewidth = 1, color = "black") +
  geom_point(data = df_subPEP3, aes(x = group, y = feature_n * scale_factor),
             inherit.aes = FALSE, size = 2.5, color = "black") +
  geom_hline(aes(yintercept = number_thres * scale_factor), data = no_thres,
             linetype = "dashed", color = "#CB181D", alpha = 0.6) +
  geom_text(aes(x = "\u2265 3 Lab(s)", y = number_thres * scale_factor, 
                label = sprintf("N = %s", format(number_thres, big.mark = ",", scientific = FALSE))),
            data = no_thres, color = "#CB181D", size = 4.5, vjust = -1) +
  scale_fill_brewer(palette = "Blues") +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  scale_y_continuous(
    name = expression(-lg(PEP)~Value),
    breaks = seq(0, 400, 100),
    expand = c(0.02, 0.01),
    sec.axis = sec_axis(~ (.) / scale_factor,
                        name = "Number of Retained Peptides",
                        breaks = seq(0, 40000, 10000))
  ) +
  labs(fill = "Filter") +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 14, face = "bold", color = "black", margin = margin(0.3, 0.3, 0.3, 0.3, "cm")),
    strip.background = element_blank(),
    axis.title.y = element_text(size = 14, color = "black", face = "bold"),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1, vjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = "cm"),
    plot.title = element_text(size = 16, color = "black", face = "bold")
  )

p_figure3c

ggsave(
  "./results/figures/supp_figure3.pdf", 
  p_figure3c, 
  width = 9.5, 
  height = 12, 
  dpi = 600
)

## Supplementary figure 6-7: D5/F7/M8/HeLa/HEK293T肽段长度/亲疏水性-----------
rm(list = ls())
gc()
passFDR0.01_pep_tables <- readRDS("./results/tables/1_qualiprop_list_PEPfdr0.01.rds")

kd_scale <- c(A =  1.8,  R = -4.5, N = -3.5, D = -3.5, C =  2.5,
              Q = -3.5,  E = -3.5, G = -0.4, H = -3.2, I =  4.5,
              L =  3.8,  K = -3.9, M =  1.9, F =  2.8, P = -1.6,
              S = -0.8,  T = -0.7, W = -0.9, Y = -1.3, V =  4.2)

length_gravy_tables <- pblapply(passFDR0.01_pep_tables, function(table_tmp) {
  sub_pep_table <- table_tmp %>%
    distinct(peptide_sequence, length) %>%
    mutate(gravy_score = sapply(peptide_sequence, function(x) {
      aa <- strsplit(x, split = "")[[1]]
      bb <- mean(kd_scale[aa], na.rm = TRUE)
      return(bb)
    }))
  return(sub_pep_table)
})

df_sub1 <- length_gravy_tables %>%
  rbindlist(., idcol = "sample") %>%
  mutate_at("sample", ~ ifelse(. %in% c("HEK293T", "HeLa"), ., paste("Quartet", .))) %>%
  filter(sample %in% c("Quartet D5", "Quartet F7", "Quartet M8", "HEK293T")) %>%
  mutate_at("sample", ~ factor(., levels = c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")))

df_sub2 <- df_sub1 %>%
  group_by(sample) %>%
  summarise(gravy_median = median(gravy_score),
            gravy_peak = {
              d <- density(gravy_score, bw = 1, na.rm = TRUE)
              d$x[which.max(d$y)]
            })

colors.sample <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745", "#E7298A", "#4D9221")
names(colors.sample) <- c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")

set.seed(2026)
p_sfigure6 <- ggplot(df_sub1, aes(x = gravy_score)) +
  geom_density(aes(fill = sample), alpha = 0.8) +
  geom_vline(data = df_sub2, mapping = aes(xintercept = gravy_peak),
             linetype = "dashed", color = "black", alpha = 0.6) +
  facet_grid(rows = vars(sample)) +
  scale_fill_manual(values = colors.sample) +
  scale_x_continuous(name = "GRAVY Score") +
  scale_y_continuous(name = "Density") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold", color = "black", margin = margin(0.3, 0.3, 0.3, 0.3, "cm")),
    strip.background = element_blank(),
    axis.title.x = element_text(size = 14, color = "black", face = "bold"),
    axis.title.y = element_text(size = 14, color = "black", face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = "cm"),
    plot.title = element_text(size = 16, color = "black", face = "bold")
  )

print(p_sfigure6)
ggsave("./results/figures/supp_figure6.pdf", p_sfigure7, width = 7, height = 9, dpi = 600, device = cairo_pdf)

df_sub3 <- df_sub1 %>%
  group_by(sample) %>%
  summarise(
    length_median = median(length),
    length_peak = {
      d <- density(length, bw = 1, na.rm = TRUE)
      d$x[which.max(d$y)]
    }
  )
set.seed(2026)
p_sfigure7 <- ggplot(df_sub1, aes(x = length)) +
  geom_density(aes(fill = sample), alpha = 0.8, bw = 1) +
  geom_vline(
    data = df_sub3,
    mapping = aes(xintercept = length_peak),
    linetype = "dashed",
    color = "black",
    alpha = 0.6
  ) +
  facet_grid(rows = vars(sample)) +
  scale_fill_manual(values = colors.sample) +
  scale_x_continuous(name = "Sequence Length") +
  scale_y_continuous(name = "Density") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold", color = "black", margin = margin(0.3, 0.3, 0.3, 0.3, "cm")),
    strip.background = element_blank(),
    axis.title.x = element_text(size = 14, color = "black", face = "bold"),
    axis.title.y = element_text(size = 14, color = "black", face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = "cm"),
    plot.title = element_text(size = 16, color = "black", face = "bold")
  )
print(p_sfigure7)
ggsave(
  "./results/figures/supp_figure7.pdf",
  p_sfigure8,
  width = 7,
  height = 9,
  dpi = 600,
  device = cairo_pdf
)

## Supplementary figure 8: D5/F7/M8/HeLa/HEK293T肽段 m/z -----------
rm(list = ls())
gc()
passFDR0.01_pep_tables <- readRDS("./results/tables/1_qualiprop_list_PEPfdr0.01.rds")

mz_tables <- pblapply(passFDR0.01_pep_tables, function(table_tmp) {
  sub_pep_table <- table_tmp %>%
    distinct(peptide_sequence, mz_ratio)
  return(sub_pep_table)
})

df_sub1 <- mz_tables %>%
  rbindlist(., idcol = "sample") %>%
  mutate_at("sample", ~ ifelse(. %in% c("HEK293T", "HeLa"), ., paste("Quartet", .))) %>%
  filter(sample %in% c("Quartet D5", "Quartet F7", "Quartet M8", "HEK293T")) %>%
  mutate_at("sample", ~ factor(., levels = c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")))

df_sub2 <- df_sub1 %>%
  group_by(sample) %>%
  summarise(
    mz_median = median(mz_ratio),
    mz_peak = {
      d <- density(mz_ratio, bw = 1, na.rm = TRUE)
      d$x[which.max(d$y)]
    }
  )

colors.sample <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745", "#E7298A", "#4D9221")
names(colors.sample) <- c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")

set.seed(2026)
p_sfigure8 <- ggplot(df_sub1, aes(x = mz_ratio)) +
  geom_density(aes(fill = sample), alpha = 0.8, bw = 20) +
  geom_vline(
    data = df_sub2,
    mapping = aes(xintercept = mz_peak),
    linetype = "dashed",
    color = "black",
    alpha = 0.6
  ) +
  facet_grid(rows = vars(sample)) +
  scale_fill_manual(values = colors.sample, name = "RM Group") +
  scale_x_continuous(name = "m/z") +
  scale_y_continuous(name = "Density") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold", color = "black", margin = margin(0.3, 0.3, 0.3, 0.3, "cm")),
    strip.background = element_blank(),
    axis.title.x = element_text(size = 14, color = "black", face = "bold"),
    axis.title.y = element_text(size = 14, color = "black", face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = "cm"),
    plot.title = element_text(size = 16, color = "black", face = "bold")
  )

print(p_sfigure8)
ggsave(
  "./results/figures/supp_figure8.pdf",
  p_sfigure9,
  width = 7,
  height = 9,
  dpi = 600,
  device = cairo_pdf
)

## supplementary figure 10a----------------------
rm(list = ls())
gc()
panel.cor <- function(x, y, digits = 2, prefix = "", ...){
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use = "pairwise.complete.obs", method = "pearson")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  test <- cor.test(x, y)
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))
  text(0.5, 0.5, txt, cex = 0.8 / strwidth(txt))
  text(.7, .9, Signif, cex = 2)
}

panel.smooth <- function(x, y, col = "black", bg = NA, pch = 18, 
                         cex = 0.4, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = col.smooth, ...)
}

labels.lab <- c("Lab-1", "Lab-2", "Lab-0", "Lab-3", "Lab-4", "Lab-5", "Lab-6", "Lab-7", "Lab-8", "Lab-9")
names(labels.lab) <- c("qinglian_bio", "phoenix", "zhejiang_university", "fudan_university",
                       "national_institute_of_methodology", "biotech_pack", "thermofisher_shanghai",
                       "omicsolution", "cas_tianjin", "academy_of_chinese_medical_sciences")

local_meta <- fread("./data/multilab/metadata_2025_10labs.csv")
all_tables <- readRDS("./data/multilab/quantdata_list_pep_2025_10labs.rds")
qualipr_tables <- readRDS("./data/multilab/qualidata_list_pep_2025_10labs.rds")

meta_quartet <- fread("./results/tables/2_outlier_madist_quartet.csv")
outliers <- meta_quartet$analysis_id[meta_quartet$is.outlier]

filtered_tables <- pblapply(1:6, function(i) {
  tmp_quant_table <- qualipr_tables[[i]] %>%
    distinct(analysis_id, peptide_sequence, protein_id) %>%
    inner_join(., all_tables[[i]], by = c("analysis_id", "peptide_sequence"), 
               relationship = "many-to-many") %>%
    filter(!analysis_id %in% outliers)
  return(tmp_quant_table)
})
names(filtered_tables) <- names(all_tables)

stat_pep_tables <- pblapply(filtered_tables, function(tmp_table) {
  stat_tmp_table <- tmp_table %>%
    ungroup() %>%
    summarise(peptide_n = length(unique(peptide_sequence)),
              protein_n = length(unique(protein_id)), .groups = "drop")
  return(stat_tmp_table)
})
stat_df0 <- rbindlist(stat_pep_tables, idcol = "sample")
stat_df0$tier <- "All quantified queries."

all_peps <- filtered_tables[[5]] %>%
  mutate_at("lab_id", ~ labels.lab[.]) %>%
  filter(lab_id != "Lab-0") %>%  # 过滤 ZJU / Lab-0
  mutate_at("value", ~ ifelse(. == 0, NA, log2(.))) %>%
  reshape2::dcast(., tube + peptide_sequence + protein_id ~ lab_id,
                  value.var = "value", fun.aggregate = median)

lab_cols <- setdiff(colnames(all_peps), c("tube", "peptide_sequence", "protein_id", "Lab-0"))
axis_lim <- range(all_peps[, lab_cols], na.rm = TRUE) * 1.1

pdf("./results/figures/supp_figure2e_1.pdf", width = 14, height = 14)
pairs(all_peps[, lab_cols], lower.panel = panel.smooth, upper.panel = panel.cor,
      ylim = axis_lim, xlim = axis_lim)
dev.off()

all_cor_tables <- pblapply(filtered_tables, function(tmp_table) {
  
  all_peps1 <- tmp_table %>%
    mutate_at("lab_id", ~ labels.lab[.]) %>%
    filter(lab_id != "Lab-0") %>%  # 过滤 ZJU / Lab-0
    mutate_at("value", ~ ifelse(. == 0, NA, log2(.))) %>%
    reshape2::dcast(., peptide_sequence + protein_id ~ lab_id + tube + injection,
                    value.var = "value", fun.aggregate = sum) %>%
    select(!contains("Lab-0"))
  
  all_peps_cor <- all_peps1[, 3:ncol(all_peps1)] %>%
    mutate_all(~ ifelse(. == 0, NA, .)) %>%
    cor_test() %>%
    filter(cor != 1, p < .05) %>%
    filter(!str_detect(var1, "^Lab-0") & !str_detect(var2, "^Lab-0")) %>%
    mutate(label = ifelse(str_extract(var1, "^Lab-\\d+") == str_extract(var2, "^Lab-\\d+"),
                          "Intra-lab",
                          "Inter-lab"))
  
  return(all_peps_cor)
})

df_cor <- all_cor_tables %>%
  rbindlist(., idcol = "Sample") %>%
  mutate_at("Sample", ~ ifelse(. %in% c("D5", "D6", "F7", "M8"), paste("Quartet", .), .)) %>%
  mutate_at("Sample", ~ factor(., levels = c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")))

fwrite(df_cor, "./results/tables/3_cor.csv")

pcc_thres <- df_cor %>%
  group_by(label, Sample) %>%
  summarise(pcc_median = median(cor), .groups = "drop")

colors.sample <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745", "#E7298A", "#4D9221")
names(colors.sample) <- c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")

p_sfigure10a <- ggplot(df_cor, aes(x = label, y = cor, fill = Sample)) +
  geom_boxplot(outliers = FALSE, width = .7, position = position_dodge(width = 1)) +
  geom_text(aes(x = label, y = pcc_median,
                label = sprintf("%.2f", pcc_median)),
            data = pcc_thres, size = 3.5, vjust = -.5, color = "black",
            position = position_dodge(width = 1)) +
  scale_y_continuous(name = "Pearson Correlation Coefficient",
                     limits = c(0, 1), breaks = seq(0, 1, .1)) +
  scale_fill_manual(values = colors.sample, name = "Sample") +
  labs(title = "Raw") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, face = "bold", color = "black"),
    strip.text = element_text(size = 14, face = "bold", color = "black", margin = margin(0.3, 0.3, 0.3, 0.3, "cm")),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, color = "black", face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 0, hjust = 0.5, vjust = 0.5),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = "cm"),
    plot.title = element_text(size = 16, color = "black", face = "bold")
  )

print(p_sfigure10a)

ggsave(
  "./results/figures/supp_figure10a.pdf",
  p_sfigure9,
  width = 8,
  height = 6,
  dpi = 600,
  device = cairo_pdf
)

## supplementary figure 10b----------------------
rm(list = ls())
gc()
panel.cor <- function(x, y, digits = 2, prefix = "", ...){
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use = "pairwise.complete.obs", method = "pearson")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  test <- cor.test(x, y)
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))
  text(0.5, 0.5, txt, cex = 0.8 / strwidth(txt))
  text(.7, .9, Signif, cex = 2)
}
## supplementary figure 9----------------------
panel.smooth <- function(x, y, col = "black", bg = NA, pch = 18, 
                         cex = 0.4, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = col.smooth, ...)
}
labels.lab <- c("Lab-1", "Lab-2", "Lab-0", "Lab-3", "Lab-4", "Lab-5", "Lab-6", "Lab-7", "Lab-8", "Lab-9")
names(labels.lab) <- c("qinglian_bio", "phoenix", "zhejiang_university", "fudan_university",
                       "national_institute_of_methodology", "biotech_pack", "thermofisher_shanghai",
                       "omicsolution", "cas_tianjin", "academy_of_chinese_medical_sciences")

local_meta <- fread("./data/multilab/metadata_2025_10labs.csv")
local_tables <- readRDS("./results/tables/2_limma_bylab_bytube.rds")

all_peps <- local_tables[[3]] %>%
  filter(adj.P.Val < .05, abs(logFC) >= 1) %>%
  mutate_at("lab_id", ~ labels.lab[.]) %>%
  filter(lab_id != "Lab-0") %>%  # 过滤 ZJU / Lab-0
  reshape2::dcast(., tube + peptide_sequence + protein_id ~ lab_id, value.var = "logFC")

lab_cols <- setdiff(colnames(all_peps), c("tube", "peptide_sequence", "protein_id", "Lab-0"))
axis_lim <- range(all_peps[, lab_cols], na.rm = TRUE) * .9

pdf("./Update/supp_figure9.pdf", width = 14, height = 14)
pairs(all_peps[, lab_cols], lower.panel = panel.smooth, upper.panel = panel.cor,
      ylim = axis_lim, xlim = axis_lim)
dev.off()

all_cor_tables <- pblapply(local_tables, function(tmp_table) {
  
  all_peps1 <- tmp_table %>%
    filter(adj.P.Val < .05, abs(logFC) >= 1) %>%
    mutate_at("lab_id", ~ labels.lab[.]) %>%
    filter(lab_id != "Lab-0") %>%  # 过滤 ZJU / Lab-0
    reshape2::dcast(., peptide_sequence + protein_id ~ lab_id + tube, value.var = "logFC") %>%
    select(!contains("Lab-0"))
  
  all_peps_cor <- all_peps1[, 3:ncol(all_peps1)] %>%
    cor_test() %>%
    filter(cor != 1, p < .05) %>%
    filter(!str_detect(var1, "^Lab-0") & !str_detect(var2, "^Lab-0")) %>%
    mutate(label = ifelse(str_extract(var1, "^Lab-\\d+") == str_extract(var2, "^Lab-\\d+"),
                          "Intra-lab",
                          "Inter-lab"))
  
  return(all_peps_cor)
})

df_cor <- all_cor_tables %>%
  rbindlist(., idcol = "Sample Pair") %>%
  mutate_at("Sample Pair", ~ factor(., levels = c("D5/D6", "F7/D6", "M8/D6", "HeLa/HEK293T")))

pcc_thres <- df_cor %>%
  group_by(label, `Sample Pair`) %>%
  summarise(pcc_median = median(cor), .groups = "drop")

colors.sample <- c("#4CC3D9", "#FFC65D", "#F16745", "#E7298A")
names(colors.sample) <- c("D5/D6", "F7/D6", "M8/D6", "HeLa/HEK293T")
p <- ggplot(df_cor, aes(x = label, y = cor, fill = `Sample Pair`)) +
  geom_boxplot(outliers = FALSE, width = .7, position = position_dodge(width = 1)) +
  geom_text(aes(x = label, y = pcc_median,
                label = sprintf("%.2f", pcc_median)),
            data = pcc_thres, size = 3.5, vjust = -.5, color = "black",
            position = position_dodge(width = 1)) +
  scale_y_continuous(name = "Pearson Correlation Coefficient",
                     limits = c(0, 1), breaks = seq(0, 1, .1)) +
  scale_fill_manual(values = colors.sample, name = "Sample Pair") +
  labs(title = "SRR-scaled") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, face = "bold", color = "black"),
    strip.text = element_text(size = 14, face = "bold", color = "black", margin = margin(0.3, 0.3, 0.3, 0.3, "cm")),
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, color = "black", face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 0, hjust = 0.5, vjust = 0.5),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = "cm"),
    plot.title = element_text(size = 16, color = "black", face = "bold")
  )
print(p)
ggsave(
  "./results/figures/supp_figure10b.pdf",
  p,
  width = 8,
  height = 5,
  dpi = 600,
  device = cairo_pdf
)