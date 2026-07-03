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


labels.lab <- c("Lab-1", "Lab-2", "Lab-0", "Lab-3", "Lab-4", "Lab-5", "Lab-6", "Lab-7", "Lab-8", "Lab-9")
names(labels.lab) <- c("QLB", "NCP", "ZJU", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")

colors.sample <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745", "#E7298A", "#4D9221")
names(colors.sample) <- c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")


## supplementary figure 1a: 作图比较鉴定数目---------
rm(list = setdiff(ls(), c("labels.lab", "colors.sample")))
gc()
all_pep <- fread("./results/tables/3_count.csv")

all_long <- all_pep %>%
  # filter(sample %in% "HeLa", tube == 1) %>%
  filter(!lab %in% "ZJU") %>%
  select(!Precursors.Identified) %>%
  mutate(source = ifelse(grepl("Lab", lab), "Public", "Local")) %>%
  reshape2::melt(., id = c(1:3, 6:7))

all_long2 <- all_long %>%
  filter(source %in% "Local") %>%
  mutate_at("sample", ~ ifelse(. %in% c("HeLa", "HEK293T"), ., paste("Quartet", .))) %>%
  mutate_at("sample", ~ factor(., levels = c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T"))) %>%
  mutate_at("lab", ~ labels.lab[.])

p <- ggplot(all_long2, aes(x = sample, y = value)) +
  stat_summary(aes(fill = sample), fun = mean, geom = "bar", width = .7,
               color = "black", linewidth = .3) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  stat_summary(aes(label = scales::label_comma(accuracy = 1)(after_stat(y)),
                   vjust = ifelse(variable %in% "Proteins.Identified", -2, -3)),
               fun = "mean", geom = "text", size = 4.5) +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 14, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_blank(),
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm")) +
  scale_y_continuous(n.breaks = 6, name = "Number of Identifications",
                     expand = expansion(mult = c(0, 0.2)),
                     labels = scales::label_comma()) +
  scale_fill_manual(values = colors.sample) +
  facet_wrap(~ variable, scales = "free_y")

ggsave(paste("./results/figures/supp_figure1a.pdf", sep = ""), p, height = 4, width = 14)


## supplementary figure 1b: 作图比较 CV ----------------------
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
  # group_by(variable, source, lab, sample, tube) %>%
  # summarise_at("cv", median) %>%
  mutate_at("sample", ~ ifelse(. %in% c("HeLa", "HEK293T"), as.character(.), paste("Quartet", .))) %>%
  mutate_at("sample", ~ factor(., levels = c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")))

cv_thres <- df_cv_test %>%
  group_by(sample, variable) %>%
  summarise(cv_median = median(cv),
            cv_mean = mean(cv),
            cv_sd = sd(cv))

p <- ggplot(df_cv_test, aes(x = sample, y = cv)) +
  geom_boxplot(aes(fill = sample), outlier.size = .1, width = .7) +
  geom_text(aes(x = sample, y = cv_median,
                label = sprintf("%.2f%%", cv_median * 100)),
            data = cv_thres, size = 3, color = "white") +
  scale_fill_manual(values = colors.sample) +
  scale_y_continuous(n.breaks = 8, name = "Coefficient of variation (CV)",
                     labels = scales::percent) +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 14, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_blank(),
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm")) +
  facet_wrap(~ variable, scales = "free_y")

ggsave(paste("./results/figures/supp_figure1b.pdf", sep = ""), p, height = 4, width = 14)


## supplementary figure 2a: 作图比较 HeLa 鉴定数目---------
rm(list = setdiff(ls(), c("labels.lab", "colors.sample", "all_long")))
gc()

all_pep <- fread("./results/tables/3_count.csv")

all_long <- all_pep %>%
  filter(sample %in% "HeLa", tube == 1, !lab %in% "ZJU") %>%
  mutate(source = ifelse(grepl("Lab", lab), "Public", "Local")) %>%
  reshape2::melt(., id = c(1:3, 7:8)) %>%
  filter(!variable %in% "Precursors.Identified")

colors.source <- c("Local" = "#A6CEE3", "Public" = "#1F78B4")

p_format <- function(p_value) {
  # 设置一个非常小的阈值
  threshold <- .0001
  # 如果P值小于阈值，则使用 "<" 符号
  if (p_value < threshold) {
    return(paste0("**** (P < 0.0001)"))
  } else if (p_value < .001) {
    # 否则，直接格式化P值并添加"P = "前缀
    return(paste0("*** (P = ", format(p_value, digits = 2), ")", sep = ""))
  } else if (p_value < .01) {
    # 否则，直接格式化P值并添加"P = "前缀
    return(paste0("** (P = ", format(p_value, digits = 2), ")", sep = ""))
  } else if (p_value < .05) {
    # 否则，直接格式化P值并添加"P = "前缀
    return(paste0("* (P = ", format(p_value, digits = 2), ")", sep = ""))
  } else {
    return(paste0("NS (P = ", format(p_value, digits = 2), ")", sep = ""))
  }
}

p <- ggplot(all_long, aes(x = source, y = value)) +
  stat_summary(aes(fill = source), fun = mean, geom = "bar", width = .7,
               color = "black", linewidth = .3) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  stat_summary(aes(label = scales::label_comma(accuracy = 1)(after_stat(y)),
                   vjust = ifelse(grepl("Protein", variable), -1.5, -4)),
               fun = "mean", geom = "text", size = 4.5) +
  # ggpubr::geom_signif(comparisons = list(c("Local", "Public")),
  #                     map_signif_level = p_format,tip_length = .05,
  #                     textsize = 4, test = "t.test") +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 14, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_blank(),
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm")) +
  scale_y_continuous(n.breaks = 6, name = "Number of Identifications",
                     expand = expansion(mult = c(0, 0.2)),
                     labels = scales::label_comma()) +
  scale_fill_manual(values = colors.source) +
  ggh4x::facet_grid2(~ variable, scales = "free_y");p

ggsave("./results/figures/supp_figure2a.pdf", p, height = 4, width = 8)


## supplementary figure 2b: 作图比较 HeLa CV ----------------------
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
  # group_by(variable, source, lab, sample, tube) %>%
  # summarise_at("cv", median) %>%
  mutate_at("variable", ~ gsub("s.Quantified", "-level", .)) %>%
  mutate_at("variable", ~ factor(., levels = c("Precursor-level", "Peptide-level", "Protein-level")))

cv_thres <- df_cv_test %>%
  group_by(source, variable) %>%
  summarise(cv_median = median(cv))

p <- ggplot(df_cv_test, aes(x = reorder(source, cv), y = cv)) +
  geom_boxplot(aes(fill = source), outlier.size = .5, width = .7) +
  geom_hline(aes(yintercept = cv_median, colour = source), cv_thres, lty = 2) +
  geom_text(aes(x = source, y = cv_median,
                label = sprintf("Median = %.2f%%", cv_median * 100)),
            data = cv_thres, size = 3, vjust = -4) +
  scale_fill_manual(values = c("Local" = "#A6CEE3", "Public" = "#1F78B4")) +
  scale_color_manual(values = c("Local" = "darkred", "Public" = "#276419")) +
  scale_y_continuous(limits = c(0, 2), n.breaks = 8, name = "Coefficient of variation (CV)",
                     labels = scales::percent) +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 14, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_blank(),
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm")) +
  ggh4x::facet_grid2( ~ variable, scales = "fixed");p

ggsave("./results/figures/supp_figure2b.pdf", p, height = 4, width = 8)


## supplementary figure 2c：定性验证-检查每个实验室的定性 ----------------------
rm(list = setdiff(ls(), c("labels.lab", "colors.sample", "all_long")))
gc()

meta_hela <- fread("./results/tables/3_outlier_madist_hela.csv")
all_tables <- readRDS("./results/tables/3_quantdata_list_pep_2025_21labs.rds")
all_tables <- all_tables %>% rbindlist %>% split(., ~ lab_id)
all_tables <- all_tables[1:20]

## 测试：1个单元管内3技术重复间的定性一致性
intra_cor_tables <- pblapply(all_tables, function(tmp_table) {
  
  test_table <- tmp_table
  test_table_cor <- test_table %>%
    reshape2::acast(., peptide_sequence ~ analysis_id,
                    fun.aggregate = length) %>%
    cor
  
  aa <- test_table_cor[1:3, 1:3]
  df_cor <- data.frame(mean_cor = median(aa[upper.tri(aa)]))
  
  return(df_cor)
})
df_cor_final <- rbindlist(intra_cor_tables, idcol = "lab_id")

## 正式检查实验室间定性
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

colors.sample <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745", "#E7298A", "#4D9221")
names(colors.sample) <- c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")
# labels.lab <- c("QLB", "NCP", "ZJU", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")
labels.lab <- c("Lab-1", "Lab-2", "Lab-0", "Lab-3", "Lab-4", "Lab-5", "Lab-6", "Lab-7", "Lab-8", "Lab-9")
names(labels.lab) <- c("qinglian_bio", "phoenix", "zhejiang_university", "fudan_university",
                       "national_institute_of_methodology", "biotech_pack", "thermofisher_shanghai",
                       "omicsolution", "cas_tianjin", "academy_of_chinese_medical_sciences")

p <- ggplot(df_inter, aes(x = reorder(lab_id, -proportion), y = proportion)) +
  geom_col(fill = "#A6CEE3", width = .7, color = "black", linewidth = .3) +
  geom_text(aes(label = scales::percent(round(proportion, digits = 4))), nudge_y = .04) +
  geom_hline(yintercept = .6, lty = 2, col = "red", size = .6) +
  scale_y_continuous(label = scales::percent, n.breaks = 10, name = "Overlapping Proportion") +
  scale_x_discrete(labels = labels.lab) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm"));p


ggsave("./results/figures/supp_figure2c.pdf", p, width = 8, height = 4, limitsize = FALSE)


## supplementary figure 2d：定量验证-检查每个实验室的CV ----------------------
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

colors.sample <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745", "#E7298A", "#4D9221")
names(colors.sample) <- c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")
# labels.lab <- c("QLB", "NCP", "ZJU", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS",
#                 paste("Lab", 1:11, sep = "_"))
labels.lab <- c("Lab-1", "Lab-2", "Lab-0", "Lab-3", "Lab-4", "Lab-5", "Lab-6", "Lab-7", "Lab-8", "Lab-9",
                paste("Lab", 1:11, sep = "_"))
names(labels.lab) <- c("qinglian_bio", "phoenix", "zhejiang_university", "fudan_university",
                       "national_institute_of_methodology", "biotech_pack", "thermofisher_shanghai",
                       "omicsolution", "cas_tianjin", "academy_of_chinese_medical_sciences",
                       paste("Lab", 1:11, sep = "_"))

df_cv <- cv_tables %>%
  rbindlist %>%
  # filter(cv < 1) %>%
  mutate_at("lab_id", ~ factor(., levels = names(labels.lab)))

cv_threshold <- df_cv %>%
  filter(grepl("Lab", lab_id)) %>%
  group_by(lab_id) %>%
  summarise_at("cv", median) %>%
  pull(cv) %>%
  max

df_cv_test <- df_cv %>% filter(!grepl("Lab", lab_id))

cv_thres <- df_cv_test %>%
  group_by(lab_id) %>%
  summarise(cv_median = median(cv))

p <- ggplot(df_cv_test, aes(x = reorder(lab_id, cv, median), y = cv)) +
  geom_boxplot(fill = "#A6CEE3", outliers = FALSE, width = .7) +
  geom_hline(yintercept = cv_threshold, lty = 2, col = "red", size = .6) +
  geom_text(aes(x = lab_id, y = cv_median,
                label = sprintf("%.2f%%", cv_median * 100)),
            data = cv_thres, size = 3.5, vjust = 1.5, color = "black") +
  annotate("text", x = "national_institute_of_methodology", y = .2,
           label = "CV = 14.85%", color = "red", size = 5) +
  scale_y_continuous(label = scales::percent, name = "Coefficient of variation (CV)",
                     limits = c(0, 1)) +
  scale_x_discrete(labels = labels.lab) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm"));p

ggsave("./results/figures/supp_figure2d.pdf",
       p, width = 8, height = 4, limitsize = FALSE)


## supplementary figure 3a：原始值验证-检查本地10家实验室间的原始值散点图 ----------------------
rm(list = ls())
gc()

panel.cor <- function(x, y, digits = 2, prefix = "", ...){
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use = "pairwise.complete.obs", method = "pearson")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  test <- cor.test(x,y)
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))
  text(0.5, 0.5, txt, cex = 0.8 / strwidth(txt))
  text(.7, .9, Signif, cex = 2)
}

panel.smooth<-function (x, y, col = "black", bg = NA, pch = 18, 
                        cex = 0.4, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = col.smooth, ...)
}

# labels.lab <- c("QLB", "NCP", "ZJU", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")
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
stat_df0$tier = "All quantified queries."

## 相关性散点图矩阵
all_peps <- filtered_tables[[5]] %>%
  mutate_at("lab_id", ~ labels.lab[.]) %>%
  mutate_at("value", ~ ifelse(. == 0,  NA, log2(.))) %>%
  reshape2::dcast(., tube + peptide_sequence + protein_id ~ lab_id,
                  value.var = "value", fun.aggregate = median) %>%
  select(!`Lab-0`)

axis_lim <- range(all_peps[, 4:12], na.rm = TRUE) * 1.1

pdf("./results/figures/supp_figure2e_1.pdf", width = 14, height = 14)

pairs(all_peps[, 4:12], lower.panel = panel.smooth, upper.panel = panel.cor,
      ylim = axis_lim, xlim = axis_lim)

dev.off()

## 相关性boxplot
all_cor_tables <- pblapply(filtered_tables, function(tmp_table) {
  
  all_peps1 <- tmp_table %>%
    mutate_at("lab_id", ~ labels.lab[.]) %>%
    mutate_at("value", ~ ifelse(. == 0,  NA, log2(.))) %>%
    reshape2::dcast(., peptide_sequence + protein_id ~ lab_id + tube + injection,
                    value.var = "value", fun.aggregate = sum) %>%
    select(!contains("Lab-0"))
  
  all_peps_cor <- all_peps1[, 3:ncol(all_peps1)] %>%
    mutate_all(~ ifelse(. == 0,  NA, .)) %>%
    cor_test %>%
    filter(cor != 1, p <.05) %>%
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
  summarise(pcc_median = median(cor))

colors.sample <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745", "#E7298A", "#4D9221")
names(colors.sample) <- c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")

p <- ggplot(df_cor, aes(x = label, y = cor, fill = Sample)) +
  # geom_violin(alpha = .7, position = position_dodge(width = 1)) +
  geom_boxplot(outliers = FALSE, width = .7, position = position_dodge(width = 1)) +
  geom_text(aes(x = label, y = pcc_median,
                label = sprintf("%.2f", pcc_median)),
            data = pcc_thres, size = 3.5, vjust = -.5, color = "black",
            position = position_dodge(width = 1)) +
  scale_y_continuous(name = "Pearson Correlation Coefficient",
                     limits = c(0, 1), breaks = seq(0, 1, .1)) +
  scale_fill_manual(values = colors.sample) +
  labs(title = "Raw") +
  theme_bw() +
  theme(legend.position = "bottom",
        title = element_text(size = 16),
        legend.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm"));p

ggsave("./results/figures/supp_figure3a.pdf",
       p, width = 8, height = 5, limitsize = FALSE)

## supplementary figure 3b：比值验证-检查本地10家实验室间的比值散点图 ----------------------
rm(list = ls())
gc()

panel.cor <- function(x, y, digits = 2, prefix = "", ...){
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use = "pairwise.complete.obs", method = "pearson")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  test <- cor.test(x,y)
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))
  text(0.5, 0.5, txt, cex = 0.8 / strwidth(txt))
  text(.7, .9, Signif, cex = 2)
}

panel.smooth<-function (x, y, col = "black", bg = NA, pch = 18, 
                        cex = 0.4, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = col.smooth, ...)
}

# labels.lab <- c("QLB", "NCP", "ZJU", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")
labels.lab <- c("Lab-1", "Lab-2", "Lab-0", "Lab-3", "Lab-4", "Lab-5", "Lab-6", "Lab-7", "Lab-8", "Lab-9")
names(labels.lab) <- c("qinglian_bio", "phoenix", "zhejiang_university", "fudan_university",
                       "national_institute_of_methodology", "biotech_pack", "thermofisher_shanghai",
                       "omicsolution", "cas_tianjin", "academy_of_chinese_medical_sciences")

local_meta <- fread("./data/multilab/metadata_2025_10labs.csv")
local_tables <- readRDS("./results/tables/2_limma_bylab_bytube.rds")

## 相关性散点图矩阵
all_peps <- local_tables[[3]] %>%
  filter(adj.P.Val < .05, abs(logFC) >= 1) %>%
  mutate_at("lab_id", ~ labels.lab[.]) %>%
  reshape2::dcast(., tube + peptide_sequence + protein_id ~ lab_id, value.var = "logFC") %>%
  select(!`Lab-0`)

axis_lim <- range(all_peps[, 4:12], na.rm = TRUE) * .9

pdf("./results/figures/supp_figure2g_1.pdf", width = 14, height = 14)

pairs(all_peps[, 4:12], lower.panel = panel.smooth, upper.panel = panel.cor,
      ylim = axis_lim, xlim = axis_lim)

dev.off()

## 相关性boxplot
all_cor_tables <- pblapply(local_tables, function(tmp_table) {
  
  all_peps1 <- tmp_table %>%
    filter(adj.P.Val < .05, abs(logFC) >= 1) %>%
    mutate_at("lab_id", ~ labels.lab[.]) %>%
    reshape2::dcast(., peptide_sequence + protein_id ~ lab_id + tube, value.var = "logFC") %>%
    select(!contains("Lab-0"))
  
  all_peps_cor <- all_peps1[, 3:ncol(all_peps1)] %>%
    cor_test %>%
    filter(cor != 1, p <.05) %>%
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
  summarise(pcc_median = median(cor))

colors.sample <- c("#4CC3D9", "#FFC65D", "#F16745", "#E7298A")
names(colors.sample) <- c("D5/D6", "F7/D6", "M8/D6", "HeLa/HEK293T")

p <- ggplot(df_cor, aes(x = label, y = cor, fill = `Sample Pair`)) +
  # geom_violin(alpha = .7, position = position_dodge(width = 1)) +
  geom_boxplot(outliers = FALSE, width = .7, position = position_dodge(width = 1)) +
  geom_text(aes(x = label, y = pcc_median,
                label = sprintf("%.2f", pcc_median)),
            data = pcc_thres, size = 3.5, vjust = -.5, color = "black",
            position = position_dodge(width = 1)) +
  scale_y_continuous(name = "Pearson Correlation Coefficient",
                     limits = c(0, 1), breaks = seq(0, 1, .1)) +
  scale_fill_manual(values = colors.sample) +
  labs(title = "SRR-scaled") +
  theme_bw() +
  theme(legend.position = "bottom",
        title = element_text(size = 16),
        legend.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm"));p

ggsave("./results/figures/supp_figure3b.pdf",
       p, width = 8, height = 4.5, limitsize = FALSE)


## 统计原始肽段/蛋白质数目 -------------------
rm(list = ls())
grouped_seq_tables <- readRDS("./results/tables/1_qualidata_list_pro_coverage_2026_10labs.rds")
grouped_pep_tables <- readRDS("./data/multilab/qualidata_list_pep_2025_10labs.rds")


## 删去ZJU ---------------------
grouped_seq_tables <- pblapply(grouped_seq_tables, function(tmp_table) {
  stat_tmp_table <- tmp_table %>%
    filter(!lab_id %in% c("zhejiang_university"))
  return(stat_tmp_table)
})
grouped_pep_tables <- pblapply(grouped_pep_tables, function(tmp_table) {
  stat_tmp_table <- tmp_table %>%
    filter(!lab_id %in% c("zhejiang_university"))
  return(stat_tmp_table)
})
stat_pep_tables <- pblapply(grouped_pep_tables, function(tmp_table) {
  stat_tmp_table <- tmp_table %>%
    ungroup() %>%
    filter(!lab_id %in% c("zhejiang_univarsity")) %>%
    summarise(peptide_n = length(unique(peptide_sequence)),
              protein_n = length(unique(protein_id)), .groups = "drop")
  return(stat_tmp_table)
})
stat_df0 <- rbindlist(stat_pep_tables, idcol = "sample")
stat_df0$tier = "All queries."

View(grouped_pep_tables[[1]] %>% filter(peptide_sequence %in% "FDSDVGEFR"))


## Figure 2b: 以D6为例展示不同实验室投票+SC阈值下的蛋白质数目 ------------------
curveSC_tables <- pblapply(grouped_seq_tables, function(table_tmp) {
  
  tmp_tables <- mclapply(c(10, 20, 30, 40, 50), function(j) {
    filtered_queries <- table_tmp %>%
      group_by(protein_id, lab_id) %>%
      summarise(pass_n = length(unique(analysis_id[coverage > j&!is.na(coverage)])),
                total_n = length(unique(analysis_id)), .groups = "drop") %>%
      mutate(coverage_label = paste("> ", j, "%", sep = ""))
    
    return(filtered_queries)
  })
  
  sub_pro_table <- tmp_tables %>%
    rbindlist %>%
    filter(total_n - pass_n == 0) %>%
    group_by(protein_id, coverage_label) %>%
    summarise(lab_n = length(unique(lab_id)), .groups = "drop")
  
  filtered_tables_i <- mclapply(1:9, function(j) {
    filtered_queries <- sub_pro_table %>%
      filter(lab_n >= j) %>%
      mutate(group = paste("\u2265", j, "Lab(s)"))
    return(filtered_queries)
  })
  df_all_i <- rbindlist(filtered_tables_i)
  
  return(df_all_i)
})
df_subSC1 <- curveSC_tables %>%
  rbindlist(., idcol = "sample") %>%
  filter(sample %in% c("Quartet D6"))

p_figure2b <- ggplot(df_subSC1) +
  geom_bar(aes(x = group, y = after_stat(count), fill = coverage_label),
           stat = "count", position = "dodge", color = "black") +
  facet_grid(cols = vars(sample)) +
  geom_hline(yintercept = 2074, linetype = "dashed", color = "#CB181D", alpha = 0.6) +
  annotate("text", x = "\u2265 3 Lab(s)", y = 2500, label = "N = 2,074",
           color = "#CB181D", size = 5) +
  scale_fill_brewer(palette = "Blues") +
  scale_x_discrete(expand = c(0.06, 0.05)) +
  scale_y_continuous(n.breaks = 10, expand = c(0.02, 0.01)) +
  labs(
    x = "Consensus Voting",
    y = "Number of Retained Proteins",
    fill = "Filter") +
  theme_bw() +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16));p_figure2b

ggsave("./results/figures/figure2b.pdf", p_figure2b, width = 10, height = 5.3)


## Supplementary figure 4: D5/F7/M8/HeLa/HEK293T展示投票+SC阈值 -------------
df_subSC2 <- curveSC_tables %>%
  rbindlist(., idcol = "sample") %>%
  filter(!sample %in% c("Quartet D6"))

no_thres <- data.frame(sample = c("Quartet D5", "Quartet F7", "Quartet M8", "HeLa", "HEK293T"),
                       number_thres = c(1806, 2032, 2052, 1745, 1875))

p_figure3b <- ggplot(df_subSC2) +
  geom_bar(aes(x = group, y = after_stat(count), fill = coverage_label),
           stat = "count", position = "dodge", color = "black") +
  facet_grid(rows = vars(sample)) +
  geom_hline(aes(yintercept = number_thres), data = no_thres,
             linetype = "dashed", color = "#CB181D", alpha = 0.6) +
  geom_text(aes(x = "\u2265 3 Lab(s)", y = number_thres, 
                label = sprintf("N = %s", format(number_thres, big.mark = ",", scientific = FALSE))),
            data = no_thres, color = "#CB181D", size = 5, vjust = -1) +
  # annotate("text", x = "\u2265 3 Lab(s)", y = 2500, label = "N = 2,074",
  #          color = "#CB181D", size = 5) +
  scale_fill_brewer(palette = "Blues") +
  scale_x_discrete(expand = c(0.06, 0.05)) +
  scale_y_continuous(n.breaks = 5, expand = c(0.02, 0.01)) +
  labs(
    x = "Consensus Voting",
    y = "Number of Retained Proteins",
    fill = "Filter") +
  theme_bw() +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16));p_figure3b

ggsave("./results/figures/supp_figure4.pdf", p_figure3b, width = 10, height = 10)


## Tier1: 至少3家实验室所有运行证据支持SC>30% ------------------------
passSC0.3_pep_tables <- pblapply(1:6, function(i) {
  
  sub_pro_table <- grouped_seq_tables[[i]] %>%
    group_by(protein_id, lab_id) %>%
    summarise(pass_n = length(unique(analysis_id[coverage > 30&!is.na(coverage)])),
              total_n = length(unique(analysis_id)), .groups = "drop") %>%
    filter(total_n - pass_n == 0) %>%
    group_by(protein_id) %>%
    summarise(lab_n = length(unique(lab_id)), .groups = "drop") %>%
    filter(lab_n >= 3) %>%
    distinct(protein_id)
  
  sub_pep_table <- grouped_pep_tables[[i]] %>%
    inner_join(., sub_pro_table, by = c("protein_id"))
  
  return(sub_pep_table)
})
names(passSC0.3_pep_tables) <- names(grouped_pep_tables)
stat_pep_tables <- pblapply(passSC0.3_pep_tables, function(tmp_table) {
  stat_tmp_table <- tmp_table %>%
    ungroup() %>%
    summarise(peptide_n = length(unique(peptide_sequence)),
              protein_n = length(unique(protein_id)))
  return(stat_tmp_table)
})
stat_df1 <- rbindlist(stat_pep_tables, idcol = "sample")
stat_df1$tier = "At least 3 labs support sequence coverage > 30%."

View(passSC0.3_pep_tables[[1]] %>% filter(peptide_sequence %in% "FDSDVGEFR"))


## Figure 2c 展示不同实验室投票+PEP阈值下的肽段-蛋白质数目 ----------------------
curvePEP_tables <- pblapply(passSC0.3_pep_tables, function(table_tmp) {
  
  sub_pro_table <- table_tmp %>%
    group_by(peptide_sequence, protein_id, lab_id) %>%
    mutate(pass_n = length(unique(analysis_id[pep < .01 &!is.na(pep)])),
           total_n = length(unique(analysis_id))) %>%
    filter(total_n - pass_n == 0) %>%
    group_by(peptide_sequence, protein_id) %>%
    summarise(lab_n = length(unique(lab_id)),
              pep_min = min(pep, na.rm = TRUE), .groups = "drop")
  
  filtered_tables_i <- mclapply(1:9, function(j) {
    filtered_queries <- sub_pro_table %>%
      filter(lab_n >= j) %>%
      mutate(group = paste("\u2265", j, "Lab(s)"))
    return(filtered_queries)
  })
  df_all_i <- rbindlist(filtered_tables_i)
  
  return(df_all_i)
})
df_subPEP1 <- curvePEP_tables %>%
  rbindlist(., idcol = "sample") %>%
  mutate_at("sample", ~ ifelse(. %in% c("D5", "D6", "F7", "M8"), paste("Quartet", .), .)) %>%
  filter(sample %in% c("Quartet D6")) %>%
  # filter(group == "\u2265 1 Lab(s)") %>%
  mutate(pep_trans = -log10(pep_min)) %>%
  mutate_at("pep_trans", ~ ifelse(is.infinite(.), 400, .))

df_subPEP2 <- df_subPEP1 %>%
  group_by(sample, group) %>%
  summarise(feature_n = length(unique(peptide_sequence)), .groups = "drop")

ymax_pep <- max(df_subPEP1$pep_trans, na.rm = TRUE)
ymax_n   <- max(df_subPEP2$feature_n, na.rm = TRUE)
scale_factor <- ymax_pep / ymax_n

colors.sample <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745", "#E7298A", "#4D9221")
names(colors.sample) <- c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")

p_figure2c <- ggplot(df_subPEP1, aes(x = group, y = pep_trans)) +
  geom_boxplot(aes(fill = group), outliers = FALSE) +
  facet_grid(cols = vars(sample)) +
  geom_line(data = df_subPEP2, aes(x = group, y = feature_n * scale_factor, group = 1),
    inherit.aes = FALSE, linewidth = 1, color = "black") +
  geom_point(data = df_subPEP2, aes(x = group, y = feature_n * scale_factor),
    inherit.aes = FALSE, size = 2.5, color = "black") +
  # geom_hline(yintercept = 23885 * scale_factor, linetype = "dashed",
  #            color = "#CB181D", alpha = 0.6) +
  annotate("text", x = "\u2265 3 Lab(s)", y = 26000 * scale_factor,
           label = "N = 23,885", color = "#CB181D", size = 5) +
  scale_fill_brewer(palette = "Blues") +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  scale_y_continuous(name = "-lg(PEP)", breaks = seq(0, 400, 50),
                     expand = c(0.02, 0.01),
                     sec.axis = sec_axis(~ (.) / scale_factor,
                                         name = "Number of Retained Peptides",
                                         breaks = seq(0, 40000, 5000))) +
  labs(x = "Consensus Voting",
       y = "-lg(PEP)",
       fill = "Filter") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16));p_figure2c

ggsave("./results/figures/figure2c.pdf", p_figure2c, width = 9.5, height = 5.3)


## Supplementary figure 5: D5/F7/M8/HeLa/HEK293T展示投票+PEP阈值 --------------
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

no_thres <- data.frame(sample = c("Quartet D5", "Quartet F7", "Quartet M8", "HeLa", "HEK293T"),
                       number_thres = c(20817, 23981, 23948, 20983, 23962))

colors.sample <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745", "#E7298A", "#4D9221")
names(colors.sample) <- c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")

p_figure3c <- ggplot(df_subPEP2, aes(x = group, y = pep_trans)) +
  geom_boxplot(aes(fill = group), outliers = FALSE) +
  facet_wrap(~ sample, ncol = 1, strip.position = "top") +
  geom_line(data = df_subPEP3, aes(x = group, y = feature_n * scale_factor, group = 1),
            inherit.aes = FALSE, linewidth = 1, color = "black") +
  geom_point(data = df_subPEP3, aes(x = group, y = feature_n * scale_factor),
             inherit.aes = FALSE, size = 2.5, color = "black") +
  geom_hline(aes(yintercept = number_thres* scale_factor), data = no_thres,
             linetype = "dashed", color = "#CB181D", alpha = 0.6) +
  geom_text(aes(x = "\u2265 3 Lab(s)", y = number_thres * scale_factor, 
                label = sprintf("N = %s", format(number_thres, big.mark = ",", scientific = FALSE))),
            data = no_thres, color = "#CB181D", size = 5, vjust = -1) +
  scale_fill_brewer(palette = "Blues") +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  scale_y_continuous(name = "-lg(PEP)", breaks = seq(0, 400, 100),
                     expand = c(0.02, 0.01),
                     sec.axis = sec_axis(~ (.) / scale_factor,
                                         name = "Number of Retained Peptides",
                                         breaks = seq(0, 40000, 10000))) +
  labs(x = "Consensus Voting", y = "-lg(PEP)", fill = "Filter") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 20),
        strip.placement = "outside",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16));p_figure3c

ggsave("./results/figures/supp_figure5.pdf", p_figure3c, width = 9.5, height = 12)


## Tier2: 至少3家实验室所有运行证据支持PEP<1% ---------------------------------
passPEP0.01_peq_tables <- pblapply(passSC0.3_pep_tables, function(tmp_table) {
  
  sub_tmp_table <- tmp_table %>%
    group_by(peptide_sequence, protein_id, lab_id) %>%
    summarise(pass_n = length(unique(analysis_id[pep < .01&!is.na(pep)])),
              total_n = length(unique(analysis_id)), .groups = "drop") %>%
    filter(total_n - pass_n == 0) %>%
    group_by(peptide_sequence, protein_id) %>%
    summarise(lab_n = length(unique(lab_id)), .groups = "drop") %>%
    filter(lab_n >= 3) %>%
    distinct(peptide_sequence, protein_id)
  
  sub_pep_table <- tmp_table %>%
    inner_join(., sub_tmp_table, by = c("peptide_sequence", "protein_id"))
  
  return(sub_pep_table)
})
stat_pep_tables <- pblapply(passPEP0.01_peq_tables, function(tmp_table) {
  stat_tmp_table <- tmp_table %>%
    ungroup() %>%
    summarise(peptide_n = length(unique(peptide_sequence)),
              protein_n = length(unique(protein_id)))
  return(stat_tmp_table)
})
stat_df2 <- rbindlist(stat_pep_tables, idcol = "sample")
stat_df2$tier = "At least 3 labs support PEP < 1%."

View(passPEP0.01_peq_tables[[1]] %>% filter(peptide_sequence %in% "FDSDVGEFR"))

## Figure 2d PEP vs Coverage ----------------------
# rm(list = setdiff(ls(), c("grouped_seq_tables", "grouped_pep_tables")))
# gc()
PEPvsSC_tables <- pblapply(1:6, function(i) {
  
  sub_pro_table <- grouped_seq_tables[[i]] %>%
    mutate(coverage_label = sapply(coverage, function(x) {
      if (x <=10) {
        y <- "0~10%"
      } else if (x <= 20 & x > 10) {
        y <- "10~20%"
      } else if (x <= 30 & x > 20) {
        y <- "20~30%"
      } else if (x <= 40 & x > 30) {
        y <- "30~40%"
      } else if (x <= 50 & x > 40) {
        y <- "40~50%"
      } else if (x <= 60 & x > 50) {
        y <- "50~60%"
      } else if (x <= 70 & x > 60) {
        y <- "60~70%"
      } else if (x <= 80 & x > 70) {
        y <- "70~80%"
      } else if (x > 80) {
        y <- "80~100%"
      }
      return(y)
    })) %>%
    distinct(protein_id, analysis_id, coverage_label)
  
  sub_pep_table <- grouped_pep_tables[[i]] %>%
    filter(!is.na(pep)) %>%
    inner_join(., sub_pro_table, by = c("analysis_id", "protein_id"))
  
  return(sub_pep_table)
})
df_sub1 <- PEPvsSC_tables[[2]] %>%
  mutate(pep_trans = -log10(pep)) %>%
  mutate_at("pep_trans", ~ ifelse(is.infinite(.), 400, .)) %>%
  mutate_at("sample", ~ ifelse(. %in% c("D5", "D6", "F7", "M8"), paste("Quartet", .), .))

colors.sample <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745", "#E7298A", "#4D9221")
names(colors.sample) <- c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")

p_figure2d <- ggplot(df_sub1, aes(x = coverage_label, y = pep_trans)) +
  geom_boxplot(aes(fill = coverage_label), outliers = FALSE) +
  facet_grid(cols = vars(sample)) +
  scale_fill_brewer(palette = "Blues", name = "Range") +
  scale_x_discrete(name = "Sequence Coverage", expand = c(0.1, 0.05)) +
  scale_y_continuous(name = "-lg(PEP)", n.breaks = 10, expand = c(0.02, 0.01)) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16));p_figure2d

ggsave("./results/figures/figure2d.pdf", p_figure2d, width = 7, height = 5)


## Supplementary figure 6: D5/F7/M8/HeLa/HEK293T展示PEP vs Coverage --------------
# rm(list = setdiff(ls(), c("PEPvsSC_tables")))
# gc()
df_sub2 <- PEPvsSC_tables[c(1, 3:6)] %>%
  rbindlist(.) %>%
  mutate_at("sample", ~ ifelse(. %in% c("D5", "D6", "F7", "M8"), paste("Quartet", .), .)) %>%
  mutate(pep_trans = -log10(pep)) %>%
  mutate_at("pep_trans", ~ ifelse(is.infinite(.), 400, .))

colors.sample <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745", "#E7298A", "#4D9221")
names(colors.sample) <- c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")

p_sfigure6 <- ggplot(df_sub2, aes(x = coverage_label, y = pep_trans)) +
  geom_boxplot(aes(fill = coverage_label), outliers = FALSE) +
  facet_grid(rows = vars(sample)) +
  scale_fill_brewer(palette = "Blues", name = "Range") +
  scale_x_discrete(name = "Sequence Coverage", expand = c(0.1, 0.05)) +
  scale_y_continuous(name = "-lg(PEP)", n.breaks = 10, expand = c(0.02, 0.01)) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16))

ggsave("./results/figures/supp_figure6.pdf", p_sfigure6, width = 9.5, height = 12)


## Tier3: 实验室间所有运行证据支持FDR校正后PEP < 1% ----------------------------
passFDR0.01_pep_tables <- pblapply(passPEP0.01_peq_tables, function(tmp_table) {
  sub_tmp_table <- tmp_table %>%
    group_by(peptide_sequence, protein_id) %>%
    mutate(fdr = p.adjust(pep, method = "BH")) %>%
    filter(fdr < .01)
  
  return(sub_tmp_table)
})
stat_pep_tables <- pblapply(passFDR0.01_pep_tables, function(tmp_table) {
  stat_tmp_table <- tmp_table %>%
    ungroup() %>%
    summarise(peptide_n = length(unique(peptide_sequence)),
              protein_n = length(unique(protein_id)))
  return(stat_tmp_table)
})
stat_df3 <- rbindlist(stat_pep_tables, idcol = "sample")
stat_df3$tier = "All evidences support FDR < 1%."

View(passFDR0.01_pep_tables[[1]] %>% filter(peptide_sequence %in% "FDSDVGEFR"))

## 合并统计结果
stat_df_filtered <- rbindlist(list(stat_df0, stat_df1, stat_df2, stat_df3))
stat_df_filtered_wide <- stat_df_filtered %>%
  mutate(tier = factor(tier, levels = unique(tier))) %>%
  mutate(`n (peptides/proteins)` = paste(comma(peptide_n), "\n(", comma(protein_n), ")", sep = "")) %>%
  reshape2::dcast(., tier ~ sample, value.var = "n (peptides/proteins)")

fwrite(stat_df_filtered_wide, "~/Desktop/tmp_定性.csv")
saveRDS(passFDR0.01_pep_tables, "./results/tables/1_qualiprop_list_PEPfdr0.01.rds")


## Figure 2e-f 肽段长度/亲疏水性 Gravy score ----------------------
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
  filter(sample %in% c("Quartet D6", "HeLa")) %>%
  mutate_at("sample", ~ factor(., levels = c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")))

df_sub2 <- df_sub1 %>%
  group_by(sample) %>%
  summarise(gravy_median = median(gravy_score),
            gravy_peak = {
              d <- density(gravy_score, bw = 1, na.rm = TRUE)
              d$x[which.max(d$y)]})

colors.sample <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745", "#E7298A", "#4D9221")
names(colors.sample) <- c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")

set.seed(2026)
p_figure2e <- ggplot(df_sub1, aes(x = gravy_score)) +
  geom_density(aes(fill = sample), alpha = .8) +
  # geom_vline(xintercept = 0, linetype = "dashed", color = "#CB181D", alpha = 0.6) +
  geom_vline(data = df_sub2, mapping = aes(xintercept = gravy_peak),
             linetype = "dashed", color = "black", alpha = 0.6) +
  facet_grid(rows = vars(sample)) +
  scale_fill_manual(values = colors.sample) +
  scale_x_continuous(name = "Gravy Score") +
  scale_y_continuous(name = "Density") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16))

ggsave("./results/figures/figure2e.pdf", p_figure2e, width = 7, height = 3)

df_sub3 <- df_sub1 %>%
  group_by(sample) %>%
  summarise(length_median = median(length),
            length_peak = {
              d <- density(length, bw = 1, na.rm = TRUE)
              d$x[which.max(d$y)]})

set.seed(2026)
p_figure2f <- ggplot(df_sub1, aes(x = length)) +
  geom_density(aes(fill = sample), alpha = .8, bw = 1) +
  # geom_vline(xintercept = 0, linetype = "dashed", color = "#CB181D", alpha = 0.6) +
  geom_vline(data = df_sub3, mapping = aes(xintercept = length_peak),
             linetype = "dashed", color = "black", alpha = 0.6) +
  facet_grid(rows = vars(sample)) +
  scale_fill_manual(values = colors.sample) +
  scale_x_continuous(name = "Sequence Length") +
  scale_y_continuous(name = "Density") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16))

ggsave("./results/figures/figure2f.pdf", p_figure2f, width = 7, height = 3)


## Supplementary figure 7-8: D5/F7/M8/HeLa/HEK293T肽段长度/亲疏水性-----------
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
              d$x[which.max(d$y)]})

colors.sample <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745", "#E7298A", "#4D9221")
names(colors.sample) <- c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")

set.seed(2026)
p_sfigure7 <- ggplot(df_sub1, aes(x = gravy_score)) +
  geom_density(aes(fill = sample), alpha = .8) +
  # geom_vline(xintercept = 0, linetype = "dashed", color = "#CB181D", alpha = 0.6) +
  geom_vline(data = df_sub2, mapping = aes(xintercept = gravy_peak),
             linetype = "dashed", color = "black", alpha = 0.6) +
  facet_grid(rows = vars(sample)) +
  scale_fill_manual(values = colors.sample) +
  scale_x_continuous(name = "Gravy Score") +
  scale_y_continuous(name = "Density") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16))

ggsave("./results/figures/supp_figure7.pdf", p_sfigure7, width = 7, height = 5)

df_sub3 <- df_sub1 %>%
  group_by(sample) %>%
  summarise(length_median = median(length),
            length_peak = {
              d <- density(length, bw = 1, na.rm = TRUE)
              d$x[which.max(d$y)]})

set.seed(2026)
p_sfigure8 <- ggplot(df_sub1, aes(x = length)) +
  geom_density(aes(fill = sample), alpha = .8, bw = 1) +
  # geom_vline(xintercept = 0, linetype = "dashed", color = "#CB181D", alpha = 0.6) +
  geom_vline(data = df_sub3, mapping = aes(xintercept = length_peak),
             linetype = "dashed", color = "black", alpha = 0.6) +
  facet_grid(rows = vars(sample)) +
  scale_fill_manual(values = colors.sample) +
  scale_x_continuous(name = "Sequence Length") +
  scale_y_continuous(name = "Density") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16))

ggsave("./results/figures/supp_figure8.pdf", p_sfigure8, width = 7, height = 5)


## Figure 2g 肽段 m/z ----------------------
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
  filter(sample %in% c("Quartet D6", "HeLa")) %>%
  mutate_at("sample", ~ factor(., levels = c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")))

df_sub2 <- df_sub1 %>%
  group_by(sample) %>%
  summarise(mz_median = median(mz_ratio),
            mz_peak = {
              d <- density(mz_ratio, bw = 1, na.rm = TRUE)
              d$x[which.max(d$y)]})

colors.sample <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745", "#E7298A", "#4D9221")
names(colors.sample) <- c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")

set.seed(2026)
p_figure2g <- ggplot(df_sub1, aes(x = mz_ratio)) +
  geom_density(aes(fill = sample), alpha = .8, bw = 20) +
  # geom_vline(xintercept = 0, linetype = "dashed", color = "#CB181D", alpha = 0.6) +
  geom_vline(data = df_sub2, mapping = aes(xintercept = mz_peak),
             linetype = "dashed", color = "black", alpha = 0.6) +
  facet_grid(rows = vars(sample)) +
  scale_fill_manual(values = colors.sample, name = "RM Group") +
  scale_x_continuous(name = "m/z") +
  scale_y_continuous(name = "Density") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16))

ggsave("./results/figures/figure2g.pdf", p_figure2g, width = 7, height = 3)


## Supplementary figure 9: D5/F7/M8/HeLa/HEK293T肽段 m/z -----------
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
  summarise(mz_median = median(mz_ratio),
            mz_peak = {
              d <- density(mz_ratio, bw = 1, na.rm = TRUE)
              d$x[which.max(d$y)]})

colors.sample <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745", "#E7298A", "#4D9221")
names(colors.sample) <- c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")

set.seed(2026)
p_sfigure9 <- ggplot(df_sub1, aes(x = mz_ratio)) +
  geom_density(aes(fill = sample), alpha = .8, bw = 20) +
  # geom_vline(xintercept = 0, linetype = "dashed", color = "#CB181D", alpha = 0.6) +
  geom_vline(data = df_sub2, mapping = aes(xintercept = mz_peak),
             linetype = "dashed", color = "black", alpha = 0.6) +
  facet_grid(rows = vars(sample)) +
  scale_fill_manual(values = colors.sample, name = "RM Group") +
  scale_x_continuous(name = "m/z") +
  scale_y_continuous(name = "Density") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16), 
        panel.spacing = unit(.5, "cm"))

ggsave("./results/figures/supp_figure9.pdf", p_sfigure9, width = 7, height = 7)


## Figure 2h 肽段-蛋白质关系 ----------------------
rm(list = ls())
gc()
passFDR0.01_pep_tables <- readRDS("./results/tables/1_qualiprop_list_PEPfdr0.01.rds")

match_tables <- pblapply(passFDR0.01_pep_tables, function(table_tmp) {
  
  sub_pep_table <- table_tmp %>%
    group_by(protein_id) %>%
    summarise(peptide_n = length(unique(peptide_sequence)))
  
  return(sub_pep_table)
})
df_sub1 <- match_tables %>%
  rbindlist(., idcol = "sample") %>%
  mutate_at("sample", ~ ifelse(. %in% c("HEK293T", "HeLa"), ., paste("Quartet", .))) %>%
  mutate(protein_label = sapply(peptide_n, function(x) {
    if (x == 1) {
      y <- "1"
    } else if (x == 2) {
      y <- "2"
    } else if (x == 3) {
      y <- "3"
    } else if (x == 4) {
      y <- "4"
    } else if (x >= 5) {
      y <- "5+"
    } 
    return(y)
  })) %>%
  mutate_at("sample", ~ factor(., levels = c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")))

df_sub2 <- df_sub1 %>%
  group_by(sample, protein_label) %>%
  summarise(protein_n = length(unique(protein_id)), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(total_n = sum(protein_n)) %>%
  mutate(percentage = protein_n/total_n)

colors.sample <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745", "#E7298A", "#4D9221")
names(colors.sample) <- c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")

set.seed(2026)
p_figure2h <- ggplot(df_sub2, aes(x = sample, y = percentage * 100)) +
  geom_bar(aes(alpha = protein_label, fill = sample), stat = "identity",
           width = .7, color = "black") +
  scale_fill_manual(values = colors.sample, name = "RM Group") +
  scale_alpha_manual(values = c(.3, .6, .8, .9, 1), name = "# Peptides") +
  guides(fill = "none") +
  scale_y_continuous(name = "Percentage of proteins (%)", breaks = seq(0, 100, 10)) +
  theme_bw() +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16))

ggsave("./results/figures/figure2h.pdf", p_figure2h, width = 10, height = 4.5)



