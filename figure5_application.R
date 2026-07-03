## Title: 以结果为导向的指标
## Author: Qiaochu Chen
## Date: Jan 15th, 2026

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
library(GGally)


## figure 5a: 作图比较 Quartet vs HeLa/HEK293T DEPs 数目------
rm(list = ls())
gc()

dda_limma_tables <- readRDS("./results/tables/dda_limma.rds")
dia_limma_tables <- readRDS("./results/tables/dia_limma.rds")

## 统计DEPs的数目
dda_dep_tables <- pblapply(dda_limma_tables, function(tmp_results) {
  
  tmp_tables <- mclapply(tmp_results, function(tmp_result) {
    
    dep_tmp_j <- tmp_result %>%
      filter(adj.P.Val<.05, abs(logFC) >= log2(2)) %>%
      group_by(group1, group2) %>%
      summarise(dep_n = length(unique(feature)))
    
    return(dep_tmp_j)
  })
  
  dep_tmp_i <- rbindlist(tmp_tables, idcol = "lab")
  
  return(dep_tmp_i)
})
names(dda_dep_tables) <- c("Precursor-level", "Peptide-level", "Protein-level")
dda_dep <- rbindlist(dda_dep_tables, idcol = "data_level")

dia_dep_tables <- pblapply(dia_limma_tables, function(tmp_results) {
  
  tmp_tables <- mclapply(tmp_results, function(tmp_result) {
    
    dep_tmp_j <- tmp_result %>%
      filter(adj.P.Val<.05, abs(logFC) >= log2(2)) %>%
      group_by(group1, group2) %>%
      summarise(dep_n = length(unique(feature)))
    
    return(dep_tmp_j)
  })
  
  dep_tmp_i <- rbindlist(tmp_tables, idcol = "lab")
  
  return(dep_tmp_i)
})
names(dia_dep_tables) <- c("Precursor-level", "Peptide-level", "Protein-level")
dia_dep <- rbindlist(dia_dep_tables, idcol = "data_level")

df_dep <- rbind(dda_dep, dia_dep) %>%
  filter(data_level %in% c("Peptide-level", "Protein-level")) %>%
  filter(!lab %in% "ZJU") %>%
  mutate(source = paste(group2, group1, sep = "/")) %>%
  # mutate_at("source", ~ ifelse(. %in% "HeLa/HEK293T", "HeLa\n/HEK293T", "Quartet")) %>%
  mutate_at("source", ~ factor(., levels = c("D5/D6", "F7/D6", "M8/D6", "HeLa/HEK293T"))) %>%
  mutate_at("data_level", ~ factor(., levels = c("Precursor-level", "Peptide-level", "Protein-level")))

dep_thres <- df_dep %>%
  group_by(source, data_level) %>%
  summarise(dep_median = median(dep_n), snr_max = max(dep_n),
            dep_mean = mean(dep_n), dep_sd = sd(dep_n))
# stat_summary(aes(label = scales::label_comma(accuracy = 1)(after_stat(y))),
#              fun = "median", geom = "text", size = 4, vjust = -1.5) +

colors.sample <- c("#4CC3D9", "#FFC65D", "#F16745", "#E7298A")
names(colors.sample) <- c("D5/D6", "F7/D6", "M8/D6", "HeLa/HEK293T")

p <- ggplot(df_dep, aes(x = source, y = dep_n)) +
  stat_summary(fun.min = min, fun.max = max, geom = "errorbar", width = .5) +
  geom_boxplot(aes(fill = source), outliers = FALSE, width = .75) +
  geom_point(aes(fill = source), shape = 21, color = "black") +
  # geom_hline(aes(yintercept = dep_median, colour = source), dep_thres, lty = 2) +
  geom_text(aes(x = source, y = dep_median, color = source,
                label = sprintf("Median = %s", scales::comma(dep_median, accuracy = 1))),
            vjust = -.2, color = "white",
            data = dep_thres, size = 3) +
  geom_text(aes(x = source, y = Inf,
                label = sprintf("Mean = %s", scales::comma(dep_mean, accuracy = 1))),
            data = dep_thres, size = 3, vjust = 3) +
  # scale_fill_manual(values = c("#8E0152", "#276419")) +
  # scale_color_manual(values = c("#8E0152", "#276419")) +
  scale_fill_manual(values = colors.sample) +
  scale_color_manual(values = colors.sample) +
  scale_y_continuous(n.breaks = 6, name = "Number of DEPs",
                     expand = expansion(mult = c(0.1, 0.2)),
                     labels = scales::label_comma()) +
  # ggh4x::facet_grid2(~data_level, scales = "free_x", space = "free") +
  facet_grid(rows = vars(data_level), scales = "free") +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 16, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm"));p

ggsave("./results/figures/figure5a.pdf", p, height = 8, width = 6.5)


## figure 5b: 作图比较 Quartet vs HeLa/HEK293T PCA ----------------------
rm(list = ls())
gc()

dda_pca_tables <- readRDS("./results/tables/dda_pca_quartet_quali_filter.rds")
dia_pca_tables <- readRDS("./results/tables/dia_pca_quartet_quali_filter.rds")

labels.lab <- c("Lab-1", "Lab-2", "Lab-0", "Lab-3", "Lab-4", "Lab-5", "Lab-6", "Lab-7", "Lab-8", "Lab-9")
names(labels.lab) <- c("QLB", "NCP", "ZJU", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")

colors.sample <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745", "#E7298A", "#4D9221")
names(colors.sample) <- c("D5", "D6", "F7", "M8", "HeLa", "HEK293T")

pca_tables_bydataset <- c(dda_pca_tables[[3]], dia_pca_tables[[3]])
names(pca_tables_bydataset) <- labels.lab[names(pca_tables_bydataset)]

p_legend <- ggplot(pca_tables_bydataset[[1]]$pcs_values, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = sample, alpha = factor(tube)), size = 5) +
  theme_bw() + theme(legend.position = "bottom") +
  scale_alpha_manual(values = c(.9, .6, .3), name = "Tube") +
  scale_color_manual(values = colors.sample, name = "Sample") +
  labs(alpha = NA)

p_pca_list <- pblapply(c(2, 1, 4, 7, 9, 5, 6, 10, 8), function(j) { ##删除ZJU
  
  sub_pca_results <- pca_tables_bydataset[[j]]
  
  p_pca <- ggplot(sub_pca_results$pcs_values, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sample, alpha = factor(tube)), size = 5) +
    scale_alpha_manual(values = c(.9, .6, .3)) +
    scale_color_manual(values = c(colors.sample)) +
    theme_bw() +
    theme(legend.position = "none",
          plot.margin = margin(.1, .1, .1, .1, "cm")) +
    labs(title = names(pca_tables_bydataset)[j],
         subtitle = sprintf("SNR = %.2f\n%s proteins, %d injections",
                            sub_pca_results$snr_results$snr,
                            scales::comma(sub_pca_results$feature_n),
                            sub_pca_results$sample_n),
         # subtitle = sprintf("SNR = %.2f, N = %d (CV < %.0f%%)",
         #                    sub_pca_results$snr_results$snr, sub_pca_results$n, 20),
         x = sprintf("PC1 (%.2f%%)", sub_pca_results$pcs_props[2, 1] * 100),
         y = sprintf("PC2 (%.2f%%)", sub_pca_results$pcs_props[2, 2] * 100))
  
  return(p_pca)
  
})

p_pca_multi <- plot_grid(plotlist = p_pca_list, ncol = 3)
p_pca_all <- plot_grid(p_pca_multi, ggpubr::get_legend(p_legend),
                       nrow = 2, rel_heights = c(1, .1))
ggsave("./results/figures/figure5b.pdf", p_pca_all,
       width = 9.5, height = 12, limitsize = FALSE)



## figure 5c: Recall---------
rm(list = ls())
gc()

all_pep <- fread("./results/tables/8_recall.csv")
all_pep <- all_pep %>% filter(!lab %in% "ZJU")

all_stat <- all_pep %>%
  filter(!lab %in% "ZJU") %>%
  group_by(lab) %>%
  summarise(recall_mean = mean(recall), 
            recall_sd = sd(recall)) %>%
  arrange(desc(recall_mean)) %>%
  mutate_at("lab", ~ factor(., levels = unique(.)))

all_pep$lab <- factor(all_pep$lab, levels = all_stat$lab)

colors.lab <- brewer.pal(10, "PiYG")
names(colors.lab) <- c("ZJU", "FDU", "CAS", "NIM", "QLB", "CMS", "NCP", "OSB", "BTP", "TFS")

labels.lab <- c("Lab-1", "Lab-2", "Lab-0", "Lab-3", "Lab-4", "Lab-5", "Lab-6", "Lab-7", "Lab-8", "Lab-9")
names(labels.lab) <- c("QLB", "NCP", "ZJU", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")

p <- ggplot(all_pep, aes(x = sample, y = recall)) +
  stat_summary(aes(fill = lab), fun = mean, geom = "bar", width = .7,
               color = "black", linewidth = .3) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .4) +
  geom_text(aes(x = "D6", y = 1,
                label = sprintf("%.3f\u00B1%.3f\n(Mean\u00B1SD)", recall_mean, recall_sd)),
            data = all_stat, size = 4.5, vjust = -.5, hjust = .35) +
  # stat_summary(aes(label = scales::percent(accuracy = .1, after_stat(y))),
  #              fun = "mean", geom = "text", size = 2.5, vjust = -1.5) +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_blank(),
        axis.title = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(.3, "cm"),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm")) +
  scale_y_continuous(breaks = seq(0, 1, .2), name = "Recall",
                     expand = expansion(mult = c(0, 0.2))) +
  # scale_fill_manual(values = colors.sample) +
  scale_fill_manual(values = colors.lab) +
  ggh4x::facet_grid2( ~ lab, scales = "free_y", labeller = as_labeller(labels.lab))

ggsave("./results/figures/figure5c.pdf", p, height = 5, width = 14)


## figure 5d: F1---------
rm(list = ls())
gc()

all_pep <- fread("./results/tables/8_f1.csv")
all_pep <- all_pep %>% filter(!lab %in% "ZJU")

all_stat <- all_pep %>%
  filter(!lab %in% "ZJU") %>%
  group_by(lab) %>%
  summarise(f1_mean = mean(f1, na.rm = TRUE), 
            f1_sd = sd(f1, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(f1_mean)) %>%
  mutate_at("lab", ~ factor(., levels = unique(.)))

all_pep <- na.omit(all_pep)
all_pep$lab <- factor(all_pep$lab, levels = all_stat$lab)

colors.lab <- brewer.pal(10, "PiYG")
names(colors.lab) <- c("CMS", "NIM", "NCP", "ZJU", "TFS", "QLB", "OSB", "BTP", "CAS", "FDU")

labels.lab <- c("Lab-1", "Lab-2", "Lab-0", "Lab-3", "Lab-4", "Lab-5", "Lab-6", "Lab-7", "Lab-8", "Lab-9")
names(labels.lab) <- c("QLB", "NCP", "ZJU", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")

p <- ggplot(all_pep, aes(x = sample, y = f1)) +
  stat_summary(aes(fill = lab), fun = mean, geom = "bar", width = .7,
               color = "black", linewidth = .3) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .4) +
  geom_text(aes(x = "F7", y = 1,
                label = sprintf("%.3f\u00B1%.3f\n(Mean\u00B1SD)", f1_mean, f1_sd)),
            data = all_stat, size = 4.5, vjust = -.5) +
  # stat_summary(aes(label = scales::percent(accuracy = .1, after_stat(y))),
  #              fun = "mean", geom = "text", size = 2.5, vjust = -1.5) +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_blank(),
        axis.title = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(.3, "cm"),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm")) +
  scale_y_continuous(breaks = seq(0, 1, .2), name = "F1 Score",
                     expand = expansion(mult = c(0, 0.2))) +
  scale_fill_manual(values = colors.lab) +
  ggh4x::facet_grid2( ~ lab, scales = "free_y", labeller = as_labeller(labels.lab))

ggsave("./results/figures/figure5d.pdf", p, height = 5, width = 14)


## Supplementary figure 12: 所有指标相关性---------
rm(list = ls())
gc()

all_metrics <- fread("./results/tables/all_metrics.csv")

## Long: 每个指标取各实验室中位值 & Scale to 1~10
all_metrics_scaled2 <- all_metrics %>% filter(!lab %in% "ZJU") %>%
  mutate_at("Delta IM", ~ ifelse(. == 0, NA, .)) %>%
  group_by(source, data_level, lab) %>%
  summarise(across(c(snr:`MS2 Accuracy`, -PC1, -PC2),
                   ~ median(abs(.), na.rm = TRUE)), .groups = "drop") %>%
  group_by(source, data_level) %>%
  mutate_at("Delta IM",
            ~ 10 - (. - min(., na.rm = TRUE)) * 9 / (max(., na.rm = TRUE) - min(., na.rm = TRUE))) %>%
  mutate_at(c("FWHM", "Calibarted Mass Error", "Delta RT", "TIC Deviation",
              "MS1 Accuracy", "MS2 Accuracy"),
            ~ 10 - (. - min(.)) * 9 / (max(.) - min(.))) %>%
  mutate_at(c("Charge_2_3_prop", "Qualification"), ~ 1 + (. - min(.)) * 9 / (max(.) - min(.))) %>%
  mutate_at("Quantification", ~ ifelse(source %in% "All", ## 代表是cv而不是f1
                                       10 - (. - min(.)) * 9 / (max(.) - min(.)),
                                       1 + (. - min(.)) * 9 / (max(.) - min(.)))) %>%
  mutate_at("snr", ~ 1 + (. - min(.)) * 9 / (max(.) - min(.))) %>%
  ungroup %>%
  select(source, data_level, lab, everything()) %>%
  mutate_at("data_level", ~ factor(., levels = c("Precursor-level", "Peptide-level", "Protein-level")))

## Long to wide
all_wide <- all_metrics_scaled2 %>%
  as.data.table(.) %>%
  # data.table::dcast(., lab + sample + tube + injection ~ source + data_level,
  data.table::dcast(., lab ~ source + data_level,
                    value.var = c("snr", "Qualification", "Quantification",
                                  "FWHM", "Calibarted Mass Error", "Delta RT",
                                  "Delta IM","TIC Deviation", "Charge_2_3_prop",
                                  "MS1 Accuracy", "MS2 Accuracy"), 
                    fun.aggregate = median, na.rm = TRUE)

all_wide[, which(sapply(all_wide, function(x) all(is.na(x)))) := NULL]

## 选择特定列
df_tmp <- all_wide %>%
  select(any_of(c("lab", 
                  "Qualification_All_Precursor-level",
                  "Quantification_All_Precursor-level",
                  "snr_High-Confidence_Protein-level",
                  "Qualification_High-Confidence_Precursor-level",
                  "Quantification_High-Confidence_Precursor-level",
                  "FWHM_High-Confidence_Precursor-level",
                  "Calibarted Mass Error_High-Confidence_Precursor-level",
                  "Delta RT_High-Confidence_Precursor-level",
                  "Delta IM_High-Confidence_Precursor-level",
                  "TIC Deviation_High-Confidence_Precursor-level",
                  "Charge_2_3_prop_High-Confidence_Precursor-level",
                  "MS1 Accuracy_All_Precursor-level",
                  "MS2 Accuracy_All_Precursor-level")))

colnames(df_tmp)[2:14] <- c("Count", "CV", "SNR", "Recall", "F1",
                            "FWHM", "Mass Error", "Delta RT", "Delta IM", "TIC",
                            "Charge", "MS1 Acc", "MS2 Acc")

# 定义支持正负染色的自定义函数
my_custom_cor <- function(data, mapping, ...) {
  
  # 1. 提取数据
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # 2. 处理缺失值
  valid_idx <- complete.cases(x, y)
  x <- x[valid_idx]
  y <- y[valid_idx]
  
  # 数据太少不计算
  if (length(x) < 3) {
    return(ggally_text("NA"))
  }
  
  # 3. 计算相关性和 P 值
  ct <- cor.test(x, y)
  r <- ct$estimate
  p <- ct$p.value
  
  # 生成星号
  stars <- ""
  if(p < 0.001) stars <- "***"
  else if(p < 0.01) stars <- "**"
  else if(p < 0.05) stars <- "*"
  
  # 4. 【核心修改】颜色逻辑
  # 默认颜色为黑色 (不显著)
  text_color <- "black"
  
  # 只有当显著时 (p < 0.05)，才根据正负变色
  if (p < 0.05) {
    if (r > 0) {
      text_color <- "#D73027"  # 正相关：红色
    } else {
      text_color <- "#4575B4"  # 负相关：蓝色
    }
  } else {
    if (r > .8) {
      text_color <- "#FDAE6B"  # 正相关：橘色
    } else if (r < -.8) {
      text_color <- "#9E9AC8"  # 负相关：紫色
    }
  }
  
  # 5. 生成标签
  lbl <- sprintf("%.3f%s", r, stars)
  
  # 6. 绘图
  ggally_text(
    label = lbl, 
    mapping = aes(),
    xP = 0.5, yP = 0.5, 
    color = text_color, 
    size = 8,
    ...
  ) + 
    theme_void() + 
    theme(panel.border = element_rect(color = "grey90", fill = NA))
}

## 绘制散点图矩阵
colors.lab <- brewer.pal(9, "PiYG")
names(colors.lab) <- c("CAS", "OSB", "NCP", "NIM", "CMS", "QLB", "FDU", "BTP", "TFS")

p <- ggpairs(df_tmp, mapping = aes(fill = lab),
             columns = 2:14, 
             # lower = list(continuous = my_scatter_matrix),
             lower = list(
               continuous = wrap("points", 
                                 size = 2,
                                 color = "transparent",
                                 shape = 21,
                                 size = 4)),
             diag = list(continuous = wrap("barDiag")),
             upper = list(continuous = my_custom_cor) 
) +
  scale_fill_manual(values = colors.lab, name = "Lab") +
  # scale_fill_brewer(palette = "PiYG") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 18, margin = unit(rep(.3, 4), "cm")),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .5, .5, .5), "cm"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 14))

ggsave("./results/figures/supp_figure12.pdf", p, height = 20, width = 20)


