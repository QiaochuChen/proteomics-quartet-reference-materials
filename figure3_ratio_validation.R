## Title: 根据HeLa标准品公共数据库数据验证DIA数据有效性。
## Author: Qiaochu Chen
## Date: Apr 30th, 2026

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


## 适配学位论文——PCA马氏距离法检查离群值+All Ratio by D6/HEK293T -----------------
grouped_tables <- readRDS("../data/multilab/quantdata_list_pep_2025_10labs.rds")
all_tables <- grouped_tables %>% rbindlist %>% split(., by = "lab_id")
rm(list = setdiff(ls(), c("all_tables")))
gc()

## 每家实验室PCA检查离群值: 马氏距离法
source("~/Desktop/utils/PCA.R")
madist_results_tables <- pblapply(all_tables, function(tmp_table) {
  
  metadata <- tmp_table %>%
    distinct(analysis_id, lab_id, sample, tube, injection, order) %>%
    tibble::column_to_rownames("analysis_id")
  
  exprdata_t <- tmp_table %>%
    reshape2::dcast(peptide_sequence ~ analysis_id, value.var = "value", fun.aggregate = sum) %>%
    select(all_of(rownames(metadata))) %>%
    mutate_all(log2) %>%
    filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
    filter(apply(., 1, function(x) sum(is.na(x)) < 0.2 * nrow(metadata))) %>%
    mutate(row_mean = apply(., 1, mean, na.rm = TRUE)) %>%
    mutate_all(~ ifelse(is.na(.), row_mean, .)) %>%
    select(!row_mean) %>%
    # na.omit %>%
    t
  
  pca_results <- calculate_pca(exprdata_t = exprdata_t, metadata = metadata,
                               center = TRUE, scale = TRUE, group = "sample")
  
  pcs_n <- which(pca_results$pcs_props[3, ] > .6)[1]
  
  sub_pca_results <- pca_results$pcs_values %>%
    select(all_of(paste("PC", 1:pcs_n, sep = "")))
  
  pcs_colmeans <- colMeans(sub_pca_results)
  pcs_cov <- cov(sub_pca_results)
  
  madist_results <- pca_results$pcs_values %>%
    select(all_of(paste("PC", 1:pcs_n, sep = ""))) %>%
    mutate(mahalanobis_dist = apply(., 1, function(x) mahalanobis(as.numeric(x), center = pcs_colmeans, pcs_cov))) %>%
    mutate(threshold_0.95 = qchisq(p = 0.95, df = pcs_n)) %>%
    mutate(order = pca_results$pcs_values$order, df = pcs_n) %>%
    mutate(is.outlier = ifelse(mahalanobis_dist > threshold_0.95, TRUE, FALSE)) %>%
    select(!starts_with("PC"))
  
  madist_results_final <- metadata %>%
    tibble::rownames_to_column("analysis_id") %>%
    left_join(., madist_results, by = c("order"))
  
  return(madist_results_final)
})
sub_meta2outlier <- rbindlist(madist_results_tables)
fwrite(sub_meta2outlier, "./results/tables/2_outlier_madist_all.csv")

## 定义ratio函数
ratio_by_ref <- function(df_quant_wide, df_meta, ref_sample = "D6", log2_transformed = FALSE) {
  
  colnames(df_quant_wide)[1] <- "feature"
  colnames(df_meta)[1] <- "library"
  
  df_corrected <- df_quant_wide %>%
    tibble::column_to_rownames("feature") %>%
    select(all_of(df_meta$library))
  
  batches <- unique(df_meta$batch)
  samples <- colnames(df_corrected)
  
  for (batch in batches) {
    match_batch <- which(df_meta$batch %in% batch)
    match_d6 <- match_batch[which(grepl(tolower(ref_sample), tolower(samples[match_batch])))]
    
    if (log2_transformed) {
      
      df_corrected <- df_corrected %>%
        as.data.frame() %>%
        mutate_at(match_batch, ~ ifelse(. == 0, NA, log2(.))) %>%
        mutate(D6mean = apply(.[match_d6], 1, mean, na.rm = TRUE)) %>%
        mutate_at(match_batch, ~ 2 ^ (. - D6mean)) %>%
        select(!D6mean)
    } else {
      
      df_corrected <- df_corrected %>%
        as.data.frame() %>%
        mutate(D6mean = apply(.[match_d6], 1, mean, na.rm = TRUE)) %>%
        mutate_at(match_batch, ~ (.) - D6mean) %>% ## 无需fot和log2，因为输入的数据就是MS log2fot、olink NPX、somaLogic
        select(!D6mean)
    }
  }
  
  df_corrected_final <- df_corrected %>%
    tibble::rownames_to_column("feature")
  
  return(df_corrected_final)
}

## 每家实验室ratio by D6
ratiobyd6_tables <- pblapply(all_tables, function(tmp_table) {
  
  meta_tmp <- tmp_table %>%
    distinct(analysis_id) %>%
    inner_join(., sub_meta2outlier, by = c("analysis_id")) %>%
    filter(is.outlier == FALSE) %>%
    # mutate(batch = ifelse(order <= 12, "phoenix_1", "phoenix_2")) %>%
    # mutate(batch = ifelse(lab_id %in% "phoenix", batch, paste(lab_id, tube, sep = "_"))) %>%
    mutate(batch = paste(lab_id, tube, sep = "_"))
  
  expr_tmp <- reshape2::dcast(tmp_table, peptide_sequence ~ analysis_id, value.var = "value", fun.aggregate = sum)
  
  expr_tmp <- expr_tmp %>% mutate_if(is.numeric, ~ ifelse(. == 0 | is.na(.) | is.nan(.), NA, .))
  
  if (unique(meta_tmp$lab_id) %in% c("phoenix", "zhejiang_university", "qinglian_bio")) {
    expr_tmp <- expr_tmp %>% mutate_if(is.numeric, ~  (.) * 10 ^ 6 / sum(., na.rm = TRUE))
  }
  
  expr_tmp <- expr_tmp %>% mutate_if(is.numeric, log2)
  
  expr_tmp_ratiobyd6 <- ratio_by_ref(expr_tmp, meta_tmp, ref_sample = "D6", log2_transformed = FALSE)
  
  tmp_table_ratiobyd6 <- expr_tmp_ratiobyd6 %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id", na.rm = TRUE) %>%
    dplyr::rename(peptide_sequence = feature) %>%
    left_join(., meta_tmp, by = "analysis_id")
  
  return(tmp_table_ratiobyd6)
})

## 合并所有数据
ratiobyd6_tables <- ratiobyd6_tables %>% rbindlist %>% split(., by = "sample")
saveRDS(ratiobyd6_tables, "./results/tables/2_quantdata_list_pep_allratiobyd6_2025_10labs.rds")

## 每家实验室ratio by HEK293T
ratiobyhek293t_tables <- pblapply(all_tables, function(tmp_table) {
  
  meta_tmp <- tmp_table %>%
    distinct(analysis_id) %>%
    inner_join(., sub_meta2outlier, by = c("analysis_id")) %>%
    filter(is.outlier == FALSE) %>%
    # mutate(batch = ifelse(order <= 12, "phoenix_1", "phoenix_2")) %>%
    # mutate(batch = ifelse(lab_id %in% "phoenix", batch, paste(lab_id, tube, sep = "_"))) %>%
    mutate(batch = paste(lab_id, tube, sep = "_"))
  
  expr_tmp <- reshape2::dcast(tmp_table, peptide_sequence ~ analysis_id, value.var = "value", fun.aggregate = sum)
  
  expr_tmp <- expr_tmp %>% mutate_if(is.numeric, ~ ifelse(. == 0 | is.na(.) | is.nan(.), NA, .))
  
  if (unique(meta_tmp$lab_id) %in% c("phoenix", "zhejiang_university", "qinglian_bio")) {
    expr_tmp <- expr_tmp %>% mutate_if(is.numeric, ~  (.) * 10 ^ 6 / sum(., na.rm = TRUE))
  }
  
  expr_tmp <- expr_tmp %>% mutate_if(is.numeric, log2)
  
  expr_tmp_ratiobyd6 <- ratio_by_ref(expr_tmp, meta_tmp, ref_sample = "HEK293T", log2_transformed = FALSE)
  
  tmp_table_ratiobyd6 <- expr_tmp_ratiobyd6 %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id", na.rm = TRUE) %>%
    dplyr::rename(peptide_sequence = feature) %>%
    left_join(., meta_tmp, by = "analysis_id")
  
  return(tmp_table_ratiobyd6)
})

## 合并所有数据
ratiobyhek293t_tables <- ratiobyhek293t_tables %>% rbindlist %>% split(., by = "sample")
saveRDS(ratiobyhek293t_tables, "./results/tables/2_quantdata_list_pep_allratiobyhek293t_2025_10labs.rds")


## 1.1: 计算Quartet合并实验室D6比例定量前后的PCA/SNR (删除ZJU) ----------------------
rm(list = ls())
gc()
source("./PCA.R")
meta_quartet <- fread("./results/tables/2_outlier_madist_all.csv")
all_tables <- readRDS("./data/multilab/quantdata_list_pep_2025_10labs.rds")
# all_tables <- readRDS("./results/tables/2_quantdata_list_pep_allratiobyd6_2025_10labs.rds")
filtered_tables_groupedbylabs <- all_tables %>% rbindlist %>% split(., by = "lab_id")
filtered_tables_groupedbylabs <- filtered_tables_groupedbylabs[c(1, 2, 4:10)]

## PCA
metadata <- meta_quartet %>%
  filter(is.outlier == FALSE) %>%
  filter(grepl("Quartet", sample)) %>%
  tibble::column_to_rownames("analysis_id")

exprdata_t <- filtered_tables_groupedbylabs %>%
  rbindlist %>%
  # reshape2::dcast(peptide_sequence ~ analysis_id, value.var = "value") %>% ## 若输入的是原始强度的数据则此步隐去
  reshape2::dcast(peptide_sequence ~ analysis_id, value.var = "value", fun.aggregate = sum) %>% ## 若输入的是Ratio by D6后的数据则此步隐去
  select(any_of(rownames(metadata))) %>%
  mutate_all(~ ifelse(. == 0, NA, .)) %>% ## 若输入的是Ratio by D6后的数据则此步隐去
  mutate_all(log2) %>% ## 若输入的是Ratio by D6后的数据则此步隐去
  filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
  filter(apply(., 1, function(x) sum(is.na(x)) < 0.2 * nrow(metadata))) %>%
  mutate(row_mean = apply(., 1, mean, na.rm = TRUE)) %>%
  mutate_all(~ ifelse(is.na(.), row_mean, .)) %>%
  select(!row_mean) %>%
  # na.omit %>%
  t

pca_results <- main_pca(exprdata_t = exprdata_t, metadata = metadata,
                        center = TRUE, scale = TRUE, group = "sample",
                        biplot = FALSE, dictGroups = metadata$sample,
                        snr = TRUE, plot = FALSE)

pca_results$feature_n <- ncol(exprdata_t)
pca_results$sample_n <- nrow(exprdata_t)

saveRDS(pca_results, "./results/tables/all_pca_quartet_combined.rds")
# saveRDS(pca_results, "./results/tables/all_pca_quartet_combined_bytube_ratiobyd6.rds")


## 1.2: 计算HeLa/HEK293T合并实验室D6比例定量前后的PCA/SNR ----------------------
rm(list = ls())
gc()
source("./PCA.R")
meta_hela <- fread("./results/tables/2_outlier_madist_all.csv")
all_tables <- readRDS("./data/multilab/quantdata_list_pep_2025_10labs.rds")
# all_tables <- readRDS("./results/tables/2_quantdata_list_pep_allratiobyd6_2025_10labs.rds")
filtered_tables_groupedbylabs <- all_tables %>% rbindlist %>% split(., by = "lab_id")
filtered_tables_groupedbylabs <- filtered_tables_groupedbylabs[c(1, 2, 4:10)]

## PCA
metadata <- meta_hela %>%
  filter(is.outlier == FALSE) %>%
  filter(sample %in% c("HeLa", "HEK293T")) %>%
  tibble::column_to_rownames("analysis_id")

exprdata_t <- filtered_tables_groupedbylabs %>%
  rbindlist %>%
  # reshape2::dcast(peptide_sequence ~ analysis_id, value.var = "value") %>% ## 若输入的是原始强度的数据则此步隐去
  reshape2::dcast(peptide_sequence ~ analysis_id, value.var = "value", fun.aggregate = sum) %>% ## 若输入的是Ratio by D6后的数据则此步隐去
  select(any_of(rownames(metadata))) %>%
  mutate_all(~ ifelse(. == 0, NA, .)) %>% ## 若输入的是Ratio by D6后的数据则此步隐去
  mutate_all(log2) %>% ## 若输入的是Ratio by D6后的数据则此步隐去
  filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
  filter(apply(., 1, function(x) sum(is.na(x)) < 0.2 * nrow(metadata))) %>%
  mutate(row_mean = apply(., 1, mean, na.rm = TRUE)) %>%
  mutate_all(~ ifelse(is.na(.), row_mean, .)) %>%
  select(!row_mean) %>%
  # na.omit %>%
  t

pca_results <- main_pca(exprdata_t = exprdata_t, metadata = metadata,
                        center = TRUE, scale = TRUE, group = "sample",
                        biplot = FALSE, dictGroups = metadata$sample,
                        snr = TRUE, plot = FALSE)

pca_results$feature_n <- ncol(exprdata_t)
pca_results$sample_n <- nrow(exprdata_t)

saveRDS(pca_results, "./results/tables/all_pca_hela_hek293t_combined.rds")
# saveRDS(pca_results, "./results/tables/all_pca_hela_hek293t_combined_bytube_ratiobyd6.rds")


## 1.3: 计算All合并实验室D6比例定量前后的PCA/SNR ----------------------
rm(list = ls())
gc()
source("./PCA.R")
meta_all <- fread("./chapter2/results/tables/2_outlier_madist_all.csv")
all_tables <- readRDS("./data/multilab/quantdata_list_pep_2025_10labs.rds")
# all_tables <- readRDS("./chapter2/results/tables/2_quantdata_list_pep_allratiobyd6_2025_10labs.rds")
filtered_tables_groupedbylabs <- all_tables %>% rbindlist %>% split(., by = "lab_id")
filtered_tables_groupedbylabs <- filtered_tables_groupedbylabs[c(1, 2, 4:10)]

## PCA
metadata <- meta_all %>%
  filter(is.outlier == FALSE) %>%
  tibble::column_to_rownames("analysis_id")

exprdata_t <- filtered_tables_groupedbylabs %>%
  rbindlist %>%
  # reshape2::dcast(peptide_sequence ~ analysis_id, value.var = "value") %>% ## 若输入的是原始强度的数据则此步隐去
  reshape2::dcast(peptide_sequence ~ analysis_id, value.var = "value", fun.aggregate = sum) %>% ## 若输入的是Ratio by D6后的数据则此步隐去
  select(any_of(rownames(metadata))) %>%
  mutate_all(~ ifelse(. == 0, NA, .)) %>% ## 若输入的是Ratio by D6后的数据则此步隐去
  mutate_all(log2) %>% ## 若输入的是Ratio by D6后的数据则此步隐去
  filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
  filter(apply(., 1, function(x) sum(is.na(x)) < 0.2 * nrow(metadata))) %>%
  mutate(row_mean = apply(., 1, mean, na.rm = TRUE)) %>%
  mutate_all(~ ifelse(is.na(.), row_mean, .)) %>%
  select(!row_mean) %>%
  # na.omit %>%
  t

pca_results <- main_pca(exprdata_t = exprdata_t, metadata = metadata,
                        center = TRUE, scale = TRUE, group = "sample",
                        biplot = FALSE, dictGroups = metadata$sample,
                        snr = TRUE, plot = FALSE)

pca_results$feature_n <- ncol(exprdata_t)
pca_results$sample_n <- nrow(exprdata_t)

saveRDS(pca_results, "./results/tables/all_pca_all_combined.rds")
# saveRDS(pca_results, "./results/tables/all_pca_all_combined_bytube_allratiobyd6.rds")


## Figure 3b: 作图检查合并实验室比例定量前后的PCA/SNR ----------------------
rm(list = ls())
gc()
colors.sample <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745", "#E7298A", "#4D9221")
names(colors.sample) <- c("Quartet D5", "Quartet D6", "Quartet F7", "Quartet M8", "HeLa", "HEK293T")
# labels.lab <- c("QLB", "NCP", "ZJU", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")
labels.lab <- c("Lab-1", "Lab-2", "Lab-0", "Lab-3", "Lab-4", "Lab-5", "Lab-6", "Lab-7", "Lab-8", "Lab-9")
names(labels.lab) <- c("qinglian_bio", "phoenix", "zhejiang_university", "fudan_university",
                       "national_institute_of_methodology", "biotech_pack", "thermofisher_shanghai",
                       "omicsolution", "cas_tianjin", "academy_of_chinese_medical_sciences")

all_files <- list.files("./results/tables", full.names = TRUE)
pca_files <- all_files[grepl("2_pca_quartet_combined", all_files)]
p_titles <- c("Ratio-by-D6", "Raw")
pca_tables <- pblapply(1:2, function(i) {
  
  pca_results <- readRDS(pca_files[i])
  
  df_pca <- pca_results$pcs_values %>%
    select(1:12) %>%
    mutate(snr = pca_results$snr_results$snr,
           feature_n = pca_results$feature_n,
           sample_n = pca_results$sample_n,
           x_title = sprintf("PC1 (%.2f%%)", pca_results$pcs_props[2, 1] * 100),
           y_title = sprintf("PC2 (%.2f%%)", pca_results$pcs_props[2, 2] * 100)) %>%
    mutate(p_title = sprintf("%s\nSNR = %.2f\nF = %s, S = %d",
                             p_titles[i], snr, feature_n, sample_n))
  
  return(df_pca)
})
df_pca_final <- pca_tables %>%
  rbindlist %>%
  filter(grepl("S = 320", p_title)) %>%
  mutate_at("sample", ~ paste("Quartet", .)) %>%
  mutate_at("p_title", ~ factor(., levels = c("Raw\nSNR = -0.05\nF = 14616, S = 320",
                                              "Ratio-by-D6\nSNR = 8.99\nF = 13254, S = 320")))

p_pca <- ggplot(df_pca_final, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = sample, shape = lab_id), size = 5) +
  scale_shape_manual(values = 1:10, labels = labels.lab) +
  scale_color_manual(values = c(colors.sample)) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        strip.background = element_blank(),
        # panel.spacing = unit(.5, "cm"),
        plot.margin = margin(.5, .5, .5, .5, "cm"),
        legend.position = "right") +
  facet_wrap( ~ p_title, ncol = 1, scales = "free") +
  labs(shape = "Lab", color = "Sample")

ggsave("./results/figures/figure3b.pdf", p_pca, width = 6, height = 10)


## Figure 3g: 整理2020年6批数据集中的比例定量值 -----------
all_tables <- readRDS("./data/stability/quantdata_list_pep_all_longterm.rds")
filtered_tables <- pblapply(all_tables, function(tmp_table) test_table <- tmp_table %>% filter(time_day <= 22))
test_tables <- filtered_tables %>% rbindlist %>% split(., by = c("lab_id"))

all_meta <- fread("./results/tables/6_outlier_madist_quartet.csv")
test_meta <- all_meta %>% filter(time_day <= 22)

## Limma: multi-lab (文章：https://doi.org/10.1038/s41467-024-47899-w)
source("./DEP.R")
limma_tables <- pblapply(test_tables, function(tmp_table) {
  
  all_sample_pairs <- list(c("Quartet D6", "Quartet D5"), c("Quartet D6", "Quartet F7"),
                           c("Quartet D6", "Quartet M8"))
  
  limma_tables_i <- mclapply(all_sample_pairs, function(sample_pair_id) {
    
    print(sample_pair_id)
    sub_wide_j <- tmp_table %>%
      filter(sample %in% sample_pair_id)
    
    if (nrow(sub_wide_j) == 0) return(NULL)
    
    sub_wide_j <- tmp_table %>%
      filter(sample %in% sample_pair_id) %>%
      reshape2::dcast(peptide_sequence ~ sample + tube + injection, value.var = "value", fun.aggregate = sum) %>%
      mutate_if(is.numeric, ~ ifelse(is.na(.), 0, .)) %>%
      tibble::column_to_rownames(c("peptide_sequence")) %>%
      select_if(is.numeric)
    
    group_j <- tmp_table %>%
      filter(sample %in% sample_pair_id) %>%
      distinct(analysis_id, sample) %>%
      pull(sample) %>%
      factor(., levels = sample_pair_id, ordered = TRUE)
    
    if (length(which(table(group_j) == 1)) > 0) return(NULL)
    
    limma_results_j <- calculate_p(exprdata = sub_wide_j, group = group_j,
                                   test_method = "limma", p_adjust = "BH",
                                   na_threshold = 0)
    
    return(limma_results_j)
  }, mc.cores = 4)
  
  limma_results_i <- limma_tables_i %>% rbindlist
  
  return(limma_results_i)
})

df_limma <- limma_tables %>%
  rbindlist(., idcol = "lab_id") %>%
  dplyr::rename(peptide_sequence = feature) %>%
  mutate(group = paste(group2, group1, sep = "/")) %>%
  mutate_at("group", ~ gsub("Quartet ", "", .)) %>%
  filter(adj.P.Val < 0.05, abs(logFC) >= log2(1.2)) %>%
  mutate(FC = 2 ^ logFC) 

ref_dt_final <- fread("./results/tables/7_quantprop_final.csv")
df_test <- df_limma %>%
  group_by(peptide_sequence, group) %>%
  summarise(FC_old = mean(FC), sd = sd(FC)) %>%
  mutate(cv = sd/FC_old) %>%
  filter(cv < .2) %>%
  inner_join(., ref_dt_final, by = c("peptide_sequence", "group"),
             relationship = "many-to-many") %>%
  dplyr::rename(FC_new = value) %>%
  select(peptide_sequence, protein_id, group, FC_new, FC_old, U)


## Figure 3g：散点图确认相关性 -------------------
colors.group <- c("#4CC3D9", "#FFC65D", "#F16745", "#999")
names(colors.group) <- c("D5/D6", "F7/D6", "M8/D6", "Reversed")

sub_test <- df_test %>%
  mutate_at(c("FC_new", "FC_old"), log2) %>%
  mutate(group_new = ifelse(FC_new * FC_old < 0, "Reversed", group))

p_cor_final <- ggplot(sub_test, aes(x = FC_old, y = FC_new)) +
  geom_point(aes(fill = group_new), size = 5, shape = 21, color = "black") +
  geom_smooth(method = "lm", se = T, linetype = 1, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "black", lty = 2) +
  ggpubr::stat_cor(size = 6) +
  # ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~`,`~")),
  #                               label.y = c(1.7, 1.7), size = 6) +
  scale_fill_manual(values = colors.group) +
  facet_grid( ~ group, scales = "fixed") +
  scale_y_continuous(limits = c(-1.8, 1.8), name = "DDA/DIA Ratios (Latest)") +
  scale_x_continuous(limits = c(-1.8, 1.8), name = "DDA/DIA Ratios (History)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        strip.background = element_blank(),
        aspect.ratio = 1,
        panel.spacing = unit(1, "cm"),
        plot.margin = unit(c(.5,.5, .5, .5), "cm"))

ggsave("./results/figures/figure3g.pdf", p_cor_final, height = 5.2, width = 14, limitsize = FALSE)



## Figure 3f: Quartet RNA 差异集 -------------------
db_uniprot <- readAAStringSet("./uniprotkb_proteome_UP000005640_2024_09_03.fasta")
db_uniprot <- data.frame(protein_sequence = as.character(db_uniprot)) %>%
  tibble::rownames_to_column("entry") %>%
  mutate(protein_id = str_extract(entry, "(?<=\\|).+(?=\\|)")) %>%
  mutate(protein_fullname = str_extract(entry, "(?<=_HUMAN ).+(?= OS\\=Homo sapiens )")) %>%
  mutate(gene_symbol = str_extract(entry, ifelse(grepl("PE=", entry),
                                                 "(?<=GN\\=).+(?= PE=)",
                                                 "(?<=GN\\=).+"))) %>%
  distinct(protein_id, protein_fullname, gene_symbol)

df_test_rna <- read_excel("./data/public/Quartet_RNA/41587_2023_1867_MOESM3_ESM.xlsx", sheet = 7, skip = 1)
# df_test_pro <- fread("./results/tables/7_quantprop_final.csv")
# df_test_combined <- df_test_rna %>%
#   # filter(!`DEG type` %in% "non-DEG") %>%
#   inner_join(., df_test_pro, by = c("Gene symbol" = "gene_symbol", "Sample pair" = "group")) %>%
#   dplyr::rename(FC_rna = FC, FC_pro = value) %>%
#   group_by(`Gene symbol`, `Sample pair`, FC_rna) %>%
#   summarise_at("FC_pro", mean)

limma_tables0 <- readRDS("./results/tables/2_limma_bylab_bytube.rds")
df_limma0 <- limma_tables0 %>%
  rbindlist(.) %>%
  mutate(group = paste(group2, group1, sep = "/")) %>%
  mutate_at("group", ~ gsub("Quartet ", "", .)) %>%
  filter(adj.P.Val < 0.05, abs(logFC) >= 1) %>%
  # group_by(peptide_sequence, protein_id, group) %>%
  mutate(FC_pep = 2 ^ (logFC)) %>%
  left_join(., db_uniprot, by = "protein_id") %>%
  filter(!(is.na(gene_symbol)|(gene_symbol %in% "")))

df_test_combined <- df_test_rna %>%
    filter(!`DEG type` %in% "non-DEG") %>%
    inner_join(., df_limma0,
               by = c("Gene symbol" = "gene_symbol", "Sample pair" = "group"),
               relationship = "many-to-many") %>%
    dplyr::rename(FC_rna = FC) %>%
    group_by(`Gene symbol`, `Sample pair`, FC_rna) %>%
    summarise_at("FC_pep", mean)


## Figure 3f: 散点图确认相关性 -------------------
colors.group <- c("#4CC3D9", "#FFC65D", "#F16745", "#999")
names(colors.group) <- c("D5/D6", "F7/D6", "M8/D6", "Reversed")

sub_test <- df_test_combined %>%
  mutate_at(3:4, log2) %>%
  filter(FC_pep * FC_rna > 0) %>%
  mutate(group_new = ifelse(FC_pep * FC_rna < 0, "Reversed", `Sample pair`))

p_cor_final <- ggplot(sub_test, aes(x = FC_rna, y = FC_pep)) +
  geom_point(aes(fill = group_new), size = 5, shape = 21, color = "black") +
  geom_smooth(method = "lm", se = T, linetype = 1, linewidth = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "black", lty = 2) +
  ggpubr::stat_cor(size = 6) +
  # ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~`,`~")),
  #                               label.y = c(8, 8), size = 6) +
  scale_fill_manual(values = colors.group) +
  facet_wrap( ~ `Sample pair`, scales = "fixed") +
  scale_y_continuous(limits = c(-8, 8), name = "DDA/DIA Ratios (Latest)") +
  scale_x_continuous(limits = c(-8, 8), name = "RNA-seq Ratios") +
  theme_bw() +
  theme(legend.position = "none",
        aspect.ratio = 1,
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        strip.background = element_blank(),
        panel.spacing = unit(1, "cm"),
        plot.margin = unit(c(.5,.5, .5, .5), "cm"))

ggsave("./results/figures/figure3f.pdf", p_cor_final, height = 5.2, width = 14, limitsize = FALSE)


## Figure 3e: 更正研制报告的数据 -------------------------
target_10peptides <- c("FDSDVGEFR", "LGVIEDHSNR",
                       "GVVDSDDLPLNVSR", "TVDNFVALATGEK",
                       "ELEEIVQPIISK", "ITPSYVAFTPEGER",
                       "FQSSHHPTDITSLDQYVER", "SGEVYTCQVEHPSVTSPLTVEWR",
                       "NEAIQAAHDAVAQEGQCR", "AAVDTYCR")

val_dt <- data.frame(
  peptide_sequence = rep(c("FDSDVGEFR", "LGVIEDHSNR",
                           "GVVDSDDLPLNVSR", "TVDNFVALATGEK",
                           "ELEEIVQPIISK", "ITPSYVAFTPEGER",
                           "FQSSHHPTDITSLDQYVER", "SGEVYTCQVEHPSVTSPLTVEWR",
                           "NEAIQAAHDAVAQEGQCR", "AAVDTYCR"), 3),
  group = rep(c("D5/D6", "F7/D6", "M8/D6"), each = 10),
  FC_validation = c(0.83, 0.46, 0.48, 0.55, 0.47, 0.41, 0.24, 0.91, 3.12, 0.87,
                    0.48, 0.93, 0.93, 1.04, 1.10, 1.09, 0.91, 0.41, 0.55, 0.37,
                    0.74, 0.54, 0.57, 0.60, 0.52, 0.47, 0.39, 0.85, 2.55, 0.98)
)

limma_tables0 <- readRDS("./results/tables/2_limma_bylab_bytube.rds")
df_limma <- limma_tables0 %>%
  rbindlist(.) %>%
  mutate(group = paste(group2, group1, sep = "/")) %>%
  mutate_at("group", ~ gsub("Quartet ", "", .)) %>%
  filter(adj.P.Val < 0.05) %>%
  mutate(FC = 2 ^ logFC) 

target_10peptides <- c("FDSDVGEFR", "LGVIEDHSNR",
                       "GVVDSDDLPLNVSR", "TVDNFVALATGEK",
                       "ELEEIVQPIISK", "ITPSYVAFTPEGER",
                       "FQSSHHPTDITSLDQYVER", "SGEVYTCQVEHPSVTSPLTVEWR",
                       "NEAIQAAHDAVAQEGQCR", "AAVDTYCR")

target_10peptides_tables <- pblapply(limma_tables0, function(tmp_table) {
  tmp_table_i <- tmp_table %>% filter(peptide_sequence %in% target_10peptides)
  return(tmp_table_i)
})
target_10peptides_df <- target_10peptides_tables %>%
  rbindlist(., idcol = "tmp") %>%
  mutate(group = paste(group2, group1, sep = "/")) %>%
  mutate_at("group", ~ gsub("Quartet ", "", .)) %>%
  group_by(peptide_sequence, protein_id, group) %>%
  summarise(FC_final = 2 ^ (mean(logFC)))

sub_ref <- target_10peptides_df %>%
  filter(!group %in% "HeLa/HEK293T") %>%
  filter(!protein_id %in% c("A0A140T9P7", "A0A7P0Z497", "D6R956")) %>%
  left_join(., val_dt, by = c("peptide_sequence", "group"))


## Figure 3e: 散点图确认相关性 -------------------
colors.group <- c("#4CC3D9", "#FFC65D", "#F16745", "#999")
names(colors.group) <- c("D5/D6", "F7/D6", "M8/D6", "Reversed")

## 手动拟合线性模型
# sub_ref <- sub_ref %>% mutate_at(c("FC_final", "FC_validation"), log2)
calibrator_tables <- split(sub_ref, sub_ref$group)
linear_tables <- pblapply(calibrator_tables, function(tmp_table) {
  
  x_pos <- min(tmp_table$FC_validation) + (max(tmp_table$FC_validation) - min(tmp_table$FC_validation)) * 0.02
  y_pos <- max(tmp_table$FC_final)
  
  model <- lm(FC_final ~ FC_validation, data = tmp_table)
  
  df_tmp <- data.frame(
    r_squared = summary(model)$r.squared,
    intercept = coef(model)[1],
    slope = coef(model)[2],
    x_axis = x_pos,
    y_axis = y_pos
  )
  return(df_tmp)
})

## 整理线性方程和R2用于绘图
df_calibrator_linear <- linear_tables %>%
  rbindlist(., idcol = "group") %>%
  mutate(label = sprintf("y = %.4f %+.4f x\nR² = %.4f", intercept, slope, r_squared))

## 绘制散点图
p_cor_final <- ggplot(sub_ref, aes(x = log2(FC_validation), y = log2(FC_final))) +
  geom_point(aes(fill = group), size = 5, shape = 21, color = "black") +
  geom_smooth(method = "lm", se = T, linetype = 1, linewidth = 0.5) +
  ggpubr::stat_cor(size = 6) +
  geom_abline(slope = 1, intercept = 0, color = "black", lty = 2) +
  # ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~`,`~")), size = 6) +
  # geom_text(data = df_calibrator_linear,
  #           aes(x = x_axis, y = y_axis, label = label),
  #           size = 6, hjust = 0, vjust = 1, inherit.aes = FALSE) +
  scale_fill_manual(values = colors.group) +
  facet_wrap( ~ group) +
  scale_y_continuous(name = "DDA/DIA Ratios (Latest)") +
  scale_x_continuous(name = "ID-MS Ratios") +
  theme_bw() +
  theme(legend.position = "none",
        aspect.ratio = 1,
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        strip.background = element_blank(),
        panel.spacing = unit(1, "cm"),
        plot.margin = unit(c(.5,.5, .5, .5), "cm"))

ggsave("./results/figures/figure3e.pdf", p_cor_final, height = 5.2, width = 14, limitsize = FALSE)

