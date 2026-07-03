## Title: 确定参考量值: 2025年10组重新搜库数据。
## Author: Qiaochu Chen
## Date: Jun 5th, 2025

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
library(scales)

rm(list = ls())
gc()


## 准备数据3——合并搜库整理所有样本-肽段序列-表达谱 ------------------
all_meta <- fread("./data/multilab/metadata_2025_10labs.csv")

## 4组DDA合并搜库(凤凰中心/浙江大学/青莲百奥/复旦大学)
dda_files <- list.files("./data/multilab/2025_grouped_dda", recursive = TRUE, full.names = TRUE)
dda_pep_files <- dda_files[grepl("peptides.txt", dda_files)]
dda_all <- pblapply(dda_pep_files, function(file_id) {
  
  dda_tmp <- fread(file_id, showProgress = FALSE)
  
  if (grepl("timstof", file_id)) {
    sample_group <- str_extract(file_id, "(?<=2025_grouped_dda\\/).+(?=\\/timstof)")
    sample_group <- ifelse(sample_group %in% c("D5", "D6", "F7", "M8"),
                           paste("Quartet-", sample_group, sep = ""),
                           sample_group)
    sample_group <- ifelse(sample_group %in% "HeLa", "Hela", sample_group)
    id_prefix <- paste("_cpfs01_projects-HDD_cfff-e44ef5cf7aa5_HDD_cqc_21112030002_data_proteomics_fudan_university_HT-SLM-ZYT--CQC-20250324-",
                       sample_group, "-", sep = "")
  } else {
    id_prefix <- "_cpfs01_projects-HDD_cfff-e44ef5cf7aa5_HDD_cqc_21112030002_data_proteomics_"
  }
  
  dda_tmp <- dda_tmp %>%
    select(Sequence, starts_with("Intensity ")) %>%
    dplyr::rename(peptide_sequence = Sequence) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id") %>%
    filter(value != 0) %>%
    mutate_at("analysis_id", ~ paste(id_prefix, str_extract(., "(?<=Intensity ).+"), sep = "")) %>%
    left_join(., all_meta, by = "analysis_id")
  
  return(dda_tmp)
})

## 6组DIA合并搜库(复旦/计量院/中科院天津所/中医科学院/赛默飞世尔/百泰派克/易算)
dia_files <- list.files("./data/multilab/2025_grouped_dia", recursive = TRUE, full.names = TRUE)
dia_pep_files <- dia_files[grepl("report.pr_matrix.tsv", dia_files)]
dia_all <- pblapply(dia_pep_files, function(file_id) {
  
  dia_tmp <- fread(file_id, showProgress = FALSE)
  
  dia_tmp <- dia_tmp %>%
    select(Stripped.Sequence, contains("cpfs01")) %>%
    dplyr::rename(peptide_sequence = Stripped.Sequence) %>%
    dplyr::rename_if(is.numeric, ~ str_extract(., "[^/]+(?=\\.)")) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id", na.rm = TRUE) %>%
    left_join(., all_meta, by = "analysis_id")
  
  return(dia_tmp)
})

## 合并所有数据
all_tables <- c(dda_all, dia_all)
grouped_tables <- all_tables %>% rbindlist %>% split(., by = "sample")
saveRDS(grouped_tables, "./data/multilab/quantdata_list_pep_2025_10labs.rds")


## 准备数据4——PCA马氏距离法检查离群值+All Ratio by D6/HEK293T -----------------
grouped_tables <- readRDS("./data/multilab/quantdata_list_pep_2025_10labs.rds")
all_tables <- grouped_tables %>% rbindlist %>% split(., by = "lab_id")
rm(list = setdiff(ls(), c("all_tables")))
gc()

## 每家实验室PCA检查离群值: 马氏距离法
source("./PCA.R")
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


## 准备数据5——PCA马氏距离法检查离群值+Quartet Ratio by D6 ------------------
grouped_tables <- readRDS("./data/multilab/quantdata_list_pep_2025_10labs.rds")
all_tables <- grouped_tables %>% rbindlist %>% split(., by = "lab_id")
rm(list = setdiff(ls(), c("all_tables")))
gc()

## 每家实验室PCA检查离群值: 马氏距离法(只对Quartet样本)
source("./PCA.R")
madist_results_tables <- pblapply(all_tables, function(tmp_table) {
 
  metadata <- tmp_table %>%
    distinct(analysis_id, lab_id, sample, tube, injection, order) %>%
    filter(sample %in% c("D5", "D6", "F7", "M8")) %>%
    tibble::column_to_rownames("analysis_id")
  
  exprdata_t <- tmp_table %>%
    reshape2::dcast(peptide_sequence ~ analysis_id, value.var = "value", fun.aggregate = sum) %>%
    select(all_of(rownames(metadata))) %>%
    mutate_all(log2) %>%
    filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
    # filter(apply(., 1, function(x) sum(is.na(x)) < 0.2 * nrow(metadata))) %>%
    # mutate_all(~ ifelse(is.na(.), median(., na.rm = TRUE), .)) %>%
    na.omit %>%
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
fwrite(sub_meta2outlier, "./results/tables/2_outlier_madist_quartet.csv")

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

## 每家实验室Quartet ratio by d6(只对Quartet样本)
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
saveRDS(ratiobyd6_tables, "./results/tables/2_quantdata_list_pep_ratiobyd6_2025_10labs.rds")


## 准备数据6——PCA马氏距离法检查离群值+HeLa Ratio by HEK293T ------------------
grouped_tables <- readRDS("./data/multilab/quantdata_list_pep_2025_10labs.rds")
sub_grouped_tables <- grouped_tables[c("HeLa", "HEK293T")]
all_tables <- sub_grouped_tables %>% rbindlist %>% split(., by = "lab_id")
rm(list = setdiff(ls(), c("all_tables")))
gc()

## 每家实验室PCA检查离群值: 马氏距离法(只对HeLa+HEK293T样本)
source("./PCA.R")
madist_results_tables <- pblapply(all_tables, function(tmp_table) {
  
  metadata <- tmp_table %>%
    distinct(analysis_id, lab_id, sample, tube, injection, order) %>%
    filter(sample %in% c("HeLa", "HEK293T")) %>%
    tibble::column_to_rownames("analysis_id")
  
  exprdata_t <- tmp_table %>%
    reshape2::dcast(peptide_sequence ~ analysis_id, value.var = "value", fun.aggregate = sum) %>%
    select(all_of(rownames(metadata))) %>%
    mutate_all(log2) %>%
    filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
    # filter(apply(., 1, function(x) sum(is.na(x)) < 0.2 * nrow(metadata))) %>%
    # mutate_all(~ ifelse(is.na(.), median(., na.rm = TRUE), .)) %>%
    na.omit %>%
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
fwrite(sub_meta2outlier, "./results/tables/2_outlier_madist_hela_hek293t.csv")

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

## 每家实验室Quartet ratio by HEK293T(只对HeLa样本)
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
  
  expr_tmp_ratiobyhek293t <- ratio_by_ref(expr_tmp, meta_tmp, ref_sample = "HEK293T", log2_transformed = FALSE)
  
  tmp_table_ratiobyhek293t <- expr_tmp_ratiobyhek293t %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id", na.rm = TRUE) %>%
    dplyr::rename(peptide_sequence = feature) %>%
    left_join(., meta_tmp, by = "analysis_id")
  
  return(tmp_table_ratiobyhek293t)
})

## 合并所有数据
ratiobyhek293t_tables <- ratiobyhek293t_tables %>% rbindlist %>% split(., by = "sample")
saveRDS(ratiobyhek293t_tables, "./results/tables/2_quantdata_list_pep_ratiobyhek293t_2025_10labs.rds")


## 准备数据7——检查每个实验室的PCA/SNR ----------------------
rm(list = ls())
gc()
source("./PCA.R")
grouped_tables <- readRDS("./data/multilab/quantdata_list_pep_2025_10labs.rds")
# grouped_tables <- readRDS("./results/tables/2_quantdata_list_pep_ratiobyd6_2025_10labs.rds")
all_tables <- grouped_tables %>% rbindlist %>% split(., by = "lab_id")

labels.lab <- c("QLB", "NCP", "ZJU", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")
names(labels.lab) <- c("qinglian_bio", "phoenix", "zhejiang_university", "fudan_university",
                "national_institute_of_methodology", "biotech_pack", "thermofisher_shanghai",
                "omicsolution", "cas_tianjin", "academy_of_chinese_medical_sciences")

all_tables <- all_tables[names(labels.lab)]

pca_tables_bydataset <- pblapply(all_tables, function(tmp_table) {
  
  metadata <- tmp_table %>%
    distinct(analysis_id, lab_id, sample, tube, injection, order) %>%
    filter(sample %in% c("D5", "D6", "F7", "M8")) %>%
    tibble::column_to_rownames("analysis_id")

  exprdata_t <- tmp_table %>%
    # reshape2::dcast(peptide_sequence ~ analysis_id, value.var = "value") %>% ## 若输入的是原始强度的数据则此步隐去
    reshape2::dcast(peptide_sequence ~ analysis_id, value.var = "value", fun.aggregate = sum) %>% ## 若输入的是Ratio by D6后的数据则此步隐去
    select(all_of(rownames(metadata))) %>%
    mutate_all(log2) %>% ## 若输入的是Ratio by D6后的数据则此步隐去
    filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
    # filter(apply(., 1, function(x) sum(is.na(x)) < 0.2 * nrow(metadata))) %>%
    # mutate_all(~ ifelse(is.na(.), median(., na.rm = TRUE), .)) %>%
    # mutate_all(~ ifelse(. == 0, NA, .)) %>%
    na.omit %>%
    t
  
  metadata <- metadata %>% filter(!sample %in% "Quartet D6")
  exprdata_t <- exprdata_t[all_of(rownames(metadata)), ]
  
  pca_results <- main_pca(exprdata_t = exprdata_t, metadata = metadata,
                          center = TRUE, scale = TRUE, group = "sample",
                          biplot = FALSE, dictGroups = metadata$sample,
                          snr = TRUE, plot = FALSE)
  
  pca_results$feature_n <- ncol(exprdata_t)
  pca_results$sample_n <- nrow(exprdata_t)
  
  return(pca_results)
})

saveRDS(pca_tables_bydataset, "results/tables/2_pca_bylab.rds")
# saveRDS(pca_tables_bydataset, "results/tables/2_pca_bylab_bytube_ratiobyd6.rds")


## 准备数据8——计算Quartet合并实验室D6比例定量前后的PCA/SNR (删除ZJU) ----------------------
rm(list = ls())
gc()
# meta_quartet <- fread("./results/tables/2_outlier_madist_quartet.csv")
meta_quartet <- fread("./results/tables/2_outlier_madist_all.csv")
# all_tables <- readRDS("./data/multilab/quantdata_list_pep_2025_10labs.rds")
all_tables <- readRDS("./results/tables/2_quantdata_list_pep_allratiobyd6_2025_10labs.rds")
filtered_tables_groupedbylabs <- all_tables %>% rbindlist %>% split(., by = "lab_id")
filtered_tables_groupedbylabs <- filtered_tables_groupedbylabs[c(1, 2, 4:10)]

## PCA
metadata <- meta_quartet %>%
  filter(is.outlier == FALSE) %>%
  filter(!lab_id %in% "zhejiang_university") %>%
  filter(sample %in% c("D5", "D6", "F7", "M8")) %>%
  tibble::column_to_rownames("analysis_id")

exprdata_t <- filtered_tables_groupedbylabs %>%
  rbindlist %>%
  reshape2::dcast(peptide_sequence ~ analysis_id, value.var = "value") %>% ## 若输入的是原始强度的数据则此步隐去
  # reshape2::dcast(peptide_sequence ~ analysis_id, value.var = "value", fun.aggregate = sum) %>% ## 若输入的是Ratio by D6后的数据则此步隐去
  select(any_of(rownames(metadata))) %>%
  # mutate_all(~ ifelse(. == 0, NA, .)) %>% ## 若输入的是Ratio by D6后的数据则此步隐去
  # mutate_all(log2) %>% ## 若输入的是Ratio by D6后的数据则此步隐去
  filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
  filter(apply(., 1, function(x) sum(is.na(x)) < 0.2 * nrow(metadata))) %>%
  mutate(row_mean = apply(., 1, mean, na.rm = TRUE)) %>%
  mutate_all(~ ifelse(is.na(.), row_mean, .)) %>%
  select(!row_mean) %>%
  # na.omit %>%
  t

dim(exprdata_t)

rm(list = setdiff(ls(), c("exprdata_t", "metadata")))
gc()

source("./PCA.R")
pca_results <- main_pca(exprdata_t = exprdata_t, metadata = metadata,
                        center = TRUE, scale = TRUE, group = "sample",
                        biplot = FALSE, dictGroups = metadata$sample,
                        snr = TRUE, plot = FALSE)

pca_results$feature_n <- ncol(exprdata_t)
pca_results$sample_n <- nrow(exprdata_t)

# saveRDS(pca_results, "./results/tables/2_pca_quartet_combined.rds")
saveRDS(pca_results, "./results/tables/2_pca_quartet_combined_bytube_ratiobyd6.rds")


## 准备数据9——limma检验计算比例定量值 ------------------
rm(list = ls())
gc()
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

rm(list = setdiff(ls(), c("filtered_tables")))
gc()

## Limma: multi-lab (文章：https://doi.org/10.1038/s41467-024-47899-w)
all_tables <- filtered_tables %>% rbindlist %>% split(., by = c("lab_id", "tube"))
rm(list = setdiff(ls(), c("all_tables")))
gc()

source("./DEP.R")
limma_tables <- pblapply(all_tables, function(tmp_table) {
  
  all_sample_pairs <- list(c("D6", "D5"), c("D6", "F7"),
                           c("D6", "M8"), c("HEK293T", "HeLa"))
  
  limma_tables_i <- mclapply(all_sample_pairs, function(sample_pair_id) {
    
    print(sample_pair_id)
    sub_wide_j <- tmp_table %>%
      filter(sample %in% sample_pair_id)
    
    if (nrow(sub_wide_j) == 0) return(NULL)
    
    sub_wide_j <- tmp_table %>%
      filter(sample %in% sample_pair_id) %>%
      reshape2::dcast(peptide_sequence + protein_id ~ sample + tube + injection, value.var = "value", fun.aggregate = sum) %>%
      mutate_if(is.numeric, ~ ifelse(is.na(.), 0, .)) %>%
      mutate(rowname = paste(peptide_sequence, protein_id)) %>%
      tibble::column_to_rownames(c("rowname")) %>%
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

limma_tables <- limma_tables %>%
  rbindlist(., idcol = "tmp") %>%
  tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
  tidyr::separate(tmp, c("lab_id", "tube"), sep = "\\.") %>%
  mutate(group = paste(group2, group1, sep = "/")) %>%
  mutate_at("group", ~ gsub("Quartet ", "", .)) %>%
  split(., ~ group)

saveRDS(limma_tables, "./results/tables/2_limma_bylab_bytube.rds")


## 根据比例量值可重复性筛选肽段标称特性定量值(剔除离群样本) ------------------
limma_tables <- readRDS("./results/quartet/tables/2_limma_bylab_bytube.rds")
rm(list = setdiff(ls(), c("limma_tables")))
gc()

qualipr_tables <- readRDS("./results/tables/1_qualiprop_list_PEPfdr0.01.rds")
qualipr_tables1 <- pblapply(qualipr_tables[c(1, 3, 6)], function(tmp_table) {
  tmp_table_i <- tmp_table %>%
    distinct(lab_id, tube, peptide_sequence, protein_id) %>%
    inner_join(., qualipr_tables[[2]], by = c("lab_id", "tube", "peptide_sequence", "protein_id"))
  return(tmp_table_i)
})
qualipr_tables2 <- pblapply(qualipr_tables[5], function(tmp_table) {
  tmp_table_i <- tmp_table %>%
    distinct(lab_id, tube, peptide_sequence, protein_id) %>%
    inner_join(., qualipr_tables[[4]], by = c("lab_id", "tube", "peptide_sequence", "protein_id"))
  return(tmp_table_i)
})
qualipr_tables <- c(qualipr_tables1[1:2], qualipr_tables2, qualipr_tables1[3])
names(qualipr_tables) <- names(limma_tables)

limma_tables <- pblapply(1:4, function(i) {
  tmp_table_i <- qualipr_tables[[i]] %>%
    distinct(lab_id, tube, peptide_sequence, protein_id) %>%
    mutate_at("tube", as.character) %>%
    inner_join(., limma_tables[[i]], c("lab_id", "tube", "peptide_sequence", "protein_id"))
  return(tmp_table_i)
})
names(limma_tables) <- names(qualipr_tables)

rm(list = setdiff(ls(), c("limma_tables")))
gc()
stat_pep_tables <- pblapply(limma_tables, function(tmp_table) {
  stat_tmp_table <- tmp_table %>%
    ungroup() %>%
    summarise(peptide_n = length(unique(peptide_sequence)),
              protein_n = length(unique(protein_id)), .groups = "drop")
  return(stat_tmp_table)
})
stat_df0 <- rbindlist(stat_pep_tables, idcol = "sample")
stat_df0$tier = "All paired qualitative properties."

## Tier1: 至少3家实验室所有管间可重复性满足CV小于20%
passCV0.2_pep_tables <- pblapply(limma_tables, function(tmp_table) {
  
  df_cv_i <- tmp_table %>%
    mutate(value = as.numeric(2 ^ logFC)) %>%
    group_by(peptide_sequence, protein_id, lab_id) %>%
    summarise(cv = sd(value) / mean(value), .groups = "drop") %>%
    na.omit
  
  sub_tmp_table <- df_cv_i %>%
    group_by(peptide_sequence, protein_id) %>%
    summarise(lab_n = length(unique(lab_id[cv < 0.2])), .groups = "drop") %>%
    filter(lab_n >= 3) %>%
    distinct(peptide_sequence, protein_id)
  
  sub_pep_table <- tmp_table %>%
    inner_join(., sub_tmp_table, by = c("peptide_sequence", "protein_id"))
  
  return(sub_pep_table)
})
stat_pep_tables <- pblapply(passCV0.2_pep_tables, function(tmp_table) {
  stat_tmp_table <- tmp_table %>%
    ungroup() %>%
    summarise(peptide_n = length(unique(peptide_sequence)),
              protein_n = length(unique(protein_id)), .groups = "drop")
  return(stat_tmp_table)
})
stat_df1 <- rbindlist(stat_pep_tables, idcol = "sample")
stat_df1$tier = "At least 3 labs support CV < 20%."

## Tier2: 至少3家实验室所有管满足p小于0.05+倍数变化至少1.2
passP0.05_pep_tables <- pblapply(passCV0.2_pep_tables, function(tmp_table) {
  
  sub_tmp_table <- sub_tmp_table <- tmp_table %>%
    group_by(peptide_sequence, lab_id) %>%
    summarise(pass_n = length(unique(tube[P.Value < .05 & abs(logFC) >= log2(1.2)])),
              total_n = length(unique(tube)), .groups = "drop") %>%
    filter(total_n - pass_n == 0) %>%
    group_by(peptide_sequence) %>%
    summarise(lab_n = length(unique(lab_id)), .groups = "drop") %>%
    filter(lab_n >= 3) %>%
    distinct(peptide_sequence)
  
  sub_pep_table <- tmp_table %>%
    inner_join(., sub_tmp_table, by = c("peptide_sequence"))
  
  return(sub_pep_table)
})
stat_pep_tables <- pblapply(passP0.05_pep_tables, function(tmp_table) {
  stat_tmp_table <- tmp_table %>%
    ungroup() %>%
    summarise(peptide_n = length(unique(peptide_sequence)),
              protein_n = length(unique(protein_id)), .groups = "drop")
  return(stat_tmp_table)
})
stat_df2 <- rbindlist(stat_pep_tables, idcol = "sample")
stat_df2$tier = "At least 3 labs support limma'test P < 0.05 and |log2FC| > 1.2."

## Tier3: 实验室间所有管支持FDR校正后P < 0.05
passFDR0.05_pep_tables <- pblapply(passP0.05_pep_tables, function(tmp_table) {
  sub_tmp_table <- tmp_table %>%
    group_by(peptide_sequence, protein_id) %>%
    mutate(fdr = p.adjust(P.Value, method = "BH")) %>%
    filter(fdr < .05)
  
  return(sub_tmp_table)
})
stat_pep_tables <- pblapply(passFDR0.05_pep_tables, function(tmp_table) {
  stat_tmp_table <- tmp_table %>%
    ungroup() %>%
    summarise(peptide_n = length(unique(peptide_sequence)),
              protein_n = length(unique(protein_id)), .groups = "drop")
  return(stat_tmp_table)
})
stat_df3 <- rbindlist(stat_pep_tables, idcol = "sample")
stat_df3$tier = "All evidences support limma'test P (FDR adjusted) < 0.05."

## 合并统计结果
stat_df_filtered <- rbindlist(list(stat_df0, stat_df1, stat_df2, stat_df3))
stat_df_filtered_wide <- stat_df_filtered %>%
  mutate(tier = factor(tier, levels = unique(tier))) %>%
  mutate(`n (peptides/proteins)` = paste(comma(peptide_n), "\n(", comma(protein_n), ")", sep = "")) %>%
  reshape2::dcast(., tier ~ sample, value.var = "n (peptides/proteins)")

fwrite(stat_df_filtered_wide, "~/Desktop/tmp_定量.csv")
saveRDS(passFDR0.05_pep_tables, "./results/tables/2_quantprop_list_Pfdr0.05.rds")

