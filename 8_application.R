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

rm(list = ls())
gc()


## 0.1 DDA: 从 evidence/peptides/proteinGroups 中统计 DEPs ------
rm(list = ls())
gc()
source("./DEP.R")

meta_all <- fread("./data/combined_meta.csv")
all_files <- list.files("./data/multilab", recursive = TRUE, full.names = TRUE)

dda_pre_files <- all_files[grepl("2025_grouped_dda.+evidence.txt", all_files)]
dda_pre_tables <- pblapply(dda_pre_files, function(file_id) {
  
  dda_tmp <- fread(file_id, showProgress = FALSE)
  dda_tmp <- dda_tmp %>%
    dplyr::rename(analysis_id = `Raw file`, value = Intensity) %>%
    mutate(feature = paste(`Modified sequence`, Charge)) %>%
    select(feature, analysis_id, value) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dda_tmp)
})
dda_pre_tables <- dda_pre_tables %>% rbindlist %>% split(., by = "lab")

dda_pep_files <- all_files[grepl("2025_grouped_dda.+peptides.txt", all_files)]
dda_pep_tables <- pblapply(dda_pep_files, function(file_id) {
  
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
    dplyr::rename(feature = Sequence) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id") %>%
    filter(value != 0) %>%
    mutate_at("analysis_id", ~ paste(id_prefix, str_extract(., "(?<=Intensity ).+"), sep = "")) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dda_tmp)
})
dda_pep_tables <- dda_pep_tables %>% rbindlist %>% split(., by = "lab")

dda_pro_files <- all_files[grepl("2025_grouped_dda.+proteinGroups.txt", all_files)]
dda_pro_tables <- pblapply(dda_pro_files, function(file_id) {
  
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
    select(`Protein names`, starts_with("LFQ intensity ")) %>%
    dplyr::rename(feature = `Protein names`) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id") %>%
    filter(value != 0) %>%
    filter(!feature %in% "") %>%
    mutate_at("analysis_id", ~ paste(id_prefix, str_extract(., "(?<=LFQ intensity ).+"), sep = "")) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dda_tmp)
})
dda_pro_tables <- dda_pro_tables %>% rbindlist %>% split(., by = "lab")

all_sample_pairs <- list(c("D6", "D5"), c("D6", "F7"),
                         c("D6", "M8"), c("HEK293T", "HeLa"))
dda_limma_tables <- pblapply(list(dda_pre_tables, dda_pep_tables, dda_pro_tables), function(tmp_tables) {
  
  limma_results_i <- pblapply(tmp_tables, function(tmp_table) {
    
    limma_tables_j <- mclapply(all_sample_pairs, function(sample_pair_id) {
      sub_wide_k <- tmp_table %>%
        filter(tube == 1) %>%
        filter(sample %in% sample_pair_id) %>%
        reshape2::dcast(feature ~ sample + tube + injection, value.var = "value", fun.aggregate = sum) %>%
        mutate_if(is.numeric, ~ ifelse(is.na(.), 0, .)) %>%
        tibble::column_to_rownames(c("feature"))
      
      group_k <- tmp_table %>%
        filter(tube == 1) %>%
        filter(sample %in% sample_pair_id) %>%
        distinct(analysis_id, sample) %>%
        pull(sample) %>%
        factor(., levels = sample_pair_id, ordered = TRUE)
      
      if (length(which(table(group_k) == 1)) > 0) return(NULL)
      
      limma_results_k <- calculate_p(exprdata = sub_wide_k, group = group_k,
                                     test_method = "limma", p_adjust = "BH",
                                     na_threshold = 0)
      
      return(limma_results_k)
    })
    
    limma_results_j <- limma_tables_j %>% rbindlist
    return(limma_results_j)
  })
  return(limma_results_i)
})

saveRDS(dda_limma_tables, "./results/tables/dda_limma.rds")


## 0.2 DIA: 从 report.parquet/report.pr_matrix/report.pg_matrix 中统计 DEPs ----
rm(list = ls())
gc()
source("./DEP.R")

meta_all <- fread("./data/combined_meta.csv")
all_files <- list.files("./data/multilab", recursive = TRUE, full.names = TRUE)

columns2keep <- c("Run", "Precursor.Id", "Ms1.Area")
dia_pre_files <- all_files[grepl("report.parquet", all_files)]
dia_pre_tables <- pblapply(dia_pre_files, function(file_id) {
  
  dia_tmp <- read_parquet(file_id, col_select = columns2keep)
  dia_tmp <- dia_tmp %>%
    dplyr::rename(analysis_id = Run, feature = Precursor.Id, value = Ms1.Area) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dia_tmp)
})
dia_pre_tables <- dia_pre_tables %>% rbindlist %>% split(., by = "lab")
gc()

dia_pep_files <- all_files[grepl("report.pr_matrix", all_files)]
dia_pep_tables <- pblapply(dia_pep_files, function(file_id) {
  
  dia_tmp <- fread(file_id)
  dia_tmp <- dia_tmp %>%
    select(Stripped.Sequence, contains("cpfs01")) %>%
    dplyr::rename(feature = Stripped.Sequence) %>%
    dplyr::rename_if(is.numeric, ~ str_extract(., "[^/]+(?=\\.)")) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id", na.rm = TRUE) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dia_tmp)
})
dia_pep_tables <- dia_pep_tables %>% rbindlist %>% split(., by = "lab")
gc()

dia_pro_files <- all_files[grepl("report.pg_matrix", all_files)]
dia_pro_tables <- pblapply(dia_pro_files, function(file_id) {
  
  dia_tmp <- fread(file_id)
  dia_tmp <- dia_tmp %>%
    select(First.Protein.Description, contains("cpfs01")) %>%
    dplyr::rename(feature = First.Protein.Description) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id", na.rm = TRUE) %>%
    mutate_at("analysis_id", ~ str_extract(., "[^/]+(?=\\.)")) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dia_tmp)
})
dia_pro_tables <- dia_pro_tables %>% rbindlist %>% split(., by = "lab")
gc()

all_sample_pairs <- list(c("D6", "D5"), c("D6", "F7"),
                         c("D6", "M8"), c("HEK293T", "HeLa"))
dia_limma_tables <- pblapply(list(dia_pre_tables, dia_pep_tables, dia_pro_tables), function(tmp_tables) {
  
  limma_results_i <- pblapply(tmp_tables, function(tmp_table) {
    
    limma_tables_j <- mclapply(all_sample_pairs, function(sample_pair_id) {
      sub_wide_k <- tmp_table %>%
        filter(tube == 1) %>%
        filter(sample %in% sample_pair_id) %>%
        reshape2::dcast(feature ~ sample + tube + injection, value.var = "value", fun.aggregate = sum) %>%
        mutate_if(is.numeric, ~ ifelse(is.na(.), 0, .)) %>%
        tibble::column_to_rownames(c("feature"))
      
      group_k <- tmp_table %>%
        filter(tube == 1) %>%
        filter(sample %in% sample_pair_id) %>%
        distinct(analysis_id, sample) %>%
        pull(sample) %>%
        factor(., levels = sample_pair_id, ordered = TRUE)
      
      if (length(which(table(group_k) == 1)) > 0) return(NULL)
      
      limma_results_k <- calculate_p(exprdata = sub_wide_k, group = group_k,
                                     test_method = "limma", p_adjust = "BH",
                                     na_threshold = 0)
      
      return(limma_results_k)
    })
    
    limma_results_j <- limma_tables_j %>% rbindlist
    return(limma_results_j)
  })
  return(limma_results_i)
})

saveRDS(dia_limma_tables, "./results/tables/dia_limma.rds")


## 1.1 DDA: 从 evidence/peptides/proteinGroups 中统计全特征 SNR 指标 ------
rm(list = ls())
gc()
source("./PCA.R")

meta_all <- fread("./data/combined_meta.csv")
all_files <- list.files("./data/multilab", recursive = TRUE, full.names = TRUE)

dda_pre_files <- all_files[grepl("2025_grouped_dda/(D5|D6|F7|M8).+evidence.txt", all_files)]
dda_pre_tables <- pblapply(dda_pre_files, function(file_id) {
  
  dda_tmp <- fread(file_id, showProgress = FALSE)
  dda_tmp <- dda_tmp %>%
    dplyr::rename(analysis_id = `Raw file`, value = Intensity) %>%
    mutate(feature = paste(`Modified sequence`, Charge)) %>%
    select(feature, analysis_id, value) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dda_tmp)
})
dda_pre_tables <- dda_pre_tables %>% rbindlist %>% split(., by = "lab")

dda_pep_files <- all_files[grepl("2025_grouped_dda/(D5|D6|F7|M8).+peptides.txt", all_files)]
dda_pep_tables <- pblapply(dda_pep_files, function(file_id) {
  
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
    dplyr::rename(feature = Sequence) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id") %>%
    filter(value != 0) %>%
    mutate_at("analysis_id", ~ paste(id_prefix, str_extract(., "(?<=Intensity ).+"), sep = "")) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dda_tmp)
})
dda_pep_tables <- dda_pep_tables %>% rbindlist %>% split(., by = "lab")

dda_pro_files <- all_files[grepl("2025_grouped_dda/(D5|D6|F7|M8).+proteinGroups.txt", all_files)]
dda_pro_tables <- pblapply(dda_pro_files, function(file_id) {
  
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
    select(`Protein names`, starts_with("LFQ intensity ")) %>%
    dplyr::rename(feature = `Protein names`) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id") %>%
    filter(value != 0) %>%
    filter(!feature %in% "") %>%
    mutate_at("analysis_id", ~ paste(id_prefix, str_extract(., "(?<=LFQ intensity ).+"), sep = "")) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dda_tmp)
})
dda_pro_tables <- dda_pro_tables %>% rbindlist %>% split(., by = "lab")

## 原始信号强度水平 + 全特征范围
dda_pca_tables <- pblapply(list(dda_pre_tables, dda_pep_tables, dda_pro_tables), function(tmp_tables) {
  
  pca_results_i <- mclapply(tmp_tables, function(tmp_table) {
    
    metadata <- tmp_table %>%
      distinct(analysis_id, lab, sample, tube, injection, order) %>%
      tibble::column_to_rownames("analysis_id")
    
    exprdata_t <- tmp_table %>%
      reshape2::dcast(feature ~ analysis_id, value.var = "value", fun.aggregate = sum) %>% 
      select(all_of(rownames(metadata))) %>%
      mutate_all(~ ifelse(. == 0, NA, log2(.))) %>%
      filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
      # filter(apply(., 1, function(x) sum(is.na(x)) < 0.2 * nrow(metadata))) %>%
      # mutate_all(~ ifelse(is.na(.), median(., na.rm = TRUE), .)) %>%
      mutate_all(~ ifelse(is.na(.), 0, .)) %>%
      # na.omit %>%
      t
    
    pca_results_j <- main_pca(exprdata_t = exprdata_t, metadata = metadata,
                              center = TRUE, scale = TRUE, group = "sample",
                              biplot = FALSE, dictGroups = metadata$sample,
                              snr = TRUE, plot = FALSE)
    
    pca_results_j$feature_n <- ncol(exprdata_t)
    pca_results_j$sample_n <- nrow(exprdata_t)
    
    return(pca_results_j)
  })
  return(pca_results_i)
})
saveRDS(dda_pca_tables, "./results/tables/dda_pca_quartet.rds")

## 比例定量水平
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
dda_pca_tables <- pblapply(list(dda_pre_tables, dda_pep_tables, dda_pro_tables), function(tmp_tables) {
  
  pca_results_i <- mclapply(tmp_tables, function(tmp_table) {
    
    df_meta_tmp <- tmp_table %>%
      distinct(analysis_id, lab, sample, tube, injection, order) %>%
      mutate(batch = paste(lab, tube))
    
    df_quant_wide_tmp <- tmp_table %>%
      reshape2::dcast(feature ~ analysis_id, value.var = "value", fun.aggregate = sum)
    
    df_ratio_wide_tmp <- ratio_by_ref(df_quant_wide_tmp, df_meta_tmp, log2_transformed = TRUE)
    
    metadata <- tmp_table %>%
      distinct(analysis_id, lab, sample, tube, injection, order) %>%
      filter(!sample %in% "D6") %>%
      tibble::column_to_rownames("analysis_id")
    
    exprdata_t <- df_ratio_wide_tmp %>% 
      select(all_of(rownames(metadata))) %>%
      mutate_all(~ ifelse(. == 0, NA, log2(.))) %>%
      filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
      # filter(apply(., 1, function(x) sum(is.na(x)) < 0.2 * nrow(metadata))) %>%
      # mutate_all(~ ifelse(is.na(.), median(., na.rm = TRUE), .)) %>%
      mutate_all(~ ifelse(is.na(.), 0, .)) %>%
      # na.omit %>%
      t
    
    pca_results_j <- main_pca(exprdata_t = exprdata_t, metadata = metadata,
                              center = TRUE, scale = TRUE, group = "sample",
                              biplot = FALSE, dictGroups = metadata$sample,
                              snr = TRUE, plot = FALSE)
    
    pca_results_j$feature_n <- ncol(exprdata_t)
    pca_results_j$sample_n <- nrow(exprdata_t)
    
    return(pca_results_j)
  })
  return(pca_results_i)
})
saveRDS(dda_pca_tables, "./results/tables/dda_pca_quartet_ratiobyd6.rds")




## 1.2 DIA: 从 report.parquet/report.pr_matrix/report.pg_matrix 中统计全特征 SNR 指标 ------
rm(list = ls())
gc()
source("./PCA.R")

meta_all <- fread("./data/combined_meta.csv")
all_files <- list.files("./data/multilab", recursive = TRUE, full.names = TRUE)

columns2keep <- c("Run", "Precursor.Id", "Ms1.Area")
dia_pre_files <- all_files[grepl("2025_grouped_dia/(D5|D6|F7|M8).+report.parquet", all_files)]
dia_pre_tables <- pblapply(dia_pre_files, function(file_id) {
  
  dia_tmp <- read_parquet(file_id, col_select = columns2keep)
  dia_tmp <- dia_tmp %>%
    dplyr::rename(analysis_id = Run, feature = Precursor.Id, value = Ms1.Area) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dia_tmp)
})
dia_pre_tables <- dia_pre_tables %>% rbindlist %>% split(., by = "lab")

dia_pep_files <- all_files[grepl("2025_grouped_dia/(D5|D6|F7|M8).+report.pr_matrix", all_files)]
dia_pep_tables <- pblapply(dia_pep_files, function(file_id) {
  
  dia_tmp <- fread(file_id)
  dia_tmp <- dia_tmp %>%
    select(Stripped.Sequence, contains("cpfs01")) %>%
    dplyr::rename(feature = Stripped.Sequence) %>%
    dplyr::rename_if(is.numeric, ~ str_extract(., "[^/]+(?=\\.)")) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id", na.rm = TRUE) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dia_tmp)
})
dia_pep_tables <- dia_pep_tables %>% rbindlist %>% split(., by = "lab")

dia_pro_files <- all_files[grepl("2025_grouped_dia/(D5|D6|F7|M8).+report.pg_matrix", all_files)]
dia_pro_tables <- pblapply(dia_pro_files, function(file_id) {
  
  dia_tmp <- fread(file_id)
  dia_tmp <- dia_tmp %>%
    select(First.Protein.Description, contains("cpfs01")) %>%
    dplyr::rename(feature = First.Protein.Description) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id", na.rm = TRUE) %>%
    mutate_at("analysis_id", ~ str_extract(., "[^/]+(?=\\.)")) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dia_tmp)
})
dia_pro_tables <- dia_pro_tables %>% rbindlist %>% split(., by = "lab")

## 原始信号强度水平
dia_pca_tables <- pblapply(list(dia_pre_tables, dia_pep_tables, dia_pro_tables), function(tmp_tables) {
  
  pca_results_i <- mclapply(tmp_tables, function(tmp_table) {
    
    metadata <- tmp_table %>%
      distinct(analysis_id, lab, sample, tube, injection, order) %>%
      tibble::column_to_rownames("analysis_id")
    
    exprdata_t <- tmp_table %>%
      reshape2::dcast(feature ~ analysis_id, value.var = "value", fun.aggregate = sum) %>% 
      select(all_of(rownames(metadata))) %>%
      mutate_all(~ ifelse(. == 0, NA, log2(.))) %>%
      filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
      # filter(apply(., 1, function(x) sum(is.na(x)) < 0.2 * nrow(metadata))) %>%
      # mutate_all(~ ifelse(is.na(.), median(., na.rm = TRUE), .)) %>%
      mutate_all(~ ifelse(is.na(.), 0, .)) %>%
      # na.omit %>%
      t
    
    pca_results_j <- main_pca(exprdata_t = exprdata_t, metadata = metadata,
                              center = TRUE, scale = TRUE, group = "sample",
                              biplot = FALSE, dictGroups = metadata$sample,
                              snr = TRUE, plot = FALSE)
    
    pca_results_j$feature_n <- ncol(exprdata_t)
    pca_results_j$sample_n <- nrow(exprdata_t)
    
    return(pca_results_j)
  })
  return(pca_results_i)
})
saveRDS(dia_pca_tables, "./results/tables/dia_pca_quartet.rds")

## 比例定量水平
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
dia_pca_tables <- pblapply(list(dia_pre_tables, dia_pep_tables, dia_pro_tables), function(tmp_tables) {
  
  pca_results_i <- mclapply(tmp_tables, function(tmp_table) {
    
    df_meta_tmp <- tmp_table %>%
      distinct(analysis_id, lab, sample, tube, injection, order) %>%
      mutate(batch = paste(lab, tube))
    
    df_quant_wide_tmp <- tmp_table %>%
      reshape2::dcast(feature ~ analysis_id, value.var = "value", fun.aggregate = sum)
    
    df_ratio_wide_tmp <- ratio_by_ref(df_quant_wide_tmp, df_meta_tmp, log2_transformed = TRUE)
    
    metadata <- tmp_table %>%
      distinct(analysis_id, lab, sample, tube, injection, order) %>%
      filter(!sample %in% "D6") %>%
      tibble::column_to_rownames("analysis_id")
    
    exprdata_t <- df_ratio_wide_tmp %>% 
      select(all_of(rownames(metadata))) %>%
      mutate_all(~ ifelse(. == 0, NA, log2(.))) %>%
      filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
      # filter(apply(., 1, function(x) sum(is.na(x)) < 0.2 * nrow(metadata))) %>%
      # mutate_all(~ ifelse(is.na(.), median(., na.rm = TRUE), .)) %>%
      mutate_all(~ ifelse(is.na(.), 0, .)) %>%
      # na.omit %>%
      t
    
    pca_results_j <- main_pca(exprdata_t = exprdata_t, metadata = metadata,
                              center = TRUE, scale = TRUE, group = "sample",
                              biplot = FALSE, dictGroups = metadata$sample,
                              snr = TRUE, plot = FALSE)
    
    pca_results_j$feature_n <- ncol(exprdata_t)
    pca_results_j$sample_n <- nrow(exprdata_t)
    
    return(pca_results_j)
  })
  return(pca_results_i)
})
saveRDS(dia_pca_tables, "./results/tables/dia_pca_quartet_ratiobyd6.rds")


## 1.3 DDA: 从 evidence/peptides/proteinGroups 中统计高置信定性特征 SNR 指标 ------
rm(list = ls())
gc()
source("./PCA.R")

ref_dt_final <- fread("./results/tables/7_qualiprop_final.csv")
ref_dt_final <- ref_dt_final %>% mutate_at("sample", ~ gsub("Quartet ", "", .))

meta_all <- fread("./data/combined_meta.csv")
all_files <- list.files("./data/multilab", recursive = TRUE, full.names = TRUE)

dda_pre_files <- all_files[grepl("2025_grouped_dda/(D5|D6|F7|M8).+evidence.txt", all_files)]
dda_pre_tables <- pblapply(dda_pre_files, function(file_id) {
  
  sample_id <- str_extract(file_id, "(?<=grouped_dda/)[^/]+")
  ref_tmp <- ref_dt_final %>% filter(sample %in% sample_id)
  
  dda_tmp <- fread(file_id, showProgress = FALSE)
  dda_tmp <- dda_tmp %>%
    dplyr::rename(analysis_id = `Raw file`, value = Intensity) %>%
    inner_join(., ref_tmp, by = c("Sequence" = "peptide_sequence",
                                  "Leading razor protein" = "protein_id")) %>%
    mutate(feature = paste(`Modified sequence`, Charge)) %>%
    select(feature, analysis_id, value) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dda_tmp)
})
dda_pre_tables <- dda_pre_tables %>% rbindlist %>% split(., by = "lab")

dda_pep_files <- all_files[grepl("2025_grouped_dda/(D5|D6|F7|M8).+peptides.txt", all_files)]
dda_pep_tables <- pblapply(dda_pep_files, function(file_id) {
  
  sample_id <- str_extract(file_id, "(?<=grouped_dda/)[^/]+")
  ref_tmp <- ref_dt_final %>% filter(sample %in% sample_id)
  
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
    inner_join(., ref_tmp, by = c("Sequence" = "peptide_sequence",
                                  "Leading razor protein" = "protein_id")) %>%
    select(Sequence, starts_with("Intensity ")) %>%
    dplyr::rename(feature = Sequence) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id") %>%
    filter(value != 0) %>%
    mutate_at("analysis_id", ~ paste(id_prefix, str_extract(., "(?<=Intensity ).+"), sep = "")) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dda_tmp)
})
dda_pep_tables <- dda_pep_tables %>% rbindlist %>% split(., by = "lab")

dda_pro_files <- all_files[grepl("2025_grouped_dda/(D5|D6|F7|M8).+proteinGroups.txt", all_files)]
dda_pro_tables <- pblapply(dda_pro_files, function(file_id) {
  
  sample_id <- str_extract(file_id, "(?<=grouped_dda/)[^/]+")
  ref_tmp <- ref_dt_final %>% filter(sample %in% sample_id)
  
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
    mutate(feature = `Protein names`) %>%
    inner_join(., ref_tmp, by = c("Protein names" = "protein_fullname"),
               relationship = "many-to-many") %>%
    select(`Protein names`, starts_with("LFQ intensity ")) %>%
    dplyr::rename(feature = `Protein names`) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id") %>%
    filter(value != 0) %>%
    filter(!feature %in% "") %>%
    mutate_at("analysis_id", ~ paste(id_prefix, str_extract(., "(?<=LFQ intensity ).+"), sep = "")) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dda_tmp)
})
dda_pro_tables <- dda_pro_tables %>% rbindlist %>% split(., by = "lab")

## 原始信号强度水平
dda_pca_tables <- pblapply(list(dda_pre_tables, dda_pep_tables, dda_pro_tables), function(tmp_tables) {
  
  pca_results_i <- mclapply(tmp_tables, function(tmp_table) {
    
    metadata <- tmp_table %>%
      distinct(analysis_id, lab, sample, tube, injection, order) %>%
      tibble::column_to_rownames("analysis_id")
    
    exprdata_t <- tmp_table %>%
      reshape2::dcast(feature ~ analysis_id, value.var = "value", fun.aggregate = sum) %>% 
      select(all_of(rownames(metadata))) %>%
      mutate_all(~ ifelse(. == 0, NA, log2(.))) %>%
      filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
      # filter(apply(., 1, function(x) sum(is.na(x)) < 0.2 * nrow(metadata))) %>%
      # mutate_all(~ ifelse(is.na(.), median(., na.rm = TRUE), .)) %>%
      mutate_all(~ ifelse(is.na(.), 0, .)) %>%
      # na.omit %>%
      t
    
    pca_results_j <- main_pca(exprdata_t = exprdata_t, metadata = metadata,
                              center = TRUE, scale = TRUE, group = "sample",
                              biplot = FALSE, dictGroups = metadata$sample,
                              snr = TRUE, plot = FALSE)
    
    pca_results_j$feature_n <- ncol(exprdata_t)
    pca_results_j$sample_n <- nrow(exprdata_t)
    
    return(pca_results_j)
  })
  return(pca_results_i)
})
saveRDS(dda_pca_tables, "./results/tables/dda_pca_quartet_quali_filter.rds")


## 1.4 DIA: 从 report.parquet/report.pr_matrix/report.pg_matrix 中统计高置信定性特征 SNR 指标 ------
rm(list = ls())
gc()
source("./PCA.R")

ref_dt_final <- fread("./results/tables/7_qualiprop_final.csv")
ref_dt_final <- ref_dt_final %>% mutate_at("sample", ~ gsub("Quartet ", "", .))

meta_all <- fread("./data/combined_meta.csv")
all_files <- list.files("./data/multilab", recursive = TRUE, full.names = TRUE)

columns2keep <- c("Run", "Precursor.Id", "Stripped.Sequence", "Protein.Group", "Ms1.Area")
dia_pre_files <- all_files[grepl("2025_grouped_dia/(D5|D6|F7|M8).+report.parquet", all_files)]
dia_pre_tables <- pblapply(dia_pre_files, function(file_id) {
  
  sample_id <- str_extract(file_id, "(?<=grouped_dia/)[^/]+")
  ref_tmp <- ref_dt_final %>% filter(sample %in% sample_id)
  
  dia_tmp <- read_parquet(file_id, col_select = columns2keep)
  dia_tmp <- dia_tmp %>%
    dplyr::rename(analysis_id = Run, feature = Precursor.Id, value = Ms1.Area) %>%
    tidyr::separate_rows(Protein.Group, sep = ";") %>%
    inner_join(., ref_tmp, by = c("Stripped.Sequence" = "peptide_sequence",
                                  "Protein.Group" = "protein_id")) %>%
    select(feature, analysis_id, value) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dia_tmp)
})
dia_pre_tables <- dia_pre_tables %>% rbindlist %>% split(., by = "lab")

dia_pep_files <- all_files[grepl("2025_grouped_dia/(D5|D6|F7|M8).+report.pr_matrix", all_files)]
dia_pep_tables <- pblapply(dia_pep_files, function(file_id) {
  
  sample_id <- str_extract(file_id, "(?<=grouped_dia/)[^/]+")
  ref_tmp <- ref_dt_final %>% filter(sample %in% sample_id)
  
  dia_tmp <- fread(file_id)
  dia_tmp <- dia_tmp %>%
    tidyr::separate_rows(Protein.Group, sep = ";") %>%
    inner_join(., ref_tmp, by = c("Stripped.Sequence" = "peptide_sequence",
                                  "First.Protein.Description" = "protein_fullname")) %>%
    select(Stripped.Sequence, contains("cpfs01")) %>%
    dplyr::rename(feature = Stripped.Sequence) %>%
    dplyr::rename_if(is.numeric, ~ str_extract(., "[^/]+(?=\\.)")) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id", na.rm = TRUE) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dia_tmp)
})
dia_pep_tables <- dia_pep_tables %>% rbindlist %>% split(., by = "lab")

dia_pro_files <- all_files[grepl("2025_grouped_dia/(D5|D6|F7|M8).+report.pg_matrix", all_files)]
dia_pro_tables <- pblapply(dia_pro_files, function(file_id) {
  
  sample_id <- str_extract(file_id, "(?<=grouped_dia/)[^/]+")
  ref_tmp <- ref_dt_final %>% filter(sample %in% sample_id)
  
  dia_tmp <- fread(file_id)
  dia_tmp <- dia_tmp %>%
    inner_join(., ref_tmp, by = c("First.Protein.Description" = "protein_fullname"),
               relationship = "many-to-many") %>%
    select(First.Protein.Description, contains("cpfs01")) %>%
    distinct %>%
    dplyr::rename(feature = First.Protein.Description) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id", na.rm = TRUE) %>%
    mutate_at("analysis_id", ~ str_extract(., "[^/]+(?=\\.)")) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dia_tmp)
})
dia_pro_tables <- dia_pro_tables %>% rbindlist %>% split(., by = "lab")

## 原始信号强度水平
dia_pca_tables <- pblapply(list(dia_pre_tables, dia_pep_tables, dia_pro_tables), function(tmp_tables) {
  
  pca_results_i <- mclapply(tmp_tables, function(tmp_table) {
    
    metadata <- tmp_table %>%
      distinct(analysis_id, lab, sample, tube, injection, order) %>%
      tibble::column_to_rownames("analysis_id")
    
    exprdata_t <- tmp_table %>%
      reshape2::dcast(feature ~ analysis_id, value.var = "value", fun.aggregate = sum) %>% 
      select(all_of(rownames(metadata))) %>%
      mutate_all(~ ifelse(. == 0, NA, log2(.))) %>%
      filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
      # filter(apply(., 1, function(x) sum(is.na(x)) < 0.2 * nrow(metadata))) %>%
      # mutate_all(~ ifelse(is.na(.), median(., na.rm = TRUE), .)) %>%
      mutate_all(~ ifelse(is.na(.), 0, .)) %>%
      # na.omit %>%
      t
    
    pca_results_j <- main_pca(exprdata_t = exprdata_t, metadata = metadata,
                              center = TRUE, scale = TRUE, group = "sample",
                              biplot = FALSE, dictGroups = metadata$sample,
                              snr = TRUE, plot = FALSE)
    
    pca_results_j$feature_n <- ncol(exprdata_t)
    pca_results_j$sample_n <- nrow(exprdata_t)
    
    return(pca_results_j)
  })
  return(pca_results_i)
})
saveRDS(dia_pca_tables, "./results/tables/dia_pca_quartet_quali_filter.rds")



## 1.5 DDA: 从 evidence/peptides/proteinGroups 中统计高置信定量特征 SNR 指标 ------
rm(list = ls())
gc()
source("./PCA.R")

ref_dt_final <- fread("./results/tables/7_quantprop_final.csv")

meta_all <- fread("./data/combined_meta.csv")
all_files <- list.files("./data/multilab", recursive = TRUE, full.names = TRUE)

dda_pre_files <- all_files[grepl("2025_grouped_dda/(D5|F7|M8).+evidence.txt", all_files)]
dda_pre_tables <- pblapply(dda_pre_files, function(file_id) {
  
  sample_id <- str_extract(file_id, "(?<=grouped_dda/)[^/]+")
  ref_tmp <- ref_dt_final %>% filter(group %in% paste(sample_id, "/D6", sep = ""))
  
  file_path <- str_extract(file_id, ".+grouped_dda/")
  file_name <- str_extract(file_id, "(?<=(D5|F7|M8)/).+")
  dda_d6 <- fread(paste(file_path, "D6/", file_name, sep = ""), showProgress = FALSE)
  dda_d6 <- dda_d6 %>%
    dplyr::rename(analysis_id = `Raw file`) %>%
    left_join(., meta_all, by = "analysis_id") %>%
    mutate(feature = paste(Sequence, `Leading razor protein`, sep = "_")) %>%
    inner_join(ref_tmp, ., by = c("peptide_sequence" = "Sequence",
                                  "protein_id" = "Leading razor protein"),
               relationship = "many-to-many") %>%
    group_by(feature, lab, tube) %>%
    summarise(d6_intensity = mean(log2(Intensity)), .groups = "drop")
  
  dda_tmp <- fread(file_id, showProgress = FALSE)
  dda_tmp <- dda_tmp %>%
    dplyr::rename(analysis_id = `Raw file`) %>%
    left_join(., meta_all, by = "analysis_id") %>%
    mutate(feature = paste(Sequence, `Leading razor protein`, sep = "_")) %>%
    inner_join(., dda_d6, by = c("feature", "lab", "tube")) %>%
    mutate(value = 2 ^ (log2(Intensity) - d6_intensity)) %>%
    mutate(feature = paste(`Modified sequence`, Charge)) %>%
    select(feature, analysis_id, value) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dda_tmp)
})
dda_pre_tables <- dda_pre_tables %>% rbindlist %>% split(., by = "lab")

dda_pep_files <- all_files[grepl("2025_grouped_dda/(D5|F7|M8).+peptides.txt", all_files)]
dda_pep_tables <- pblapply(dda_pep_files, function(file_id) {
  
  if (grepl("timstof", file_id)) {
    sample_group <- str_extract(file_id, "(?<=2025_grouped_dda\\/).+(?=\\/timstof)")
    sample_group <- ifelse(sample_group %in% c("D5", "D6", "F7", "M8"),
                           paste("Quartet-", sample_group, sep = ""),
                           sample_group)
    sample_group <- ifelse(sample_group %in% "HeLa", "Hela", sample_group)
    id_prefix <- paste("_cpfs01_projects-HDD_cfff-e44ef5cf7aa5_HDD_cqc_21112030002_data_proteomics_fudan_university_HT-SLM-ZYT--CQC-20250324-",
                       sample_group, "-", sep = "")
    id_prefix_d6 <- paste("_cpfs01_projects-HDD_cfff-e44ef5cf7aa5_HDD_cqc_21112030002_data_proteomics_fudan_university_HT-SLM-ZYT--CQC-20250324-",
                          "Quartet-D6", "-", sep = "")
  } else {
    id_prefix <- "_cpfs01_projects-HDD_cfff-e44ef5cf7aa5_HDD_cqc_21112030002_data_proteomics_"
    id_prefix_d6 <- id_prefix
  }
  
  sample_id <- str_extract(file_id, "(?<=grouped_dda/)[^/]+")
  ref_tmp <- ref_dt_final %>% filter(group %in% paste(sample_id, "/D6", sep = ""))
  
  file_path <- str_extract(file_id, ".+grouped_dda/")
  file_name <- str_extract(file_id, "(?<=(D5|F7|M8)/).+")
  dda_d6 <- fread(paste(file_path, "D6/", file_name, sep = ""), showProgress = FALSE)
  dda_d6 <- dda_d6 %>%
    inner_join(., ref_tmp, by = c("Sequence" = "peptide_sequence",
                                  "Leading razor protein" = "protein_id")) %>%
    select(Sequence, starts_with("Intensity ")) %>%
    dplyr::rename(feature = Sequence) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id") %>%
    filter(value != 0) %>%
    mutate_at("analysis_id", ~ paste(id_prefix_d6, str_extract(., "(?<=Intensity ).+"), sep = "")) %>%
    left_join(., meta_all, by = "analysis_id") %>%
    group_by(feature, lab, tube) %>%
    summarise(d6_intensity = mean(log2(value)), .groups = "drop")
  
  dda_tmp <- fread(file_id, showProgress = FALSE)
  dda_tmp <- dda_tmp %>%
    select(Sequence, starts_with("Intensity ")) %>%
    dplyr::rename(feature = Sequence) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id", value.name = "Intensity") %>%
    filter(Intensity != 0) %>%
    mutate_at("analysis_id", ~ paste(id_prefix, str_extract(., "(?<=Intensity ).+"), sep = "")) %>%
    left_join(., meta_all, by = "analysis_id") %>%
    inner_join(., dda_d6, by = c("feature", "lab", "tube")) %>%
    mutate(value = 2 ^ (log2(Intensity) - d6_intensity)) %>%
    select(feature, analysis_id, value) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dda_tmp)
})
dda_pep_tables <- dda_pep_tables %>% rbindlist %>% split(., by = "lab")

dda_pro_files <- all_files[grepl("2025_grouped_dda/(D5|F7|M8).+proteinGroups.txt", all_files)]
dda_pro_tables <- pblapply(dda_pro_files, function(file_id) {
  
  if (grepl("timstof", file_id)) {
    sample_group <- str_extract(file_id, "(?<=2025_grouped_dda\\/).+(?=\\/timstof)")
    sample_group <- ifelse(sample_group %in% c("D5", "D6", "F7", "M8"),
                           paste("Quartet-", sample_group, sep = ""),
                           sample_group)
    sample_group <- ifelse(sample_group %in% "HeLa", "Hela", sample_group)
    id_prefix <- paste("_cpfs01_projects-HDD_cfff-e44ef5cf7aa5_HDD_cqc_21112030002_data_proteomics_fudan_university_HT-SLM-ZYT--CQC-20250324-",
                       sample_group, "-", sep = "")
    id_prefix_d6 <- paste("_cpfs01_projects-HDD_cfff-e44ef5cf7aa5_HDD_cqc_21112030002_data_proteomics_fudan_university_HT-SLM-ZYT--CQC-20250324-",
                          "Quartet-D6", "-", sep = "")
  } else {
    id_prefix <- "_cpfs01_projects-HDD_cfff-e44ef5cf7aa5_HDD_cqc_21112030002_data_proteomics_"
    id_prefix_d6 <- id_prefix
  }
  
  sample_id <- str_extract(file_id, "(?<=grouped_dda/)[^/]+")
  ref_tmp <- ref_dt_final %>% filter(group %in% paste(sample_id, "/D6", sep = ""))
  
  file_path <- str_extract(file_id, ".+grouped_dda/")
  file_name <- str_extract(file_id, "(?<=(D5|F7|M8)/).+")
  dda_d6 <- fread(paste(file_path, "D6/", file_name, sep = ""), showProgress = FALSE)
  dda_d6 <- dda_d6 %>%
    inner_join(., ref_tmp, by = c("Protein names" = "protein_fullname"),
               relationship = "many-to-many") %>%
    select(`Protein names`, starts_with("LFQ intensity ")) %>%
    dplyr::rename(feature = `Protein names`) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id") %>%
    filter(value != 0) %>%
    mutate_at("analysis_id", ~ paste(id_prefix_d6, str_extract(., "(?<=LFQ intensity ).+"), sep = "")) %>%
    left_join(., meta_all, by = "analysis_id") %>%
    group_by(feature, lab, tube) %>%
    summarise(d6_intensity = mean(log2(value)), .groups = "drop")
  
  dda_tmp <- fread(file_id, showProgress = FALSE)
  
  dda_tmp <- dda_tmp %>%
    select(`Protein names`, starts_with("LFQ intensity ")) %>%
    dplyr::rename(feature = `Protein names`) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id", value.name = "Intensity") %>%
    filter(Intensity != 0) %>%
    mutate_at("analysis_id", ~ paste(id_prefix, str_extract(., "(?<=LFQ intensity ).+"), sep = "")) %>%
    left_join(., meta_all, by = "analysis_id") %>%
    inner_join(., dda_d6, by = c("feature", "lab", "tube")) %>%
    mutate(value = 2 ^ (log2(Intensity) - d6_intensity)) %>%
    select(feature, analysis_id, value) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dda_tmp)
})
dda_pro_tables <- dda_pro_tables %>% rbindlist %>% split(., by = "lab")

## 比例定量水平
dda_pca_tables <- pblapply(list(dda_pre_tables, dda_pep_tables, dda_pro_tables), function(tmp_tables) {
  
  pca_results_i <- mclapply(tmp_tables, function(tmp_table) {
    
    metadata <- tmp_table %>%
      distinct(analysis_id, lab, sample, tube, injection, order) %>%
      tibble::column_to_rownames("analysis_id")
    
    exprdata_t <- tmp_table %>%
      reshape2::dcast(feature ~ analysis_id, value.var = "value", fun.aggregate = sum) %>% 
      select(all_of(rownames(metadata))) %>%
      mutate_all(~ ifelse(. == 0, NA, log2(.))) %>%
      filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
      # filter(apply(., 1, function(x) sum(is.na(x)) < 0.2 * nrow(metadata))) %>%
      # mutate_all(~ ifelse(is.na(.), median(., na.rm = TRUE), .)) %>%
      mutate_all(~ ifelse(is.na(.), 0, .)) %>%
      # na.omit %>%
      t
    
    pca_results_j <- main_pca(exprdata_t = exprdata_t, metadata = metadata,
                              center = TRUE, scale = TRUE, group = "sample",
                              biplot = FALSE, dictGroups = metadata$sample,
                              snr = TRUE, plot = FALSE)
    
    pca_results_j$feature_n <- ncol(exprdata_t)
    pca_results_j$sample_n <- nrow(exprdata_t)
    
    return(pca_results_j)
  })
  return(pca_results_i)
})
saveRDS(dda_pca_tables, "./results/tables/dda_pca_quartet_quant_filter.rds")


## 1.6 DIA: 从 report.parquet/report.pr_matrix/report.pg_matrix 中统计定量特征 SNR 指标 ------
rm(list = ls())
gc()
source("./PCA.R")

ref_dt_final <- fread("./results/tables/7_quantprop_final.csv")

meta_all <- fread("./data/combined_meta.csv")
all_files <- list.files("./data/multilab", recursive = TRUE, full.names = TRUE)

columns2keep <- c("Run", "Precursor.Id", "Stripped.Sequence", "Protein.Group", "Ms1.Area")
dia_pre_files <- all_files[grepl("2025_grouped_dia/(D5|F7|M8).+report.parquet", all_files)]
dia_pre_tables <- pblapply(dia_pre_files, function(file_id) {
  
  sample_id <- str_extract(file_id, "(?<=grouped_dia/)[^/]+")
  ref_tmp <- ref_dt_final %>% filter(group %in% paste(sample_id, "/D6", sep = ""))
  
  file_path <- str_extract(file_id, ".+grouped_dia/")
  file_name <- str_extract(file_id, "(?<=(D5|F7|M8)/).+")
  dia_d6 <- read_parquet(paste(file_path, "D6/", file_name, sep = ""), col_select = columns2keep)
  dia_d6 <- dia_d6 %>%
    dplyr::rename(analysis_id = Run, feature = Precursor.Id) %>%
    tidyr::separate_rows(Protein.Group, sep = ";") %>%
    inner_join(., ref_tmp, by = c("Stripped.Sequence" = "peptide_sequence",
                                  "Protein.Group" = "protein_id")) %>%
    left_join(., meta_all, by = "analysis_id") %>%
    group_by(feature, lab, tube) %>%
    summarise(d6_intensity = mean(log2(Ms1.Area)), .groups = "drop")
  
  dia_tmp <- read_parquet(file_id, col_select = columns2keep)
  dia_tmp <- dia_tmp %>%
    dplyr::rename(analysis_id = Run, feature = Precursor.Id, intensity = Ms1.Area) %>%
    tidyr::separate_rows(Protein.Group, sep = ";") %>%
    left_join(., meta_all, by = "analysis_id") %>%
    inner_join(., dia_d6, by = c("feature", "lab", "tube")) %>%
    mutate(value = 2 ^ (log2(intensity) - d6_intensity)) %>%
    select(feature, analysis_id, value) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dia_tmp)
})
dia_pre_tables <- dia_pre_tables %>% rbindlist %>% split(., by = "lab")

dia_pep_files <- all_files[grepl("2025_grouped_dia/(D5|F7|M8).+report.pr_matrix", all_files)]
dia_pep_tables <- pblapply(dia_pep_files, function(file_id) {
  
  sample_id <- str_extract(file_id, "(?<=grouped_dia/)[^/]+")
  ref_tmp <- ref_dt_final %>% filter(group %in% paste(sample_id, "/D6", sep = ""))
  
  file_path <- str_extract(file_id, ".+grouped_dia/")
  file_name <- str_extract(file_id, "(?<=(D5|F7|M8)/).+")
  dia_d6 <- fread(paste(file_path, "D6/", file_name, sep = ""))
  dia_d6 <- dia_d6 %>%
    inner_join(., ref_tmp, by = c("Stripped.Sequence" = "peptide_sequence",
                                  "First.Protein.Description" = "protein_fullname")) %>%
    select(Stripped.Sequence, contains("cpfs01")) %>%
    dplyr::rename(feature = Stripped.Sequence) %>%
    dplyr::rename_if(is.numeric, ~ str_extract(., "[^/]+(?=\\.)")) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id", na.rm = TRUE) %>%
    left_join(., meta_all, by = "analysis_id") %>%
    group_by(feature, lab, tube) %>%
    summarise(d6_intensity = mean(log2(value)), .groups = "drop")
  
  dia_tmp <- fread(file_id)
  dia_tmp <- dia_tmp %>%
    select(Stripped.Sequence, contains("cpfs01")) %>%
    dplyr::rename(feature = Stripped.Sequence) %>%
    dplyr::rename_if(is.numeric, ~ str_extract(., "[^/]+(?=\\.)")) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id", value.name = "intensity", na.rm = TRUE) %>%
    left_join(., meta_all, by = "analysis_id") %>%
    inner_join(., dia_d6, by = c("feature", "lab", "tube")) %>%
    mutate(value = 2 ^ (log2(intensity) - d6_intensity)) %>%
    select(feature, analysis_id, value) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dia_tmp)
})
dia_pep_tables <- dia_pep_tables %>% rbindlist %>% split(., by = "lab")

dia_pro_files <- all_files[grepl("2025_grouped_dia/(D5|F7|M8).+report.pg_matrix", all_files)]
dia_pro_tables <- pblapply(dia_pro_files, function(file_id) {
  
  sample_id <- str_extract(file_id, "(?<=grouped_dia/)[^/]+")
  ref_tmp <- ref_dt_final %>% filter(group %in% paste(sample_id, "/D6", sep = ""))
  
  file_path <- str_extract(file_id, ".+grouped_dia/")
  file_name <- str_extract(file_id, "(?<=(D5|F7|M8)/).+")
  dia_d6 <- fread(paste(file_path, "D6/", file_name, sep = ""))
  dia_d6 <- dia_d6 %>%
    inner_join(., ref_tmp, by = c("First.Protein.Description" = "protein_fullname"),
               relationship = "many-to-many") %>%
    select(First.Protein.Description, contains("cpfs01")) %>%
    dplyr::rename(feature = First.Protein.Description) %>%
    dplyr::rename_if(is.numeric, ~ str_extract(., "[^/]+(?=\\.)")) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id", na.rm = TRUE) %>%
    left_join(., meta_all, by = "analysis_id") %>%
    group_by(feature, lab, tube) %>%
    summarise(d6_intensity = mean(log2(value)), .groups = "drop")
  
  dia_tmp <- fread(file_id)
  dia_tmp <- dia_tmp %>%
    select(First.Protein.Description, contains("cpfs01")) %>%
    distinct %>%
    dplyr::rename(feature = First.Protein.Description) %>%
    dplyr::rename_if(is.numeric, ~ str_extract(., "[^/]+(?=\\.)")) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id", value.name = "intensity", na.rm = TRUE) %>%
    left_join(., meta_all, by = "analysis_id") %>%
    inner_join(., dia_d6, by = c("feature", "lab", "tube")) %>%
    mutate(value = 2 ^ (log2(intensity) - d6_intensity)) %>%
    select(feature, analysis_id, value) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dia_tmp)
})
dia_pro_tables <- dia_pro_tables %>% rbindlist %>% split(., by = "lab")

## 原始信号强度水平
dia_pca_tables <- pblapply(list(dia_pre_tables, dia_pep_tables, dia_pro_tables), function(tmp_tables) {
  
  pca_results_i <- mclapply(tmp_tables, function(tmp_table) {
    
    metadata <- tmp_table %>%
      distinct(analysis_id, lab, sample, tube, injection, order) %>%
      tibble::column_to_rownames("analysis_id")
    
    exprdata_t <- tmp_table %>%
      reshape2::dcast(feature ~ analysis_id, value.var = "value", fun.aggregate = sum) %>% 
      select(all_of(rownames(metadata))) %>%
      mutate_all(~ ifelse(. == 0, NA, log2(.))) %>%
      filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
      # filter(apply(., 1, function(x) sum(is.na(x)) < 0.2 * nrow(metadata))) %>%
      # mutate_all(~ ifelse(is.na(.), median(., na.rm = TRUE), .)) %>%
      mutate_all(~ ifelse(is.na(.), 0, .)) %>%
      # na.omit %>%
      t
    
    pca_results_j <- main_pca(exprdata_t = exprdata_t, metadata = metadata,
                              center = TRUE, scale = TRUE, group = "sample",
                              biplot = FALSE, dictGroups = metadata$sample,
                              snr = TRUE, plot = FALSE)
    
    pca_results_j$feature_n <- ncol(exprdata_t)
    pca_results_j$sample_n <- nrow(exprdata_t)
    
    return(pca_results_j)
  })
  return(pca_results_i)
})
saveRDS(dia_pca_tables, "./results/tables/dia_pca_quartet_quant_filter.rds")



## 2 鉴定 vs Recall: 从report.parquet/evidence中计算Recall (满足3针中至少2针检出) ------
rm(list = ls())
gc()

## 读取鉴定数目结果
all_count <- fread("./results/tables/3_count.csv")

## 计算Recall结果
ref_dt_final <- fread("./results/tables/7_qualiprop_final.csv")
ref_dt_final <- ref_dt_final %>% mutate_at("sample", ~ gsub("Quartet ", "", .))
ref_dt <- ref_dt_final %>%
  mutate(feature = paste(peptide_sequence, protein_id, sep = "_")) %>%
  distinct(sample, feature)

meta_all <- fread("./data/combined_meta.csv")
all_files <- list.files("./data/multilab", recursive = TRUE, full.names = TRUE)

dda_pep_files <- all_files[grepl("(D5|D6|F7|M8).+evidence.txt", all_files)]
dda_pep_tables <- pblapply(dda_pep_files, function(file_id) {
  
  dda_tmp <- fread(file_id, showProgress = FALSE)
  dda_tmp <- dda_tmp %>%
    dplyr::rename(analysis_id = `Raw file`) %>%
    left_join(., meta_all, by = "analysis_id")
  
  sample_id <- str_extract(file_id, "(?<=grouped_dda/)[^/]+")
  ref_tmp <- ref_dt %>% filter(sample %in% sample_id)
  
  dda_with_ref <- dda_tmp %>%
    mutate(feature = paste(Sequence, `Leading razor protein`, sep = "_")) %>%
    reshape2::dcast(., feature + lab + sample + tube ~ injection,
                    value.var = "feature", fun.aggregate = length) %>%
    filter(apply(., 1, function(x) sum(is.na(x)) <= 1)) %>%
    distinct(feature, lab, sample, tube) %>%
    reshape2::dcast(., feature ~ lab + sample + tube,
                    value.var = "feature", length, fill = 0) %>%
    left_join(ref_tmp, ., by = "feature") %>%
    select_if(is.numeric) %>%
    mutate_all(~ ifelse(is.na(.), 0, .))
  
  dda_recall <- data.frame(recall = colSums(dda_with_ref) / nrow(dda_with_ref)) %>%
    tibble::rownames_to_column("tmp_id") %>%
    tidyr::separate(tmp_id, c("lab", "sample", "tube")) %>%
    mutate_at("tube", as.integer)
  
  return(dda_recall)
})
dda_pep <- dda_pep_tables %>%
  rbindlist %>%
  na.omit %>%
  mutate(mode = "DDA")

columns2keep <- c("Run", "Precursor.Id", "Stripped.Sequence", "Protein.Group",
                  "Proteotypic", "Protein.Q.Value",
                  "Q.Value", "Global.Q.Value", "Lib.Q.Value",
                  "PG.Q.Value", "Global.PG.Q.Value", "Lib.PG.Q.Value",
                  "GG.Q.Value", "Global.PG.Q.Value")
dia_pep_files <- all_files[grepl("2025_grouped_dia/(D5|D6|F7|M8).+report.parquet", all_files)]
dia_pep_tables <- pblapply(dia_pep_files, function(file_id) {
  
  dia_tmp <- read_parquet(file_id, col_select = columns2keep)
  dia_tmp <- dia_tmp %>%
    dplyr::rename(analysis_id = Run) %>%
    left_join(., meta_all, by = "analysis_id")
  
  sample_id <- str_extract(file_id, "(?<=grouped_dia/)[^/]+")
  ref_tmp <- ref_dt %>% filter(sample %in% sample_id)
  
  dia_with_ref <- dia_tmp %>%
    tidyr::separate_rows(Protein.Group, sep = ";") %>%
    mutate(feature = paste(Stripped.Sequence, Protein.Group, sep = "_")) %>%
    reshape2::dcast(., feature + lab + sample + tube ~ injection,
                    value.var = "feature", fun.aggregate = length) %>%
    filter(apply(., 1, function(x) sum(is.na(x)) <= 1)) %>%
    distinct(feature, lab, sample, tube) %>%
    reshape2::dcast(., feature ~ lab + sample + tube,
                    value.var = "feature", length, fill = 0) %>%
    left_join(ref_tmp, ., by = "feature") %>%
    select_if(is.numeric) %>%
    mutate_all(~ ifelse(is.na(.), 0, .))
  
  dia_recall <- data.frame(recall = colSums(dia_with_ref) / nrow(dia_with_ref)) %>%
    tibble::rownames_to_column("tmp_id") %>%
    tidyr::separate(tmp_id, c("lab", "sample", "tube")) %>%
    mutate_at("tube", as.integer)
  
  return(dia_recall)
})
dia_pep <- dia_pep_tables %>%
  rbindlist %>%
  mutate(mode = "DIA")

## 合并鉴定数目+Recall
all_recall <- rbind(dda_pep, dia_pep)
all_count_recall <- inner_join(all_count, all_recall, by = c("mode", "lab", "sample", "tube"))
fwrite(all_recall, "./results/tables/8_recall.csv")
fwrite(all_count_recall, "./results/tables/all_count_recall.csv")


## 3 CV vs F1: 从report.parquet/evidence中计算F1 (满足3针中至少2针检出) ------
rm(list = ls())
gc()

## 读取第二章CV结果
dda_cv_tables <- readRDS("./results/tables/dda_cv.rds")
dia_cv_tables <- readRDS("./results/tables/dia_cv.rds")

df_cv1 <- c(dia_cv_tables[[1]], dda_cv_tables[[1]]) %>%
  rbindlist %>%
  filter(sample %in% c("D5", "D6", "F7", "M8")) %>%
  group_by(mode, lab, sample, tube) %>%
  summarise_at("cv", median) %>%
  mutate(variable = "Precursor-level") 

df_cv2 <- c(dia_cv_tables[[2]], dda_cv_tables[[2]]) %>%
  rbindlist %>%
  filter(sample %in% c("D5", "D6", "F7", "M8")) %>%
  group_by(mode, lab, sample, tube) %>%
  summarise_at("cv", median) %>%
  mutate(variable = "Peptide-level") 

df_cv3 <- c(dia_cv_tables[[3]], dda_cv_tables[[3]]) %>%
  rbindlist %>%
  filter(sample %in% c("D5", "D6", "F7", "M8")) %>%
  group_by(mode, lab, sample, tube) %>%
  summarise_at("cv", median) %>%
  mutate(variable = "Protein-level") 

df_cv_test <- df_cv1 %>%
  rbind(., df_cv2, df_cv3) %>%
  reshape2::dcast(., mode + lab + sample + tube ~ variable, value.var = "cv")

## 计算F1结果 
ref_dt_final <- fread("./results/tables/7_quantprop_final.csv")
ref_dt <- ref_dt_final %>%
  mutate(feature = paste(peptide_sequence, protein_id, sep = "_")) %>%
  mutate(log2_floor = log2(value * (1 - U)),
         log2_ceiling = log2(value * (1 + U))) %>%
  distinct(group, feature, log2_floor, log2_ceiling)

meta_all <- fread("./data/combined_meta.csv")
all_files <- list.files("./data/multilab", recursive = TRUE, full.names = TRUE)

dda_pep_files <- all_files[grepl("(D5|F7|M8).+evidence.txt", all_files)]
dda_pep_tables <- pblapply(dda_pep_files, function(file_id) {
  
  file_path <- str_extract(file_id, ".+grouped_dda/")
  file_name <- str_extract(file_id, "(?<=(D5|F7|M8)/).+")
  dda_d6 <- fread(paste(file_path, "D6/", file_name, sep = ""), showProgress = FALSE)
  dda_d6 <- dda_d6 %>%
    dplyr::rename(analysis_id = `Raw file`) %>%
    left_join(., meta_all, by = "analysis_id") %>%
    mutate(feature = paste(Sequence, `Leading razor protein`, sep = "_")) %>%
    filter(feature %in% ref_dt$feature) %>%
    group_by(feature, lab, tube) %>%
    summarise(d6_intensity = mean(log2(Intensity)), .groups = "drop")
  
  sample_id <- str_extract(file_id, "(?<=grouped_dda/)[^/]+")
  ref_tmp <- ref_dt %>% filter(group %in% paste(sample_id, "/D6", sep = ""))
  
  dda_tmp <- fread(file_id, showProgress = FALSE)
  dda_tmp <- dda_tmp %>%
    dplyr::rename(analysis_id = `Raw file`) %>%
    left_join(., meta_all, by = "analysis_id") %>%
    mutate(feature = paste(Sequence, `Leading razor protein`, sep = "_")) %>%
    inner_join(., dda_d6, by = c("feature", "lab", "tube")) %>%
    mutate(log2_test = log2(Intensity) - d6_intensity) %>%
    reshape2::dcast(., feature ~ lab + tube, fun.aggregate = mean,
                    value.var = "log2_test")
  
  dda_with_ref <- dda_tmp %>%
    left_join(ref_tmp, ., by = "feature") %>%
    reshape2::melt(., id = 1:4, value.name = "log2_test") %>%
    tidyr::separate(variable, c("lab", "tube"), sep = "_") %>%
    mutate(label = ifelse(log2_test > log2_floor & log2_test < log2_ceiling, "TP", "FP")) %>%
    mutate_at("label", ~ ifelse(is.na(log2_test), "FN", .))
  
  return(dda_with_ref)
})
dda_pep <- dda_pep_tables %>%
  rbindlist %>%
  mutate(mode = "DDA")

columns2keep <- c("Run", "Stripped.Sequence", "Protein.Group", "Precursor.Normalised")
dia_pep_files <- all_files[grepl("2025_grouped_dia/(D5|F7|M8).+report.parquet", all_files)]
dia_pep_tables <- pblapply(dia_pep_files, function(file_id) {
  
  file_path <- str_extract(file_id, ".+grouped_dia/")
  file_name <- str_extract(file_id, "(?<=(D5|F7|M8)/).+")
  dia_d6 <- read_parquet(paste(file_path, "D6/", file_name, sep = ""), col_select = columns2keep)
  dia_d6 <- dia_d6 %>%
    dplyr::rename(analysis_id = Run) %>%
    left_join(., meta_all, by = "analysis_id") %>%
    tidyr::separate_rows(Protein.Group, sep = ";") %>%
    mutate(feature = paste(Stripped.Sequence, Protein.Group, sep = "_")) %>%
    filter(feature %in% ref_dt$feature) %>%
    group_by(feature, lab, tube) %>%
    summarise(d6_intensity = mean(log2(Precursor.Normalised)), .groups = "drop")
  
  sample_id <- str_extract(file_id, "(?<=grouped_dia/)[^/]+")
  ref_tmp <- ref_dt %>% filter(group %in% paste(sample_id, "/D6", sep = ""))
  
  dia_tmp <- read_parquet(file_id, col_select = columns2keep)
  dia_tmp <- dia_tmp %>%
    dplyr::rename(analysis_id = Run) %>%
    left_join(., meta_all, by = "analysis_id") %>%
    tidyr::separate_rows(Protein.Group, sep = ";") %>%
    mutate(feature = paste(Stripped.Sequence, Protein.Group, sep = "_")) %>%
    inner_join(., dia_d6, by = c("feature", "lab", "tube")) %>%
    mutate(log2_test = log2(Precursor.Normalised) - d6_intensity) %>%
    reshape2::dcast(., feature ~ lab + tube, fun.aggregate = mean,
                    value.var = "log2_test")
  
  dia_with_ref <- dia_tmp %>%
    left_join(ref_tmp, ., by = "feature") %>%
    reshape2::melt(., id = 1:4, value.name = "log2_test") %>%
    tidyr::separate(variable, c("lab", "tube"), sep = "_") %>%
    mutate(label = ifelse(log2_test > log2_floor & log2_test < log2_ceiling, "TP", "FP")) %>%
    mutate_at("label", ~ ifelse(is.na(log2_test), "FN", .))
  
  return(dia_with_ref)
})
dia_pep <- dia_pep_tables %>%
  rbindlist %>%
  mutate(mode = "DIA")

## 合并CV+F1
all_with_ref <- rbind(dda_pep, dia_pep)
all_f1 <- all_with_ref %>%
  mutate(sample = gsub("/D6", "", group)) %>%
  reshape2::dcast(., mode + lab + sample + tube ~ label, value.var = "label",
                  fun.aggregate = length) %>%
  mutate_at("tube", as.integer) %>%
  mutate_if(is.integer, as.numeric) %>%
  mutate(recall = TP / (TP + FN),
         precision = TP / (TP + FP),
         csi = TP / (TP + FP + FN),
         fdr = FP/(FP+TP),
         f1 = 2 * precision * recall / (precision + recall))
fwrite(all_f1, "./results/tables/8_f1.csv")

all_cv_f1 <- left_join(df_cv_test, all_f1, by = c("mode", "lab", "sample", "tube"))
fwrite(all_cv_f1, "./results/tables/all_cv_f1.csv")


## 4.1: 提取高置信定性特征所有指标 ------
ref_dt_final <- fread("./results/tables/7_qualiprop_final.csv")
ref_dt_final <- ref_dt_final %>% mutate_at("sample", ~ gsub("Quartet ", "", .))

meta_all <- fread("./data/combined_meta.csv")
all_files <- list.files("./data/multilab", recursive = TRUE, full.names = TRUE)

dda_pre_files <- all_files[grepl("2025_grouped_dda/(D5|D6|F7|M8).+evidence.txt", all_files)]
dda_pre_tables <- pblapply(dda_pre_files, function(file_id) {
  
  sample_id <- str_extract(file_id, "(?<=grouped_dda/)[^/]+")
  ref_tmp <- ref_dt_final %>% filter(sample %in% sample_id)
  
  dda_tmp <- fread(file_id, showProgress = FALSE)
  dda_tmp <- dda_tmp %>%
    dplyr::rename(analysis_id = `Raw file`,
                  peptide_sequence = Sequence,
                  protein_id = `Leading razor protein`) %>%
    inner_join(., ref_tmp, by = c("peptide_sequence", "protein_id")) %>%
    mutate(FWHM = `Retention length` / 1.7,
           TIC = Intensity,
           `MS1 Accuracy` = `Mass error [ppm]`,
           `Delta RT` = `Calibrated retention time` - `Retention time`,
           `Delta IM` = NA)
  
  if (grepl("timstof", file_id)) {
    dda_tmp <- dda_tmp %>% mutate(`Delta IM` = `1/K0` - `Calibrated 1/K0`)
  }
  
  dda_final <- dda_tmp %>%
    select(peptide_sequence, protein_id, analysis_id,
           FWHM, TIC, Charge, `MS1 Accuracy`, `Delta RT`, `Delta IM`) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dda_final)
})
dda_pre_tables <- dda_pre_tables %>% rbindlist %>% split(., by = "lab")

columns2keep <- c("Run", "Stripped.Sequence", "Protein.Group",
                  "Precursor.Charge", "Ms1.Area",
                  "Precursor.Mz", "Ms1.Apex.Mz.Delta",
                  "FWHM", "RT", "Predicted.RT",
                  "IM", "Predicted.IM")
dia_pre_files <- all_files[grepl("2025_grouped_dia/(D5|D6|F7|M8).+report.parquet", all_files)]
dia_pre_tables <- pblapply(dia_pre_files, function(file_id) {
  
  sample_id <- str_extract(file_id, "(?<=grouped_dia/)[^/]+")
  ref_tmp <- ref_dt_final %>% filter(sample %in% sample_id)
  
  dia_tmp <- read_parquet(file_id, col_select = columns2keep)
  dia_tmp <- dia_tmp %>%
    dplyr::rename(analysis_id = Run,
                  peptide_sequence = Stripped.Sequence,
                  protein_id = Protein.Group) %>%
    tidyr::separate_rows(protein_id, sep = ";") %>%
    inner_join(., ref_tmp, by = c("peptide_sequence", "protein_id")) %>%
    mutate(TIC = Ms1.Area,
           `MS1 Accuracy` = Ms1.Apex.Mz.Delta*10^6/Precursor.Mz,
           `Delta RT` = RT - Predicted.RT,
           `Delta IM` = IM - Predicted.IM,
           Charge = Precursor.Charge) %>%
    select(peptide_sequence, protein_id, analysis_id,
           FWHM, TIC, Charge, `MS1 Accuracy`, `Delta RT`, `Delta IM`) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dia_tmp)
})
dia_pre_tables <- dia_pre_tables %>% rbindlist %>% split(., by = "lab")

all_pre_tables <- c(dda_pre_tables, dia_pre_tables)
saveRDS(all_pre_tables, "./results/tables/all_lc_ms_metrics_quali_filter.rds")


## 4.2: 提取全特征所有指标 ------
meta_all <- fread("./data/combined_meta.csv")
all_files <- list.files("./data/multilab", recursive = TRUE, full.names = TRUE)

dda_pre_files <- all_files[grepl("2025_grouped_dda/(D5|D6|F7|M8).+evidence.txt", all_files)]
dda_pre_tables <- pblapply(dda_pre_files, function(file_id) {
  
  dda_tmp <- fread(file_id, showProgress = FALSE)
  dda_tmp <- dda_tmp %>%
    dplyr::rename(analysis_id = `Raw file`,
                  peptide_sequence = Sequence,
                  protein_id = `Leading razor protein`) %>%
    mutate(FWHM = `Retention length` / 1.7,
           TIC = Intensity,
           `MS1 Accuracy` = `Mass error [ppm]`,
           `Delta RT` = `Calibrated retention time` - `Retention time`,
           `Delta IM` = NA)
  
  if (grepl("timstof", file_id)) {
    dda_tmp <- dda_tmp %>% mutate(`Delta IM` = `1/K0` - `Calibrated 1/K0`)
  }
  
  dda_final <- dda_tmp %>%
    select(peptide_sequence, protein_id, analysis_id,
           FWHM, TIC, Charge, `MS1 Accuracy`, `Delta RT`, `Delta IM`) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dda_final)
})
dda_pre_tables <- dda_pre_tables %>% rbindlist %>% split(., by = "lab")

columns2keep <- c("Run", "Stripped.Sequence", "Protein.Group",
                  "Precursor.Charge", "Ms1.Area",
                  "Precursor.Mz", "Ms1.Apex.Mz.Delta",
                  "FWHM", "RT", "Predicted.RT",
                  "IM", "Predicted.IM")
dia_pre_files <- all_files[grepl("2025_grouped_dia/(D5|D6|F7|M8).+report.parquet", all_files)]
dia_pre_tables <- pblapply(dia_pre_files, function(file_id) {
  
  dia_tmp <- read_parquet(file_id, col_select = columns2keep)
  dia_tmp <- dia_tmp %>%
    dplyr::rename(analysis_id = Run,
                  peptide_sequence = Stripped.Sequence,
                  protein_id = Protein.Group) %>%
    tidyr::separate_rows(protein_id, sep = ";") %>%
    mutate(TIC = Ms1.Area,
           `MS1 Accuracy` = Ms1.Apex.Mz.Delta*10^6/Precursor.Mz,
           `Delta RT` = RT - Predicted.RT,
           `Delta IM` = IM - Predicted.IM,
           Charge = Precursor.Charge) %>%
    select(peptide_sequence, protein_id, analysis_id,
           FWHM, TIC, Charge, `MS1 Accuracy`, `Delta RT`, `Delta IM`) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dia_tmp)
})
dia_pre_tables <- dia_pre_tables %>% rbindlist %>% split(., by = "lab")

all_pre_tables <- c(dda_pre_tables, dia_pre_tables)
saveRDS(all_pre_tables, "./results/tables/all_lc_ms_metrics.rds")


## 4.3: 提取全特征原始质量精度中位值 ------
meta_all <- fread("./data/combined_meta.csv")
all_files <- list.files("./data/multilab", recursive = TRUE, full.names = TRUE)

dda_pre_files1 <- all_files[grepl("2025_grouped_dda/(D5|D6|F7|M8).+evidence.txt", all_files)]
dda_pre_tables1 <- pblapply(dda_pre_files1, function(file_id) {
  
  dda_tmp <- fread(file_id, showProgress = FALSE)
  dda_tmp <- dda_tmp %>% dplyr::rename(analysis_id = `Raw file`)
  
  dda_final <- dda_tmp %>%
    group_by(analysis_id) %>%
    summarise(value = median(abs(`Uncalibrated mass error [ppm]`), na.rm = TRUE)) %>%
    mutate(variable = "MS1 Accuracy") %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dda_final)
})

dda_pre_files2 <- all_files[grepl("2025_grouped_dda/(D5|D6|F7|M8).+msms.txt", all_files)]
dda_pre_tables2 <- pblapply(dda_pre_files2, function(file_id) {
  
  dda_tmp <- fread(file_id, showProgress = FALSE)
  dda_tmp <- dda_tmp %>%
    dplyr::rename(analysis_id = `Raw file`,
                  peptide_sequence = Sequence) %>%
    filter(!is.na(`Mass deviations [ppm]`) & `Mass deviations [ppm]` != "")
  
  dda_tmp[, `MS2 Mass deviations [ppm]` := sapply(`Mass deviations [ppm]`, function(x) {
    vals <- as.numeric(unlist(strsplit(x, ";", fixed = TRUE)))
    return(median(abs(vals), na.rm = TRUE))
  })]
  
  dda_final <- dda_tmp %>%
    group_by(analysis_id) %>%
    summarise(value = median(abs(`MS2 Mass deviations [ppm]`), na.rm = TRUE)) %>%
    mutate(variable = "MS2 Accuracy") %>%
    left_join(., meta_all, by = "analysis_id")
  
  rm(dda_tmp)
  gc()
  
  return(dda_final)
})

dda_pre <- c(dda_pre_tables1, dda_pre_tables2) %>%
  rbindlist %>%
  reshape2::dcast(., analysis_id ~ variable, value.var = "value") %>%
  left_join(., meta_all, by = "analysis_id")

dia_pre_files <- all_files[grepl("2025_grouped_dia/(D5|D6|F7|M8).+report.stats", all_files)]
dia_pre_tables <- pblapply(dia_pre_files, function(file_id) {
  
  dia_tmp <- fread(file_id)
  
  dia_tmp <- dia_tmp %>%
    mutate(analysis_id = str_extract(File.Name, "[^/]+(?=\\.)")) %>%
    filter(!grepl("fudan_university", analysis_id))  %>%
    dplyr::rename(`MS1 Accuracy` = Median.Mass.Acc.MS1,
                  `MS2 Accuracy` = Median.Mass.Acc.MS2) %>%
    select(analysis_id, `MS1 Accuracy`, `MS2 Accuracy`) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dia_tmp)
})

all_pre <- rbindlist(c(list(dda_pre), dia_pre_tables))
fwrite(all_pre, "./results/tables/all_ms_accuracy.csv")


## 5.1: 作图: FWHM动态变化 ------
rm(list = ls())
gc()

all_files <- list.files("./results/tables", recursive = TRUE, full.names = TRUE)
all_lcms_files <- all_files[grepl("lc_ms", all_files)]

labels.source <- c("High-Confidence")
names(labels.source) <- c("_quali_filter")
lcms_tables <- pblapply(all_lcms_files, function(file_id) {
  
  source_id <- labels.source[str_extract(file_id, "(?<=lc_ms_metrics).*(?=\\.rds)")]
  source_id <- ifelse(is.na(source_id), "All", source_id)
  
  lcms_tables_i <- readRDS(file_id)
  qc_tables_i <- pblapply(lcms_tables_i, function(tmp_table) {
    
    qc_df_i <- tmp_table %>%
      group_by(lab, sample, tube, injection, order) %>%
      summarise_at("FWHM", median, na.rm = TRUE)
    
    return(qc_df_i)
  })
  lcms_df <- rbindlist(qc_tables_i) %>% mutate(source = source_id)
  
  return(lcms_df)
})

df_lcms <- lcms_tables %>% rbindlist %>% arrange(lab, order)
rm(list = setdiff(ls(), c("df_lcms")))
gc()

labels.order <- 1:36
names(labels.order) <- as.character(unique(df_lcms$order))

df_tmp <- df_lcms %>%
  mutate(order_code = labels.order[as.character(order)]) %>%
  arrange(source, desc(FWHM))

colors.lab <- brewer.pal(10, "PiYG")
names(colors.lab) <- unique(df_tmp$lab)

p <- ggplot(df_tmp, aes(y = FWHM)) +
  stat_summary(aes(x = order_code, color = lab), fun = median, geom = "line", linewidth = 1) +
  stat_summary(aes(x = order_code, fill = lab), fun = median, geom = "point",
               size = 5, shape = 21, color = "black") +
  scale_y_continuous(n.breaks = 10, name = "FWHM (min)") +
  scale_x_continuous(breaks = 1:36, name = "Injection order") +
  scale_fill_manual(values = colors.lab, name = "Lab") +
  scale_color_manual(values = colors.lab, name = "Lab") +
  facet_wrap(~ source, ncol = 1) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm"))

ggsave("./chapter4_4_7.pdf", p, height = 10, width = 14)


## 5.2: 作图: FWHM vs 流速 ------
rm(list = setdiff(ls(), c("df_lcms")))
gc()

labels.lc <- c("EASY-nLC 1200", "EASY-nLC 1200", "EASY-nLC 1200", "nanoElute",
               "EASY-nLC 1200", "Vanquish Neo", "Vanquish Neo", "Vanquish Neo",
               "nanoElute", "Waters M class")
names(labels.lc) <- c("ZJU", "QLB", "NCP", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")

labels.flowrate <- c(450, 300, 600, 300, 600, 600, 800, 300, 600, 350)
names(labels.flowrate) <- c("ZJU", "QLB", "NCP", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")

labels.class <- c("Mid", "Low", "High", "Low", "High", "High", "High", "Low", "High", "Mid")
names(labels.class) <- c("ZJU", "QLB", "NCP", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")

df_tmp <- df_lcms %>%
  mutate(lc = labels.lc[lab]) %>%
  mutate(flowrate = labels.flowrate[lab]) %>%
  mutate_at("flowrate", as.character) %>%
  mutate(class = labels.class[lab]) %>%
  mutate_at("class", ~ factor(., levels = c("High", "Mid", "Low")))

fwhm_thres <- df_tmp %>%
  group_by(lc, source) %>%
  summarise(rt_median = median(FWHM), rt_max = max(FWHM),
            rt_mean = mean(FWHM), rt_sd = sd(FWHM),
            .groups = "drop") %>%
  arrange(desc(rt_median), lc)

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

colors.lc <- brewer.pal(4, "PiYG")
names(colors.lc) <- unique(fwhm_thres$lc)

p <- ggplot(df_tmp, aes(x = lc, y = FWHM)) +
  # geom_hline(aes(yintercept = rt_median), fwhm_thres, lty = 2) +
  # stat_summary(aes(group = source), fun = mean, geom = "line", color = "#999") +
  # stat_summary(aes(group = source), fun.data = "mean_se", geom = "errorbar", width = .2) +
  # stat_summary(aes(fill = source), fun = mean, geom = "point", color = "black", shape = 21, size = 5) +
  stat_summary(fun.min = min, fun.max = max, geom = "errorbar", width = .5) +
  geom_boxplot(aes(fill = lc), outliers = FALSE, width = .8) +
  geom_text(aes(x = lc, y = rt_median,
                label = sprintf("Median = %.2f", rt_median)),
            vjust = -.5, color = "black",
            data = fwhm_thres, size = 4.5) +
  ggpubr::geom_signif(comparisons = list(c("EASY-nLC 1200", "nanoElute"),
                                         c("EASY-nLC 1200", "Vanquish Neo"),
                                         c("EASY-nLC 1200", "Waters M class")),
                      map_signif_level = p_format,
                      vjust = -.5,
                      y_position = .3,
                      step_increase = .2,
                      tip_length = .05,
                      textsize = 6,
                      test = "t.test") +
  scale_y_continuous(n.breaks = 10,
                     expand = expansion(mult = c(.05, 0.15)),
                     name = "FWHM (min)") +
  scale_fill_manual(values = colors.lc)+
  # scale_color_manual(values = c("black", "white")) +
  facet_wrap(~ source, ncol = 2, scales = "free_x") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm"))

ggsave("./chapter4_4_8.pdf", p, height = 6.5, width = 14)


## 6.1: 作图: Delta RT动态变化 ------
rm(list = ls())
gc()

all_files <- list.files("./results/tables", recursive = TRUE, full.names = TRUE)
all_lcms_files <- all_files[grepl("lc_ms", all_files)]

labels.source <- c("High-Confidence")
names(labels.source) <- c("_quali_filter")
lcms_tables <- pblapply(all_lcms_files, function(file_id) {
  
  source_id <- labels.source[str_extract(file_id, "(?<=lc_ms_metrics).*(?=\\.rds)")]
  source_id <- ifelse(is.na(source_id), "All", source_id)
  
  lcms_tables_i <- readRDS(file_id)
  qc_tables_i <- pblapply(lcms_tables_i, function(tmp_table) {
    
    qc_df_i <- tmp_table %>%
      group_by(lab, sample, tube, injection, order) %>%
      summarise_at("Delta RT", ~ median(abs(.), na.rm = TRUE))
    
    return(qc_df_i)
  })
  lcms_df <- rbindlist(qc_tables_i) %>% mutate(source = source_id)
  
  return(lcms_df)
})

df_lcms <- lcms_tables %>% rbindlist %>% arrange(lab, order)
rm(list = setdiff(ls(), c("df_lcms")))
gc()

labels.order <- 1:36
names(labels.order) <- as.character(unique(df_lcms$order))

df_tmp <- df_lcms %>%
  mutate(order_code = labels.order[as.character(order)]) %>%
  arrange(source, desc(`Delta RT`))

colors.lab <- brewer.pal(10, "PiYG")
names(colors.lab) <- unique(df_tmp$lab)

p <- ggplot(df_tmp, aes(y = `Delta RT`)) +
  stat_summary(aes(x = order_code, color = lab), fun = median, geom = "line", linewidth = 1) +
  stat_summary(aes(x = order_code, fill = lab), fun = median, geom = "point",
               size = 5, shape = 21, color = "black") +
  scale_y_continuous(n.breaks = 10, name = "Delta RT (min)") +
  scale_x_continuous(breaks = 1:36, name = "Injection order") +
  scale_fill_manual(values = colors.lab, name = "Lab") +
  scale_color_manual(values = colors.lab, name = "Lab") +
  facet_wrap(~ source, ncol = 1) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm"))

ggsave("./chapter4_4_9.pdf", p, height = 10, width = 14)


## 6.2: 作图: Delta RT vs 色谱仪型号 ------
rm(list = setdiff(ls(), c("df_lcms")))
gc()

labels.lc <- c("EASY-nLC 1200", "EASY-nLC 1200", "EASY-nLC 1200", "nanoElute",
               "EASY-nLC 1200", "Vanquish Neo", "Vanquish Neo", "Vanquish Neo",
               "nanoElute", "Waters M class")
names(labels.lc) <- c("ZJU", "QLB", "NCP", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")

df_tmp <- df_lcms %>%
  mutate(class = labels.lc[lab])

rt_thres <- df_tmp %>%
  group_by(class, source) %>%
  summarise(rt_median = median(`Delta RT`), rt_max = max(`Delta RT`),
            rt_mean = mean(`Delta RT`), rt_sd = sd(`Delta RT`),
            .groups = "drop") %>%
  arrange(desc(rt_median), lc)

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

colors.lc <- brewer.pal(4, "PiYG")
names(colors.lc) <- unique(rt_thres$lc)

p <- ggplot(df_tmp, aes(x = class, y = `Delta RT`)) +
  stat_summary(fun.min = min, fun.max = max, geom = "errorbar", width = .5) +
  geom_boxplot(aes(fill = class), outliers = FALSE, width = .8) +
  geom_text(aes(x = class, y = rt_median,
                label = sprintf("Median = %.2f", rt_median)),
            vjust = -.5, color = "black",
            data = rt_thres, size = 4.5) +
  ggpubr::geom_signif(comparisons = list(c("EASY-nLC 1200", "nanoElute"),
                                         c("EASY-nLC 1200", "Vanquish Neo"),
                                         c("EASY-nLC 1200", "Waters M class")),
                      map_signif_level = p_format,
                      vjust = -.5,
                      y_position = 20,
                      step_increase = .2,
                      tip_length = .05,
                      textsize = 6,
                      test = "t.test") +
  scale_y_continuous(n.breaks = 10,
                     expand = expansion(mult = c(.05, 0.15)),
                     name = "Delta RT (min)") +
  scale_fill_manual(values = colors.lc)+
  facet_wrap(~ source, ncol = 2) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm"))

ggsave("./chapter4_4_10.pdf", p, height = 6.5, width = 14)


## 7.1: 作图: 电荷分布动态变化 ------
rm(list = ls())
gc()

all_files <- list.files("./results/tables", recursive = TRUE, full.names = TRUE)
all_lcms_files <- all_files[grepl("lc_ms", all_files)]

labels.source <- c("High-Confidence")
names(labels.source) <- c("_quali_filter")
lcms_tables <- pblapply(all_lcms_files, function(file_id) {
  
  source_id <- labels.source[str_extract(file_id, "(?<=lc_ms_metrics).*(?=\\.rds)")]
  source_id <- ifelse(is.na(source_id), "All", source_id)
  
  lcms_tables_i <- readRDS(file_id)
  charge_tables_i <- pblapply(lcms_tables_i, function(tmp_table) {
    
    charge_df_i <- tmp_table %>%
      reshape2::dcast(., lab + sample + tube + injection + order ~ Charge,
                      value.var = "peptide_sequence",
                      fun.aggregate = length)
    
    prop_df_i <- charge_df_i %>%
      mutate(Total = rowSums(pick(matches("^\\d+$")), na.rm = TRUE)) %>%
      mutate(Charge_1_prop = `1` / Total,
             Charge_2_prop = `2` / Total,
             Charge_3_prop = `3` / Total,
             Charge_4_prop = 4 / Total) %>%
      select(lab, sample, tube, injection, order,
             Charge_1_prop, Charge_2_prop, Charge_3_prop, Charge_4_prop)
    
    return(prop_df_i)
  })
  lcms_df <- rbindlist(charge_tables_i) %>% mutate(source = source_id)
  
  return(lcms_df)
})

df_lcms <- lcms_tables %>% rbindlist %>% arrange(lab, order)
rm(list = setdiff(ls(), c("df_lcms")))
gc()

labels.order <- 1:36
names(labels.order) <- as.character(unique(df_lcms$order))

df_tmp <- df_lcms %>%
  mutate(order_code = labels.order[as.character(order)]) %>%
  arrange(source, Charge_2_prop + Charge_3_prop)

colors.lab <- brewer.pal(10, "PiYG")
names(colors.lab) <- unique(df_tmp$lab)

p <- ggplot(df_tmp, aes(y = Charge_2_prop + Charge_3_prop)) +
  stat_summary(aes(x = order_code, color = lab), fun = median, geom = "line", linewidth = 1) +
  stat_summary(aes(x = order_code, fill = lab), fun = median, geom = "point",
               size = 5, shape = 21, color = "black") +
  scale_y_continuous(n.breaks = 10,
                     name = "Proportion of 2+,3+ Precursors",
                     labels = scales::percent) +
  scale_x_continuous(breaks = 1:36, name = "Injection order") +
  scale_fill_manual(values = colors.lab, name = "Lab") +
  scale_color_manual(values = colors.lab, name = "Lab") +
  facet_wrap(~ source, ncol = 1) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm"))

ggsave("./chapter4_4_11.pdf", p, height = 10, width = 14)


## 7.2: 作图: 电荷分布 vs 质谱仪型号 ------
rm(list = setdiff(ls(), c("df_lcms")))
gc()

labels.ms <- c("Orbitrap Exploris 480", "Q Extractive HF-X", "Orbitrap Fusion",
               "timsTOF HT", "Orbitrap Fusion Lumos",	"Orbitrap Exploris 480",
               "Orbitrap Astral",	"timsTOF HT",	"timsTOF Pro2",	"ZenoTOF 7600")
names(labels.ms) <- c("ZJU", "QLB", "NCP", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")

labels.class <- c("Orbitrap", "Orbitrap", "Orbitrap", "timsTOF", "Orbitrap",
                  "Orbitrap", "Orbitrap",	"timsTOF",	"timsTOF",	"ZenoTOF")
names(labels.class) <- c("ZJU", "QLB", "NCP", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")

df_tmp <- df_lcms %>%
  mutate(ms = labels.ms[lab]) %>%
  mutate(class = labels.class[lab])

rt_thres <- df_tmp %>%
  group_by(class, source) %>%
  summarise(rt_median = median(Charge_2_prop + Charge_3_prop),
            rt_max = max(Charge_2_prop + Charge_3_prop),
            rt_mean = mean(Charge_2_prop + Charge_3_prop),
            rt_sd = sd(Charge_2_prop + Charge_3_prop),
            .groups = "drop") %>%
  arrange(source, rt_median, class)

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

colors.ms <- brewer.pal(11, "PiYG")[c(2, 4, 10)]
names(colors.ms) <- unique(rt_thres$class)

df_tmp$class <- factor(df_tmp$class, levels = unique(rt_thres$class))

p <- ggplot(df_tmp, aes(x = class, y = Charge_2_prop + Charge_3_prop)) +
  stat_summary(fun.min = min, fun.max = max, geom = "errorbar", width = .5) +
  geom_boxplot(aes(fill = class), outliers = FALSE, width = .9) +
  geom_text(aes(x = class, y = rt_median,
                label = sprintf("Median = %.2f", rt_median)),
            vjust = -.5, color = "black",
            data = rt_thres, size = 4.5) +
  ggpubr::geom_signif(comparisons = list(c("Orbitrap", "ZenoTOF"),
                                         c("timsTOF", "ZenoTOF"),
                                         c("Orbitrap", "timsTOF")),
                      map_signif_level = p_format,
                      vjust = -.5,
                      # y_position = 1.01,
                      step_increase = .2,
                      tip_length = .05,
                      textsize = 6,
                      test = "t.test") +
  scale_y_continuous(breaks = seq(.9, 1, .02),
                     expand = expansion(mult = c(.05, 0.15)),
                     name = "Proportion",
                     labels = scales::percent) +
  scale_fill_manual(values = colors.ms)+
  facet_wrap(~ source, ncol = 2) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 20),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm"))

ggsave("./chapter4_4_12.pdf", p, height = 5, width = 14)


## 8.1: 作图: TIC 动态变化 ------
rm(list = ls())
gc()

all_files <- list.files("./results/tables", recursive = TRUE, full.names = TRUE)
all_lcms_files <- all_files[grepl("lc_ms", all_files)]

labels.source <- c("High-Confidence")
names(labels.source) <- c("_quali_filter")
lcms_tables <- pblapply(all_lcms_files, function(file_id) {
  
  source_id <- labels.source[str_extract(file_id, "(?<=lc_ms_metrics).*(?=\\.rds)")]
  source_id <- ifelse(is.na(source_id), "All", source_id)
  
  lcms_tables_i <- readRDS(file_id)
  qc_tables_i <- pblapply(lcms_tables_i, function(tmp_table) {
    
    qc_df_i <- tmp_table %>%
      group_by(lab, sample, tube, injection, order) %>%
      summarise_at("TIC", median, na.rm = TRUE)
    
    return(qc_df_i)
  })
  lcms_df <- rbindlist(qc_tables_i) %>% mutate(source = source_id)
  
  return(lcms_df)
})

df_lcms <- lcms_tables %>% rbindlist %>% arrange(lab, order)
rm(list = setdiff(ls(), c("df_lcms")))
gc()

labels.order <- 1:36
names(labels.order) <- as.character(unique(df_lcms$order))

df_tmp <- df_lcms %>%
  mutate(order_code = labels.order[as.character(order)])

rt_thres <- df_tmp %>%
  group_by(lab, source) %>%
  summarise(rt_median = median(TIC),
            rt_max = max(TIC),
            rt_mean = mean(TIC),
            rt_sd = sd(TIC),
            .groups = "drop") %>%
  arrange(source, desc(rt_sd))

colors.lab <- brewer.pal(10, "PiYG")
names(colors.lab) <- unique(rt_thres$lab)

p <- ggplot(df_tmp, aes(y = TIC)) +
  stat_summary(aes(x = order_code, color = lab), fun = median, geom = "line", linewidth = 1) +
  stat_summary(aes(x = order_code, fill = lab), fun = median, geom = "point",
               size = 5, shape = 21, color = "black") +
  scale_y_continuous(n.breaks = 10, name = "TIC") +
  scale_x_continuous(breaks = 1:36, name = "Injection order") +
  scale_fill_manual(values = colors.lab, name = "Lab") +
  scale_color_manual(values = colors.lab, name = "Lab") +
  facet_wrap(~ source, ncol = 1) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm"))

ggsave("./chapter4_4_13.pdf", p, height = 10, width = 14)


## 8.2: 作图: TIC vs 质谱仪型号 ------
rm(list = setdiff(ls(), c("df_lcms")))
gc()

labels.ms <- c("Orbitrap Exploris 480", "Q Extractive HF-X", "Orbitrap Fusion",
               "timsTOF HT", "Orbitrap Fusion Lumos",	"Orbitrap Exploris 480",
               "Orbitrap Astral",	"timsTOF HT",	"timsTOF Pro2",	"ZenoTOF 7600")
names(labels.ms) <- c("ZJU", "QLB", "NCP", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")

labels.volume <- c(1, 1, .5, 1, 1, .5, .5, .4, .214, .1)
names(labels.volume) <- c("ZJU", "QLB", "NCP", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")

labels.conc <- c(500, 500, 250, 500, 100, 500, 500, 500, 71.3, 50)
names(labels.volume) <- c("ZJU", "QLB", "NCP", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")

labels.class1 <- c("High", "High", "Low", "High", "Low",
                   "High", "High", "High", "Low", "Low")
names(labels.class1) <- c("ZJU", "QLB", "NCP", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")

labels.class2 <- c("Orbitrap", "Orbitrap", "Orbitrap", "TOF", "Orbitrap",
                   "Orbitrap", "Orbitrap",	"TOF",	"TOF",	"TOF")
names(labels.class2) <- c("ZJU", "QLB", "NCP", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")

df_tmp <- df_lcms %>%
  group_by(source, lab, sample, tube) %>%
  summarise_at("TIC", ~ sd(log10(.))) %>%
  mutate(class1 = labels.class1[lab]) %>%
  mutate_at("class1", ~ factor(., levels = c("High", "Mid", "Low"))) %>%
  mutate(class2 = labels.class2[lab])

rt_thres <- df_tmp %>%
  group_by(class1, class2, source) %>%
  summarise(rt_median = median(TIC),
            rt_max = max(TIC),
            rt_mean = mean(TIC),
            rt_sd = sd(TIC),
            .groups = "drop") %>%
  arrange(source, desc(rt_median))

colors.ms <- brewer.pal(11, "PiYG")[c(2, 10)]
names(colors.ms) <- unique(rt_thres$class1)

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

p <- ggplot(df_tmp, aes(x = class1, y = TIC)) +
  stat_summary(fun.min = min, fun.max = max, geom = "errorbar", width = .5) +
  geom_boxplot(aes(fill = class1), outliers = FALSE, width = .85) +
  geom_text(aes(x = class1, y = rt_median, color = class1,
                label = sprintf("Median = %.2f", rt_median)),
            vjust = -.5, data = rt_thres, size = 4.5) +
  ggpubr::geom_signif(comparisons = list(c("High", "Low")),
                      map_signif_level = p_format,
                      vjust = -.5,
                      y_position = .26,
                      step_increase = .2,
                      tip_length = .05,
                      textsize = 6,
                      test = "t.test") +
  scale_y_continuous(n.breaks = 10,
                     expand = expansion(mult = c(.05, 0.2)),
                     name = "TIC Deviation") +
  scale_fill_manual(values = colors.ms)+
  scale_color_manual(values = c("black", "white"))+
  ggh4x::facet_grid2( ~ class2 + source) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 20),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm"))

ggsave("./chapter4_4_14.pdf", p, height = 5, width = 14)


## 9.1: 作图: 校正后MS1质量精度动态变化 ------
rm(list = ls())
gc()

all_files <- list.files("./results/tables", recursive = TRUE, full.names = TRUE)
all_lcms_files <- all_files[grepl("lc_ms", all_files)]

labels.source <- c("High-Confidence")
names(labels.source) <- c("_quali_filter")
lcms_tables <- pblapply(all_lcms_files, function(file_id) {
  
  source_id <- labels.source[str_extract(file_id, "(?<=lc_ms_metrics).*(?=\\.rds)")]
  source_id <- ifelse(is.na(source_id), "All", source_id)
  
  lcms_tables_i <- readRDS(file_id)
  qc_tables_i <- pblapply(lcms_tables_i, function(tmp_table) {
    
    qc_df_i <- tmp_table %>%
      group_by(lab, sample, tube, injection, order) %>%
      summarise_at("MS1 Accuracy", ~ median(abs(.), na.rm = TRUE))
    
    return(qc_df_i)
  })
  lcms_df <- rbindlist(qc_tables_i) %>% mutate(source = source_id)
  
  return(lcms_df)
})

df_lcms <- lcms_tables %>% rbindlist %>% arrange(lab, order)
rm(list = setdiff(ls(), c("df_lcms")))
gc()

labels.order <- 1:36
names(labels.order) <- as.character(unique(df_lcms$order))

df_tmp <- df_lcms %>%
  mutate(order_code = labels.order[as.character(order)]) %>%
  arrange(source, desc(`MS1 Accuracy`))

colors.lab <- brewer.pal(10, "PiYG")
names(colors.lab) <- unique(df_tmp$lab)

p <- ggplot(df_tmp, aes(y = `MS1 Accuracy`)) +
  stat_summary(aes(x = order_code, color = lab), fun = median, geom = "line", linewidth = 1) +
  stat_summary(aes(x = order_code, fill = lab), fun = median, geom = "point",
               size = 5, shape = 21, color = "black") +
  scale_y_continuous(n.breaks = 10, name = "Mass Error (ppm)") +
  scale_x_continuous(breaks = 1:36, name = "Injection order") +
  scale_fill_manual(values = colors.lab, name = "Lab") +
  scale_color_manual(values = colors.lab, name = "Lab") +
  facet_wrap(~ source, ncol = 1) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm"))

ggsave("./chapter4_4_15.pdf", p, height = 10, width = 14)


## 9.2: 作图: 校正后质量精度 全特征 vs 高置信特征 ------
rm(list = setdiff(ls(), c("df_lcms", "colors.lab")))
gc()

labels.lab <- 1:10
names(labels.lab) <- c("TFS", "BTP", "FDU", "QLB", "CMS", "NIM", "NCP", "OSB", "ZJU", "CAS")

df_tmp <- df_lcms %>%
  mutate(lab_code = labels.lab[lab])

rt_thres <- df_tmp %>%
  group_by(source) %>%
  summarise(rt_median = median(`MS1 Accuracy`),
            rt_max = max(`MS1 Accuracy`),
            rt_mean = mean(`MS1 Accuracy`),
            rt_sd = sd(`MS1 Accuracy`),
            .groups = "drop")

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

p1 <- ggplot(df_tmp, aes(x = lab_code, y = `MS1 Accuracy`)) +
  geom_hline(aes(yintercept = rt_median), rt_thres, lty = 2) +
  stat_summary(aes(group = source), fun = mean, geom = "line", color = "#999") +
  stat_summary(aes(group = source), fun.data = "mean_se", geom = "errorbar", width = .2) +
  stat_summary(aes(fill = source), fun = mean, geom = "point", color = "black", shape = 21, size = 5) +
  geom_text(aes(x = 1.1, y = rt_median,
                vjust = ifelse(source %in% "All", -.3, 1.2),
                label = sprintf("%.2f\u00B1%.2f\n(Mean\u00B1SD)", rt_mean, rt_sd)),
            data = rt_thres, color = "black", size = 4.5) +
  scale_fill_manual(values = c("#74C476", "#006D2C")) +
  scale_x_continuous(breaks = 1:10, labels = names(labels.lab)) +
  scale_y_continuous(n.breaks = 8, name = "Mass Error (ppm)",
                     expand = expansion(mult = c(0.05, 0.15))) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 15),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .3, .5, .5), units = "cm"));p1

p2 <- ggplot(df_tmp, aes(x = source, y = `MS1 Accuracy`)) +
  stat_summary(fun.min = min, fun.max = max, geom = "errorbar", width = .5) +
  geom_hline(aes(yintercept = rt_median), rt_thres, lty = 2) +
  geom_boxplot(aes(fill = source), outliers = FALSE, width = .7) +
  geom_text(aes(x = source, y = rt_median, color = source,
                label = sprintf("Median = %.2f", rt_median)),
            vjust = -.5,
            data = rt_thres, size = 4.5) +
  ggpubr::geom_signif(comparisons = list(c("All", "High-Confidence")),
                      map_signif_level = p_format,
                      vjust = -.5,
                      step_increase = .2,
                      tip_length = .05,
                      textsize = 4,
                      test = "t.test") +
  scale_fill_manual(values = c("#74C476", "#006D2C")) +
  scale_color_manual(values = c("black", "white")) +
  scale_y_continuous(n.breaks = 8,
                     expand = expansion(mult = c(0.05, 0.15))) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 20),
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .5, .5, 0), units = "cm"));p2

p <- plot_grid(p1, p2, ncol = 2, rel_widths = c(1, .5))

ggsave("./results/figures/chapter4_test_mass_error.pdf", p, height = 4, width = 14)


## 10.1: 作图: 原始质量精度动态变化 ------
rm(list = setdiff(ls(), c("colors.lab")))
gc()

df_lcms <- fread("./results/tables/all_ms_accuracy.csv")
df_lcms <- df_lcms %>%
  reshape2::melt(., id = c(1, 4:12), variable.name = "source",
                 value.name = "MS Accuracy") %>%
  filter(!lab %in% "NCP") %>%
  arrange(lab, order)

labels.order <- 1:36
names(labels.order) <- as.character(unique(df_lcms$order))

df_tmp <- df_lcms %>%
  mutate(order_code = labels.order[as.character(order)])

p <- ggplot(df_tmp, aes(y = `MS Accuracy`)) +
  stat_summary(aes(x = order_code, color = lab), fun = median, geom = "line", linewidth = 1) +
  stat_summary(aes(x = order_code, fill = lab), fun = median, geom = "point",
               size = 5, shape = 21, color = "black") +
  scale_y_continuous(n.breaks = 5, name = "Mass Accuracy") +
  scale_x_continuous(breaks = 1:36, name = "Injection order") +
  # scale_fill_manual(values = colors.sample, name = "Sample") +
  scale_fill_manual(values = colors.lab, name = "Lab")+
  scale_color_manual(values = colors.lab, name = "Lab")+
  ggh4x::facet_wrap2(source ~., scales = "free_y", ncol = 1) +
  theme_bw() +
  theme(legend.position = "right",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm"))

ggsave("./chapter4_test_ms1_ms2_acc.pdf", p, height = 10, width = 14)


## 10.2: 作图: 质量精度 MS1 vs MS2 ------
rm(list = setdiff(ls(), c("colors.lab")))
gc()

df_lcms <- fread("./results/tables/all_ms_accuracy.csv")
df_lcms <- df_lcms %>%
  reshape2::melt(., id = c(1, 4:12), variable.name = "source",
                 value.name = "MS Accuracy") %>%
  arrange(lab, order)

labels.class <- c("Orbitrap", "Orbitrap", "Orbitrap", "timsTOF", "Orbitrap",
                  "Orbitrap", "Orbitrap",	"timsTOF",	"timsTOF",	"ZenoTOF")
names(labels.class) <- c("ZJU", "QLB", "NCP", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")

df_tmp <- df_lcms %>%
  mutate(class = labels.class[lab]) %>%
  mutate_at("source", ~ gsub(" Accuracy", "", .))

rt_thres <- df_tmp %>%
  group_by(lab, source) %>%
  summarise(rt_median = median(`MS Accuracy`),
            rt_max = max(`MS Accuracy`),
            rt_mean = mean(`MS Accuracy`),
            rt_sd = sd(`MS Accuracy`),
            .groups = "drop") %>%
  arrange(source, desc(rt_median))

colors.lab <- brewer.pal(10, "PiYG")
names(colors.lab) <- rt_thres$lab[11:20]

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

df_tmp1 <- df_tmp %>%
  filter(!lab %in% "NCP") %>%
  mutate_at("lab", ~ factor(., levels = rt_thres$lab[12:20]))

rt_thres1 <- rt_thres %>%
  filter(!lab %in% "NCP") %>%
  mutate_at("lab", ~ factor(., levels = rt_thres$lab[12:20]))

p1 <- ggplot(df_tmp1, aes(x = lab, y = `MS Accuracy`)) +
  stat_summary(fun.min = min, fun.max = max, geom = "errorbar", width = .5) +
  geom_boxplot(aes(fill = lab), outliers = FALSE, width = .7) +
  geom_text(aes(x = lab, y = 10,
                label = sprintf("%.2f\u00B1%.2f\n(Mean\u00B1SD)", rt_mean, rt_sd)),
            vjust = -.5, color = "black",
            data = rt_thres1, size = 4.5) +
  scale_y_continuous(n.breaks = 10,
                     expand = expansion(mult = c(.05, 0.25)),
                     name = "Mass Accuracy (ppm)") +
  scale_fill_manual(values = colors.lab)+
  facet_wrap(~ source, ncol = 1) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 16),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm"))

p2 <- ggplot(df_tmp[df_tmp$lab %in% "NCP", ], aes(x = lab, y = `MS Accuracy`)) +
  stat_summary(fun.min = min, fun.max = max, geom = "errorbar", width = .5) +
  geom_boxplot(aes(fill = lab), outliers = FALSE, width = .7) +
  geom_text(aes(x = lab, y = 120,
                label = sprintf("%.2f\u00B1%.2f\n(Mean\u00B1SD)", rt_mean, rt_sd)),
            vjust = -.5, color = "black",
            data = rt_thres[rt_thres$lab %in% "NCP", ], size = 4.5) +
  scale_x_discrete(label = ~ gsub(" Accuracy", "", .)) +
  scale_y_continuous(n.breaks = 10,
                     expand = expansion(mult = c(.05, 0.25)),
                     name = "Mass Accuracy (ppm)") +
  scale_fill_manual(values = colors.lab)+
  facet_wrap(~ source, ncol = 1) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 16),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(.5, "cm"),
        plot.margin = unit(c(.5, .5, .5, 0), units = "cm"))

p <- plot_grid(p1, p2, nrow = 1, rel_widths = c(1, .17))

ggsave("./chapter4_4_16.pdf", p, height = 10, width = 14)


## 11: 作图: IM动态变化 ------
rm(list = ls())
gc()

all_files <- list.files("./results/tables", recursive = TRUE, full.names = TRUE)
all_lcms_files <- all_files[grepl("lc_ms", all_files)]

labels.source <- c("High-Confidence")
names(labels.source) <- c("_quali_filter")
lcms_tables <- pblapply(all_lcms_files, function(file_id) {
  
  source_id <- labels.source[str_extract(file_id, "(?<=lc_ms_metrics).*(?=\\.rds)")]
  source_id <- ifelse(is.na(source_id), "All", source_id)
  
  lcms_tables_i <- readRDS(file_id)
  qc_tables_i <- pblapply(lcms_tables_i, function(tmp_table) {
    
    qc_df_i <- tmp_table %>%
      group_by(lab, sample, tube, injection, order) %>%
      summarise_at("Delta IM", ~ median(abs(.), na.rm = TRUE))
    
    return(qc_df_i)
  })
  lcms_df <- rbindlist(qc_tables_i) %>% mutate(source = source_id)
  
  return(lcms_df)
})

df_lcms <- lcms_tables %>% rbindlist %>% arrange(lab, order)
rm(list = setdiff(ls(), c("df_lcms")))
gc()

labels.order <- 1:36
names(labels.order) <- as.character(unique(df_lcms$order))

df_tmp <- df_lcms %>%
  mutate(order_code = labels.order[as.character(order)]) %>%
  filter(!(is.na(`Delta IM`)|`Delta IM`==0)) %>%
  arrange(source, desc(`Delta IM`)) 

colors.lab <- brewer.pal(10, "PiYG")[c(2, 5, 10)]
names(colors.lab) <- unique(df_tmp$lab)

p <- ggplot(df_tmp, aes(y = `Delta IM`)) +
  stat_summary(aes(x = order_code, color = lab), fun = median, geom = "line", linewidth = 1) +
  stat_summary(aes(x = order_code, fill = lab), fun = median, geom = "point",
               size = 5, shape = 21, color = "black") +
  scale_y_continuous(n.breaks = 10, name = "Delta IM") +
  scale_x_continuous(breaks = 1:36, name = "Injection order") +
  scale_fill_manual(values = colors.lab, name = "Lab") +
  scale_color_manual(values = colors.lab, name = "Lab") +
  facet_wrap(~ source, ncol = 1) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        strip.background = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        plot.margin = unit(c(.5, .5, .5, .5), units = "cm"))

ggsave("./chapter4_4_17.pdf", p, height = 10, width = 14)


## 整理附录4.1 --------------------
all_metrics <- fread("./results/tables/all_metrics.csv")

all_stat <- all_metrics %>%
  select(any_of(c("lab", "data_level", "source", "snr", "Qualification", "Quantification"))) %>%
  group_by(lab, data_level, source) %>%
  summarise_if(is.numeric, ~ sprintf("%.3f\u00B1%.3f",
                                     mean(., na.rm = TRUE),
                                     sd(., na.rm = TRUE)),
               .groups = "drop")

fwrite(all_stat, "~/Desktop/tmp.csv")


## 整理附录4.2 --------------------
all_metrics <- fread("./results/tables/all_metrics.csv")

all_stat <- all_metrics %>%
  filter(source %in% "All", data_level %in% "Precursor-level") %>%
  select(any_of(c("lab", "FWHM", "Delta RT",
                  "Charge_2_3_prop", "TIC Deviation", "Calibarted Mass Error",
                  "MS1 Accuracy", "MS2 Accuracy", "Delta IM"))) %>%
  mutate_at("Delta IM", ~ ifelse(. == 0, NA, .)) %>%
  group_by(lab) %>%
  summarise_if(is.numeric, ~ sprintf("%.3f\u00B1%.3f",
                                     mean(., na.rm = TRUE),
                                     sd(., na.rm = TRUE)),
               .groups = "drop")

fwrite(all_stat, "~/Desktop/tmp.csv")


## 所有指标相关性矩阵 ---------
rm(list = ls())
gc()

## Count & Recall
all_count_recall <- fread("./results/tables/all_count_recall.csv")
all_count_recall <- all_count_recall %>%
  select(lab, sample, tube, recall, contains("s.Identified")) %>%
  dplyr::rename_at(5:7, ~ gsub("s.Identified", "-level", ., fixed = TRUE)) %>%
  reshape2::melt(., id = 1:4, variable.name = "data_level", value.name = "count") %>%
  mutate_at("recall", ~ ifelse(data_level %in% "Precursor-level", ., NA)) %>%
  reshape2::melt(., id = c(1:3, 5), variable.name = "source", value.name = "Qualification") %>%
  mutate_at("source", ~ ifelse(. %in% "recall", "High-Confidence", "All")) %>%
  na.omit

## CV & F1
all_cv_f1 <- fread("./results/tables/all_cv_f1.csv")
all_cv_f1 <- all_cv_f1 %>%
  select(lab, sample, tube, f1, contains("level")) %>%
  reshape2::melt(., id = 1:4, variable.name = "data_level", value.name = "cv") %>%
  mutate_at("f1", ~ ifelse(data_level %in% "Precursor-level", ., NA)) %>%
  reshape2::melt(., id = c(1:3, 5), variable.name = "source", value.name = "Quantification") %>%
  mutate_at("source", ~ ifelse(. %in% "f1", "High-Confidence", "All")) %>%
  na.omit

## PCA & SNR
all_files <- list.files("./results/tables", recursive = TRUE, full.names = TRUE)
all_pca_files <- all_files[grepl("a_pca_quartet", all_files)]
labels.source <- c("High-Confidence (Raw)", "High-Confidence (Ratio)", "All (Ratio)")
names(labels.source) <- c("_quali_filter", "_quant_filter", "_ratiobyd6")
pca_tables <- pblapply(all_pca_files, function(file_id) {
  
  source_id <- labels.source[str_extract(file_id, "(?<=quartet).*(?=\\.rds)")]
  source_id <- ifelse(is.na(source_id), "All (Raw)", source_id)
  
  pca_tables_i <- readRDS(file_id)
  pca_tables_i <- mclapply(pca_tables_i, function(tmp_results) {
    
    tmp_tables <- mclapply(tmp_results, function(tmp_result) {
      
      pca_tmp_j <- tmp_result$pcs_values %>%
        select(lab, sample, tube, injection, PC1, PC2) %>%
        mutate(snr = tmp_result$snr_results$snr)
      
      return(pca_tmp_j)
    })
    
    pca_tmp_i <- rbindlist(tmp_tables)
    
    return(pca_tmp_i)
  })
  names(pca_tables_i) <- c("Precursor-level", "Peptide-level", "Protein-level")
  pca_df <- rbindlist(pca_tables_i, idcol = "data_level") %>% mutate(source = source_id)
  
  return(pca_df)
})
df_pca <- pca_tables %>%
  rbindlist %>%
  filter(source %in% c("All (Raw)", "High-Confidence (Raw)")) %>%
  mutate(source = str_extract(source, "All|High-Confidence")) %>%
  select(source, lab, sample, tube, injection, data_level, snr, contains("PC"))

## LC-MS metrics
all_files <- list.files("./results/tables", recursive = TRUE, full.names = TRUE)
all_lcms_files <- all_files[grepl("lc_ms", all_files)]
labels.source <- c("High-Confidence")
names(labels.source) <- c("_quali_filter")
lcms_tables <- pblapply(all_lcms_files, function(file_id) {
  
  source_id <- labels.source[str_extract(file_id, "(?<=lc_ms_metrics).*(?=\\.rds)")]
  source_id <- ifelse(is.na(source_id), "All", source_id)
  
  lcms_tables_i <- readRDS(file_id)
  qc_tables_i <- pblapply(lcms_tables_i, function(tmp_table) {
    
    qc_df_i_1 <- tmp_table %>%
      group_by(lab, sample, tube, injection) %>%
      summarise_at(c("FWHM", "MS1 Accuracy", "Delta RT", "Delta IM"), ~ median(abs(.), na.rm = TRUE))
    
    qc_df_i_2 <- tmp_table %>%
      group_by(lab, sample, tube) %>%
      summarise(`TIC Deviation` = sd(TIC, na.rm = TRUE) / mean(TIC, na.rm = TRUE),
                .groups = "drop")
    
    charge_df_i <- tmp_table %>%
      reshape2::dcast(., lab + sample + tube + injection ~ Charge,
                      value.var = "peptide_sequence",
                      fun.aggregate = length) %>%
      mutate(Total = rowSums(pick(matches("^\\d+$")), na.rm = TRUE)) %>%
      mutate(Charge_1_prop = `1` / Total,
             Charge_2_prop = `2` / Total,
             Charge_3_prop = `3` / Total,
             Charge_4_prop = 4 / Total) %>%
      mutate(Charge_2_3_prop = Charge_2_prop + Charge_3_prop) %>%
      select(lab, sample, tube, injection, Charge_2_3_prop)
    
    qc_df_i <- qc_df_i_1 %>%
      full_join(., qc_df_i_2, by = c("lab", "sample", "tube")) %>%
      full_join(., charge_df_i, by = c("lab", "sample", "tube", "injection")) 
    
    return(qc_df_i)
  })
  lcms_df <- rbindlist(qc_tables_i) %>% mutate(source = source_id)
  
  return(lcms_df)
})
df_lcms <- lcms_tables %>%
  rbindlist %>%
  dplyr::rename(`Calibarted Mass Error` = `MS1 Accuracy`) %>%
  mutate(data_level = "Precursor-level")

## Mass Accuracy
df_mass_acc <- fread("./results/tables/all_ms_accuracy.csv")
df_mass_acc <- df_mass_acc %>%
  select(lab, sample, tube, injection, `MS1 Accuracy`, `MS2 Accuracy`) %>%
  mutate(source = "All", data_level = "Precursor-level")

## Long: Combined all
all_metrics <- df_pca %>%
  full_join(., df_lcms, by = c("source", "lab", "sample", "tube", "injection", "data_level")) %>%
  full_join(., all_count_recall, by = c("lab", "sample", "tube", "source", "data_level")) %>%
  full_join(., all_cv_f1, by = c("lab", "sample", "tube", "source", "data_level")) %>%
  full_join(., df_mass_acc, by = c("lab", "sample", "tube", "injection", "source", "data_level"))

# fwrite(all_metrics, "./results/tables/all_metrics.csv")

## Long: Scale to 1~10
all_metrics_scaled <- all_metrics %>%
  mutate_at("Delta IM", ~ ifelse(. == 0, NA, .)) %>%
  group_by(source, data_level, lab) %>%
  mutate_at("Delta IM",
            ~ 10 - (. - min(., na.rm = TRUE)) * 9 / (max(., na.rm = TRUE) - min(., na.rm = TRUE))) %>%
  mutate_at(c("FWHM", "Calibarted Mass Error", "Delta RT", "TIC Deviation",
              "MS1 Accuracy", "MS2 Accuracy"),
            ~ 10 - (. - min(.)) * 9 / (max(.) - min(.))) %>%
  mutate_at(c("PC1", "PC2", "Charge_2_3_prop", "Qualification"),
            ~ 1 + (. - min(.)) * 9 / (max(.) - min(.))) %>%
  mutate_at("Quantification", ~ ifelse(source %in% "All", ## 代表是cv而不是f1
                                       10 - (. - min(.)) * 9 / (max(.) - min(.)),
                                       1 + (. - min(.)) * 9 / (max(.) - min(.)))) %>%
  ungroup %>%
  mutate_at("snr", ~ 1 + (. - min(.)) * 9 / (max(.) - min(.))) %>%
  select(source, data_level, lab, sample, tube, injection, everything()) %>%
  mutate_at("data_level", ~ factor(., levels = c("Precursor-level", "Peptide-level", "Protein-level")))

fwrite(all_metrics_scaled, "./results/tables/all_metrics_scaled.csv")

## Long: 每个指标取各实验室中位值 & Scale to 1~10
all_metrics_scaled2 <- all_metrics %>%
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

fwrite(all_metrics_scaled2, "./results/tables/all_metrics_scaled_median_bylab.csv")
