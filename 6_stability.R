## Title: Short- and Long-term stability.
## Author: Qiaochu Chen
## Date: Jun 7th, 2026

library(readxl)
library(data.table)
library(rstatix)
library(pbapply)
library(parallel)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(lubridate)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)


## 准备数据10--所有稳定性数据参考量值ratiobyd6 ----------------
all_tables <- readRDS("./data/stability/quantdata_list_pep_all_longterm.rds")
allbylab_tables <- all_tables %>%
  rbindlist %>%
  mutate(dataset = paste(lab_id, time_day)) %>%
  split(., ~ dataset)

rm(list = setdiff(ls(), c("all_tables", "allbylab_tables")))
gc()

## 每家实验室PCA检查离群值: 马氏距离法(只对Quartet样本)
source("./PCA.R")
madist_results_tables <- pblapply(allbylab_tables, function(tmp_table) {
  
  metadata <- tmp_table %>%
    select(!value & !peptide_sequence) %>%
    distinct() %>%
    filter(grepl("Quartet", sample)) %>%
    tibble::column_to_rownames("analysis_id")
  
  exprdata_t <- tmp_table %>%
    distinct(peptide_sequence, analysis_id, value) %>%
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
fwrite(sub_meta2outlier, "./results/tables/6_outlier_madist_quartet.csv")

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
ratiobyd6_tables <- pblapply(allbylab_tables, function(tmp_table) {
  
  meta_tmp <- tmp_table %>%
    distinct(analysis_id) %>%
    inner_join(., sub_meta2outlier, by = c("analysis_id")) %>%
    filter(!is.outlier) %>%
    mutate(batch = ifelse(time_day < 1608, paste(lab_id, time_day, sep = "_"),
                          paste(lab_id, tube, sep = "_")))
  
  expr_tmp <- reshape2::dcast(tmp_table, peptide_sequence ~ analysis_id, value.var = "value", fun.aggregate = sum)
  
  expr_tmp <- expr_tmp %>% mutate_if(is.numeric, ~ ifelse(. == 0 | is.na(.) | is.nan(.), NA, .))
  
  if (unique(meta_tmp$lab_id) %in% c("phoenix", "zhejiang_university", "qinglian_bio") | grepl("DDA", unique(meta_tmp$lab_id))) {
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
saveRDS(ratiobyd6_tables, "./results/tables/6_quantdata_list_pep_ratiobyd6_all_longterm.rds")


## 标称特性的稳定性检验: CFFF计算 ------------------------------
# all_tables <- readRDS("./data/stability/quantdata_list_pep_all_longterm.rds")
rm(list = setdiff(ls(), c("all_tables")))
gc()

## 输入数据
passPEPfdr0.01_pep_tables <- readRDS("./results/tables/1_qualiprop_list_PEPfdr0.01.rds")
passPEPfdr0.01_pep_tables <- passPEPfdr0.01_pep_tables[c(1:3, 6)]

filtered_tables <- pblapply(1:4, function(i) {
  
  qualiprop_tables_i <- passPEPfdr0.01_pep_tables[[i]] %>%
    distinct(peptide_sequence, protein_id) %>%
    inner_join(., all_tables[[i]], by =c("peptide_sequence"),
               relationship = "many-to-many")
  
  return(qualiprop_tables_i)
})
names(filtered_tables) <- names(all_tables)

rm(list = setdiff(ls(), c("filtered_tables")))
gc()

all_meta <- fread("./results/tables/6_outlier_madist_quartet.csv")

## Fisher's Exact Test/Chi-square test: 参考ISO 33406:2024-5.5.3.4.1
fisher_tables <- pblapply(filtered_tables, function(tmp_table_i) {
  
  print(unique(tmp_table_i$sample))
  all_meta_i <- all_meta %>% filter(sample %in% unique(tmp_table_i$sample))
  tmp_table_i <- tmp_table_i %>% mutate(feature = paste(peptide_sequence, protein_id))
  
  all_features <- unique(tmp_table_i$feature)
  
  fisher_tables_i <- mclapply(all_features, function(feature_id) {
    
    tryCatch({sub_tmp_table_j <- tmp_table_i %>%
      filter(feature %in% feature_id) %>%
      ungroup() %>%
      select(!value) %>%
      distinct(feature, analysis_id) %>%
      inner_join(all_meta_i, ., by = "analysis_id") %>%
      filter(!is.outlier)
    
    ## 短期稳定性
    sub_table_j1 <- sub_tmp_table_j %>%
      filter(lab_id %in% c("NIM_DDA", "HUP_DIA")) %>%
      filter(time_day >= 1290) %>%
      mutate(time = ifelse(lab_id %in% "HUP_DIA", 1320 - time_day, time_day - 1551)) %>%
      filter(time != 30) %>%
      mutate_at("temperature", ~ ifelse(. %in% "-80°C", "-20°C_-78.5°C", .)) %>%
      tidyr::separate_rows(temperature, sep = "_") %>%
      mutate(dataset_group = paste(lab_id, temperature))
    
    all_dataset_groups <- unique(sub_table_j1$dataset_group)
    
    fisher_tables_j1 <- mclapply(all_dataset_groups, function(dataset_group_id) {
      
      sub_table_k <- sub_table_j1 %>%
        filter(dataset_group %in% dataset_group_id) %>%
        reshape2::acast(., feature ~ lab_id + time_day,
                        value.var = "feature",
                        fun.aggregate = length)
      
      sub_table_k_theoretical <- all_meta_i %>%
        mutate_at("temperature", ~ ifelse(time_day == 1320, "-20°C_-78.5°C", .)) %>%
        tidyr::separate_rows(temperature, sep = "_") %>%
        mutate(time = ifelse(lab_id %in% "HUP_DIA", 1320 - time_day, time_day - 1551)) %>%
        filter(time != 30) %>%
        mutate(dataset_group = paste(lab_id, temperature)) %>%
        filter(dataset_group %in% dataset_group_id) %>%
        group_by(lab_id, time_day) %>%
        summarise(rep = length(order), feature = feature_id, .groups = "drop") 
      
      theoretical_reps <- sub_table_k_theoretical %>%
        reshape2::acast(., feature ~ lab_id + time_day, value.var = "rep")
      
      if (nrow(sub_table_k) == 1 & rownames(sub_table_k)[1] %in% "NA") return(NULL)
      
      if (ncol(sub_table_k) < 4 & grepl("HUP_DIA", dataset_group_id)) return(NULL)
      
      if (ncol(sub_table_k) < 5 & grepl("NIM_DDA", dataset_group_id)) return(NULL)
      
      if (nrow(sub_table_k) == 1 & !rownames(sub_table_k)[1] %in% "NA") sub_table_k <- rbind(sub_table_k, theoretical_reps-sub_table_k)
       
      aa <- fisher.test(sub_table_k)
      
      all_times <- sub_table_k_theoretical %>%
        mutate(time = ifelse(lab_id %in% "HUP_DIA", 1320 - time_day, time_day - 1551)) %>%
        pull(time) %>%
        sort
      
      fisher_table_k <- data.frame(feature = feature_id,
                                   p = aa$p.value,
                                   n_timepoints = ncol(sub_table_k),
                                   period = paste0(all_times, collapse = ", ")) %>%
        mutate(dataset_group = dataset_group_id) %>%
        tidyr::separate(dataset_group, c("category", "temperature"), sep = " ")
      
      return(fisher_table_k)
      
    }, mc.cores = 4)
    
    fisher_df_j1 <- rbindlist(fisher_tables_j1, use.names = TRUE)
    
    if (nrow(fisher_df_j1) == 0) fisher_df_j1 <- NULL
    
    if (!is.null(fisher_df_j1)) {
      fisher_df_j1 <- fisher_df_j1 %>%
        mutate_at("category", ~ sapply(., function(x) {
          if (as.character(x) %in% "NIM_DDA") {
            a <- "Short-term-Use"
          } else if (as.character(x) %in% "HUP_DIA") {
            a <- "Short-term-Transport"
          }
        }))
    }
    
    ## 长期稳定性
    sub_table_j2 <- sub_tmp_table_j %>%
      reshape2::acast(., feature ~ lab_id + time_month,
                      value.var = "feature",
                      fun.aggregate = length)
    
    sub_table_j2_theoretical <- all_meta_i %>%
      group_by(lab_id, time_month) %>%
      summarise(rep = length(order), .groups = "drop") %>%
      mutate(dataset_group = paste(lab_id, time_month, sep = "_"), feature = feature_id) %>%
      filter(dataset_group %in% colnames(sub_table_j2)) 
    
    theoretical_reps <- sub_table_j2_theoretical %>%
      reshape2::acast(., feature ~ lab_id + time_month, value.var = "rep")
    
    if (nrow(sub_table_j2) == 1 & rownames(sub_table_j2)[1] %in% "NA") return(fisher_df_j1)
    
    if (ncol(sub_table_j2) < 12) return(fisher_df_j1)
    
    if (!"fudan_university_58.03" %in% colnames(sub_table_j2)) return(fisher_df_j1)
    
    if (!"NVG_DDA_0" %in% colnames(sub_table_j2)) return(fisher_df_j1)
    
    if (nrow(sub_table_j2) == 1 & !rownames(sub_table_j2)[1] %in% "NA") sub_table_j2 <- rbind(sub_table_j2, theoretical_reps-sub_table_j2)
    
    if (sum(rowSums(sub_table_j2) == 0) > 0) {
      aa <- fisher.test(sub_table_j2)
    } else {
      aa <- fisher.test(sub_table_j2, simulate.p.value = TRUE)
    }
    
    all_times <- sub_table_j2_theoretical %>%
      mutate(time = time_month) %>%
      pull(time) %>%
      sort
    
    fisher_df_j2 <- data.frame(feature = feature_id,
                               p = aa$p.value,
                               n_timepoints = ncol(sub_table_j2),
                               period = paste0(all_times, collapse = ", ")) %>%
      mutate(category = "Long-term", temperature = "-80°C")
    
    fisher_df_j <- rbind(fisher_df_j1, fisher_df_j2)
    
    return(fisher_df_j)}, error = function(e) {
      message("Error message: ", e$message)
      stop("Stopping execution due to error at feature_id = ", feature_id)})
    
  }, mc.cores = 12)
  
  fisher_df_i <- rbindlist(fisher_tables_i)
  
  return(fisher_df_i)
})

bb <- fisher_tables %>% rbindlist %>% mutate(is.stab = ifelse(p > .05, TRUE, FALSE))
nrow(bb[is.stab == TRUE, ])

## 保存结果
saveRDS(fisher_tables, "./results/tables/6_qualiprop_stabtest_fisher.rds")


## 参考量值的稳定性检验 ------------------------------
rm(list = ls())
gc()

## 输入数据
all_meta <- fread("./results/tables/6_outlier_madist_quartet.csv")

passPfdr0.05_pep_tables <- readRDS("./results/tables/2_quantprop_list_Pfdr0.05.rds")
passPfdr0.05_pep_tables <- passPfdr0.05_pep_tables[c(1:2, 4)]

ratiobyd6_tables <- readRDS("./results/tables/6_quantdata_list_pep_ratiobyd6_all_longterm.rds")
ratiobyd6_tables <- ratiobyd6_tables[c(1, 3:4)]
names(ratiobyd6_tables) <- names(passPfdr0.05_pep_tables)[1:3]

all_tables <- ratiobyd6_tables
quantprop_tables <- pblapply(1:length(all_tables), function(i) {
  
  quantprop_tables_i <- passPfdr0.05_pep_tables[[i]] %>%
    distinct(peptide_sequence, protein_id) %>%
    inner_join(., all_tables[[i]], by =c("peptide_sequence"),
               relationship = "many-to-many")
  
  return(quantprop_tables_i)
})
names(quantprop_tables) <- names(all_tables)

rm(list = setdiff(ls(), c("quantprop_tables")))
gc()

## 线性模型: 参考ISO 33405:2024-8.5.2.3和JJF 1343-2022-6.4.3
mlr_tables <- pblapply(quantprop_tables, function(tmp_table_i) {
  
  tmp_table_i <- tmp_table_i %>%
    mutate(feature = paste(peptide_sequence, protein_id)) %>%
    mutate_at("value", ~ 2 ^ (.)) %>%
    # filter(time_month >= 41) %>% ##原计划以第41月为0点监测长期稳定性，被弃用
    # mutate_at("time_month", ~ . - 41) %>%
    filter(!is.infinite(value))
  
  all_features <- unique(tmp_table_i$feature)
  
  mlr_tables_i <- mclapply(all_features, function(feature_id) {
    
    tryCatch({
    sub_expr_tmp <- tmp_table_i %>%
      filter(feature %in% feature_id) %>%
      ungroup()
    
    ## 短期稳定性
    sub_table_j1 <- sub_expr_tmp %>%
      filter(lab_id %in% c("NIM_DDA", "HUP_DIA")) %>%
      filter(time_day >= 1290) %>%
      mutate_at("temperature", ~ ifelse(. %in% "-80°C", "-20°C_-78.5°C", .)) %>%
      tidyr::separate_rows(temperature, sep = "_") %>%
      mutate(dataset_group = paste(lab_id, temperature)) %>%
      mutate(time = ifelse(lab_id %in% "HUP_DIA", 1320 - time_day, time_day - 1551))
    
    all_dataset_groups <- unique(sub_table_j1$dataset_group)
    
    mlr_tables_j1 <- mclapply(all_dataset_groups, function(dataset_group_id) {
      
      sub_table_k <- sub_table_j1 %>%
        filter(dataset_group %in% dataset_group_id)
      
      if (nrow(sub_table_k) == 0) return(NULL)
      
      if (length(unique(sub_table_k$time)) < 3) return(NULL)
      
      all_times <- sort(unique(sub_table_k$time))
      
      if (length(all_times) >= 3) {
        
        all_periods <- list(all_times[1:3])
        if(length(all_times) >= 4) {
          for (i in 4:length(all_times)) {
            all_periods <- c(all_periods, list(all_times[1:i]))
          }
        }
        
        stab_tables_k <- mclapply(all_periods, function(period_id) {
          
          sub_expr_k <- sub_table_k %>%
            filter(time %in% period_id)
          
          lm_l <- lm(formula = value ~ time, data = sub_expr_k)
          
          sub_lm_k <- sub_expr_k %>%
            mutate(value_predicted = lm_l$fitted.values) %>%
            summarise(feature = feature_id,
                      dataset_group = dataset_group_id,
                      period = paste0(period_id, collapse = ", "),
                      b0 = lm_l$coefficients[1],
                      b1 = lm_l$coefficients[2],
                      n = length(period_id),
                      s = sqrt(sum((value - value_predicted) ^ 2) / (n - 2)),
                      s_b1 = s / sqrt(sum((time - mean(time)) ^ 2)),
                      `s_b1'` = s_b1 * qt(0.95, n - 2),
                      pass = ifelse(abs(b1) < `s_b1'`, "Yes", "No")) %>%
            tidyr::separate(dataset_group, c("category", "temperature"), sep = " ")
          
          return(sub_lm_k)
          
        })
        
        sub_stab_k <- rbindlist(stab_tables_k) %>%
          mutate_at("category", ~ sapply(., function(x) {
            if (as.character(x) %in% "NIM_DDA") {
              a <- "Short-term-Use"
            } else if (as.character(x) %in% "HUP_DIA") {
              a <- "Short-term-Transport"
            }
          }))
        
      } else {
        sub_stab_k <- NULL
      }
      
      return(sub_stab_k)
    })
    
    mlr_df_j1 <- rbindlist(mlr_tables_j1)
    
    ## 长期稳定性
    sub_table_j2 <- sub_expr_tmp %>%
      mutate(time = time_month)
    
    if (nrow(sub_table_j2) == 0) sub_stab_k <- NULL
    
    if (length(unique(sub_table_j2$time)) < 3) sub_stab_k <- NULL
    
    all_times <- sort(unique(sub_table_j2$time))
    
    if (length(all_times) >= 12) {
      
      all_periods <- list(all_times[1:12])
      if(length(all_times) >= 13) {
        for (i in 13:length(all_times)) {
          all_periods <- c(all_periods, list(all_times[1:i]))
        }
      }
      
      stab_tables_k <- mclapply(all_periods, function(period_id) {
        
        # print(period_id)
        sub_expr_k <- sub_table_j2 %>%
          filter(time %in% period_id)
        
        lm_l <- lm(formula = value ~ time, data = sub_expr_k)
        
        sub_lm_k <- sub_expr_k %>%
          mutate(value_predicted = lm_l$fitted.values) %>%
          summarise(feature = feature_id,
                    period = paste0(period_id, collapse = ", "),
                    b0 = lm_l$coefficients[1],
                    b1 = lm_l$coefficients[2],
                    n = length(period_id),
                    s = sqrt(sum((value - value_predicted) ^ 2) / (n - 2)),
                    s_b1 = s / sqrt(sum((time - mean(time)) ^ 2)),
                    `s_b1'` = s_b1 * qt(0.95, n - 2),
                    pass = ifelse(abs(b1) < `s_b1'`, "Yes", "No"))
        
        return(sub_lm_k)
        
      })
      
      sub_stab_k <- stab_tables_k %>%
        rbindlist(.)  %>%
        mutate(category = "Long-term", temperature = "-80°C")
      
    } else {
      sub_stab_k <- NULL
    }
    
    mlr_df_j2 <- sub_stab_k
    
    mlr_df_j <- rbind(mlr_df_j1, mlr_df_j2, use.names=TRUE)
    
    return(mlr_df_j)}, error = function(e) {
      message("Error message: ", e$message)
      stop("Stopping execution due to error at feature_id = ", feature_id)
    })
  })
  
  stab_results_i <- rbindlist(mlr_tables_i, use.names=TRUE) %>%
    tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ")
  
  return(stab_results_i)
})

cc <- mlr_tables %>% rbindlist
nrow(cc[pass %in% "Yes", ])

saveRDS(mlr_tables, "./results/tables/6_quantprop_stabtest_linear.rds")


## test1 -----------------------
stab_tables <- readRDS("./results/tables/6_quantprop_stabtest_linear.rds")
# stab_tables <- mlr_tables
rm(list = setdiff(ls(), c("stab_tables")))
gc()

all_meta <- fread("./results/tables/6_outlier_madist_quartet.csv")
norm_tables <- readRDS("./results/tables/4_quantprop_normtest_shapirowilk.rds")
norm_tables <- norm_tables[c(1:2, 4)]
df_uchar <- norm_tables %>% rbindlist(., idcol = "group")

stat_stab1 <- stab_tables %>%
  rbindlist(., idcol = "group", use.names=TRUE) %>%
  # filter(category %in% "Short-term-Transport", grepl("14", period), !grepl("30", period)) %>%
  filter(category %in% "Short-term-Transport", period %in% ("0, 3, 7, 14")) %>%
  mutate_at(7:12, as.numeric) %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  distinct(peptide_sequence, group, temperature, pass) %>%
  reshape2::dcast(., group + temperature ~ pass, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(n = Yes + No) %>%
  mutate(prop = paste(round(Yes/n * 100, 2), "%"))

stat_stab1 <- stab_tables %>%
  rbindlist(., idcol = "group", use.names=TRUE) %>%
  # filter(category %in% "Short-term-Transport", grepl("14", period), !grepl("30", period)) %>%
  filter(category %in% "Short-term-Transport", period %in% ("0, 3, 7, 14")) %>%
  mutate_at(7:12, as.numeric) %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  distinct(protein_id, group, temperature, pass) %>%
  reshape2::dcast(., group + temperature ~ pass, value.var = "protein_id",
                  fun.aggregate = length) %>%
  mutate(n = Yes + No) %>%
  mutate(prop = paste(round(Yes/n * 100, 2), "%"))

fwrite(stat_stab1, "~/Desktop/tmp.csv")

stat_stab2 <- stab_tables %>%
  rbindlist(., idcol = "group", use.names=TRUE) %>%
  # filter(category %in% "Short-term-Use", grepl("7", period), !grepl("30", period)) %>%
  filter(category %in% "Short-term-Use", period %in% ("0, 1, 2, 4, 7")) %>%
  mutate_at(7:12, as.numeric) %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  distinct(peptide_sequence, group, temperature, pass) %>%
  reshape2::dcast(., group + temperature ~ pass, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(n = Yes + No) %>%
  mutate(prop = paste(round(Yes/n * 100, 2), "%"))

stat_stab2 <- stab_tables %>%
  rbindlist(., idcol = "group", use.names=TRUE) %>%
  # filter(category %in% "Short-term-Use", grepl("7", period), !grepl("30", period)) %>%
  filter(category %in% "Short-term-Use", period %in% ("0, 1, 2, 4, 7")) %>%
  mutate_at(7:12, as.numeric) %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  distinct(protein_id, group, temperature, pass) %>%
  reshape2::dcast(., group + temperature ~ pass, value.var = "protein_id",
                  fun.aggregate = length) %>%
  mutate(n = Yes + No) %>%
  mutate(prop = paste(round(Yes/n * 100, 2), "%"))

fwrite(stat_stab2, "~/Desktop/tmp.csv")

stat_stab3 <- stab_tables %>%
  rbindlist(., idcol = "group", use.names=TRUE) %>%
  filter(category %in% "Long-term", grepl("0", period) & grepl("58.03", period)) %>%
  mutate_at(7:12, as.numeric) %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  distinct(peptide_sequence, group, temperature, pass) %>%
  reshape2::dcast(., group ~ pass, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(n = Yes + No) %>%
  mutate(prop = paste(round(Yes/n * 100, 2), "%"))

stat_stab3 <- stab_tables %>%
  rbindlist(., idcol = "group", use.names=TRUE) %>%
  filter(category %in% "Long-term", grepl("0", period) & grepl("58.03", period)) %>%
  mutate_at(7:12, as.numeric) %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  distinct(protein_id, group, temperature, pass) %>%
  reshape2::dcast(., group ~ pass, value.var = "protein_id",
                  fun.aggregate = length) %>%
  mutate(n = Yes + No) %>%
  mutate(prop = paste(round(Yes/n * 100, 2), "%"))

fwrite(stat_stab3, "~/Desktop/tmp.csv")


## test2 -----------------------
stab_tables <- readRDS("./results/tables/6_qualiprop_stabtest_fisher.rds")
rm(list = setdiff(ls(), c("stab_tables")))
gc()

bb <- stab_tables$`Quartet D5`

aa <- stab_tables %>%
  rbindlist(., idcol = "group") %>%
  mutate(pass = ifelse((p > .05) & (!is.na(p)), "Yes", "No")) %>%
  tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
  distinct(peptide_sequence, protein_id, group, pass) %>%
  reshape2::dcast(., peptide_sequence + protein_id + group ~ pass, value.var = "pass") %>%
  na.omit ##不含缺失值说明不存在同一肽段yes/no的矛盾，否则考虑是否同一肽段匹配到了不同蛋白质，因而计算的P不同。

stat_stab1 <- stab_tables %>%
  rbindlist(., idcol = "group") %>%
  filter(category %in% "Short-term-Transport") %>%
  tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
  mutate_at("p", ~ ifelse(is.na(.), 0, .)) %>%
  group_by(peptide_sequence, group, temperature) %>%
  summarise_at("p", max) %>%
  mutate(pass = ifelse((p > .05) & (!is.na(p)), "Yes", "No")) %>%
  reshape2::dcast(., group + temperature ~ pass, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(n = Yes + No) %>%
  mutate(prop = paste(round(Yes/n * 100, 2), "%"))

stat_stab1 <- stab_tables %>%
  rbindlist(., idcol = "group") %>%
  filter(category %in% "Short-term-Transport") %>%
  tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
  mutate_at("p", ~ ifelse(is.na(.), 0, .)) %>%
  group_by(protein_id, group, temperature) %>%
  summarise_at("p", max) %>%
  mutate(pass = ifelse((p > .05) & (!is.na(p)), "Yes", "No")) %>%
  reshape2::dcast(., group + temperature ~ pass, value.var = "protein_id",
                  fun.aggregate = length) %>%
  mutate(n = Yes + No) %>%
  mutate(prop = paste(round(Yes/n * 100, 2), "%"))

fwrite(stat_stab1, "~/Desktop/tmp.csv")

stat_stab2 <- stab_tables %>%
  rbindlist(., idcol = "group") %>%
  filter(category %in% "Short-term-Use") %>%
  tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
  mutate_at("p", ~ ifelse(is.na(.), 0, .)) %>%
  group_by(peptide_sequence, group, temperature) %>%
  summarise_at("p", max) %>%
  mutate(pass = ifelse((p > .05) & (!is.na(p)), "Yes", "No")) %>%
  reshape2::dcast(., group + temperature ~ pass, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(n = Yes + No) %>%
  mutate(prop = paste(round(Yes/n * 100, 2), "%"))

stat_stab2 <- stab_tables %>%
  rbindlist(., idcol = "group") %>%
  filter(category %in% "Short-term-Use") %>%
  tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
  mutate_at("p", ~ ifelse(is.na(.), 0, .)) %>%
  group_by(protein_id, group, temperature) %>%
  summarise_at("p", max) %>%
  mutate(pass = ifelse((p > .05) & (!is.na(p)), "Yes", "No")) %>%
  reshape2::dcast(., group + temperature ~ pass, value.var = "protein_id",
                  fun.aggregate = length) %>%
  mutate(n = Yes + No) %>%
  mutate(prop = paste(round(Yes/n * 100, 2), "%"))

fwrite(stat_stab2, "~/Desktop/tmp.csv")

stat_stab3 <- stab_tables %>%
  rbindlist(., idcol = "group") %>%
  filter(category %in% "Long-term") %>%
  tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
  mutate_at("p", ~ ifelse(is.na(.), 0, .)) %>%
  group_by(peptide_sequence, group, temperature) %>%
  summarise_at("p", max) %>%
  mutate(pass = ifelse((p > .05) & (!is.na(p)), "Yes", "No")) %>%
  reshape2::dcast(., group ~ pass, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(n = Yes + No) %>%
  mutate(prop = paste(round(Yes/n * 100, 2), "%"))

stat_stab3 <- stab_tables %>%
  rbindlist(., idcol = "group") %>%
  filter(category %in% "Long-term") %>%
  tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
  mutate_at("p", ~ ifelse(is.na(.), 0, .)) %>%
  group_by(protein_id, group, temperature) %>%
  summarise_at("p", max) %>%
  mutate(pass = ifelse((p > .05) & (!is.na(p)), "Yes", "No")) %>%
  reshape2::dcast(., group ~ pass, value.var = "protein_id",
                  fun.aggregate = length) %>%
  mutate(n = Yes + No) %>%
  mutate(prop = paste(round(Yes/n * 100, 2), "%"))

fwrite(stat_stab3, "~/Desktop/tmp.csv")


## 定性作图 -----------------------
all_meta <- fread("./results/tables/6_outlier_madist_quartet.csv")
all_tables <- readRDS("../data/stability/quantdata_list_pep_all_longterm.rds")
stab_tables <- readRDS("./results/tables/6_qualiprop_stabtest_fisher.rds")
rm(list = setdiff(ls(), c("all_meta", "all_tables", "stab_tables")))
gc()

theoretical_list <- all_meta %>%
  group_by(sample, lab_id, time_day) %>%
  summarise(theoretical_freq = length(order), .groups = "drop")

nopass_list <- stab_tables %>%
  rbindlist(., idcol = "sample") %>%
  filter(category %in% "Long-term") %>%
  filter(n_timepoints == 25) %>%
  mutate(pass = ifelse(p > 0.05, "Yes", "No")) %>%
  tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
  distinct(sample, peptide_sequence, protein_id, pass)

nopass_tables <- pblapply(all_tables, function(tmp_table) {
  
  tmp_nopass_table <- tmp_table %>%
    inner_join(., nopass_list, by = c("sample","peptide_sequence"),
               relationship = "many-to-many") %>%
    distinct(sample, peptide_sequence, protein_id, time_day, lab_id, tube, injection, pass) %>%
    reshape2::dcast(.,
                    sample + peptide_sequence + protein_id + time_day + lab_id + pass ~.,
                    value.var = "peptide_sequence", fun.aggregate = length) %>%
    plyr::rename(c("." = "freq"))
  
  return(tmp_nopass_table)
})

df_nopass_wide <- nopass_tables %>%
  rbindlist %>%
  full_join(theoretical_list, ., by = c("sample","lab_id", "time_day"),
            relationship = "many-to-many") %>%
  mutate(missed_prop = (theoretical_freq - freq) / theoretical_freq) %>%
  reshape2::acast(., pass + sample + peptide_sequence + protein_id ~ time_day + lab_id,
                  value.var = "missed_prop", fill = 1)


## 漏检比例热图
days <- as.numeric(str_extract(colnames(df_nopass_wide), "^[0-9\\.]+"))
months <- round(days/30, digits = 2)

group_list <- ifelse(str_detect(rownames(df_nopass_wide), "Yes"), "Pass", "Fail")
names(group_list) <- rownames(df_nopass_wide)

group.colors = c("Pass" = "#276419", "Fail" = "#8E0152")
cell.colors = colorRamp2(range(df_nopass_wide), c("white", "#999"))
time.colors <- colorRampPalette(c("#DEEBF7", "#08519C"))(length(months))
names(time.colors) <- months

row_ha = rowAnnotation(Group = anno_block(gp = gpar(fill = c("#8E0152", "#276419"), col = NA),
                                          labels = c("Fail", "Pass"),
                                          labels_gp = gpar(col = "white", fontsize = 10, fontface = "bold"),
                                          width = unit(5, "mm")))

column_ha = HeatmapAnnotation(Month = months,
                              col = list(Month = time.colors),
                              annotation_name_side = "right")

ht <- Heatmap(df_nopass_wide,
              name = "Missing Prop",
              col = cell.colors,
              row_split = group_list,
              row_gap = unit(3, "mm"),
              row_title = NULL,
              cluster_columns = FALSE,
              top_annotation = column_ha,
              left_annotation = row_ha,
              show_column_names = FALSE,
              cluster_rows = TRUE,
              show_row_dend = FALSE,
              show_row_names = FALSE,
              row_dend_reorder = TRUE,
              use_raster = FALSE,
              raster_quality = 5,
              border = FALSE)

pdf("./results/figures/supp_figurexx.pdf", width = 10, height = 8)
draw(ht)
dev.off()
