## Title: Homogeneity.
## Author: Qiaochu Chen
## Date: Jun 6th, 2026

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


## 标称特性的均匀性检验: 输入数据 ---------------------
passPEPfdr0.01_pep_tables <- readRDS("./results/tables/1_qualiprop_list_PEPfdr0.01.rds")
all_meta <- fread("./data/multilab/metadata_2025_10labs.csv")


## 标称特性的均匀性检验(跨9家实验室): CFFF计算 ------------------------------
rm(list = setdiff(ls(), c("passPEPfdr0.01_pep_tables", "all_meta")))
gc()

## Fisher's Exact Test: 参考ISO 33406:2024-5.5.3.4.1
fisher_tables <- pblapply(passPEPfdr0.01_pep_tables, function(tmp_table_i) {
  
  all_meta_i <- all_meta %>%
    filter(!lab_id %in% "zhejiang_university") %>%
    filter(sample %in% unique(tmp_table_i$sample))
  
  tmp_table_i <- tmp_table_i %>%
    filter(!lab_id %in% "zhejiang_university") %>%
    mutate(feature = paste(peptide_sequence, protein_id))
  
  all_features <- unique(tmp_table_i$feature)
  
  fisher_tables_i <- mclapply(all_features, function(feature_id) {
    tryCatch({
      sub_tmp_table_j <- tmp_table_i %>%
      filter(feature %in% feature_id) %>%
      ungroup() %>%
      distinct(analysis_id, peptide_sequence) %>%
      inner_join(all_meta_i, ., by = "analysis_id") %>%
      reshape2::acast(., peptide_sequence ~ lab_id + tube,
                      value.var = "peptide_sequence",
                      fun.aggregate = length)
    
    if (nrow(sub_tmp_table_j) == 0) return(NULL)
      
    if (ncol(sub_tmp_table_j) < 11) return(NULL)
    
    if (nrow(sub_tmp_table_j) == 1) sub_tmp_table_j <- rbind(sub_tmp_table_j, 3-sub_tmp_table_j)
    
    aa <- fisher.test(sub_tmp_table_j)
    
    fisher_table_j <- data.frame(feature = feature_id,
                                 `fisher'p` = aa$p.value,
                                 n_units = ncol(sub_tmp_table_j),
                                 n_replicates = 3) %>%
      tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ")
    
    return(fisher_table_j)}, error = function(e) {
      message("Error message: ", e$message)
      stop("Stopping execution due to error at feature_id = ", feature_id)})
  }, mc.cores = 8)
  
  fisher_df_i <- rbindlist(fisher_tables_i)
  
  return(fisher_df_i)
})

bb <- fisher_tables %>% rbindlist %>% mutate(is.homo = ifelse(fisher.p > .05, TRUE, FALSE))
nrow(bb[is.homo == TRUE, ])

## 保存结果
saveRDS(fisher_tables, "./results/tables/5_qualiprop_homotest_fisher.rds")


## 标称特性的均匀性检验(跨DDA3家实验室): CFFF计算 ------------------------------
rm(list = setdiff(ls(), c("passPEPfdr0.01_pep_tables", "all_meta")))
gc()

## Fisher's Exact Test: 参考ISO 33406:2024-5.5.3.4.1
fisher_tables <- pblapply(passPEPfdr0.01_pep_tables, function(tmp_table_i) {
  
  all_meta_i <- all_meta %>%
    filter(!lab_id %in% "zhejiang_university") %>%
    filter(mode %in% "DDA") %>%
    filter(sample %in% unique(tmp_table_i$sample))
  
  tmp_table_i <- tmp_table_i %>%
    filter(!lab_id %in% "zhejiang_university") %>%
    filter(mode %in% "DDA") %>%
    mutate(feature = paste(peptide_sequence, protein_id))
  
  all_features <- unique(tmp_table_i$feature)
  
  fisher_tables_i <- mclapply(all_features, function(feature_id) {
    tryCatch({
      sub_tmp_table_j <- tmp_table_i %>%
        filter(feature %in% feature_id) %>%
        ungroup() %>%
        distinct(analysis_id, peptide_sequence) %>%
        inner_join(all_meta_i, ., by = "analysis_id") %>%
        reshape2::acast(., peptide_sequence ~ lab_id + tube,
                        value.var = "peptide_sequence",
                        fun.aggregate = length)
      
      if (nrow(sub_tmp_table_j) == 0) return(NULL)
      
      if (ncol(sub_tmp_table_j) < 9) return(NULL) #因为3家实验室最多9单元管
      
      if (nrow(sub_tmp_table_j) == 1) sub_tmp_table_j <- rbind(sub_tmp_table_j, 3-sub_tmp_table_j)
      
      aa <- fisher.test(sub_tmp_table_j)
      
      fisher_table_j <- data.frame(feature = feature_id,
                                   `fisher'p` = aa$p.value,
                                   n_units = ncol(sub_tmp_table_j),
                                   n_replicates = 3) %>%
        tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ")
      
      return(fisher_table_j)}, error = function(e) {
        message("Error message: ", e$message)
        stop("Stopping execution due to error at feature_id = ", feature_id)})
  }, mc.cores = 8)
  
  fisher_df_i <- rbindlist(fisher_tables_i)
  
  return(fisher_df_i)
})

bb <- fisher_tables %>% rbindlist %>% mutate(is.homo = ifelse(fisher.p > .05, TRUE, FALSE))
nrow(bb[is.homo == TRUE, ])

## 保存结果
saveRDS(fisher_tables, "./results/tables/5_qualiprop_homotest_fisher_dda.rds")


## 标称特性的均匀性检验(跨DIA6家实验室): CFFF计算 ------------------------------
rm(list = setdiff(ls(), c("passPEPfdr0.01_pep_tables", "all_meta")))
gc()

## Fisher's Exact Test: 参考ISO 33406:2024-5.5.3.4.1
fisher_tables <- pblapply(passPEPfdr0.01_pep_tables, function(tmp_table_i) {
  
  all_meta_i <- all_meta %>%
    filter(!lab_id %in% "zhejiang_university") %>%
    filter(mode %in% "DIA") %>%
    filter(sample %in% unique(tmp_table_i$sample))
  
  tmp_table_i <- tmp_table_i %>%
    filter(!lab_id %in% "zhejiang_university") %>%
    filter(mode %in% "DIA") %>%
    mutate(feature = paste(peptide_sequence, protein_id))
  
  all_features <- unique(tmp_table_i$feature)
  
  fisher_tables_i <- mclapply(all_features, function(feature_id) {
    tryCatch({
      sub_tmp_table_j <- tmp_table_i %>%
        filter(feature %in% feature_id) %>%
        ungroup() %>%
        distinct(analysis_id, peptide_sequence) %>%
        inner_join(all_meta_i, ., by = "analysis_id") %>%
        reshape2::acast(., peptide_sequence ~ lab_id + tube,
                        value.var = "peptide_sequence",
                        fun.aggregate = length)
      
      if (nrow(sub_tmp_table_j) == 0) return(NULL)
      
      if (ncol(sub_tmp_table_j) < 9) return(NULL) ##与DDA保持一致
      
      if (nrow(sub_tmp_table_j) == 1) sub_tmp_table_j <- rbind(sub_tmp_table_j, 3-sub_tmp_table_j)
      
      aa <- fisher.test(sub_tmp_table_j)
      
      fisher_table_j <- data.frame(feature = feature_id,
                                   `fisher'p` = aa$p.value,
                                   n_units = ncol(sub_tmp_table_j),
                                   n_replicates = 3) %>%
        tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ")
      
      return(fisher_table_j)}, error = function(e) {
        message("Error message: ", e$message)
        stop("Stopping execution due to error at feature_id = ", feature_id)})
  }, mc.cores = 8)
  
  fisher_df_i <- rbindlist(fisher_tables_i)
  
  return(fisher_df_i)
})

bb <- fisher_tables %>% rbindlist %>% mutate(is.homo = ifelse(fisher.p > .05, TRUE, FALSE))
nrow(bb[is.homo == TRUE, ])

## 保存结果
saveRDS(fisher_tables, "./results/tables/5_qualiprop_homotest_fisher_dia.rds")


## 标称特性的均匀性检验(跨Orbitrap5家实验室): CFFF计算 ------------------------------
rm(list = setdiff(ls(), c("passPEPfdr0.01_pep_tables", "all_meta")))
gc()

## Fisher's Exact Test: 参考ISO 33406:2024-5.5.3.4.1
fisher_tables <- pblapply(passPEPfdr0.01_pep_tables, function(tmp_table_i) {
  
  all_meta_i <- all_meta %>%
    filter(!lab_id %in% "zhejiang_university") %>%
    filter(lab %in% c("NCP", "QLB", "BTP", "TFS", "NIM")) %>%
    filter(sample %in% unique(tmp_table_i$sample))
  
  tmp_table_i <- tmp_table_i %>%
    filter(!lab_id %in% "zhejiang_university") %>%
    filter(lab %in% c("NCP", "QLB", "BTP", "TFS", "NIM")) %>%
    mutate(feature = paste(peptide_sequence, protein_id))
  
  all_features <- unique(tmp_table_i$feature)
  
  fisher_tables_i <- mclapply(all_features, function(feature_id) {
    tryCatch({
      sub_tmp_table_j <- tmp_table_i %>%
        filter(feature %in% feature_id) %>%
        ungroup() %>%
        distinct(analysis_id, peptide_sequence) %>%
        inner_join(all_meta_i, ., by = "analysis_id") %>%
        reshape2::acast(., peptide_sequence ~ lab_id + tube,
                        value.var = "peptide_sequence",
                        fun.aggregate = length)
      
      if (nrow(sub_tmp_table_j) == 0) return(NULL)
      
      if (ncol(sub_tmp_table_j) < 11) return(NULL)
      
      if (nrow(sub_tmp_table_j) == 1) sub_tmp_table_j <- rbind(sub_tmp_table_j, 3-sub_tmp_table_j)
      
      aa <- fisher.test(sub_tmp_table_j)
      
      fisher_table_j <- data.frame(feature = feature_id,
                                   `fisher'p` = aa$p.value,
                                   n_units = ncol(sub_tmp_table_j),
                                   n_replicates = 3) %>%
        tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ")
      
      return(fisher_table_j)}, error = function(e) {
        message("Error message: ", e$message)
        stop("Stopping execution due to error at feature_id = ", feature_id)})
  }, mc.cores = 8)
  
  fisher_df_i <- rbindlist(fisher_tables_i)
  
  return(fisher_df_i)
})

bb <- fisher_tables %>% rbindlist %>% mutate(is.homo = ifelse(fisher.p > .05, TRUE, FALSE))
nrow(bb[is.homo == TRUE, ])

## 保存结果
saveRDS(fisher_tables, "./results/tables/5_qualiprop_homotest_fisher_orbitrap.rds")


## 标称特性的均匀性检验(跨TOF4家实验室): CFFF计算 ------------------------------
rm(list = setdiff(ls(), c("passPEPfdr0.01_pep_tables", "all_meta")))
gc()

## Fisher's Exact Test: 参考ISO 33406:2024-5.5.3.4.1
fisher_tables <- pblapply(passPEPfdr0.01_pep_tables, function(tmp_table_i) {
  
  all_meta_i <- all_meta %>%
    filter(!lab_id %in% "zhejiang_university") %>%
    filter(lab %in% c("FDU", "OSB", "CAS", "CMS")) %>%
    filter(sample %in% unique(tmp_table_i$sample))
  
  tmp_table_i <- tmp_table_i %>%
    filter(!lab_id %in% "zhejiang_university") %>%
    filter(lab %in% c("FDU", "OSB", "CAS", "CMS")) %>%
    mutate(feature = paste(peptide_sequence, protein_id))
  
  all_features <- unique(tmp_table_i$feature)
  
  fisher_tables_i <- mclapply(all_features, function(feature_id) {
    tryCatch({
      sub_tmp_table_j <- tmp_table_i %>%
        filter(feature %in% feature_id) %>%
        ungroup() %>%
        distinct(analysis_id, peptide_sequence) %>%
        inner_join(all_meta_i, ., by = "analysis_id") %>%
        reshape2::acast(., peptide_sequence ~ lab_id + tube,
                        value.var = "peptide_sequence",
                        fun.aggregate = length)
      
      if (nrow(sub_tmp_table_j) == 0) return(NULL)
      
      if (ncol(sub_tmp_table_j) < 11) return(NULL)
      
      if (nrow(sub_tmp_table_j) == 1) sub_tmp_table_j <- rbind(sub_tmp_table_j, 3-sub_tmp_table_j)
      
      aa <- fisher.test(sub_tmp_table_j)
      
      fisher_table_j <- data.frame(feature = feature_id,
                                   `fisher'p` = aa$p.value,
                                   n_units = ncol(sub_tmp_table_j),
                                   n_replicates = 3) %>%
        tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ")
      
      return(fisher_table_j)}, error = function(e) {
        message("Error message: ", e$message)
        stop("Stopping execution due to error at feature_id = ", feature_id)})
  }, mc.cores = 8)
  
  fisher_df_i <- rbindlist(fisher_tables_i)
  
  return(fisher_df_i)
})

bb <- fisher_tables %>% rbindlist %>% mutate(is.homo = ifelse(fisher.p > .05, TRUE, FALSE))
nrow(bb[is.homo == TRUE, ])

## 保存结果
saveRDS(fisher_tables, "./results/tables/5_qualiprop_homotest_fisher_tof.rds")


## 参考量值的均匀性检验: 输入数据 --------------------------
all_meta <- fread("./data/multilab/metadata_2025_10labs.csv")

passPfdr0.05_pep_tables <- readRDS("./results/tables/2_quantprop_list_Pfdr0.05.rds")
passPfdr0.05_pep_tables <- passPfdr0.05_pep_tables[c(1:2, 4, 3)]

ratiobyd6_tables <- readRDS("./results/tables/2_quantdata_list_pep_ratiobyd6_2025_10labs.rds")
ratiobyd6_tables <- ratiobyd6_tables[c(1, 3:4)]
names(ratiobyd6_tables) <- names(passPfdr0.05_pep_tables)[1:3]

ratiobyhek293t_tables <- readRDS("./results/tables/2_quantdata_list_pep_ratiobyhek293t_2025_10labs.rds")
ratiobyhek293t_tables <- ratiobyhek293t_tables[1]
names(ratiobyhek293t_tables) <- names(passPfdr0.05_pep_tables)[4]

all_tables <- c(ratiobyd6_tables, ratiobyhek293t_tables)
quantprop_tables <- pblapply(1:4, function(i) {
  
  quantprop_tables_i <- passPfdr0.05_pep_tables[[i]] %>%
    distinct(peptide_sequence, protein_id, lab_id, tube) %>%
    mutate_at("tube", as.numeric) %>%
    inner_join(., all_tables[[i]], by =c("lab_id", "tube", "peptide_sequence"),
               relationship = "many-to-many")
  
  return(quantprop_tables_i)
})
names(quantprop_tables) <- names(all_tables)


## 参考量值的均匀性检验(跨9家实验室): ------------------------------
rm(list = setdiff(ls(), c("quantprop_tables")))

## 限制性极大似然估计: 参考ISO 33405:2024-7.7.6和JJF 1343-2022-5.9.7
reml_tables <- pblapply(quantprop_tables, function(tmp_table_i) {
  
  tmp_table_i <- tmp_table_i %>%
    mutate(feature = paste(peptide_sequence, protein_id))
  
  all_features <- unique(tmp_table_i$feature)
  
  reml_tables_i <- mclapply(all_features, function(feature_id) {
    # print(feature_id)
    sub_expr_tmp <- tmp_table_i %>%
      filter(feature %in% feature_id) %>%
      ungroup()
    
    if (nrow(sub_expr_tmp) == 0) return(NULL)
    
    if (length(unique(sub_expr_tmp$lab_id)) == 1) return(NULL)
    
    n_total <- nrow(sub_expr_tmp)
    n_unit <- sub_expr_tmp %>% distinct(lab_id, tube) %>% nrow
    
    if (n_unit < 11) return(NULL)
    
    ## 混合线性模型估计方差
    lm_results1 <- lmer(value ~ 1 + (1 | lab_id) + (1 | lab_id:tube), data = sub_expr_tmp, REML = TRUE)
    lm_results2 <- lmer(value ~ 1 + (1 | lab_id), data = sub_expr_tmp, REML = TRUE)
    anova_results <- anova(lm_results1, lm_results2, refit = FALSE) ## 禁用refitting model(s) with ML (instead of REML)
    var_results1 <- summary(lm_results1)
    
    ## 单元间与单元内不确定度
    homo_results_j <- data.frame(feature = feature_id) %>%
      tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
      mutate(`Chisq'p` = anova_results$`Pr(>Chisq)`[2],
             method = "REML") %>%
      mutate(is.homogeneous = ifelse(`Chisq'p` > .05, TRUE, FALSE)) %>%
      mutate(s_bu = sqrt(var_results1$varcor$`lab_id:tube`[1]),
             s_wu = var_results1$sigma,
             u_bu = ifelse(s_bu > s_wu,
                           sqrt((s_bu ^ 2 - s_wu ^ 2) / 3),
                           sqrt(s_wu ^ 2 / 3) * (2 / (n_total - n_unit)) ^ 0.25),
             u_wu = sqrt(s_wu ^ 2 / (n_total - n_unit)),
             u_homo = sqrt(u_bu ^ 2 + u_wu ^ 2))

    return(homo_results_j)
  }, mc.cores = 4)
  
  homo_results_i <- rbindlist(reml_tables_i)
  
  return(homo_results_i)
})

saveRDS(reml_tables, "./results/tables/5_quantprop_homotest_reml.rds")


## 参考量值的均匀性检验(跨DDA3家实验室): ------------------------------
rm(list = setdiff(ls(), c("quantprop_tables")))

## 限制性极大似然估计: 参考ISO 33405:2024-7.7.6和JJF 1343-2022-5.9.7
reml_tables <- pblapply(quantprop_tables, function(tmp_table_i) {
  
  tmp_table_i <- tmp_table_i %>%
    filter(lab_id %in% c("qinglian_bio", "phoenix", "fudan_university")) %>%
    mutate(feature = paste(peptide_sequence, protein_id))
  
  all_features <- unique(tmp_table_i$feature)
  
  reml_tables_i <- mclapply(all_features, function(feature_id) {
    # print(feature_id)
    sub_expr_tmp <- tmp_table_i %>%
      filter(feature %in% feature_id) %>%
      ungroup()
    
    if (nrow(sub_expr_tmp) == 0) return(NULL)
    
    if (length(unique(sub_expr_tmp$lab_id)) == 1) return(NULL)
    
    n_total <- nrow(sub_expr_tmp)
    n_unit <- sub_expr_tmp %>% distinct(lab_id, tube) %>% nrow
    
    if (n_unit < 4) return(NULL) ## 3家实验室最多9单元管;2026-06-22修改按照11/27比例计算
    
    ## 混合线性模型估计方差
    lm_results1 <- lmer(value ~ 1 + (1 | lab_id) + (1 | lab_id:tube), data = sub_expr_tmp, REML = TRUE)
    lm_results2 <- lmer(value ~ 1 + (1 | lab_id), data = sub_expr_tmp, REML = TRUE)
    anova_results <- anova(lm_results1, lm_results2, refit = FALSE) ## 禁用refitting model(s) with ML (instead of REML)
    var_results1 <- summary(lm_results1)
    
    ## 单元间与单元内不确定度
    homo_results_j <- data.frame(feature = feature_id) %>%
      tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
      mutate(`Chisq'p` = anova_results$`Pr(>Chisq)`[2],
             method = "REML") %>%
      mutate(is.homogeneous = ifelse(`Chisq'p` > .05, TRUE, FALSE)) %>%
      mutate(s_bu = sqrt(var_results1$varcor$`lab_id:tube`[1]),
             s_wu = var_results1$sigma,
             u_bu = ifelse(s_bu > s_wu,
                           sqrt((s_bu ^ 2 - s_wu ^ 2) / 3),
                           sqrt(s_wu ^ 2 / 3) * (2 / (n_total - n_unit)) ^ 0.25),
             u_wu = sqrt(s_wu ^ 2 / (n_total - n_unit)),
             u_homo = sqrt(u_bu ^ 2 + u_wu ^ 2))
    
    return(homo_results_j)
  }, mc.cores = 4)
  
  homo_results_i <- rbindlist(reml_tables_i)
  
  return(homo_results_i)
})

saveRDS(reml_tables, "./results/tables/5_quantprop_homotest_reml_dda.rds")


## 参考量值的均匀性检验(跨DIA6家实验室): ------------------------------
rm(list = setdiff(ls(), c("quantprop_tables")))

## 限制性极大似然估计: 参考ISO 33405:2024-7.7.6和JJF 1343-2022-5.9.7
reml_tables <- pblapply(quantprop_tables, function(tmp_table_i) {
  
  tmp_table_i <- tmp_table_i %>%
    filter(lab_id %in% c("thermofisher_shanghai", "omicsolution", "cas_tianjin",
                         "biotech_pack", "academy_of_chinese_medical_sciences",
                         "national_institute_of_methodology")) %>%
    mutate(feature = paste(peptide_sequence, protein_id))
  
  all_features <- unique(tmp_table_i$feature)
  
  reml_tables_i <- mclapply(all_features, function(feature_id) {
    # print(feature_id)
    sub_expr_tmp <- tmp_table_i %>%
      filter(feature %in% feature_id) %>%
      ungroup()
    
    if (nrow(sub_expr_tmp) == 0) return(NULL)
    
    if (length(unique(sub_expr_tmp$lab_id)) == 1) return(NULL)
    
    n_total <- nrow(sub_expr_tmp)
    n_unit <- sub_expr_tmp %>% distinct(lab_id, tube) %>% nrow
    
    if (n_unit < 7) return(NULL) ## 2026-06-22修改按照11/27比例计算
    
    ## 混合线性模型估计方差
    lm_results1 <- lmer(value ~ 1 + (1 | lab_id) + (1 | lab_id:tube), data = sub_expr_tmp, REML = TRUE)
    lm_results2 <- lmer(value ~ 1 + (1 | lab_id), data = sub_expr_tmp, REML = TRUE)
    anova_results <- anova(lm_results1, lm_results2, refit = FALSE) ## 禁用refitting model(s) with ML (instead of REML)
    var_results1 <- summary(lm_results1)
    
    ## 单元间与单元内不确定度
    homo_results_j <- data.frame(feature = feature_id) %>%
      tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
      mutate(`Chisq'p` = anova_results$`Pr(>Chisq)`[2],
             method = "REML") %>%
      mutate(is.homogeneous = ifelse(`Chisq'p` > .05, TRUE, FALSE)) %>%
      mutate(s_bu = sqrt(var_results1$varcor$`lab_id:tube`[1]),
             s_wu = var_results1$sigma,
             u_bu = ifelse(s_bu > s_wu,
                           sqrt((s_bu ^ 2 - s_wu ^ 2) / 3),
                           sqrt(s_wu ^ 2 / 3) * (2 / (n_total - n_unit)) ^ 0.25),
             u_wu = sqrt(s_wu ^ 2 / (n_total - n_unit)),
             u_homo = sqrt(u_bu ^ 2 + u_wu ^ 2))
    
    return(homo_results_j)
  }, mc.cores = 4)
  
  homo_results_i <- rbindlist(reml_tables_i)
  
  return(homo_results_i)
})

saveRDS(reml_tables, "./results/tables/5_quantprop_homotest_reml_dia.rds")


## 参考量值的均匀性检验(跨Orbitrap5家实验室): ------------------------------
rm(list = setdiff(ls(), c("quantprop_tables")))

## 限制性极大似然估计: 参考ISO 33405:2024-7.7.6和JJF 1343-2022-5.9.7
reml_tables <- pblapply(quantprop_tables, function(tmp_table_i) {
  
  tmp_table_i <- tmp_table_i %>%
    filter(lab_id %in% c("thermofisher_shanghai", "biotech_pack", "qinglian_bio",
                         "phoenix", "national_institute_of_methodology")) %>%
    mutate(feature = paste(peptide_sequence, protein_id))
  
  all_features <- unique(tmp_table_i$feature)
  
  reml_tables_i <- mclapply(all_features, function(feature_id) {
    # print(feature_id)
    sub_expr_tmp <- tmp_table_i %>%
      filter(feature %in% feature_id) %>%
      ungroup()
    
    if (nrow(sub_expr_tmp) == 0) return(NULL)
    
    if (length(unique(sub_expr_tmp$lab_id)) == 1) return(NULL)
    
    n_total <- nrow(sub_expr_tmp)
    n_unit <- sub_expr_tmp %>% distinct(lab_id, tube) %>% nrow
    
    if (n_unit < 6) return(NULL) ## 2026-06-22修改按照11/27比例计算
    
    ## 混合线性模型估计方差
    lm_results1 <- lmer(value ~ 1 + (1 | lab_id) + (1 | lab_id:tube), data = sub_expr_tmp, REML = TRUE)
    lm_results2 <- lmer(value ~ 1 + (1 | lab_id), data = sub_expr_tmp, REML = TRUE)
    anova_results <- anova(lm_results1, lm_results2, refit = FALSE) ## 禁用refitting model(s) with ML (instead of REML)
    var_results1 <- summary(lm_results1)
    
    ## 单元间与单元内不确定度
    homo_results_j <- data.frame(feature = feature_id) %>%
      tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
      mutate(`Chisq'p` = anova_results$`Pr(>Chisq)`[2],
             method = "REML") %>%
      mutate(is.homogeneous = ifelse(`Chisq'p` > .05, TRUE, FALSE)) %>%
      mutate(s_bu = sqrt(var_results1$varcor$`lab_id:tube`[1]),
             s_wu = var_results1$sigma,
             u_bu = ifelse(s_bu > s_wu,
                           sqrt((s_bu ^ 2 - s_wu ^ 2) / 3),
                           sqrt(s_wu ^ 2 / 3) * (2 / (n_total - n_unit)) ^ 0.25),
             u_wu = sqrt(s_wu ^ 2 / (n_total - n_unit)),
             u_homo = sqrt(u_bu ^ 2 + u_wu ^ 2))
    
    return(homo_results_j)
  }, mc.cores = 4)
  
  homo_results_i <- rbindlist(reml_tables_i)
  
  return(homo_results_i)
})

saveRDS(reml_tables, "./results/tables/5_quantprop_homotest_reml_orbitrap.rds")


## 参考量值的均匀性检验(跨TOF4家实验室): ------------------------------
rm(list = setdiff(ls(), c("quantprop_tables")))

## 限制性极大似然估计: 参考ISO 33405:2024-7.7.6和JJF 1343-2022-5.9.7
reml_tables <- pblapply(quantprop_tables, function(tmp_table_i) {
  
  tmp_table_i <- tmp_table_i %>%
    filter(lab_id %in% c("omicsolution", "cas_tianjin", "fudan_university",
                         "academy_of_chinese_medical_sciences")) %>%
    mutate(feature = paste(peptide_sequence, protein_id))
  
  all_features <- unique(tmp_table_i$feature)
  
  reml_tables_i <- mclapply(all_features, function(feature_id) {
    # print(feature_id)
    sub_expr_tmp <- tmp_table_i %>%
      filter(feature %in% feature_id) %>%
      ungroup()
    
    if (nrow(sub_expr_tmp) == 0) return(NULL)
    
    if (length(unique(sub_expr_tmp$lab_id)) == 1) return(NULL)
    
    n_total <- nrow(sub_expr_tmp)
    n_unit <- sub_expr_tmp %>% distinct(lab_id, tube) %>% nrow
    
    if (n_unit < 5) return(NULL) ## 2026-06-22修改按照11/27比例计算
    
    ## 混合线性模型估计方差
    lm_results1 <- lmer(value ~ 1 + (1 | lab_id) + (1 | lab_id:tube), data = sub_expr_tmp, REML = TRUE)
    lm_results2 <- lmer(value ~ 1 + (1 | lab_id), data = sub_expr_tmp, REML = TRUE)
    anova_results <- anova(lm_results1, lm_results2, refit = FALSE) ## 禁用refitting model(s) with ML (instead of REML)
    var_results1 <- summary(lm_results1)
    
    ## 单元间与单元内不确定度
    homo_results_j <- data.frame(feature = feature_id) %>%
      tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
      mutate(`Chisq'p` = anova_results$`Pr(>Chisq)`[2],
             method = "REML") %>%
      mutate(is.homogeneous = ifelse(`Chisq'p` > .05, TRUE, FALSE)) %>%
      mutate(s_bu = sqrt(var_results1$varcor$`lab_id:tube`[1]),
             s_wu = var_results1$sigma,
             u_bu = ifelse(s_bu > s_wu,
                           sqrt((s_bu ^ 2 - s_wu ^ 2) / 3),
                           sqrt(s_wu ^ 2 / 3) * (2 / (n_total - n_unit)) ^ 0.25),
             u_wu = sqrt(s_wu ^ 2 / (n_total - n_unit)),
             u_homo = sqrt(u_bu ^ 2 + u_wu ^ 2))
    
    return(homo_results_j)
  }, mc.cores = 4)
  
  homo_results_i <- rbindlist(reml_tables_i)
  
  return(homo_results_i)
})

saveRDS(reml_tables, "./results/tables/5_quantprop_homotest_reml_tof.rds")


## test1 -----------------------
homo_tables <- readRDS("./results/tables/5_quantprop_homotest_reml_tof.rds")
rm(list = setdiff(ls(), c("homo_tables")))
gc()

aa <- homo_tables %>%
  rbindlist(., idcol = "group") %>%
  dplyr::rename(p = `Chisq'p`) %>%
  mutate(pass = ifelse((p > .05) & (!is.na(p)), "Yes", "No")) %>%
  distinct(peptide_sequence, protein_id, group, pass) %>%
  reshape2::dcast(., peptide_sequence + protein_id + group ~ pass, value.var = "pass") %>%
  na.omit ##不含缺失值说明不存在同一肽段yes/no的矛盾，否则考虑是否同一肽段匹配到了不同蛋白质，因而计算的P不同。

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
  mutate(prop = paste(round(Yes/n * 100, 2), "%"))

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
  mutate(prop = paste(round(Yes/n * 100, 2), "%"))


## test2 -----------------------
passPEPfdr0.01_pep_tables <- readRDS("./results/tables/1_qualiprop_list_PEPfdr0.01.rds")
homo_tables <- readRDS("./results/tables/5_qualiprop_homotest_fisher_tof.rds")

diff_peptides <- pblapply(1:6, function(i) {
  
  df_qualiprop <- passPEPfdr0.01_pep_tables[[i]]
  df_homo <- homo_tables[[i]]
  
  setdiff(df_qualiprop$peptide_sequence, df_homo$peptide_sequence)
  
  df_qualiprop_tmp <- df_qualiprop %>%
    filter(peptide_sequence %in% "KLFVGGLK")
})

rm(list = setdiff(ls(), c("homo_tables", "passPEPfdr0.01_pep_tables")))
gc()

aa <- homo_tables %>%
  rbindlist(., idcol = "group") %>%
  dplyr::rename(p = `fisher.p`) %>%
  mutate(pass = ifelse((p > .05) & (!is.na(p)), "Yes", "No")) %>%
  distinct(peptide_sequence, protein_id, group, pass) %>%
  reshape2::dcast(., peptide_sequence + protein_id + group ~ pass, value.var = "pass") %>%
  na.omit ##不含缺失值说明不存在同一肽段yes/no的矛盾，否则考虑是否同一肽段匹配到了不同蛋白质，因而计算的P不同。

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
  mutate(prop = paste(round(Yes/n * 100, 2), "%"))

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
  mutate(prop = paste(round(Yes/n * 100, 2), "%"))

fwrite(stat_sw, "~/Desktop/tmp.csv")


## 作图 -----------------------
homo_tables <- readRDS("./results/tables/5_quantprop_homotest_reml.rds")

ratioby293t_tables <- readRDS("./results/tables/2_quantdata_list_pep_ratiobyhek293t_2025_10labs.rds")
ratioby293t_tables <- ratioby293t_tables[1]
names(ratioby293t_tables) <- names(homo_tables)[4]

ratiobyd6_tables <- readRDS("./results/tables/2_quantdata_list_pep_ratiobyd6_2025_10labs.rds")
ratiobyd6_tables <- ratiobyd6_tables[c(1, 3:4)]
names(ratiobyd6_tables) <- names(homo_tables)[1:3]

df_homo <- homo_tables[1:3] %>%
  rbindlist(., idcol = "group")

## 挑选组织特异性低（所有细胞均表达）的蛋白质: P29966/P02545/Q9UHB9/P42704
## 挑选淋巴组织高表达的蛋白质: P50453/P80723/P20592/B0I1T2
## 挑选第7部分认定值表得到的值: AQLLELPYAR-P50453
## 学位论文修改；根据最终定值结果挑选特征
df_final <- fread("./results/tables/7_quantprop_final.csv")
vip_peptides <- df_final %>%
  reshape2::dcast(., peptide_sequence + protein_id ~ group, value.var = "value") %>%
  na.omit

vip_peptides <- vip_peptides[1:3, ]

df_ratio_vip <- ratiobyd6_tables %>%
  rbindlist(., idcol = "group") %>%
  filter(!lab_id %in% "zhejiang_university") %>%
  inner_join(., vip_peptides, by = "peptide_sequence", relationship = "many-to-many")

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
  mutate_at("group", ~ factor(., levels = c("HeLa/HEK293T", "D5/D6", "F7/D6", "M8/D6")))

df_homo_tmp <- df_homo %>%
  inner_join(., vip_peptides, by = c("peptide_sequence", "protein_id"))

colors.class <- c("#54278F", "#9E9AC8")
names(colors.class) <- c("Inter-vial", "Intra-vial")
colors.group <- c("#4CC3D9", "#FFC65D", "#F16745", "#E7298A")
names(colors.group) <- c("D5/D6", "F7/D6", "M8/D6", "HeLa/HEK293T")

set.seed(1)
p <- ggplot() +
  geom_boxplot(aes(fill = group, alpha = class, x = class, y = value),
               data = df_ratio_vip_final, outliers = FALSE) +
  geom_jitter(aes(fill = group, alpha = class, x = class, y = value),
              data = df_ratio_vip_final, width = 0.3, shape = 21, color = "black") +
  geom_text(aes(x = "Inter-vial", y = Inf, label = sprintf("Chisq'p (REML) = %.2f", `Chisq'p`)),
            data = df_homo_tmp, size = 4, vjust = 1.5) +
  facet_grid(peptide_sequence + protein_id ~ group, scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "none",
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 16)) +
  scale_y_continuous(name = "Residual (log2 transformed)",
                     expand = expansion(mult = c(0.05, 0.15))) +
  scale_alpha_manual(values = c(1, .6)) +
  scale_fill_manual(values = colors.group);p

ggsave("./chapter3_3_2.pdf", p, height = 6, width = 12)


