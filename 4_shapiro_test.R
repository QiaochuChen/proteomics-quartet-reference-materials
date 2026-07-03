## Title: normality.
## Author: Qiaochu Chen
## Date: Jun 6th, 2026

library(readxl)
library(data.table)
library(rstatix)
library(pbapply)
library(parallel)
library(dplyr)
library(stringr)

rm(list = ls())
gc()


## 读取肽段参考量值数据 ------------------------------------
passPfdr0.05_pep_tables <- readRDS("./results/tables/2_quantprop_list_Pfdr0.05.rds")
meta <- fread("./results/tables/2_outlier_madist_quartet.csv")
meta <- meta %>% filter(!is.outlier)

fwrite(meta %>% filter(is.outlier), "./meta_tmp.csv")


## Shapiro-Wilk tests: --------------------
norm_tables <- pblapply(passPfdr0.05_pep_tables, function(tmp_table) {
  
  sub_limma <- tmp_table %>%
    reshape2::dcast(., peptide_sequence + protein_id ~ lab_id + tube, value.var = "logFC")
  
  sub_expr_j <- sub_limma %>%
    filter(apply(., 1, function(x) sum(!is.na(x[-(1:3)])) >= 3)) %>%## 注意是因为前三列是特征和样本对信息
    select(peptide_sequence, protein_id, where(is.numeric))
  
  n_features <- nrow(sub_expr_j)
  
  norm_tables_j <- mclapply(1:n_features, function(k) {
    
    peptide_id <- sub_expr_j$peptide_sequence[k]
    protein_name <- sub_expr_j$protein_id[k]
    
    n_total <- sub_expr_j[k, ] %>%
      select_if(is.numeric) %>%
      t %>%
      na.omit %>%
      length
    
    sub_expr_k <- sub_expr_j[k, ] %>%
      reshape2::melt(., id = 1:2, variable.name = "id", na.rm = TRUE)
    
    # ## Box-Cox
    # sub_expr_k <- sub_expr_k %>%
    #   mutate_at("value", ~ 2 ^ (.))
    # 
    # results <- MASS::boxcox(sub_expr_k$value ~ 1)
    # lambda <- results$x[which.max(results$y)]
    # 
    # sub_expr_k <- sub_expr_k %>%
    #   mutate_at("value", ~ ((.) ^ lambda - 1) / lambda)
    
    ## remove outliers: MAD
    sub_expr_k <- sub_expr_k %>%
      mutate_at("value", ~ 2 ^ (.)) %>%
      mutate(median = median(value)) %>%
      mutate(deviation = abs(value - median)) %>%
      mutate(mad = median(deviation)) %>%
      # mutate_at("value", ~ ifelse(deviation <= 2.5 * mad, ., median + 2.5 * mad)) %>%
      filter(deviation <= 2.5 * mad) %>%
      select(peptide_sequence, id, value, everything())
    
    n_used <- nrow(sub_expr_k)
    
    if (n_used >= 3) {
      sub_sw_k <- sub_expr_k %>%
        shapiro_test(value) %>%
        select(!variable) %>%
        mutate(peptide_sequence = peptide_id, 
               protein_id = protein_name) %>%
        mutate(n_total = n_total, n_used = n_used)
      
    } else {
      sub_sw_k <- data.frame(statistic = NA, p = NA) %>%
        mutate(peptide_sequence = peptide_id, 
               protein_id = protein_name) %>%
        mutate(n_total = n_total, n_used = n_used)
      
    }
    
    y <- sub_expr_k$value
    
    sub_sw_k <- sub_sw_k %>%
      mutate(FC_mean = mean(y)) %>%
      mutate(FC_median = median(y)) %>%
      mutate(uchar = sqrt(sum((y - mean(y)) ^ 2)/(n_used * (n_used - 1)))) ## 根据JJF 1343-2022 7.9表3确定er=1

    return(sub_sw_k)
    
  }, mc.cores = 4)
  
  df_sw_j <- rbindlist(norm_tables_j)

  return(df_sw_j)
  
})

saveRDS(norm_tables, "./results/tables/4_quantprop_normtest_shapirowilk.rds")


## test1 -----------------------
norm_tables <- readRDS("./results/quartet/tables/4_quantprop_normtest_shapirowilk.rds")
rm(list = setdiff(ls(), c("norm_tables")))
gc()

aa <- norm_tables %>%
  rbindlist(., idcol = "group") %>%
  mutate(pass = ifelse((p > .05) & (!is.na(p)), "Yes", "No")) %>%
  distinct(peptide_sequence, protein_id, group, pass) %>%
  reshape2::dcast(., peptide_sequence + protein_id + group ~ pass, value.var = "pass") %>%
  na.omit ##不含缺失值说明不存在同一肽段yes/no的矛盾，否则考虑是否同一肽段匹配到了不同蛋白质，因而计算的P不同。

stat_missing <- norm_tables %>%
  rbindlist(., idcol = "group") %>%
  mutate(n_miss = n_total - n_used) %>%
  reshape2::dcast(., peptide_sequence + protein_id ~ group, value.var = "n_miss")

stat_sw <- norm_tables %>%
  rbindlist(., idcol = "group") %>%
  mutate_at("p", ~ ifelse(is.na(.), 0, .)) %>%
  group_by(peptide_sequence, group) %>%
  summarise_at("p", max) %>%
  mutate(pass = ifelse((p > .05) & (!is.na(p)), "Yes", "No")) %>%
  reshape2::dcast(., group ~ pass, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(n = Yes + No) %>%
  mutate(prop = paste(round(Yes/n * 100, 2), "%"))

fwrite(stat_sw, "~/Desktop/tmp.csv")

stat_sw <- norm_tables %>%
  rbindlist(., idcol = "group") %>%
  mutate_at("p", ~ ifelse(is.na(.), 0, .)) %>%
  group_by(protein_id, group) %>%
  summarise_at("p", max) %>%
  mutate(pass = ifelse((p > .05) & (!is.na(p)), "Yes", "No")) %>%
  reshape2::dcast(., group ~ pass, value.var = "protein_id",
                  fun.aggregate = length) %>%
  mutate(n = Yes + No) %>%
  mutate(prop = paste(round(Yes/n * 100, 2), "%"))

fwrite(stat_sw, "~/Desktop/tmp.csv")


