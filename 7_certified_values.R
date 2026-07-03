## Title: Certified values in the reference data set.
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
library(Biostrings)

rm(list = ls())
options(digits = 4)

db_uniprot <- readAAStringSet("./uniprotkb_proteome_UP000005640_2024_09_03.fasta")
db_uniprot <- data.frame(protein_sequence = as.character(db_uniprot)) %>%
  tibble::rownames_to_column("entry") %>%
  mutate(protein_id = str_extract(entry, "(?<=\\|).+(?=\\|)")) %>%
  mutate(protein_fullname = str_extract(entry, "(?<=_HUMAN ).+(?= OS\\=Homo sapiens )")) %>%
  mutate(gene_symbol = str_extract(entry, ifelse(grepl("PE=", entry),
                                                 "(?<=GN\\=).+(?= PE=)",
                                                 "(?<=GN\\=).+"))) %>%
  distinct(protein_id, protein_fullname, gene_symbol)

## 肽段标称特性的不确定度 ---------------------
qualiprop_tables <- readRDS("./results/tables/1_qualiprop_list_PEPfdr0.01.rds")
qualiprop_tables <- qualiprop_tables[c(1:3, 6)]
stat_qualiprops <- qualiprop_tables %>% rbindlist %>% pull(peptide_sequence) %>% unique %>% length
stat_qualiprops <- qualiprop_tables %>% rbindlist %>% pull(protein_id) %>% unique %>% length

homo_tables <- readRDS("./results/tables/5_qualiprop_homotest_fisher.rds")
homo_tables <- homo_tables[c(1:3, 6)]
stab_tables <- readRDS("./results/tables/6_qualiprop_stabtest_fisher.rds")

ucert_tables <- pblapply(1:4, function(i) {
  
  sub_quantprop <- qualiprop_tables[[i]] %>%
    tidyr::separate_rows(gene_symbol, sep = ";") %>%
    group_by(peptide_sequence, protein_id) %>%
    summarise(fdr = median(fdr), .groups = "drop")
  
  sub_homo <- homo_tables[[i]] %>%
    filter(fisher.p > .05) %>%
    dplyr::rename(p_homo = fisher.p) %>%
    inner_join(., sub_quantprop, by = c("peptide_sequence", "protein_id"))
  
  sub_stab1 <- stab_tables[[i]] %>%
    tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
    filter(category %in% "Short-term-Transport", temperature %in% "-78.5°C", period %in% ("0, 3, 7, 14")) %>%
    filter(p > .05) %>%
    dplyr::rename(p_stst = p, period_stst = period) %>%
    select(peptide_sequence, protein_id, p_stst, period_stst) %>%
    inner_join(., sub_homo, by = c("peptide_sequence", "protein_id"))
  
  sub_stab2 <- stab_tables[[i]] %>%
    tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
    filter(category %in% "Short-term-Use", temperature %in% "4°C", period %in% ("0, 1, 2, 4, 7")) %>%
    filter(p > .05) %>%
    dplyr::rename(p_stsu = p, period_stsu = period) %>%
    select(peptide_sequence, protein_id, p_stsu, period_stsu) %>%
    inner_join(., sub_stab1, by = c("peptide_sequence", "protein_id"))
  
  sub_stab3 <- stab_tables[[i]] %>%
    tidyr::separate(feature, c("peptide_sequence", "protein_id"), sep = " ") %>%
    filter(category %in% "Long-term") %>%
    filter(p > .05) %>%
    dplyr::rename(p_lts = p, period_lts = period) %>%
    select(peptide_sequence, protein_id, p_lts, period_lts) %>%
    inner_join(., sub_stab2, by = c("peptide_sequence", "protein_id"))
  
  sub_final <- sub_stab3 %>%
    select(peptide_sequence, protein_id, fdr, everything())
  
  stat_all1 <- data.frame(tier4 = length(unique(sub_homo$peptide_sequence)),
                          tier5 = length(unique(sub_stab1$peptide_sequence)),
                          tier6 = length(unique(sub_stab2$peptide_sequence)),
                          tier7 = length(unique(sub_stab3$peptide_sequence)))
  stat_all2 <- data.frame(tier4 = length(unique(sub_homo$protein_id)),
                          tier5 = length(unique(sub_stab1$protein_id)),
                          tier6 = length(unique(sub_stab2$protein_id)),
                          tier7 = length(unique(sub_stab3$protein_id)))
  
  return(sub_final)
})
names(ucert_tables) <- names(qualiprop_tables)

ucert_df <- ucert_tables %>%
  rbindlist(., idcol = "sample") %>%
  left_join(., db_uniprot, by = c("protein_id")) %>%
  select(sample, peptide_sequence, protein_id, protein_fullname, gene_symbol, everything())

fwrite(ucert_df, "./results/tables/7_qualiprop_final.csv")


## 肽段参考量值的不确定度 ---------------------
norm_tables <- readRDS("./results/tables/4_quantprop_normtest_shapirowilk.rds")
norm_tables <- norm_tables[c(1:2, 4)]
stat_quantprops <- norm_tables %>% rbindlist %>% pull(peptide_sequence) %>% unique %>% length
stat_quantprops <- norm_tables %>% rbindlist %>% pull(protein_id) %>% unique %>% length

homo_tables <- readRDS("./results/tables/5_quantprop_homotest_reml.rds")
homo_tables <- homo_tables[1:3]
stab_tables <- readRDS("./results/tables/6_quantprop_stabtest_linear.rds")

df_uchar <- norm_tables %>%
  rbindlist(., idcol = "group")

## 定值引入的不确定度
stat_u <- df_uchar %>%
  filter(p > .05) %>%
  mutate(u = uchar / FC_mean) %>%
  group_by(peptide_sequence, group) %>%
  summarise_at("u", min) %>%
  # filter(u < 0.3) %>%
  mutate_at("u", ~ . * 100) %>%
  mutate(class = ifelse(u <= 20, "0~20%", "20%~")) %>%
  reshape2::dcast(., group ~ class, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(n = apply(., 1, function(x) sum(as.numeric(x[2:3]))))

stat_u <- df_uchar %>%
  filter(p > .05) %>%
  mutate(u = uchar / FC_mean) %>%
  group_by(protein_id, group) %>%
  summarise_at("u", min) %>%
  # filter(u < 0.3) %>%
  mutate_at("u", ~ . * 100) %>%
  mutate(class = ifelse(u <= 20, "0~20%", "20%~")) %>%
  reshape2::dcast(., group ~ class, value.var = "protein_id",
                  fun.aggregate = length) %>%
  mutate(n = apply(., 1, function(x) sum(as.numeric(x[2:3]))))

stat_u[, 2:4] <- stat_u[, 2:4]/stat_u[, 4]

## 单元间不均匀性引入的不确定度
stat_u <- homo_tables %>%
  rbindlist(., idcol = "group") %>%
  filter(is.homogeneous) %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  mutate(u = u_bu / FC_mean) %>%
  group_by(peptide_sequence, group) %>%
  summarise_at("u", min) %>%
  # filter(u < 0.3) %>%
  mutate_at("u", ~ . * 100) %>%
  mutate(class = ifelse(u <= 20, "0~20%", "20%~")) %>%
  reshape2::dcast(., group ~ class, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(n = apply(., 1, function(x) sum(as.numeric(x[2:3]))))

stat_u <- homo_tables %>%
  rbindlist(., idcol = "group") %>%
  filter(is.homogeneous) %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  mutate(u = u_bu / FC_mean) %>%
  group_by(protein_id, group) %>%
  summarise_at("u", min) %>%
  # filter(u < 0.3) %>%
  mutate_at("u", ~ . * 100) %>%
  mutate(class = ifelse(u <= 20, "0~20%", "20%~")) %>%
  reshape2::dcast(., group ~ class, value.var = "protein_id",
                  fun.aggregate = length) %>%
  mutate(n = apply(., 1, function(x) sum(as.numeric(x[2:3]))))

stat_u[, 2:4] <- stat_u[, 2:4]/stat_u[, 4]

## 单元内不均匀性引入的不确定度
stat_u <- homo_tables %>%
  rbindlist(., idcol = "group") %>%
  filter(is.homogeneous) %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  mutate(u = u_wu / FC_mean) %>%
  group_by(peptide_sequence, group) %>%
  summarise_at("u", min) %>%
  # filter(u < 0.3) %>%
  mutate_at("u", ~ . * 100) %>%
  mutate(class = ifelse(u <= 20, "0~20%", "20%~")) %>%
  reshape2::dcast(., group ~ class, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(n = apply(., 1, function(x) sum(as.numeric(x[2:3]))))

stat_u <- homo_tables %>%
  rbindlist(., idcol = "group") %>%
  filter(is.homogeneous) %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  mutate(u = u_wu / FC_mean) %>%
  group_by(protein_id, group) %>%
  summarise_at("u", min) %>%
  # filter(u < 0.3) %>%
  mutate_at("u", ~ . * 100) %>%
  mutate(class = ifelse(u <= 20, "0~20%", "20%~")) %>%
  reshape2::dcast(., group ~ class, value.var = "protein_id",
                  fun.aggregate = length) %>%
  mutate(n = apply(., 1, function(x) sum(as.numeric(x[2:3]))))

stat_u[, 2:4] <- stat_u[, 2:4]/stat_u[, 4]

## 不均匀性引入的不确定度
stat_u <- homo_tables %>%
  rbindlist(., idcol = "group") %>%
  filter(is.homogeneous) %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  mutate(u = u_homo / FC_mean) %>%
  group_by(peptide_sequence, group) %>%
  summarise_at("u", min) %>%
  # filter(u < 0.3) %>%
  mutate_at("u", ~ . * 100) %>%
  mutate(class = ifelse(u <= 20, "0~20%", "20%~")) %>%
  reshape2::dcast(., group ~ class, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(n = apply(., 1, function(x) sum(as.numeric(x[2:3]))))

stat_u <- homo_tables %>%
  rbindlist(., idcol = "group") %>%
  filter(is.homogeneous) %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  mutate(u = u_homo / FC_mean) %>%
  group_by(protein_id, group) %>%
  summarise_at("u", min) %>%
  # filter(u < 0.3) %>%
  mutate_at("u", ~ . * 100) %>%
  mutate(class = ifelse(u <= 20, "0~20%", "20%~")) %>%
  reshape2::dcast(., group ~ class, value.var = "protein_id",
                  fun.aggregate = length) %>%
  mutate(n = apply(., 1, function(x) sum(as.numeric(x[2:3]))))

stat_u[, 2:4] <- stat_u[, 2:4]/stat_u[, 4]

## 运输不稳定性引入的不确定度
stat_u <- stab_tables %>%
  rbindlist(., idcol = "group", use.names=TRUE) %>%
  filter(category %in% "Short-term-Transport", temperature %in% "-78.5°C", period %in% ("0, 3, 7, 14")) %>%
  mutate_at(7:12, as.numeric) %>%
  mutate(u_sts1 = s_b1 * 0.47) %>%
  filter(pass %in% "Yes") %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  mutate(u = u_sts1 / FC_mean) %>%
  group_by(peptide_sequence, group) %>%
  summarise_at("u", min) %>%
  # filter(u < 0.3) %>%
  mutate_at("u", ~ . * 100) %>%
  mutate(class = ifelse(u <= 20, "0~20%", "20%~")) %>%
  reshape2::dcast(., group ~ class, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(n = apply(., 1, function(x) sum(as.numeric(x[2:3]))))

stat_u <- stab_tables %>%
  rbindlist(., idcol = "group", use.names=TRUE) %>%
  filter(category %in% "Short-term-Transport", temperature %in% "-78.5°C", period %in% ("0, 3, 7, 14")) %>%
  mutate_at(7:12, as.numeric) %>%
  mutate(u_sts1 = s_b1 * 0.47) %>%
  filter(pass %in% "Yes") %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  mutate(u = u_sts1 / FC_mean) %>%
  group_by(protein_id, group, pass) %>%
  summarise_at("u", min) %>%
  # filter(u < 0.3) %>%
  mutate_at("u", ~ . * 100) %>%
  mutate(class = ifelse(u <= 20, "0~20%", "20%~")) %>%
  reshape2::dcast(., group ~ class, value.var = "protein_id",
                  fun.aggregate = length) %>%
  mutate(n = apply(., 1, function(x) sum(as.numeric(x[2:3]))))

stat_u[, 2:4] <- stat_u[, 2:4]/stat_u[, 4]

## 使用不稳定性引入的不确定度
stat_u <- stab_tables %>%
  rbindlist(., idcol = "group", use.names=TRUE) %>%
  filter(category %in% "Short-term-Use", temperature %in% "4°C", period %in% ("0, 1, 2, 4, 7")) %>%
  mutate_at(7:12, as.numeric) %>%
  mutate(u_sts2 = s_b1 * 0.23) %>%
  filter(pass %in% "Yes") %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  mutate(u = u_sts2 / FC_mean) %>%
  group_by(peptide_sequence, group) %>%
  summarise_at("u", min) %>%
  # filter(u < 0.3) %>%
  mutate_at("u", ~ . * 100) %>%
  mutate(class = ifelse(u <= 20, "0~20%", "20%~")) %>%
  reshape2::dcast(., group ~ class, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(n = apply(., 1, function(x) sum(as.numeric(x[2:3]))))

stat_u <- stab_tables %>%
  rbindlist(., idcol = "group", use.names=TRUE) %>%
  filter(category %in% "Short-term-Use", temperature %in% "4°C", period %in% ("0, 1, 2, 4, 7")) %>%
  mutate_at(7:12, as.numeric) %>%
  mutate(u_sts2 = s_b1 * 0.23) %>%
  filter(pass %in% "Yes") %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  mutate(u = u_sts2 / FC_mean) %>%
  group_by(protein_id, group) %>%
  summarise_at("u", min) %>%
  # filter(u < 0.3) %>%
  mutate_at("u", ~ . * 100) %>%
  mutate(class = ifelse(u <= 20, "0~20%", "20%~")) %>%
  reshape2::dcast(., group ~ class, value.var = "protein_id",
                  fun.aggregate = length) %>%
  mutate(n = apply(., 1, function(x) sum(as.numeric(x[2:3]))))

stat_u[, 2:4] <- stat_u[, 2:4]/stat_u[, 4]

## 长期不稳定性引入的不确定度
stat_u <- stab_tables %>%
  rbindlist(., idcol = "group", use.names=TRUE) %>%
  filter(category %in% "Long-term", grepl("0", period) & grepl("58.03", period)) %>%
  mutate_at(7:12, as.numeric) %>%
  mutate(u_lts = s_b1 * 12) %>%
  filter(pass %in% "Yes") %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  mutate(u = u_lts / FC_mean) %>%
  group_by(peptide_sequence, group) %>%
  summarise_at("u", min) %>%
  # filter(u < 0.3) %>%
  mutate_at("u", ~ . * 100) %>%
  mutate(class = ifelse(u <= 20, "0~20%", "20%~")) %>%
  reshape2::dcast(., group ~ class, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(n = apply(., 1, function(x) sum(as.numeric(x[2:3]))))

stat_u <- stab_tables %>%
  rbindlist(., idcol = "group", use.names=TRUE) %>%
  filter(category %in% "Long-term", grepl("0", period) & grepl("58.03", period)) %>%
  mutate_at(7:12, as.numeric) %>%
  mutate(u_lts = s_b1 * 12) %>%
  filter(pass %in% "Yes") %>%
  inner_join(., df_uchar, by = c("group", "peptide_sequence", "protein_id")) %>%
  mutate(u = u_lts / FC_mean) %>%
  group_by(protein_id, group) %>%
  summarise_at("u", min) %>%
  # filter(u < 0.3) %>%
  mutate_at("u", ~ . * 100) %>%
  mutate(class = ifelse(u <= 20, "0~20%", "20%~")) %>%
  reshape2::dcast(., group ~ class, value.var = "protein_id",
                  fun.aggregate = length) %>%
  mutate(n = apply(., 1, function(x) sum(as.numeric(x[2:3]))))

stat_u[, 2:4] <- stat_u[, 2:4]/stat_u[, 4]

## 不稳定性引入的不确定度
ustab_tables <- pblapply(stab_tables, function(tmp_table) {
  
  sub_stab1 <- tmp_table %>%
    # filter(category %in% "Short-term-Transport", grepl("0", period) & grepl("14", period) & !grepl("30", period)) %>%
    filter(category %in% "Short-term-Transport", period %in% ("0, 3, 7, 14")) %>%
    filter(temperature %in% "-78.5°C") %>%
    filter(pass %in% "Yes") %>%
    mutate_at(6:11, as.numeric) %>%
    mutate(u_sts1 = s_b1 * 0.47) %>%
    select(peptide_sequence, protein_id, u_sts1)
  
  sub_stab2 <- tmp_table %>%
    # filter(category %in% "Short-term-Use", grepl("0", period) & grepl("7", period)) %>%
    filter(category %in% "Short-term-Use", period %in% ("0, 1, 2, 4, 7")) %>%
    filter(temperature %in% "4°C") %>%
    filter(pass %in% "Yes") %>%
    mutate_at(6:11, as.numeric) %>%
    mutate(u_sts2 = s_b1 * 0.23) %>%
    select(peptide_sequence, protein_id, u_sts2)
  
  sub_stab3 <- tmp_table %>%
    filter(category %in% "Long-term", grepl("0", period) & grepl("58.03", period)) %>%
    filter(pass %in% "Yes") %>%
    mutate_at(6:11, as.numeric) %>%
    mutate(u_lts = s_b1 * 12) %>%
    select(peptide_sequence, protein_id, u_lts)
  
  sub_all <- sub_stab1 %>%
    inner_join(., sub_stab2, by = c("peptide_sequence", "protein_id")) %>%
    inner_join(., sub_stab3, by = c("peptide_sequence", "protein_id")) %>%
    inner_join(., df_uchar, by = c("peptide_sequence", "protein_id")) %>%
    mutate(u_stab = sqrt(u_sts1 ^ 2 + u_sts2 ^ 2 + u_lts ^ 2))
  
  return(sub_all)
})
stat_u <- ustab_tables %>%
  rbindlist(., use.names=TRUE) %>%
  mutate(u = u_stab / FC_mean) %>%
  group_by(peptide_sequence, group) %>%
  summarise_at("u", min) %>%
  # filter(u < 0.3) %>%
  mutate_at("u", ~ . * 100) %>%
  mutate(class = ifelse(u <= 20, "0~20%", "20%~")) %>%
  reshape2::dcast(., group ~ class, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(n = apply(., 1, function(x) sum(as.numeric(x[2:3]))))

stat_u <- ustab_tables %>%
  rbindlist(., use.names=TRUE) %>%
  mutate(u = u_stab / FC_mean) %>%
  group_by(protein_id, group) %>%
  summarise_at("u", min) %>%
  # filter(u < 0.3) %>%
  mutate_at("u", ~ . * 100) %>%
  mutate(class = ifelse(u <= 20, "0~20%", "20%~")) %>%
  reshape2::dcast(., group ~ class, value.var = "protein_id",
                  fun.aggregate = length) %>%
  mutate(n = apply(., 1, function(x) sum(as.numeric(x[2:3]))))

stat_u[, 2:4] <- stat_u[, 2:4]/stat_u[, 4]

## 合成不确定度
ucert_tables <- pblapply(1:3, function(i) {
  
  sub_uchar <- norm_tables[[i]] %>%
    filter(p > .05) %>%
    # mutate(value = ifelse(p > .05, FC_mean, FC_median)) %>%## ISO 33405:2024-A.2.5.3
    mutate(value = FC_mean, u_char = uchar) %>%## ISO 33405:2024-A.2.5.3
    select(peptide_sequence, protein_id, value, u_char)
  
  sub_homo <- homo_tables[[i]] %>%
    filter(is.homogeneous == TRUE) %>%
    inner_join(., sub_uchar, by = c("peptide_sequence", "protein_id")) %>%
    select(peptide_sequence, protein_id, value, u_char, u_bu, u_wu, u_homo)
  
  sub_stab <- ustab_tables[[i]] %>%
    inner_join(., sub_homo, by = c("peptide_sequence", "protein_id")) %>%
    select(peptide_sequence, protein_id, value, u_char, u_bu, u_wu, u_homo, u_sts1, u_sts2, u_lts, u_stab)
  
  stat_all1 <- data.frame(tier4 = length(unique(sub_uchar$peptide_sequence)),
                          tier5 = length(unique(sub_homo$peptide_sequence)),
                          tier6 = length(unique(sub_stab$peptide_sequence)))
  stat_all2 <- data.frame(tier4 = length(unique(sub_uchar$protein_id)),
                          tier5 = length(unique(sub_homo$protein_id)),
                          tier6 = length(unique(sub_stab$protein_id)))
                                                   
  sub_all <- sub_stab %>%
    mutate_at(c("u_bu", "u_wu", "u_sts1", "u_sts2", "u_lts"), ~ ifelse(. < u_char / 3, 0, .)) %>%
    mutate(u_comb = sqrt(u_char ^ 2 + u_bu ^ 2 + u_wu ^ 2 + u_sts1 ^ 2 + u_sts2 ^ 2 + u_lts ^ 2))
  
  return(sub_all)
})
names(ucert_tables) <- names(stab_tables)

## 新增B类不确定度
u_b_all <- c("D5/D6" = 0.0331, "F7/D6" = 0.0306, "M8/D6" = 0.0320)

ucert_df <- ucert_tables %>%
  rbindlist(., idcol = "group") %>%
  mutate(u_b = u_b_all[group]) %>%
  mutate(u_comb2 = sqrt(u_comb ^ 2 + u_b ^ 2)) %>%
  mutate(U2 = u_comb2 * 2) %>%
  mutate(U = u_comb * 2)

ucert_df_final <- ucert_df %>%
  mutate_at(5:17, ~ . / value) %>%
  filter(U < .2) %>%
  distinct() %>%
  # mutate_at(5:14, ~ sprintf("%.2f%%", . * 100)) %>%
  left_join(., db_uniprot, by = c("protein_id")) %>%
  select(group, peptide_sequence, protein_id, protein_fullname, gene_symbol, everything())

ucert_df_10peps <- ucert_df_final %>%
  filter(peptide_sequence %in% c("FDSDVGEFR", "LGVIEDHSNR",
                                 "GVVDSDDLPLNVSR", "TVDNFVALATGEK",
                                 "ELEEIVQPIISK", "ITPSYVAFTPEGER",
                                 "FQSSHHPTDITSLDQYVER", "SGEVYTCQVEHPSVTSPLTVEWR",
                                 "NEAIQAAHDAVAQEGQCR", "AAVDTYCR"))

fwrite(ucert_df_10peps, "~/Desktop/tmp.csv")
# fwrite(ucert_df_final, "./results/tables/7_quantprop_final.csv")

stat_ucomb <- ucert_df %>%
  mutate_at(5:17, ~ . / value) %>%
  group_by(peptide_sequence, group) %>%
  summarise_at("u_comb", min) %>%
  # filter(u_comb < 0.3) %>%
  mutate_at("u_comb", ~ . * 100) %>%
  mutate(class = ifelse(u_comb <= 20, "0~20%", "20%~")) %>%
  reshape2::dcast(., group ~ class, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(n = apply(., 1, function(x) sum(as.numeric(x[2:3]))))

fwrite(stat_ucomb, "~/Desktop/tmp.csv")

stat_ucomb[, 2:4] <- stat_ucomb[, 2:4]/stat_ucomb[, 4]

stat_ucert <- ucert_df %>%
  mutate_at(5:17, ~ . / value) %>%
  group_by(peptide_sequence, group) %>%
  summarise_at(c("U", "U2"), min) %>%
  filter(U < 0.2) %>%
  mutate_at("U2", ~ . * 100) %>%
  mutate(class = sapply(U2, function(x) {
    if (x <= 15) {
      return("0~15%")
    } else if (x >15 & x <= 20) {
      return("15%~20%")
    } else {
      return("20%~")
    }
  })) %>%
  reshape2::dcast(., group ~ class, value.var = "peptide_sequence",
                  fun.aggregate = length) %>%
  mutate(n = apply(., 1, function(x) sum(as.numeric(x[2:4]))))

stat_ucert[, 2:5] <- stat_ucert[, 2:5]/stat_ucert[, 5]

stat_ucert <- ucert_df %>%
  mutate_at(5:17, ~ . / value) %>%
  group_by(protein_id, group) %>%
  summarise_at(c("U", "U2"), min) %>%
  filter(U < 0.2) %>%
  mutate_at("U2", ~ . * 100) %>%
  mutate(class = sapply(U2, function(x) {
    if (x <= 15) {
      return("0~15%")
    } else if (x >15 & x <= 20) {
      return("15%~20%")
    } else {
      return("20%~")
    }
  })) %>%
  reshape2::dcast(., group ~ class, value.var = "protein_id",
                  fun.aggregate = length) %>%
  mutate(n = apply(., 1, function(x) sum(as.numeric(x[2:4]))))

stat_ucert[, 2:5] <- stat_ucert[, 2:5]/stat_ucert[, 5]

