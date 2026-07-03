## Title: 数据清理，根据HeLa标准品公共数据库数据验证DIA数据有效性。
## Author: Qiaochu Chen
## Date: Jun 10th, 2026

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


## 统计鉴定数目: DDA + DIA (满足3针中至少2针检出) ----------------
rm(list = ls())
gc()

meta_all <- fread("./data/combined_meta.csv")
all_files <- list.files("./data", recursive = TRUE, full.names = TRUE)

dda_pep_files <- all_files[grepl("evidence.txt", all_files)]
dda_pep_tables <- pblapply(dda_pep_files, function(file_id) {
  
  dda_tmp <- fread(file_id, showProgress = FALSE)
  dda_tmp <- dda_tmp %>%
    dplyr::rename(analysis_id = `Raw file`) %>%
    left_join(., meta_all, by = "analysis_id")
  
  dda_tmp1 <- dda_tmp %>%
    filter(!grepl("CON__", `Leading razor protein`)) %>%
    mutate(Precursor.Id = paste(`Modified sequence`, Charge)) %>%
    reshape2::dcast(., Precursor.Id + lab + sample + tube ~ injection,
                    value.var = "Precursor.Id", fun.aggregate = length) %>%
    filter(apply(., 1, function(x) sum(is.na(x)) <= 1)) %>%
    group_by(lab, sample, tube) %>%
    summarise(Precursors.Identified = length(unique(Precursor.Id)),
              .groups = "drop")
  
  dda_tmp2 <- dda_tmp %>%
    filter(!grepl("CON__", `Leading razor protein`)) %>%
    filter(PEP < .01) %>%
    distinct(Sequence, lab, sample, tube, injection) %>%
    reshape2::dcast(., Sequence + lab + sample + tube ~ injection,
                    value.var = "Sequence") %>%
    filter(apply(., 1, function(x) sum(is.na(x)) <= 1)) %>%
    group_by(lab, sample, tube) %>%
    summarise(Peptides.Identified = length(unique(Sequence)),
              .groups = "drop")
  
  dda_tmp3 <- dda_tmp %>%
    filter(!grepl("CON__", `Leading razor protein`)) %>%
    filter(PEP < .01) %>%
    distinct(`Leading razor protein`, lab, sample, tube, injection) %>%
    reshape2::dcast(., `Leading razor protein` + lab + sample + tube ~ injection,
                    value.var = "Leading razor protein") %>%
    filter(apply(., 1, function(x) sum(is.na(x)) <= 1)) %>%
    group_by(lab, sample, tube) %>%
    summarise(Proteins.Identified = length(unique(`Leading razor protein`)),
              .groups = "drop")
  
  dda_stats_final <- dda_tmp1 %>%
    full_join(., dda_tmp2, by = c("lab", "sample", "tube")) %>%
    full_join(., dda_tmp3, by = c("lab", "sample", "tube"))
  
  return(dda_stats_final)
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
dia_pep_files <- all_files[grepl("report.parquet", all_files)]
dia_pep_files <- dia_pep_files[!grepl("crap_irt.+report.parquet", dia_pep_files)]
dia_pep_tables <- pblapply(dia_pep_files, function(file_id) {
  
  dia_tmp <- read_parquet(file_id, col_select = columns2keep)
  dia_tmp <- dia_tmp %>%
    dplyr::rename(analysis_id = Run) %>%
    left_join(., meta_all, by = "analysis_id")
  
  dia_tmp1 <- dia_tmp %>%
    reshape2::dcast(., Precursor.Id + lab + sample + tube ~ injection,
                    value.var = "Precursor.Id") %>%
    filter(apply(., 1, function(x) sum(is.na(x)) <= 1)) %>%
    group_by(lab, sample, tube) %>%
    summarise(Precursors.Identified = length(unique(Precursor.Id)),
              .groups = "drop")
  
  dia_tmp2 <- dia_tmp %>%
    filter(Q.Value <= .01, Global.Q.Value <= .01, Lib.Q.Value<= .01,
           Proteotypic == 1) %>%
    distinct(Stripped.Sequence, lab, sample, tube, injection) %>%
    reshape2::dcast(., Stripped.Sequence + lab + sample + tube ~ injection,
                    value.var = "Stripped.Sequence") %>%
    filter(apply(., 1, function(x) sum(is.na(x)) <= 1)) %>%
    group_by(lab, sample, tube) %>%
    summarise(Peptides.Identified = length(unique(Stripped.Sequence)),
              .groups = "drop")
  
  dia_tmp3 <- dia_tmp %>%
    filter(Protein.Q.Value <= .01,
           Q.Value <= .01, Global.Q.Value <= .01, Lib.Q.Value <= .01,
           PG.Q.Value <= 0.01, Global.PG.Q.Value <= 0.01, Lib.PG.Q.Value <= .01,
           GG.Q.Value <= 0.01, Global.PG.Q.Value <= .01,
           Proteotypic == 1) %>%
    distinct(Protein.Group, lab, sample, tube, injection) %>%
    reshape2::dcast(., Protein.Group + lab + sample + tube ~ injection,
                    value.var = "Protein.Group") %>%
    filter(apply(., 1, function(x) sum(is.na(x)) <= 1)) %>%
    group_by(lab, sample, tube) %>%
    summarise(Proteins.Identified = length(unique(Protein.Group)),
              .groups = "drop")
  
  dia_stats_final <- dia_tmp1 %>%
    full_join(., dia_tmp2, by = c("lab", "sample", "tube")) %>%
    full_join(., dia_tmp3, by = c("lab", "sample", "tube"))
  
  return(dia_stats_final)
})
dia_pep <- dia_pep_tables %>%
  rbindlist %>%
  mutate(mode = "DIA")

all_count <- rbind(dda_pep, dia_pep)

fwrite(all_count, "./results/tables/3_count.csv")


## 统计CV指标: DDA ------------------------------
rm(list = ls())
gc()

meta_all <- fread("./data/combined_meta.csv")
all_files <- list.files("./data", recursive = TRUE, full.names = TRUE)

## 本地数据集
dda_pre_files <- all_files[grepl("2025_grouped_dda.+evidence.txt", all_files)]
dda_pre_tables1 <- pblapply(dda_pre_files, function(file_id) {
  
  dda_tmp <- fread(file_id, showProgress = FALSE)
  dda_tmp <- dda_tmp %>%
    dplyr::rename(analysis_id = `Raw file`) %>%
    mutate(Precursors.Quantified = paste(`Modified sequence`, Charge)) %>%
    left_join(., meta_all, by = "analysis_id")
  
  dda_cv <- dda_tmp %>%
    group_by(Precursors.Quantified, lab, sample, tube, mode) %>%
    summarise(cv = sd(Intensity) / mean(Intensity),
              .groups = "drop") %>%
    na.omit
  
  return(dda_cv)
})

dda_pep_files <- all_files[grepl("2025_grouped_dda.+peptides.txt", all_files)]
dda_pep_tables1 <- pblapply(dda_pep_files, function(file_id) {
  
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
    dplyr::rename(Peptides.Quantified = Sequence) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id") %>%
    filter(value != 0) %>%
    mutate_at("analysis_id", ~ paste(id_prefix, str_extract(., "(?<=Intensity ).+"), sep = "")) %>%
    left_join(., meta_all, by = "analysis_id")
  
  dda_cv <- dda_tmp %>%
    group_by(Peptides.Quantified, lab, sample, tube, mode) %>%
    summarise(cv = sd(value) / mean(value), .groups = "drop") %>%
    na.omit
  
  return(dda_cv)
})

dda_pro_files <- all_files[grepl("2025_grouped_dda.+proteinGroups.txt", all_files)]
dda_pro_tables1 <- pblapply(dda_pro_files, function(file_id) {
  
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
    dplyr::rename(Proteins.Quantified = `Protein names`) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id") %>%
    filter(value != 0) %>%
    filter(!Proteins.Quantified %in% "") %>%
    mutate_at("analysis_id", ~ paste(id_prefix, str_extract(., "(?<=LFQ intensity ).+"), sep = "")) %>%
    left_join(., meta_all, by = "analysis_id")
  
  dda_cv <- dda_tmp %>%
    group_by(Proteins.Quantified, lab, sample, tube, mode) %>%
    summarise(cv = sd(value) / mean(value), .groups = "drop") %>%
    na.omit
  
  return(dda_cv)
})

## 公共数据集
dda_pre_files <- all_files[grepl("PRIDE_PXD042233.+evidence.txt", all_files)]
dda_pre_tables2 <- pblapply(dda_pre_files, function(file_id) {
  
  dda_tmp <- fread(file_id, showProgress = FALSE)
  dda_tmp <- dda_tmp %>%
    dplyr::rename(analysis_id = `Raw file`) %>%
    mutate(Precursors.Quantified = paste(`Modified sequence`, Charge)) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dda_tmp)
})
dda_pre_tables2 <- dda_pre_tables2 %>% rbindlist %>% split(., by = "lab")
dda_pre_tables2 <- pblapply(dda_pre_tables2, function(tmp_table) {
  
  dda_cv <- tmp_table %>%
    group_by(Precursors.Quantified, lab, sample, tube, mode) %>%
    summarise(cv = sd(Intensity) / mean(Intensity),
              .groups = "drop") %>%
    na.omit
  
  return(dda_cv)
})
dda_pre_tables2 <- dda_pre_tables2 %>% rbindlist %>% list(.)

dda_pep_files <- all_files[grepl("PRIDE_PXD042233.+peptides.txt", all_files)]
dda_pep_tables2 <- pblapply(dda_pep_files, function(file_id) {
  
  dda_tmp <- fread(file_id, showProgress = FALSE)
  
  evidence_file_id <- gsub("peptides.txt", "evidence.txt", file_id)
  evidence_tmp <- fread(evidence_file_id, showProgress = FALSE)
  
  dda_tmp <- dda_tmp %>%
    select(Sequence, Intensity) %>%
    dplyr::rename(Peptides.Quantified = Sequence, value = Intensity) %>%
    filter(value != 0) %>%
    mutate(analysis_id = unique(evidence_tmp$`Raw file`)) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dda_tmp)
})
dda_pep_tables2 <- dda_pep_tables2 %>% rbindlist %>% split(., by = "lab")
dda_pep_tables2 <- pblapply(dda_pep_tables2, function(tmp_table) {
  
  dda_cv <- tmp_table %>%
    group_by(Peptides.Quantified, lab, sample, tube, mode) %>%
    summarise(cv = sd(value) / mean(value), .groups = "drop") %>%
    na.omit
  
  return(dda_cv)
})
dda_pep_tables2 <- dda_pep_tables2 %>% rbindlist %>% list(.)

dda_pro_files <- all_files[grepl("PRIDE_PXD042233.+proteinGroups.txt", all_files)]
dda_pro_tables2 <- pblapply(dda_pro_files, function(file_id) {
  
  dda_tmp <- fread(file_id, showProgress = FALSE)
  
  evidence_file_id <- gsub("proteinGroups.txt", "evidence.txt", file_id)
  evidence_tmp <- fread(evidence_file_id, showProgress = FALSE)
  
  dda_tmp <- dda_tmp %>%
    select(`Protein names`, Intensity) %>%
    dplyr::rename(Proteins.Quantified = `Protein names`, value = Intensity) %>%
    filter(value != 0) %>%
    mutate(analysis_id = unique(evidence_tmp$`Raw file`)) %>%
    left_join(., meta_all, by = "analysis_id")
  
  return(dda_tmp)
})
dda_pro_tables2 <- dda_pro_tables2 %>% rbindlist %>% split(., by = "lab")
dda_pro_tables2 <- pblapply(dda_pro_tables2, function(tmp_table) {
  
  dda_cv <- tmp_table %>%
    group_by(Proteins.Quantified, lab, sample, tube, mode) %>%
    summarise(cv = sd(value) / mean(value), .groups = "drop") %>%
    na.omit
  
  return(dda_cv)
})
dda_pro_tables2 <- dda_pro_tables2 %>% rbindlist %>% list(.)

## 合并
dda_pre_tables <- c(dda_pre_tables1, dda_pre_tables2)
dda_pep_tables <- c(dda_pep_tables1, dda_pep_tables2)
dda_pro_tables <- c(dda_pro_tables1, dda_pro_tables2)

dda_cv_tables <- list(dda_pre_tables, dda_pep_tables, dda_pro_tables)
saveRDS(dda_cv_tables, "./results/tables/3_dda_cv.rds")

## 统计CV指标: DIA ------
rm(list = ls())
gc()

meta_all <- fread("./data/combined_meta.csv")
all_files <- list.files("./data", recursive = TRUE, full.names = TRUE)
all_files <- all_files[!grepl("crap_irt.+report.parquet", all_files)]

columns2keep <- c("Run", "Precursor.Id", "Ms1.Area")
dia_pre_files <- all_files[grepl("report.parquet", all_files)]
dia_pre_tables <- pblapply(dia_pre_files, function(file_id) {
  
  dia_tmp <- read_parquet(file_id, col_select = columns2keep)
  dia_tmp <- dia_tmp %>%
    dplyr::rename(analysis_id = Run, Precursors.Quantified = Precursor.Id) %>%
    left_join(., meta_all, by = "analysis_id")
  
  dia_cv <- dia_tmp %>%
    group_by(Precursors.Quantified, lab, sample, tube, mode) %>%
    summarise(cv = sd(Ms1.Area) / mean(Ms1.Area),
              .groups = "drop") %>%
    na.omit
  
  return(dia_cv)
})

dia_pep_files <- all_files[grepl("report.pr_matrix", all_files)]
dia_pep_tables <- pblapply(dia_pep_files, function(file_id) {
  
  dia_tmp <- fread(file_id)
  dia_tmp <- dia_tmp %>%
    select(Stripped.Sequence, contains("cpfs01")) %>%
    dplyr::rename(Peptides.Quantified = Stripped.Sequence) %>%
    dplyr::rename_if(is.numeric, ~ str_extract(., "[^/]+(?=\\.)")) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id", na.rm = TRUE) %>%
    left_join(., meta_all, by = "analysis_id")
  
  dia_cv <- dia_tmp %>%
    group_by(Peptides.Quantified, lab, sample, tube, mode) %>%
    summarise(cv = sd(value) / mean(value), .groups = "drop") %>%
    na.omit
  
  return(dia_cv)
})

dia_pro_files <- all_files[grepl("report.pg_matrix", all_files)]
dia_pro_tables <- pblapply(dia_pro_files, function(file_id) {
  
  dia_tmp <- fread(file_id)
  dia_tmp <- dia_tmp %>%
    select(First.Protein.Description, contains("cpfs01")) %>%
    dplyr::rename(Proteins.Quantified = First.Protein.Description) %>%
    dplyr::rename_if(is.numeric, ~ str_extract(., "[^/]+(?=\\.)")) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id", na.rm = TRUE) %>%
    left_join(., meta_all, by = "analysis_id")
  
  dia_cv <- dia_tmp %>%
    group_by(Proteins.Quantified, lab, sample, tube, mode) %>%
    summarise(cv = sd(value) / mean(value), .groups = "drop") %>%
    na.omit
  
  return(dia_cv)
})

dia_cv_tables <- list(dia_pre_tables, dia_pep_tables, dia_pro_tables)

saveRDS(dia_cv_tables, "./results/tables/3_dia_cv.rds")


## 准备数据1——整理hela公共数据库样本-肽段序列-蛋白质ID-基因Symbol等信息表 ------
rm(list = ls())

dia_files <- list.files("data/public/MassIVE_MSV000084976/HeLa", recursive = TRUE, full.names = TRUE)
dia_pep_files <- dia_files[grepl("report.parquet", dia_files)]
dia_pep_tables <- pblapply(dia_pep_files, function(file_id) {
  
  columns2keep <- c("Run", "Stripped.Sequence", "Protein.Group", "Genes",
                    "Precursor.Charge", "Precursor.Mz","RT", "PEP", "Precursor.Quantity")
  
  dia_tmp <- read_parquet(file_id, col_select = columns2keep)
  
  dia_tmp_pep <- dia_tmp %>%
    mutate(length = nchar(Stripped.Sequence)) %>%
    tidyr::separate_rows(Protein.Group, sep = ";") %>%
    dplyr::rename(analysis_id = Run, peptide_sequence = Stripped.Sequence,
                  protein_id = Protein.Group, gene_symbol = Genes,
                  mz_ratio = Precursor.Mz, charge = Precursor.Charge,
                  retention_time = RT, intensity = Precursor.Quantity) %>%
    dplyr::rename_all(tolower)
  
  return(dia_tmp_pep)
  
})

## 整理元信息
rm(list = setdiff(ls(), c("dia_pep_tables")))

all_meta <- dia_pep_tables %>%
  rbindlist %>%
  distinct(analysis_id) %>%
  mutate(lab_id = str_extract(analysis_id, "(?i).+(?=-(Day_1))")) %>%
  mutate(sample = "HeLa", tube = 1) %>%
  mutate(set_order = as.numeric(gsub("Lab_", "", lab_id))) %>%
  mutate(injection = str_extract(analysis_id, "(?<=Rep_).+")) %>%
  arrange(set_order, injection, sample) %>%
  select(!set_order) %>%
  mutate(order = rep(1:3, 11))

fwrite(all_meta, "./data/public/metadata_2025_11labs.csv")

## 合并(及时释放内存)
all_pep_tables2meta <- pblapply(dia_pep_tables, function(tmp_table) {
  tmp_table2meta <- tmp_table %>% left_join(., all_meta, by = "analysis_id")
  return(tmp_table2meta)
})

rm(list = setdiff(ls(), c("all_pep_tables2meta")))
gc()

grouped_pep_tables <- all_pep_tables2meta %>% rbindlist %>% split(., by = "sample")

saveRDS(grouped_pep_tables, "./data/public/qualidata_list_pep_2025_11labs.rds")


## 准备数据2——计算肽序列覆盖率: CFFF ------------------
rm(list = ls())
dia_files <- list.files("/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/cqc_21112030002/data_proteomics/eth_zurich_pmid33067419_hela/tsv_2025-02-24", recursive = TRUE, full.names = TRUE)
# dia_files <- list.files("data/public", recursive = TRUE, full.names = TRUE)
dia_pep_files <- dia_files[grepl("report.parquet", dia_files)]

dia_sequence_tables <- pblapply(dia_pep_files, function(file_id) {
  
  columns2keep <- c("Run", "Stripped.Sequence", "Protein.Group")
  
  dia_tmp <- read_parquet(file_id, col_select = columns2keep)
  
  dia_tmp_sequence <- dia_tmp %>%
    tidyr::separate_rows(Protein.Group, sep = ";") %>%
    dplyr::rename(analysis_id = Run, peptide_sequence = Stripped.Sequence,
                  protein_id = Protein.Group)
  
  return(dia_tmp_sequence)
  
})

print("dia_sequence_tables done.")

## 准备uniprot数据库
print("step1: preparing db_uniprot...")

rm(list = setdiff(ls(), c("dia_sequence_tables")))

meta_alllabs <- fread("/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/cqc_21112030002/data_proteomics/eth_zurich_pmid33067419_hela/metadata_2025_11labs.csv")
# meta_alllabs <- fread("./data/quartet/public/metadata_2025_11labs.csv")
db_uniprot <- readAAStringSet("/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/cqc_21112030002/data_proteomics/uniprotkb_proteome_UP000005640_2024_09_03.fasta")
# db_uniprot <- readAAStringSet("~/Desktop/utils/databases/uniprot/uniprotkb_proteome_UP000005640_2024_09_03.fasta")
db_uniprot <- data.frame(protein_sequence = as.character(db_uniprot)) %>%
  tibble::rownames_to_column("entry") %>%
  mutate(protein_id = str_extract(entry, "(?<=\\|).+(?=\\|)"))

print("db_uniprot done.")

## 手动计算肽序列覆盖率
print("step2: calculating sequence_coverage_tables..")

sequence_coverage_tables <- mclapply(1, function(i) {
  
  tryCatch({
    dia_sequence_match_i <- dia_sequence_tables[[i]] %>%
      dplyr::left_join(., db_uniprot, by = "protein_id") %>%
      mutate(peptide_len = nchar(peptide_sequence)) %>%
      mutate(start_pos = apply(., 1, function(x) {
        peptide_seq <- as.character(x[2])
        protein_seq <- as.character(x[5])
        start_pos <- gregexpr(peptide_seq, protein_seq)[[1]][1]
        return(start_pos)
      })) %>%
      mutate(end_pos = start_pos + peptide_len - 1)
    
    all_analysis_ids <- unique(dia_sequence_match_i$analysis_id)
    
    sequence_coverage_tables_i <- mclapply(all_analysis_ids, function(analysis_tmp_id) {
      
      tryCatch({
        dia_sequence_match_j <- dia_sequence_match_i %>%
          filter(analysis_id %in% analysis_tmp_id)
        
        all_protein_ids <- unique(dia_sequence_match_j$protein_id)
        
        ## 手动计算
        sequence_coverage_tables_j <- mclapply(all_protein_ids, function(protein_tmp_id) {
          
          tryCatch({# 逐个蛋白质计算
            positions <- dia_sequence_match_j %>%
              filter(protein_id %in% protein_tmp_id) %>%
              arrange(start_pos)
            
            # 检查数据无误
            if (nrow(positions) == 0|is.na(positions$entry[1])) return(NULL)
            
            # 按起始位置排序
            positions <- positions[order(positions$start_pos), ]
            
            # 合并重叠区间
            merged <- list()
            current <- c(positions$start_pos[1], positions$end_pos[1])
            for (i in 1:nrow(positions)) {
              if (positions$start_pos[i] <= current[2] + 1) {
                current[2] <- max(current[2], positions$end_pos[i])
              } else {
                merged <- c(merged, list(current))
                current <- c(positions$start_pos[i], positions$end_pos[i])
              }
            }
            merged <- c(merged, list(current))
            
            # 计算覆盖率
            dia_sequence_coverage_k <- positions %>%
              distinct(analysis_id, protein_id, protein_sequence) %>%
              mutate(merged_len = sum(sapply(merged, function(x) x[2] - x[1] + 1))) %>%
              mutate(protein_len = nchar(protein_sequence)) %>%
              mutate(coverage = merged_len / protein_len * 100)
            
            return(dia_sequence_coverage_k)}, error = function(e) {
              message("Error message: ", e$message)
              stop("Stopping execution due to error at protein_id = ", protein_tmp_id)
              
            })
        }, mc.cores = 6)
        
        dia_sequence_coverage_j <- rbindlist(sequence_coverage_tables_j)
        return(dia_sequence_coverage_j)
      }, error = function(e) {
        message("Error message: ", e$message)
        stop("Stopping execution due to error at analysis_id = ", analysis_tmp_id)
        
      }, mc.cores = 3)
    })
    
    dia_sequence_coverage_i <- sequence_coverage_tables_i %>%
      rbindlist(.) %>%
      inner_join(meta_alllabs, ., by = "analysis_id")
    
    return(dia_sequence_coverage_i)
    print(paste("step2-", i, " done.", sep = ""))
  }, error = function(e) {
    message("Error message: ", e$message)
    stop("Stopping execution due to error at i = ", i)
    
  })
  
}, mc.cores = 3)

print("sequence_coverage_tables done.")

## 储存数据
print("step3: writing all_sequence_coverage...")

all_sequence_coverage <- rbindlist(sequence_coverage_tables)
write.csv(all_sequence_coverage, "/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/cqc_21112030002/all_sequence_coverage_publicdata.csv", row.names = FALSE)

grouped_seq_tables <- all_sequence_coverage %>% split(., by = "sample")
saveRDS(grouped_seq_tables, "/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/cqc_21112030002/qualidata_list_pro_coverage_2025_11labs.rds")

## 本地文件: ./results/tables/3_qualidata_list_pro_coverage_hela_pmid33067419.rds
print("all_sequence_coverage done.")


## 准备数据3——整理hela公共数据库样本-肽段序列-表达谱 ------------------
rm(list = ls())
all_meta <- fread("./data/public/metadata_2025_11labs.csv")
dia_files <- list.files("./data/public", recursive = TRUE, full.names = TRUE)
dia_pep_files <- dia_files[grepl("report.pr_matrix.tsv", dia_files)]
dia_all <- pblapply(dia_pep_files, function(file_id) {
  
  dia_tmp <- fread(file_id, showProgress = FALSE)
  
  dia_i <- dia_tmp %>%
    select(Stripped.Sequence, starts_with("/cpfs01")) %>%
    plyr::rename(c("Stripped.Sequence" = "peptide_sequence")) %>%
    rename_if(is.numeric, ~ str_extract(., "(?<=mzML_2025-02-24\\/).+(?=(\\.mzML)|(\\.raw)|(\\.d))")) %>%
    reshape2::melt(., id = 1, variable.name = "analysis_id", na.rm = TRUE)
  
  return(dia_i)
})

## 合并所有数据
all_tables <- dia_all
all_tables2meta <- pblapply(all_tables, function(tmp_table) {
  tmp_table2meta <- tmp_table %>% left_join(., all_meta, by = "analysis_id")
  return(tmp_table2meta)
})
grouped_tables <- all_tables2meta %>% rbindlist %>% split(., by = "sample")
saveRDS(grouped_tables, "./data/public/quantdata_list_pep_2025_11labs.rds")


