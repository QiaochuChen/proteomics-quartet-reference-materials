## Title: 确定标称特性: 2025年10组重新搜库数据。
## Author: Qiaochu Chen
## Date: Jun 5th, 2026

library(arrow)
library(data.table)
library(dplyr)
library(pbapply)
library(stringr)
library(Biostrings)
library(parallel)
library(scales)
library(ggplot2)
library(RColorBrewer)
library(cowplot)


## 整理 metadata: DDA公共数据集 -----------------------
dda_files <- list.files("./data/public/PRIDE_PXD042233/HeLa", recursive = TRUE, full.names = TRUE)
dda_pep_files <- dda_files[grepl("evidence.txt", dda_files)]
dda_meta_tables <- pblapply(dda_pep_files, function(file_id) {
  
  dda_tmp <- fread(file_id, showProgress = FALSE)
  
  meta_tmp <- dda_tmp %>%
    dplyr::rename(analysis_id = `Raw file`) %>%
    distinct(analysis_id)
  
  return(meta_tmp)
})

all_meta_public1 <- dda_meta_tables %>%
  rbindlist %>%
  mutate(lab_id = rep(c("Lab_A", "Lab_B", "Lab_C", "Lab_D", "Lab_E", "lab_F"),
                      c(2, 2, 2, 3, 3, 3))) %>%
  mutate(sample = "HeLa", tube = 1) %>%
  mutate(injection = c(rep(1:2, 3), rep(1:3, 3))) %>%
  mutate(order = injection) %>%
  mutate(lab = lab_id) %>%
  mutate(mode = "DDA") %>%
  mutate(lc = "Unknown") %>%
  mutate(ms = rep(c("Q Extractive HF-X", "Oribitrap Fusion Lumos", "Oribitrap Exploris 480"),
                  c(6, 3, 6)))

fwrite(all_meta_public1, "./data/public/PRIDE_PXD042233/metadata_PXD042233.csv")


## 整理 metadata: DIA公共数据集 --------------------
columns2keep <- c("Run")
# labels.lab <- paste("Lab", 12:22, sep = "_")
labels.lab <- paste("Lab", 1:11, sep = "_")
names(labels.lab) <- paste("Lab", 1:11, sep = "_")

labels.lc <- c("Ultimate 3000 RSLC", "Easy-nLC 1200", "Ultimate 3000 RSLC",
               "Ultimate 3000 RSLC", "Easy-nLC 1200", "Ultimate 3000 RSLC",
               "Easy-nLC 1200", "Ultimate 3000 RSLC", "Ultimate 3000 RSLC",
               "Easy-nLC 1200", "Ultimate 3000 RSLC")
names(labels.lc) <- paste("Lab", 1:11, sep = "_")

labels.ms <- rep("Q Extractive HF", 11)
names(labels.ms) <- paste("Lab", 1:11, sep = "_")

labels.mode <- rep("DIA", 11)
names(labels.mode) <- paste("Lab", 1:11, sep = "_")

dia_tmp <- read_parquet("./data/public/MassIVE_MSV000084976/HeLa/report.parquet", col_select = columns2keep)

all_meta_public2 <- dia_tmp %>%
  dplyr::rename(analysis_id = Run) %>%
  distinct(analysis_id) %>%
  mutate(lab_id = str_extract(analysis_id, "(?i).+(?=-(Day_1))")) %>%
  mutate(sample = "HeLa", tube = 1) %>%
  mutate(set_order = as.numeric(gsub("Lab_", "", lab_id))) %>%
  mutate(injection = str_extract(analysis_id, "(?<=Rep_).+")) %>%
  arrange(set_order, injection, sample) %>%
  select(!set_order) %>%
  mutate(order = rep(1:3, 11)) %>%
  mutate(lab = labels.lab[lab_id]) %>%
  mutate(mode = labels.mode[lab_id]) %>%
  mutate(lc = labels.lc[lab_id]) %>%
  mutate(ms = labels.ms[lab_id])

fwrite(all_meta_public2, "./data/public/MassIVE_MSV000084976/metadata_2020_11labs_MSV000084976.csv")


## 整理 metadata: 本地数据集 -------------------
columns2keep <- c("Run")
labels.lab <- c("QLB", "NCP", "ZJU", "FDU", "NIM", "BTP", "TFS", "OSB", "CAS", "CMS")
# labels.lab <- paste("Lab", 1:10, sep = "_")
names(labels.lab) <- c("qinglian_bio", "phoenix", "zhejiang_university", "fudan_university",
                       "national_institute_of_methodology", "biotech_pack", "thermofisher_shanghai",
                       "omicsolution", "cas_tianjin", "academy_of_chinese_medical_sciences")

labels.ms <- c("Q Extractive HF-X", "Oribitrap Fusion", "Orbitrap Exploris 480",
               "timsTOF HT", "Orbitrap Fusion Lumos", "Oribitrap Exploris 480",
               "Orbitrap Astral", "timsTOF HT", "timsTOF Pro2", "ZenoTOF 7600")
names(labels.ms) <- names(labels.lab)

labels.lc <- c("EASY-nLC 1200", "EASY-nLC 1200", "EASY-nLC 1200", "nanoElute2",
               "EASY-nLC 1200", "Vanquish Neo", "Vanquish Neo", "Vanquish Neo",
               "nanoElute", "M-Class")
names(labels.lc) <- names(labels.lab)

labels.mode <- rep(c("DDA", "DIA"), c(4, 6))
names(labels.mode) <- names(labels.lab)

dda_files <- list.files("./data/multilab/2025_grouped_dda", recursive = TRUE, full.names = TRUE)
dda_pep_files <- dda_files[grepl("evidence.txt", dda_files)]
dda_meta_tables <- pblapply(dda_pep_files, function(file_id) {
  
  dda_tmp <- fread(file_id, showProgress = FALSE)
  
  meta_tmp <- dda_tmp %>%
    dplyr::rename(analysis_id = `Raw file`) %>%
    distinct(analysis_id)
  
  return(meta_tmp)
})

dia_files <- list.files("./data/multilab/2025_grouped_dia", recursive = TRUE, full.names = TRUE)
dia_pep_files <- dia_files[grepl("report.parquet$", dia_files)]
dia_meta_tables <- pblapply(dia_pep_files, function(file_id) {
  
  dia_tmp <- read_parquet(file_id, col_select = columns2keep)
  
  meta_tmp <- dia_tmp %>%
    dplyr::rename(analysis_id = Run) %>%
    distinct(analysis_id)
  
  return(meta_tmp)
})

all_meta_local <- c(dda_meta_tables, dia_meta_tables) %>%
  rbindlist %>%
  distinct(analysis_id) %>%
  arrange(analysis_id) %>%
  mutate(lab_id = str_extract(analysis_id, "(?i)(?<=data_proteomics_).+?(?=(_|-)(mzML|HT|LUMOS2|Quartet|HeLa|HEK293T))")) %>%
  mutate(run_id = str_extract(analysis_id, "(?i)(D5|D6|F7|M8|F8|HEK293T|HeLa)(_|-)[1-3](_|-)[1-3]")) %>%
  tidyr::separate(run_id, c("sample", "tube", "injection"), sep = "-|_") %>%
  mutate_at("sample", ~ gsub("Hela", "HeLa", .)) %>%
  mutate_at("sample", ~ gsub("F8", "M8", .)) %>%
  mutate_at("sample", ~ factor(., levels = c("HeLa", "HEK293T", "D5", "D6", "F7", "M8"))) %>%
  arrange(lab_id, tube, injection, sample) %>%
  mutate(order = c(rep(1:54, 9), 1:38)) %>%
  mutate(lab = labels.lab[lab_id]) %>%
  mutate(mode = labels.mode[lab_id]) %>%
  mutate(lc = labels.lc[lab_id]) %>%
  mutate(ms = labels.ms[lab_id])

fwrite(all_meta_local, "./data/multilab/metadata_2025_10labs.csv")


## 整理 metadata: 合并 ------------------
all_meta <- rbind(all_meta_local, all_meta_public1, all_meta_public2)
fwrite(all_meta, "./data/combined_meta.csv")


## 准备数据1——整理所有样本-肽段序列-蛋白质ID-基因Symbol等信息表 ------------------
rm(list = ls())
gc()

## DDA
dda_files <- list.files("data/multilab/2025_grouped_dda", recursive = TRUE, full.names = TRUE)
dda_pep_files <- dda_files[grepl("evidence.txt", dda_files)]
dda_pep_tables <- pblapply(dda_pep_files, function(file_id) {
  
  dda_tmp <- fread(file_id, showProgress = FALSE)
  
  dda_tmp_pep <- dda_tmp %>%
    select(`Raw file`, Sequence, `Leading razor protein`, `Gene names`,
           Charge, `m/z`, `Retention time`, PEP, Intensity, Length) %>%
    dplyr::rename(analysis_id = `Raw file`, peptide_sequence = Sequence,
                  protein_id = `Leading razor protein`, gene_symbol = `Gene names`,
                  mz_ratio = `m/z`, retention_time = `Retention time`) %>%
    dplyr::rename_all(tolower)
  
  return(dda_tmp_pep)
  
})

## DIA
dia_files <- list.files("data/multilab/2025_grouped_dia", recursive = TRUE, full.names = TRUE)
dia_pep_files <- dia_files[grepl("report.parquet", dia_files)]
dia_pep_tables <- NULL
for (file_id in dia_pep_files) {
  
  print(file_id)
  columns2keep <- c("Run", "Stripped.Sequence", "Protein.Group", "Genes",
                    "Precursor.Charge", "Precursor.Mz","RT", "PEP", "Precursor.Quantity")
  
  dia_tmp <- read_parquet(file_id, col_select = columns2keep)
  
  dia_tmp_pep <- dia_tmp %>%
    dplyr::mutate(length = nchar(Stripped.Sequence)) %>%
    tidyr::separate_rows(Protein.Group, sep = ";") %>%
    dplyr::rename(analysis_id = Run, peptide_sequence = Stripped.Sequence,
                  protein_id = Protein.Group, gene_symbol = Genes,
                  mz_ratio = Precursor.Mz, charge = Precursor.Charge,
                  retention_time = RT, intensity = Precursor.Quantity) %>%
    dplyr::rename_all(tolower)
  
  dia_pep_tables <- c(dia_pep_tables, list(dia_tmp_pep))
}

## 合并(及时释放内存)
all_pep_tables <- c(dda_pep_tables, dia_pep_tables)
rm(list = setdiff(ls(), c("all_pep_tables")))
gc()

all_meta <- fread("./data/multilab/metadata_2025_10labs.csv")
all_pep_tables2meta <- pblapply(all_pep_tables, function(tmp_table) {
  tmp_table2meta <- tmp_table %>% left_join(., all_meta, by = "analysis_id")
  return(tmp_table2meta)
})

rm(list = setdiff(ls(), c("all_pep_tables2meta")))
gc()

grouped_pep_tables <- all_pep_tables2meta %>% rbindlist %>% split(., by = "sample")

saveRDS(grouped_pep_tables, "./data/multilab/qualidata_list_pep_2025_10labs.rds")


## 准备数据2——计算肽序列覆盖率: CFFF ------------------
rm(list = ls())
gc()

## 读取所有样本-肽段序列-蛋白质ID一对一表
print("step1: preparing dda_sequence_tables...")

grouped_pep_tables <- readRDS("./data/multilab/qualidata_list_pep_2025_10labs.rds")
all_sequence_tables <- pblapply(grouped_pep_tables, function(tmp_table) {
  
  tmp_sequence_table <- tmp_table %>%
    distinct(analysis_id, peptide_sequence, protein_id)
  
  return(tmp_sequence_table)
})
print("all_sequence_tables done.")

## 准备uniprot数据库
rm(list = setdiff(ls(), c("all_sequence_tables")))
gc()

print("step2: preparing db_uniprot...")

all_meta <- fread("/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/cqc_21112030002/data_proteomics/metadata_2025_10labs.csv")
# all_meta <- fread("./data/multilab/metadata_2025_10labs.csv")
db_uniprot <- readAAStringSet("/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/cqc_21112030002/data_proteomics/uniprotkb_proteome_UP000005640_2024_09_03.fasta")
# db_uniprot <- readAAStringSet("./uniprotkb_proteome_UP000005640_2024_09_03.fasta")
db_uniprot <- data.frame(protein_sequence = as.character(db_uniprot)) %>%
  tibble::rownames_to_column("entry") %>%
  mutate(protein_id = str_extract(entry, "(?<=\\|).+(?=\\|)"))

print("db_uniprot done.")

## 手动计算蛋白质序列覆盖率
print("step3: calculating sequence_coverage_tables..")

sequence_coverage_tables <- mclapply(all_sequence_tables, function(tmp_table) {
  
  tryCatch({
    sequence_match_i <- tmp_table %>%
      # filter(protein_id %in% c("P55011", "Q86U42-2")) %>% ## test only
      dplyr::left_join(., db_uniprot, by = "protein_id") %>%
      mutate(peptide_len = nchar(peptide_sequence)) %>%
      mutate(start_pos = apply(., 1, function(x) {
        peptide_seq <- as.character(x[2])
        protein_seq <- as.character(x[5])
        start_pos <- gregexpr(peptide_seq, protein_seq)[[1]][1]
        return(start_pos)
      })) %>%
      mutate(end_pos = start_pos + peptide_len - 1)
    
    all_analysis_ids <- unique(sequence_match_i$analysis_id)
    
    sequence_coverage_tables_i <- mclapply(all_analysis_ids, function(analysis_tmp_id) {
      
      tryCatch({
        sequence_match_j <- sequence_match_i %>%
          filter(analysis_id %in% analysis_tmp_id)
        
        all_protein_ids <- unique(sequence_match_j$protein_id)
        
        ## 手动计算
        sequence_coverage_tables_j <- mclapply(all_protein_ids, function(protein_tmp_id) {
          
          tryCatch({# 逐个蛋白质计算
            positions <- sequence_match_j %>%
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
            sequence_coverage_k <- positions %>%
              distinct(analysis_id, protein_id, protein_sequence) %>%
              mutate(merged_len = sum(sapply(merged, function(x) x[2] - x[1] + 1))) %>%
              mutate(protein_len = nchar(protein_sequence)) %>%
              mutate(coverage = merged_len / protein_len * 100)
            
            return(sequence_coverage_k)}, error = function(e) {
              message("Error message: ", e$message)
              stop("Stopping execution due to error at protein_id = ", protein_tmp_id)
              
            })
        }, mc.cores = 6)
        
        sequence_coverage_j <- rbindlist(sequence_coverage_tables_j)
        return(sequence_coverage_j)
      }, error = function(e) {
        message("Error message: ", e$message)
        stop("Stopping execution due to error at analysis_id = ", analysis_tmp_id)
        
      })
    }, mc.cores = 3)
    
    sequence_coverage_i <- sequence_coverage_tables_i %>%
      rbindlist(.) %>%
      inner_join(all_meta, ., by = "analysis_id")
    
    return(sequence_coverage_i)
    print(paste("step4-", i, " done.", sep = ""))
  }, error = function(e) {
    message("Error message: ", e$message)
    stop("Stopping execution due to error at i = ", i)
    
  })
  
}, mc.cores = 3)

print("sequence_coverage_tables done.")

## 储存数据
print("step4: writing sequence_coverage_tables...")

saveRDS(sequence_coverage_tables, "/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/cqc_21112030002/qualidata_list_pro_coverage_2026_10labs.rds")

print("sequence_coverage_tables saved.")


## 根据后验错误概率PEP及其BH校正后FDR筛选肽段 ------------------
## 统计原始肽段/蛋白质数目
rm(list = ls())
grouped_seq_tables <- readRDS("./results/tables/qualidata_list_pro_coverage_2026_10labs.rds")
grouped_pep_tables <- readRDS("./data/multilab/qualidata_list_pep_2025_10labs.rds")

## 删去ZJU
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

## Tier1: 至少3家实验室所有运行证据支持SC>30%
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

## Tier2: 至少3家实验室所有运行证据支持PEP<1%
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

## Tier3: 实验室间所有运行证据支持FDR校正后PEP < 1%
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

