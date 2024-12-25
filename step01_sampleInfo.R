require(readr)
require(tidyverse)
require(dplyr)
require(caret)
require(edgeR)
require(ggplot2)
require(pROC)
require(reshape2)
require(pROC)
require(readxl)
require(VennDiagram)
require(sva)
require(cluster)

source('/storeData/home/xucongmin/R_functions_bulkRNA.R')

raw_data = readRDS('/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/raw_data.rds')
sample_info <- readxl::read_xlsx("/storeData/project/19-preRD_isomir_CRC/code_figures/CRC_clinical.xlsx",sheet = "汇总")
names(sample_info)[1] <- "sample"
sample_info = as.data.frame(sample_info)
rownames(sample_info) = sample_info$sample
sample_info = sample_info[colnames(raw_data),]

#将数据里的Meta替换成英文
sample_info$性别 <- gsub("男","Male",sample_info$性别)
sample_info$性别 <- gsub("女","Female",sample_info$性别)
sample_info$建模分组 <- gsub("normal","Normal",sample_info$建模分组)
sample_info$建模分组 <- gsub("benign lesions|NAA","NAA+BL",sample_info$建模分组)
sample_info$建模分组 <- gsub("CRC-stageIII|CRC-stageIV","CRC III/IV",sample_info$建模分组)
sample_info$建模分组 <- gsub("CRC-stageI|CRC-stageII","CRC I/II",sample_info$建模分组)

#新建Type列，区分Control和Cancer
sample_info$Type <- sample_info$建模分组
sample_info$Type <- gsub("Normal|NAA+BL","Control",sample_info$Type) # For unknown reason, NAA+BL cannot be replaced via gsub
sample_info$Type[sample_info$Type=='NAA+BL'] = 'Control'
sample_info$Type <- gsub("CRC I/II|CRC III/IV","Cancer",sample_info$Type)

#sample processing date
sample_info$year <- sub("-.*", "", sample_info$miRNA样本处理时间)

saveRDS(sample_info, '/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/sample_info.rds')

