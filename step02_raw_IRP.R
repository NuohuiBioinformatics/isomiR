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


data_file = "/storeData/project/19-preRD_isomir_CRC/code_figures/10_new_analysis_output/merge_isomiR/merge_isomiR_833_ReadsCount.tsv"
raw_data <- readr::read_tsv(data_file)
raw_data = as.data.frame(raw_data)
rownames(raw_data) = raw_data$isomiRs  
raw_data = raw_data[, -1]

#colnames(raw_data) =  gsub("R03|\\.isomiRs-detection", "", colnames(raw_data))

colnames(raw_data) =  gsub("R03", "", colnames(raw_data))
colnames(raw_data) =  sub("-.*", "", colnames(raw_data))

# keep only those features with more than 5 reads in at least 10 samples, to reduce computation burdon
keep = rowSums(raw_data > 5) >= 10
cpm = calculate_cpm(raw_data)
irp = calculate_cpf(raw_data)

raw_data = raw_data[keep,]
cpm = cpm[keep,]
irp = irp[keep,]

saveRDS(raw_data, '/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/raw_data.rds')
saveRDS(cpm, '/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/cpm.rds')
saveRDS(irp, '/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/irp.rds')

