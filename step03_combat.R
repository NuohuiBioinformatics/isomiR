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
cpm = readRDS('/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/cpm.rds')
irp = readRDS('/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/irp.rds')
sample_info = readRDS('/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/sample_info.rds')

x <- raw_data
batch <- factor(sample_info$year)
#design <- model.matrix(~as.factor(Type), data = sample_info)  # Define the design matrix
raw_data_combat <- ComBat(dat = x, batch = batch, mod = NULL, par.prior = TRUE)
saveRDS(raw_data_combat, file = '/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/raw_data_combat.rds')

cpm = readRDS('/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/cpm.rds')
sample_info = readRDS('/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/sample_info.rds')
x <- cpm
batch <- factor(sample_info$year)
#design <- model.matrix(~as.factor(Type), data = sample_info)  # Define the design matrix
cpm_combat <- ComBat(dat = x, batch = batch, mod = NULL, par.prior = TRUE)

saveRDS(cpm_combat, file = '/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/cpm_combat.rds')

irp = readRDS('/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/irp.rds')
sample_info = readRDS('/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/sample_info.rds')
x <- irp
batch <- factor(sample_info$year)
#design <- model.matrix(~as.factor(Type), data = sample_info)  # Define the design matrix
irp_combat <- ComBat(dat = x, batch = batch, mod = NULL, par.prior = TRUE)

saveRDS(irp_combat, file = '/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/irp_combat.rds')
