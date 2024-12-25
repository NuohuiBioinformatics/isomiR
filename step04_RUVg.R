source('R_functions.R')
source('load_lib.R')

setwd('/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/')

raw_data = readRDS('/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/raw_data.rds')
cpm = readRDS('/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/cpm.rds')
irp = readRDS('/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/irp.rds')
sample_info = readRDS('/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/sample_info.rds')

### show the distribution of house keeping genes selected from 1000 samples
hk.genes = c("hsa-miR-16-5p_delT_0_0", 
"hsa-miR-15b-5p_0_delACA_0", 
"hsa-miR-451a_delAAA_delT_0", 
"hsa-miR-19a-3p_0_delGA_0", 
"hsa-miR-15a-5p_delTAG_delG_0", 
"hsa-miR-32-5p_0_delCA_0", 
"hsa-let-7i-5p_0_delT_0", 
"hsa-miR-27a-3p_0_0_21CA", 
"hsa-miR-29c-3p_0_delA_0", 
"hsa-miR-140-3p_delT_insAC_0")

# pdf("hkgenes_raw_count_boxplot.pdf", width = 7, height = 5)
# df_t <- t(raw_data[hk.genes, ])
# boxplot(df_t, 
#         main = "Raw count",
#         xlab = "",
#         ylab = "",
#         las = 2, # Rotate x-axis labels for readability
#         col = "lightblue", # Color for the boxes
#         border = "blue")
# cv_values <- apply(df_t, 2, function(x) {  sd(x) / mean(x) * 100})
# text(x = 1:length(cv_values), y = apply(df_t, 2, max) * 1.2, 
#      labels = paste0("CV: ", round(cv_values, 1), "%"), 
#      cex = 0.7, col = "red")

# dev.off()

pdf("hkgenes_raw_count_log_boxplot.pdf", width = 7, height = 7)
df_t <- t(log(raw_data[hk.genes, ]+0.1))
par(mar = c(12, 4, 4, 2))  # Increase the bottom margin to accommodate x-axis labels
boxplot(df_t, 
        main = "Raw count",
        xlab = "",
        ylab = "",
        las = 2, # Rotate x-axis labels for readability
        col = "lightblue", # Color for the boxes
        border = "blue")
cv_values <- apply(df_t, 2, function(x) {  sd(x) / mean(x) * 100})
text(x = 1:length(cv_values), y = apply(df_t, 2, median)+1 , 
     labels = paste0("CV: ", round(cv_values, 1), "%"), 
     cex = 0.7, col = "red")
write.csv(df_t[, 8], file = 'test.csv')

dev.off()

# pdf("hkgenes_cpm_log_boxplot.pdf", width = 7, height = 7)
# df_t <- t(log(cpm[hk.genes, ]+0.1))
# par(mar = c(12, 4, 4, 2))  # Increase the bottom margin to accommodate x-axis labels
# boxplot(df_t, 
#         main = "Raw count",
#         xlab = "",
#         ylab = "",
#         las = 2, # Rotate x-axis labels for readability
#         col = "lightblue", # Color for the boxes
#         border = "blue")

# cv_values <- apply(df_t, 2, function(x) {  sd(x) / mean(x) * 100})
# text(x = 1:length(cv_values), y = apply(df_t, 2, median)+1 , 
#      labels = paste0("CV: ", round(cv_values, 1), "%"), 
#      cex = 0.7, col = "red")

# dev.off()

### ----------------------------------- RUVg noormalization --------------------------
library(RUVSeq)
# Create a DGEList object from the raw count data
dge <- DGEList(counts = raw_data)
# Create a SeqExpressionSet object
set <- newSeqExpressionSet(as.matrix(raw_data), 
                           phenoData = data.frame(row.names = colnames(raw_data)))
# Apply RUVg with the housekeeping genes
set_ruvg <- RUVg(set, hk.genes, k = 1)  # k = 1 means the first factor of unwanted variation will be removed
# You can increase k if you believe more unwanted variation needs to be accounted for.
# Get the normalized counts after RUVg
adjusted_counts <- normCounts(set_ruvg)

saveRDS(adjusted_counts, file = 'raw_count_RUVg.rds')



