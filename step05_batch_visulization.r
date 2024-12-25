source('R_functions.R')
source('load_lib.R')

setwd('/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/')

raw_data = readRDS('/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/raw_data.rds')
cpm = readRDS('/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/cpm.rds')
irp = readRDS('/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/irp.rds')
sample_info = readRDS('/storeData/project/19-preRD_isomir_CRC/code_figures/2_isomiR_marker筛选及对marker的分析/sample_info.rds')
cpm_combat = readRDS('cpm_combat.rds')

adjusted_counts = readRDS('raw_count_RUVg.rds')
##------------------- Data batch visulization---------------------------------
# CPM -------------
pca_result <- prcomp(t(cpm), scale = TRUE)
pca_data <- as.data.frame(pca_result$x)
merged_data <- cbind(pca_data, sample_info)
merged_data_cancer_normal = subset(merged_data, Type == 'Cancer'|Type == 'Control')

pdf("PCA_dotplt_CPM.pdf", width = 7, height = 7)
# Create a PCA plot with colored dots based on sample types
ggplot(data = merged_data_cancer_normal, aes(x = PC1, y = PC2, color = year)) +
  geom_point() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2")
dev.off()

asw_cpm_year = average_silhouette_width_batch(pca_data[, 1:10], sample_info$year, metric = 'cosine')
asw_cpm_type = average_silhouette_width_batch(pca_data[, 1:10], as.factor(sample_info$Type), metric = 'cosine')
asw_cpm_chip = average_silhouette_width_batch(pca_data[, 1:10], as.factor(sample_info$miRNA测序芯片号), metric = 'cosine')
#print(paste('asw_CPM_year:', asw_cpm_year, '; asw_CPM_type:', asw_cpm_type))

# RUVg----------------
pca_result <- prcomp(t(adjusted_counts), scale = TRUE)
pca_data <- as.data.frame(pca_result$x)
merged_data <- cbind(pca_data, sample_info)
merged_data_cancer_normal = subset(merged_data, Type == 'Cancer'|Type == 'Control')

pdf("PCA_dotplt_RUVg.pdf", width = 7, height = 7)
# Create a PCA plot with colored dots based on sample types
ggplot(data = merged_data_cancer_normal, aes(x = PC1, y = PC2, color = year)) +
  geom_point() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2")
dev.off()

asw_RUVg_year = average_silhouette_width_batch(pca_data[, 1:10], sample_info$year, metric = 'cosine')
asw_RUVg_type = average_silhouette_width_batch(pca_data[, 1:10], as.factor(sample_info$Type), metric = 'cosine')
asw_RUVg_chip = average_silhouette_width_batch(pca_data[, 1:10], as.factor(sample_info$miRNA测序芯片号), metric = 'cosine')
#print(paste('asw_RUVg_year:', asw_RUVg_year, '; asw_RUVg_type:', asw_RUVg_type))

# raw_data---------------
pca_result <- prcomp(t(raw_data), scale = TRUE)
pca_data <- as.data.frame(pca_result$x)
merged_data <- cbind(pca_data, sample_info)
merged_data_cancer_normal = subset(merged_data, Type == 'Cancer'|Type == 'Control')

pdf("PCA_dotplt_RawCount.pdf", width = 7, height = 7)
# Create a PCA plot with colored dots based on sample types
ggplot(data = merged_data_cancer_normal, aes(x = PC1, y = PC2, color = year)) +
  geom_point() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2")
dev.off()

asw_rawCount_year = average_silhouette_width_batch(pca_data[, 1:10], sample_info$year, metric = 'cosine')
asw_rawCount_type = average_silhouette_width_batch(pca_data[, 1:10], as.factor(sample_info$Type), metric = 'cosine')
asw_rawcount_chip = average_silhouette_width_batch(pca_data[, 1:10], as.factor(sample_info$miRNA测序芯片号), metric = 'cosine')
#print(paste('asw_rawCount_year:', asw_rawCount_year, '; asw_rawCount_type:', asw_rawCount_type))

# CPM_combat------------------
pca_result <- prcomp(t(cpm_combat), scale = TRUE)
pca_data <- as.data.frame(pca_result$x)
merged_data <- cbind(pca_data, sample_info)
merged_data_cancer_normal = subset(merged_data, Type == 'Cancer'|Type == 'Control')

pdf("PCA_dotplt_CPM_combat.pdf", width = 7, height = 7)
# Create a PCA plot with colored dots based on sample types
ggplot(data = merged_data_cancer_normal, aes(x = PC1, y = PC2, color = year)) +
  geom_point() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2")
dev.off()

asw_cpm_combat_year = average_silhouette_width_batch(pca_data[, 1:10], sample_info$year, metric = 'cosine')
asw_cpm_combat_type = average_silhouette_width_batch(pca_data[, 1:10], as.factor(sample_info$Type), metric = 'cosine')
asw_cpm_combat_chip = average_silhouette_width_batch(pca_data[, 1:10], as.factor(sample_info$miRNA测序芯片号), metric = 'cosine')
#print(paste('asw_cpm_combat_year:', asw_cpm_combat_year, '; asw_cpm_combat_type:', asw_cpm_combat_type))

asw.df = data.frame(cpm = c(asw_cpm_year, asw_cpm_type, asw_cpm_chip), 
raw_count = c(asw_rawCount_year, asw_rawCount_type, asw_rawcount_chip),
RUVg = c(asw_RUVg_year, asw_RUVg_type, asw_RUVg_chip),
cpm_combat = c(asw_cpm_combat_year, asw_cpm_combat_type, asw_cpm_combat_chip))

row.names(asw.df) = c('Processing year','Sample group', 'Different Run')

print(asw.df)

write.csv(asw.df, file = 'asw_df.csv')