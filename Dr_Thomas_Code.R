#Code for snRNA-seq analysis of thalamic and cortical organoids

#load required packages
library(Seurat)
library(voxhunt)
library(dplyr)
library(monocle3)
library(tidyr)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(NeuronChat)
library(CellChat)
library(ggalluvial)
library(SeuratObject)
library(ComplexHeatmap)
library(circlize)

#load customized graphing functions for NeuronChat analysis
source('NeuronChat_graphing_functions.R')

########################################
#hThO cluster IDs
########################################

#load hThO integrated data
snTha <- readRDS('snTha_Integrated_obj.rds')

#Assign names to clusters
new.cluster.ids <- c("ExN1", "ExN2", "ExN3", "PT/ZLI/rTh", "ExN4", 
                     "Glia", 
                     "Radial Glia", 
                     "Cycling Progenitors"
)

names(new.cluster.ids) <- levels(snTha)
snTha <- RenameIdents(snTha, new.cluster.ids)

#graph UMAP plot with new cluster IDs
DimPlot(snTha, reduction = "umap") +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_color_brewer(palette = 'Set2') +
  theme(legend.position = 'bottom',
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line= element_blank(),
        legend.text = element_text(size = 6))
ggsave('snTha_clusters.pdf', height = 10, width = 10, units = 'cm')

########################################
#VoxHunt Analysis
#Mapping ExN clusters onto e15 mouse brain
########################################

# Point VoxHunt to directory containing ABA expression data
load_aba_data('/mnt/storage1/Kristen/voxhunt/voxhunt_data/voxhunt_data')


#cTh: clusters containing thalamic excitatory neurons (caudal thalamus)
cTh <- subset(snTha, idents = c('ExN1', 'ExN2', 'ExN3', 'ExN4'))

#follow VoxHunt vignette
regional_markers <- structure_markers('E15') %>%
  group_by(group) %>%
  top_n(10, auc) %>% 
  {unique(.$gene)}
  
vox_map <- voxel_map(
  cTh, 
  stage = 'E15', 
  genes_use = regional_markers
)

#graph results: sagittal plane of e15 mouse brain
plot_map(vox_map, groups = 'none') +
  scale_fill_viridis_c(option = 'magma', direction = -1) +
  labs(fill = 'e15 mouse correlation (r)') +
  theme(legend.position = 'top')
ggsave('cTh_ABA_e15_corr.pdf', height = 6, width = 8, units = 'cm')

########################################
#VoxHunt Analysis
#Mapping PT/ZLI/rTh cluster onto E15 mouse brain
########################################

other_neurons <- subset(snTha, idents = 'PT/ZLI/rTh')

#follow VoxHunt vignette
vox_map2 <- voxel_map(
  other_neurons, 
  stage = 'E15', 
  genes_use = regional_markers
)

plot_map(vox_map2, groups = 'none') +
  scale_fill_viridis_c(option = 'magma', direction = -1) +
  labs(fill = 'e15 mouse correlation (r)') +
  theme(legend.position = 'top')
ggsave('PT.ZLI.rTh_ABA_e15_corr.pdf', height = 6, width = 8, units = 'cm')

########################################
#VoxHunt Analysis
#Mapping all hThO clusters onto BrainSpan data
########################################

#order clusters for graphing
levels(snTha) <- c("Cycling Progenitors" , 
                   "Radial Glia",
                   "Glia" ,
                   "PT/ZLI/rTh" ,
                   "ExN1",
                   "ExN2",
                   "ExN3",
                   "ExN4" 
)

data('brainspan') #load brainspan data

#follow VoxHunt vignette
ref_map <- brainspan_map(
  snTha,
  stages = 13:24, #include subjects ages 13-24 pcw
  #group_name = 'cluster',
  genes_use = regional_markers,
  pseudobulk_groups = T
)

#graph results
plot_structure_similarity(ref_map, annotation_level = 'structure_name', scale = F)
ggsave('VoxHunt_BrainSpan_snTha.pdf', height = 15, width = 10, units = 'cm')


########################################
#Pseudotime analysis of hThO clusters using Monocle3
#########################################

snTha.cds <- SeuratWrappers::as.cell_data_set(snTha)
snTha.cds <- cluster_cells(snTha.cds, reduction_method = 'UMAP')
snTha.cds <- learn_graph(snTha.cds, use_partition = TRUE)
snTha.cds <- order_cells(snTha.cds, reduction_method = 'UMAP')
#2 nodes were then chosen based on neural precursor gene expression

#add pseudotime to seurat object
snTha <- AddMetaData(
  object = snTha,
  metadata = snTha.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime"
)

#graph results
FeaturePlot(snTha, features = 'pseudotime') +
  scale_color_viridis_c(direction = -1) +
  labs(color = 'pseudotime') +
  ggtitle(NULL) +
  theme(legend.position = 'bottom',
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line= element_blank(),
        legend.text = element_text(size = 6))
ggsave('snTha_pseudotime.pdf', height = 10, width = 10, units = 'cm')

########################################
#hCO cluster IDs
########################################

#load hCO integrated data
snCBO <- readRDS('snCBO_Integrated_obj.rds')

#Assign names to clusters
new.cluster.ids <- c("Migrating ExN", "Maturing ExN", "Subplate/DL ExN", "DL ExN1", 
                     "Intermediate Progenitors", "Cycling Progenitors", 
                     "Un.ExN1",
                     "Choroid Plexus", 
                     "UL ExN", 
                     "Glia", 
                     "DL ExN2", "Un.ExN2")

names(new.cluster.ids) <- levels(snCBO)
snCBO <- RenameIdents(snCBO, new.cluster.ids)

#reorder clusters
levels(snCBO) <- c("Cycling Progenitors" , 
                   "Intermediate Progenitors",
                   "Migrating ExN" ,
                   "Maturing ExN" ,
                   "Subplate/DL ExN",
                   "DL ExN1",
                   "DL ExN2",
                   "UL ExN",
                   "Un.ExN1" ,
                   "Un.ExN2",
                   "Glia",
                   "Choroid Plexus" 
                   )

#graph UMAP plot with new cluster IDs
DimPlot(snCBO, reduction = "umap") +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_color_brewer(palette = 'Set3') +
  theme(legend.position = 'bottom',
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line= element_blank(),
        legend.text = element_text(size = 6))
ggsave('snCBO_clusters.pdf', height = 10, width = 12, units = 'cm')

########################################
#VoxHunt Analysis
#Mapping clusters onto E15 mouse brain
########################################

#follow VoxHunt vignette
vox_map3 <- voxel_map(
  snCBO, 
  stage = 'E15',  
  genes_use = regional_markers
)

#graph results: sagittal plane of e15 mouse brain
plot_map(vox_map3, groups = 'none') +
  scale_fill_viridis_c(option = 'magma', direction = -1) +
  labs(fill = 'e15 mouse correlation (r)') +
  theme(legend.position = 'top')
ggsave('hCO_ABA_e15_corr.pdf', height = 6, width = 8, units = 'cm')

########################################
#VoxHunt Analysis
#Mapping all hCO clusters onto BrainSpan data
########################################

#follow VoxHunt vignette
ref_map2 <- brainspan_map(
  snCBO,
  stages = 13:24,
  genes_use = regional_markers,
  pseudobulk_groups = T
)

#graph results
plot_structure_similarity(ref_map2, annotation_level = 'structure_name', scale = F)
ggsave('VoxHunt_BrainSpan_snCBO.pdf', height = 15, width = 10, units = 'cm')

########################################
#Pseudotime analysis of hCO clusters using Monocle3
#########################################

snCBO.cds <- SeuratWrappers::as.cell_data_set(snCBO)
snCBO.cds <- cluster_cells(snCBO.cds, reduction_method = 'UMAP')
snCBO.cds <- learn_graph(snCBO.cds, use_partition = TRUE)
snCBO.cds <- order_cells(snCBO.cds, reduction_method = 'UMAP')
#nodes were chosen based on neural precursor gene expression

#add pseudotime to seurat object
snCBO <- AddMetaData(
  object = snCBO,
  metadata = snCBO.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime"
)

#select cell trajectory of interest for Figure 2
main.trajectory <- subset(snCBO, idents = c("Cycling Progenitors" , 
                                              "Intermediate Progenitors",
                                              "Migrating ExN" ,
                                              "Maturing ExN" ,
                                              "Subplate/DL ExN",
                                              "DL ExN1",
                                              "DL ExN2",
                                              "UL ExN"))

#graph UMAP plot with pseudotime
FeaturePlot(main.trajectory, features = 'pseudotime') +
  scale_color_viridis_c(direction = -1) +
  xlim(-10,6.5) +
  labs(color = 'pseudotime') +
  ggtitle(NULL) +
  theme(legend.position = 'bottom',
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line= element_blank(),
        legend.text = element_text(size = 6))
ggsave('snCBO_pseudotime.pdf', height = 10, width = 10, units = 'cm')

########################################
#NeuronChat Analysis
#########################################

#Get normalized count data and meta data for NeuronChat analysis
#metadata
meta <- snCBO@meta.data
meta$active.ident <- snCBO@active.ident
meta <- meta[ ,c(9,26,43)]
meta <- meta %>%
  mutate(group = ifelse(grepl('Un.ExN', meta$active.ident), 'Ctx Un.ExN',
                        ifelse(grepl('ExN', meta$active.ident), 'Ctx ExN',
                        ifelse(grepl('Progenitors', meta$active.ident), 'Ctx Progenitors', 'Ctx Other'))))				
rownames(meta) <- paste(rownames(meta), '_C', sep = '')

meta_Tha <- snTha@meta.data
meta_Tha$active.ident <- snTha@active.ident
meta_Tha <- meta_Tha[ ,c(9,26,43)]
meta_Tha <- meta_Tha %>%
  mutate(group = ifelse(grepl('ExN', meta_Tha$active.ident), 'cTh',
                        ifelse(meta_Tha$active.ident == 'PT/ZLI/rTh', 'PT/ZLI/rTh', 'Th Other')))
rownames(meta_Tha) <- paste(rownames(meta_Tha), '_T', sep = '')
meta <- rbind(meta, meta_Tha)
rm(meta_Tha)

meta <- meta %>%
  mutate(active.ident.new = case_when(
    active.ident == 'Cycling Progenitors' & Organoid_protocol == 'cortex' ~ 'CO Cycling Progenitors',
    active.ident == 'Cycling Progenitors' & Organoid_protocol == 'thalamus'~ 'ThO Cycling Progenitors',
    active.ident == 'Glia' & Organoid_protocol == 'cortex'~ 'CO Glia',
    active.ident == 'Glia' & Organoid_protocol == 'thalamus'~ 'ThO Glia',
    TRUE ~ active.ident))

meta$active.ident.new <- factor(meta$active.ident.new, levels = c("CO Cycling Progenitors" , 
                   "Intermediate Progenitors",
                   "Migrating ExN" ,
                   "Maturing ExN" ,
                   "Subplate/DL ExN",
                   "DL ExN1",
                   "DL ExN2",
                   "UL ExN",
                   "Un.ExN1" ,
                   "Un.ExN2",
                   "CO Glia",
                   "Choroid Plexus", 
                   "ThO Cycling Progenitors" , 
                   "Radial Glia",
                   "ThO Glia" ,
                   "PT/ZLI/rTh" ,
                   "ExN1",
                   "ExN2",
                   "ExN3",
                   "ExN4" 
                   )
)

#normalized count data
target_df_CBO <- GetAssayData(object = snCBO, slot = "data")
colnames(target_df_CBO) <- paste(colnames(target_df_CBO), '_C', sep = '')
target_df_Tha <- GetAssayData(object = snTha, slot = "data")
colnames(target_df_Tha) <- paste(colnames(target_df_Tha), '_T', sep = '')
target_df <- RowMergeSparseMatrices(target_df_CBO, target_df_Tha) #use function from Seurat Object to perform quickly
rm(target_df_CBO, target_df_Tha)

#create group vector
group <- as.vector(meta$Organoid_protocol)
names(group) <- meta$active.ident.new
group <- group[unique(names(group))]

#create NeuronChat object and run NeuronChat
TCA_NCO <- createNeuronChat(target_df, DB = 'human', group.by = meta$active.ident.new)
TCA_NCO <- run_NeuronChat(TCA_NCO,M=100)
net_aggregated_TCA <- net_aggregation(TCA_NCO@net,method = 'weight', cut_off = 0.05)

#focus on thalamocortical interactions
net_aggregated_TCA2 <- net_aggregation(TCA_NCO@net,method = 'weight_threshold', cut_off = 0.8)
net_aggregated_TCA2 <- ifelse(net_aggregated_TCA2 < 0.6, 0, net_aggregated_TCA2)

pdf('net_cTh_PT.pdf', height = 4, width = 4)
netVisual_circle_neuron(net_aggregated_TCA2, 
                        sources.use= unique(subset(meta, group %in% c('cTh', 'PT/ZLI/rTh'))$active.ident.new), 
                        group=group, 
                        vertex.label.cex = 0.5, 
                        arrow.size = 1,
                        vertex.size.max = 8,
                        edge.width.max = 5)
dev.off()

pdf(file = 'heatmap_TC.pdf', height = 3.5, width = 8)
heatmap_aggregated2(TCA_NCO, method = 'weight', group = group,
                    sender.names= c('ExN1', 'ExN2', 'ExN3', 'ExN4', 'PT/ZLI/rTh'))
dev.off()

#focus on corticothalamic interactions
net_aggregated_TCA3 <- net_aggregation(TCA_NCO@net,method = 'weight_threshold', cut_off = 0.8)
net_aggregated_TCA3 <- ifelse(net_aggregated_TCA3 < 0.5, 0, net_aggregated_TCA3)

pdf('net_Ctx_ExN.pdf', height = 4, width = 4)
netVisual_circle_neuron(net_aggregated_TCA3, 
                        sources.use= unique(subset(meta, group == 'Ctx ExN')$active.ident.new), 
                        group=group, 
                        vertex.label.cex = 0.5, 
                        arrow.size = 1,
                        vertex.size.max = 8,
                        edge.width.max = 5)
dev.off()

pdf(file = 'heatmap_CT.pdf', height = 3.5, width = 8)
heatmap_aggregated3(TCA_NCO, method = 'weight', group = group,
                    sender.names= c('UL ExN', 'DL ExN2','DL ExN1', 'Subplate/DL ExN','Maturing ExN', 'Migrating ExN'))
dev.off()

#create heatmap focusing on glutamatergic,etc., interactions using custom function heatmap_aggregated_subset
#use for subsetted x@net objects

glu_TC <- TCA_NCO@net[grepl('Glu_', names(TCA_NCO@net))]
pdf(file = 'heatmap_glu.pdf', height = 6, width = 8)
heatmap_aggregated_subset(glu_TC, method = 'weight', group = group,
                          sender.names = c('ExN1', 'ExN2', 'ExN3', 'ExN4', 'PT/ZLI/rTh',  'UL ExN', 'DL ExN2','DL ExN1', 'Subplate/DL ExN','Maturing ExN', 'Migrating ExN'))
dev.off()


NRXN_TC <- TCA_NCO@net[grepl('NRXN', names(TCA_NCO@net))]
pdf(file = 'heatmap_nrxn.pdf', height = 6, width = 8)
heatmap_aggregated_subset(NRXN_TC, method = 'weight', group = group,
                          sender.names = c('ExN1', 'ExN2', 'ExN3', 'ExN4', 'PT/ZLI/rTh',  'UL ExN', 'DL ExN2','DL ExN1', 'Subplate/DL ExN','Maturing ExN', 'Migrating ExN'))
dev.off()
