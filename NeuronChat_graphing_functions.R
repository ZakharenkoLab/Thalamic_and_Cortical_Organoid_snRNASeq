#modify color schemes for heatmap functions from NeuronChat

#################### For Heatmap plot
#' Heatmap plot with ligand abundance and target abundance for single ligand-target pair
#' @param object a NeuronChat object
#' @param interaction_name the ligand-target interaction used to plot
#' @param sender.names sending cell groups used for plot
#' @param receiver.names receiving cell groups used for plot
#' @param group a vector indicating which cell types (cell subclass, e.g., L2/3 IT) belong to which big groups (cell class, e.g., Glutamatergic)
#' @import circlize
#' @import ComplexHeatmap
#' @importFrom CellChat scPalette
#' @return
#' @export
heatmap_single2 <- function(object, interaction_name,group=NULL, sender.names=NULL,receiver.names=NULL){
  if(is.null(sender.names)){sender.names=rownames(object@net[[1]])}
  if(is.null(receiver.names)){receiver.names=colnames(object@net[[1]])}
  col_map = colorRamp2(c(0,max(object@net[[interaction_name]][sender.names,receiver.names])/2, max(object@net[[interaction_name]][sender.names,receiver.names])), c("blue", "white", "red"))
  column_ha = HeatmapAnnotation(target = anno_barplot(object@target.abundance[receiver.names,interaction_name],border = F,height = unit(2,'cm'),gp = gpar(fill = "lightskyblue1")),annotation_name_side='left',annotation_name_rot = 0,annotation_label = 'target \n abundance')
  row_ha = rowAnnotation(ligand = anno_barplot(object@ligand.abundance[sender.names,interaction_name],border = F,width = unit(2,'cm'),gp = gpar(fill = "lightskyblue1")),annotation_label = 'ligand \n abundance')
  if(is.null(group)){
    ComplexHeatmap::Heatmap(object@net[[interaction_name]][sender.names,receiver.names], name = "Commu. \n Prob.", top_annotation = column_ha, right_annotation = row_ha,
                            cluster_rows = F,cluster_columns = F,column_names_rot=45,row_names_side = 'left',col = col_map,
                            column_title = 'Receiver',row_title = 'Sender',column_title_side = 'bottom',
                            heatmap_legend_param = list(color_bar='continuous'))} else {
                              left_Annotation = rowAnnotation(foo = anno_block(gp = gpar(fill =scPalette(length(unique(group[sender.names])))), labels = sort(unique(group[sender.names]))))
                              bottom_Annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill =scPalette(length(unique(group[receiver.names])))), labels = sort(unique(group[receiver.names]))))
                              ComplexHeatmap::Heatmap(object@net[[interaction_name]][sender.names,receiver.names], name = "Commu. \n Prob.",
                                                      top_annotation = column_ha, right_annotation = row_ha,
                                                      left_annotation =left_Annotation,bottom_annotation = bottom_Annotation,
                                                      split=as.character(group[sender.names]),column_split = as.character(group[receiver.names]),
                                                      cluster_rows = F,cluster_columns = F,column_names_rot=45,row_names_side = 'left',col = col_map,
                                                      column_title = 'Receiver',row_title = 'Sender',column_title_side = 'bottom',
                                                      heatmap_legend_param = list(color_bar='continuous'))
                            }
}


#' Heatmap plot for the aggregated communication
#' @param object a NeuronChat object
#' @param interaction_use ligand-target interaction indexes  used for aggregation ('all' means use all interation pairs)
#' @param group a vector indicating which cell types (cell subclass, e.g., L2/3 IT) belong to which big groups (cell class, e.g., Glutamatergic)
#' @param method method used for aggregation; see also function `net_aggregation`
#' @param cut_off threshold used for aggregation; see also function `net_aggregation`
#' @param sender.names sending cell groups used for plot
#' @param receiver.names receiving cell groups used for plot
#' @import circlize
#' @import ComplexHeatmap
#' @importFrom CellChat scPalette
#' @return
#' @export
heatmap_aggregated2 <- function(object, method=c('weight','count','weighted_count','weighted_count2','weight_threshold'),cut_off=0.05, interaction_use='all',group=NULL, sender.names=NULL,receiver.names=NULL){
  method <- match.arg(method)
  if(is.null(sender.names)){sender.names=rownames(object@net[[1]])}
  if(is.null(receiver.names)){receiver.names=colnames(object@net[[1]])}
  if(interaction_use=='all') {
    net_aggregated <- net_aggregation(object@net,method=method,cut_off=cut_off)} else {net_aggregated <- net_aggregation(object@net[interaction_use],method=method,cut_off=cut_off)}
  #col_map = colorRamp2(c(0,max(net_aggregated[sender.names,receiver.names])/2, max(net_aggregated[sender.names,receiver.names])), c("blue", "white", "red"))
  library(RColorBrewer); col_map = brewer.pal(8,"YlGnBu")
  if(is.null(group)){
    ComplexHeatmap::Heatmap(net_aggregated[sender.names,receiver.names], name = method,
                            cluster_rows = F,cluster_columns = F,column_names_rot=45,row_names_side = 'left',col = col_map,
                            column_title = 'Receiver',row_title = 'Sender',column_title_side = 'bottom',
                            heatmap_legend_param = list(color_bar='continuous'))} else {
                              left_Annotation = rowAnnotation(foo = anno_block(gp = gpar(fill ='magenta'), labels = sort(unique(group[sender.names]))))
                              bottom_Annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c('green', 'magenta')), labels = sort(unique(group[receiver.names]))))
                              ComplexHeatmap::Heatmap(net_aggregated[sender.names,receiver.names], name = "Commu. \n Prob.",
                                                      left_annotation =left_Annotation,bottom_annotation = bottom_Annotation,
                                                      cluster_rows = F,cluster_columns = F,column_names_rot=45,row_names_side = 'left',col = col_map,
                                                      split=as.character(group[sender.names]),column_split = as.character(group[receiver.names]),
                                                      column_title = 'Receiver',row_title = 'Sender',column_title_side = 'bottom',
                                                      heatmap_legend_param = list(color_bar='continuous'))
                            }
}

heatmap_aggregated3 <- function(object, method=c('weight','count','weighted_count','weighted_count2','weight_threshold'),cut_off=0.05, interaction_use='all',group=NULL, sender.names=NULL,receiver.names=NULL){
  method <- match.arg(method)
  if(is.null(sender.names)){sender.names=rownames(object@net[[1]])}
  if(is.null(receiver.names)){receiver.names=colnames(object@net[[1]])}
  if(interaction_use=='all') {
    net_aggregated <- net_aggregation(object@net,method=method,cut_off=cut_off)} else {net_aggregated <- net_aggregation(object@net[interaction_use],method=method,cut_off=cut_off)}
  #col_map = colorRamp2(c(0,max(net_aggregated[sender.names,receiver.names])/2, max(net_aggregated[sender.names,receiver.names])), c("blue", "white", "red"))
  library(RColorBrewer); col_map = brewer.pal(8,"YlGnBu")
  if(is.null(group)){
    ComplexHeatmap::Heatmap(net_aggregated[sender.names,receiver.names], name = method,
                            cluster_rows = F,cluster_columns = F,column_names_rot=45,row_names_side = 'left',col = col_map,
                            column_title = 'Receiver',row_title = 'Sender',column_title_side = 'bottom',
                            heatmap_legend_param = list(color_bar='continuous'))} else {
                              left_Annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c('green', 'magenta')), labels = sort(unique(group[sender.names]))))
                              bottom_Annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c('green', 'magenta')), labels = sort(unique(group[receiver.names]))))
                              ComplexHeatmap::Heatmap(net_aggregated[sender.names,receiver.names], name = "Commu. \n Prob.",
                                                      left_annotation =left_Annotation,bottom_annotation = bottom_Annotation,
                                                      cluster_rows = F,cluster_columns = F,column_names_rot=45,row_names_side = 'left',col = col_map,
                                                      split=as.character(group[sender.names]),column_split = as.character(group[receiver.names]),
                                                      column_title = 'Receiver',row_title = 'Sender',column_title_side = 'bottom',
                                                      heatmap_legend_param = list(color_bar='continuous'))
                            }
}

heatmap_aggregated_subset <- function(object, method=c('weight','count','weighted_count','weighted_count2','weight_threshold'),cut_off=0.05, interaction_use='all',group=NULL, sender.names=NULL,receiver.names=NULL){
  method <- match.arg(method)
  if(is.null(sender.names)){sender.names=rownames(object[[1]])}
  if(is.null(receiver.names)){receiver.names=colnames(object[[1]])}
  if(interaction_use=='all') {
    net_aggregated <- net_aggregation(object,method=method,cut_off=cut_off)} else {net_aggregated <- net_aggregation(object[interaction_use],method=method,cut_off=cut_off)}
  col_map = brewer.pal(8,"YlGnBu")
  #col_map = colorRamp2(seq(0,0.1, length = 9), brewer.pal(9, 'YlGnBu'))
  if(is.null(group)){
    ComplexHeatmap::Heatmap(net_aggregated[sender.names,receiver.names], name = method,
                            cluster_rows = F,cluster_columns = F,column_names_rot=45,row_names_side = 'left',col = col_map,
                            column_title = 'Receiver',row_title = 'Sender',column_title_side = 'bottom',
                            heatmap_legend_param = list(color_bar='continuous'))} else {
                              left_Annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c('green', 'magenta')), labels = sort(unique(group[sender.names]))))
                              bottom_Annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c('green', 'magenta')), labels = sort(unique(group[receiver.names]))))
                              ComplexHeatmap::Heatmap(net_aggregated[sender.names,receiver.names], name = "Commu. \n Prob.",
                                                      left_annotation =left_Annotation,bottom_annotation = bottom_Annotation,
                                                      cluster_rows = F,cluster_columns = F,column_names_rot=45,row_names_side = 'left',col = col_map,
                                                      split=as.character(group[sender.names]),column_split = as.character(group[receiver.names]),
                                                      column_title = 'Receiver',row_title = 'Sender',column_title_side = 'bottom',
                                                      heatmap_legend_param = list(color_bar='continuous'))
                            }
}

