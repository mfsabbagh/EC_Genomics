#cellranger_utils.R

lookupGeneId<-function(eset,gene_names){
  res <- rownames(fData(eset))[fData(eset)$gene_short_name %in% gene_names]
  res <- c(res,rownames(fData(eset))[rownames(fData(eset)) %in% gene_names])
  res <- unique(res)
  res
}

lookupGeneName<-function(eset,gene_id){
  res <- fData(eset[gene_id,])$gene_short_name
  res <- unique(res)
  res
}

my_visualize_umi_counts<-function (gbm, projection, limits = c(0, 10), marker_size = 0.1)
{
    gene_values <- log(colSums(Biobase::exprs(gbm)), base = 10)
    gene_values[gene_values < limits[1]] <- limits[1]
    gene_values[gene_values > limits[2]] <- limits[2]
    projection_names <- colnames(projection)
    colnames(projection) <- c("Component.1", "Component.2")
    proj_gene <- data.frame(cbind(projection, gene_values))
    p <- ggplot(proj_gene, aes(Component.1, Component.2)) + geom_point(aes(colour = gene_values),
        size = marker_size) + scale_colour_gradient(low = "blue",
        high = "red", name = "log10") + labs(x = projection_names[1],
        y = projection_names[2]) + ggtitle("Total number of UMIs") +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    return(p)
}

my_visualize_gene_markers<-function (gbm, gene_probes, projection, limits = c(0, 10), marker_size = 0.1,
    title = NULL)
{
    gbm_trunc <- trunc_gbm_by_genes(gbm, gene_probes)
    gene_values <- t(as.matrix(Biobase::exprs(gbm_trunc)))
    gene_values[gene_values < limits[1]] <- limits[1]
    gene_values[gene_values > limits[2]] <- limits[2]
    colnames(gene_values) <- gene_probes
    projection_names <- colnames(projection)
    colnames(projection) <- c("Component.1", "Component.2")
    proj_gene <- data.frame(cbind(projection, gene_values))
    proj_gene <-cbind(proj_gene,pData(gbm))
    print(head(proj_gene))
    proj_gene_melt <- melt(proj_gene, id.vars = c("Component.1",
        "Component.2","barcode","genotype"))
    print(head(proj_gene_melt))
    p <- ggplot(proj_gene_melt, aes(Component.1, Component.2)) +
        geom_point(aes(colour = value,shape = genotype), size = marker_size) +
        facet_wrap(~variable) + scale_colour_gradient(low = "grey",
        high = "red", name = "val") + labs(x = projection_names[1],
        y = projection_names[2])
    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }
    p <- p + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    return(p)
}

visualize_clusters2<-function (cluster_result, projection, colour = NULL, alpha = 1,
    marker_size = 0.1, title = NULL)
{
    projection_names <- colnames(projection)
    colnames(projection) <- c("Component.1", "Component.2")
    proj_clu <- data.frame(cbind(projection, cluster_result))
    proj_clu_melt <- melt(proj_clu, id.vars = c("Component.1",
        "Component.2"))
    print(head(proj_clu_melt))
    p <- ggplot(proj_clu_melt, aes(Component.1, Component.2)) +
        geom_point(aes(colour = value), size = marker_size,
            alpha = alpha) + guides(col = guide_legend(title = "ID",
        override.aes = list(size = 3))) + facet_wrap(~variable) +
        labs(x = projection_names[1], y = projection_names[2]) +
        theme_bw()
    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }
    if ((!is.null(colour))) {
        names(colour) <- 1:length(colour)
        p <- p + scale_color_manual(values = colour)
    }
    if (is.vector(cluster_result)) {
        p <- p + theme(plot.title = element_text(hjust = 0.5),
            strip.text.x = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    }
    else {
        p <- p + theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    }
    return(p)
}


myTSNEPlotAlpha<-function(cds,markers=NULL,logMode=T,color_by="color",shape_by=NULL,scaled=FALSE,cell_size=2){
  tmp<-pData(cds)
  if(!is.null(markers)){
    genes<-as.matrix(Biobase::exprs(cds[rownames(fData(cds)) %in% lookupGeneId(cds,markers)]))
    if(logMode){
      genes<-log10(genes+1)
    }
    geneMeans<-rowMax(genes)
    #print(geneMeans)
    if(scaled){
      genes<-genes/geneMeans
    }
    genes<-t(genes)
    #print(genes)
    genes<-melt(genes)
    colnames(genes)<-c("cell_id","gene_id","value")
    genes<-merge(genes,fData(cds),by.x="gene_id",by.y="id",all.x=TRUE,sort=FALSE)
    #print(head(genes))
    tmp<-merge(tmp,genes,by.x=0,by.y="cell_id")
    #print(head(tmp))
    p<-ggplot(tmp,aes(x=tSNE1_pos,y=tSNE2_pos))
    if(is.null(shape_by)){
      p + geom_point(aes_string(color=color_by,alpha="value"),stroke=0,size=cell_size) + facet_wrap('gene_short_name')+ theme_bw() + scale_color_brewer(palette="Set1") + monocle:::monocle_theme_opts() + scale_alpha(range=c(0.05,1)) 
    }else{
      p + geom_point(aes_string(color=color_by,alpha="value",stroke=0,shape=shape_by),size=cell_size) + facet_wrap('gene_short_name')+ theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts() + scale_alpha(range=c(0.05,1)) 
    }
  }else{
    p<-ggplot(tmp,aes(x=tSNE1_pos,y=tSNE2_pos))
    if(is.null(shape_by)){
      p + geom_point(aes_string(color=color_by),size=cell_size) + theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts() 
    }else{
      p + geom_point(aes_string(color=color_by,shape=shape_by),size=cell_size) + theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts() 
    }
  }
}



myTSNEPlotAlpha2<-function(cds,markers=NULL,logMode=T,color_by="color",shape_by=NULL,scaled=FALSE,cell_size=2){
  tmp<-pData(cds)
  if(!is.null(markers)){
    genes<-as.matrix(Biobase::exprs(cds[rownames(fData(cds)) %in% lookupGeneId(cds,markers)]))
    if(logMode){
      genes<-log10(genes+1)
    }
    geneMeans<-rowMax(genes)
    #print(geneMeans)
    if(scaled){
      genes<-genes/geneMeans
    }
    genes<-t(genes)
    #print(genes)
    genes<-melt(genes)
    colnames(genes)<-c("cell_id","gene_id","value")
    genes<-merge(genes,fData(cds),by.x="gene_id",by.y="id",all.x=TRUE,sort=FALSE)
    #print(head(genes))
    tmp<-merge(tmp,genes,by.x=0,by.y="cell_id")
    #print(head(tmp))
    p<-ggplot(tmp,aes(x=RtSNE1_pos,y=RtSNE2_pos))
    if(is.null(shape_by)){
      p + geom_point(aes_string(color=color_by,alpha="value"),stroke=0,size=cell_size) + facet_wrap('gene_short_name')+ theme_bw() + scale_color_brewer(palette="Set1") + monocle:::monocle_theme_opts() + scale_alpha(range=c(0.05,1)) 
    }else{
      p + geom_point(aes_string(color=color_by,alpha="value",stroke=0,shape=shape_by),size=cell_size) + facet_wrap('gene_short_name')+ theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts() + scale_alpha(range=c(0.05,1)) 
    }
  }else{
    p<-ggplot(tmp,aes(x=RtSNE1_pos,y=RtSNE2_pos))
    if(is.null(shape_by)){
      p + geom_point(aes_string(color=color_by),size=cell_size) + theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts() 
    }else{
      p + geom_point(aes_string(color=color_by,shape=shape_by),size=cell_size) + theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts() 
    }
  }
}

###New custom functions by JH

#visualize_clusters3

visualize_clusters3<-function (cluster_result, legend_title = "ID", projection, colour = NULL, alpha = 1,
                               marker_size = 0.1, title = NULL, axis_line_size = 0.5, axis_tick_size = 0.5)
{
  projection_names <- colnames(projection)
  colnames(projection) <- c("Component.1", "Component.2")
  proj_clu <- data.frame(cbind(projection, cluster_result))
  proj_clu_melt <- melt(proj_clu, id.vars = c("Component.1",
                                              "Component.2"))
  print(head(proj_clu_melt))
  p <- ggplot(proj_clu_melt, aes(Component.1, Component.2)) +
    geom_point(aes(colour = value), size = marker_size,
               alpha = alpha) + guides(col = guide_legend(title = legend_title,
                                                          override.aes = list(size = 3), keywidth = 0.1, keyheight = 0.3, default.unit = "inch")) + facet_wrap(~variable) +
    labs(x = projection_names[1], y = projection_names[2]) + theme_classic() 
  #theme_bw()
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  if ((!is.null(colour))) {
    names(colour) <- 1:length(colour)
    p <- p + scale_color_manual(values = colour)
  }
  if (is.vector(cluster_result)) {
    p <- p + theme(axis.line = element_line(colour = 'black', size = axis_line_size),
                   axis.ticks = element_line(colour = "black", size = axis_tick_size),
                   plot.title = element_text(hjust = 0.5),strip.background = element_blank(),
                   strip.text.x = element_blank(), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
  }
  else {
    p <- p + theme(axis.line = element_line(colour = 'black', size = axis_line_size),
                   axis.ticks = element_line(colour = "black", size = axis_tick_size),
                   plot.title = element_text(hjust = 0.5),strip.background = element_blank(),
                   strip.text.x = element_blank(),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }
  return(p)
}

###



#visualize_gene_markers2

visualize_gene_markers2 <-function (gbm, gene_probes, projection, limits = c(0, 0.5), low_col = "lightblue", high_col = "darkblue", marker_size = 0.1, 
                                    title = NULL, axis_line_size = 1.0, axis_tick_size = 1.0, panel_border_size = 1.0) 
{
  gbm_trunc <- trunc_gbm_by_genes(gbm, gene_probes)
  gene_values <- t(as.matrix(Biobase::exprs(gbm_trunc)))
  gene_values[gene_values < limits[1]] <- limits[1]
  gene_values[gene_values > limits[2]] <- limits[2]
  colnames(gene_values) <- gene_probes
  projection_names <- colnames(projection)
  colnames(projection) <- c("Component.1", "Component.2")
  proj_gene <- data.frame(cbind(projection, gene_values))
  proj_gene_melt <- melt(proj_gene, id.vars = c("Component.1", "Component.2"))
  colnames(proj_gene_melt) <- c("Component.1", "Component.2", "gene_probe", "value")
  p <- ggplot(proj_gene_melt, aes(Component.1, Component.2)) + 
    geom_point(aes(colour = value), size = marker_size) + 
    facet_wrap(~gene_probe,ncol=7, scales = "free") + scale_colour_gradient(low = low_col, 
                                                  high = high_col, name = "value") + theme(strip.background = element_rect(colour = "black", fill = "white"))
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  p <- p + theme_bw() + theme(panel.border = element_rect(linetype = "solid", size = panel_border_size, fill = NA),
                              plot.title = element_text(hjust = 0.5),
                              axis.line.x = element_line(colour = 'black', size = axis_line_size),
                              axis.line.y = element_line(colour = 'black', size = axis_line_size),
                              axis.ticks = element_line(colour = "black", size = axis_tick_size),
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                              strip.background = element_rect(colour = "black", fill = "white"),
                              legend.key.height = unit(0.85, "cm"))
                      
  return(p)
}


####





#plot_complex_cell_trajectory_markers

plot_complex_cell_trajectory_markers <-  function (cds, x = 1, y = 2, root_states = NULL, color_by = "State", 
                                            show_tree = TRUE, show_backbone = TRUE, backbone_color = "black", 
                                            markers = NULL, limits = c(0,10), low_col = "lightblue", high_col = "darkblue", nrow = 1, show_cell_names = FALSE, cell_size = 1.5, 
                                            alpha = alpha, cell_link_size = 0.75, branch_point_size = 2.0, cell_name_size = 2, legend_text_size = 12, show_branch_points = TRUE, 
                                            theta = 0, legend.position = "right", axis_line_size = 0.5, axis_tick_size = 0.5) 
{
  gene_short_name <- NA
  sample_name <- NA
  data_dim_1 <- NA
  data_dim_2 <- NA
  lib_info_with_pseudo <- pData(cds)
  if (is.null(cds@dim_reduce_type)) {
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
  }
  if (cds@dim_reduce_type == "ICA") {
    reduced_dim_coords <- reducedDimS(cds)
  }
  else if (cds@dim_reduce_type %in% c("SimplePPT", "DDRTree", 
                                      "SGL-tree")) {
    reduced_dim_coords <- reducedDimK(cds)
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  }
  else {
    stop("Error: unrecognized dimensionality reduction method.")
  }
  if (is.null(reduced_dim_coords)) {
    stop("You must first call reduceDimension() before using this function")
  }
  dp_mst <- minSpanningTree(cds)
  if (is.null(root_states)) {
    if (is.null(lib_info_with_pseudo$Pseudotime)) {
      root_cell <- row.names(lib_info_with_pseudo)[degree(dp_mst) == 
                                                     1][1]
    }
    else root_cell <- row.names(subset(lib_info_with_pseudo, 
                                       Pseudotime == 0))
    if (cds@dim_reduce_type != "ICA") 
      root_cell <- igraph::V(dp_mst)$name[cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[root_cell, 
                                                                                                ]]
  }
  else {
    candidate_root_cells <- row.names(subset(pData(cds), 
                                             State %in% root_states))
    if (cds@dim_reduce_type == "ICA") {
      root_cell <- candidate_root_cells[which(degree(dp_mst, 
                                                     candidate_root_cells) == 1)]
    }
    else {
      Y_candidate_root_cells <- igraph::V(dp_mst)$name[cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[candidate_root_cells, 
                                                                                                             ]]
      root_cell <- Y_candidate_root_cells[which(degree(dp_mst, 
                                                       Y_candidate_root_cells) == 1)]
    }
  }
  tree_coords <- igraph::layout_as_tree(dp_mst, root = root_cell)
  ica_space_df <- data.frame(tree_coords)
  row.names(ica_space_df) <- colnames(reduced_dim_coords)
  colnames(ica_space_df) <- c("prin_graph_dim_1", "prin_graph_dim_2")
  ica_space_df$sample_name <- row.names(ica_space_df)
  if (is.null(dp_mst)) {
    stop("You must first call orderCells() before using this function")
  }
  edge_list <- as.data.frame(igraph::get.edgelist(dp_mst))
  colnames(edge_list) <- c("source", "target")
  edge_df <- merge(ica_space_df, edge_list, by.x = "sample_name", 
                   by.y = "source", all = TRUE)
  edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = "source_prin_graph_dim_1", 
                                     prin_graph_dim_2 = "source_prin_graph_dim_2"))
  edge_df <- merge(edge_df, ica_space_df[, c("sample_name", 
                                             "prin_graph_dim_1", "prin_graph_dim_2")], by.x = "target", 
                   by.y = "sample_name", all = TRUE)
  edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = "target_prin_graph_dim_1", 
                                     prin_graph_dim_2 = "target_prin_graph_dim_2"))
  if (cds@dim_reduce_type == "ICA") {
    S_matrix <- tree_coords[, ]
  }
  else if (cds@dim_reduce_type %in% c("DDRTree", "SimplePPT", 
                                      "SGL-tree")) {
    S_matrix <- tree_coords[closest_vertex, ]
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  }
  data_df <- data.frame(S_matrix)
  row.names(data_df) <- colnames(reducedDimS(cds))
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- row.names(data_df)
  data_df <- merge(data_df, lib_info_with_pseudo, by.x = "sample_name", 
                   by.y = "row.names")
  return_rotation_mat <- function(theta) {
    theta <- theta/180 * pi
    matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 
           nrow = 2)
  }
  tmp <- return_rotation_mat(theta) %*% t(as.matrix(data_df[, 
                                                            c(2, 3)]))
  data_df$data_dim_1 <- tmp[1, ]
  data_df$data_dim_2 <- tmp[2, ]
  tmp <- return_rotation_mat(theta = theta) %*% t(as.matrix(edge_df[, 
                                                                    c("source_prin_graph_dim_1", "source_prin_graph_dim_2")]))
  edge_df$source_prin_graph_dim_1 <- tmp[1, ]
  edge_df$source_prin_graph_dim_2 <- tmp[2, ]
  tmp <- return_rotation_mat(theta) %*% t(as.matrix(edge_df[, 
                                                            c("target_prin_graph_dim_1", "target_prin_graph_dim_2")]))
  edge_df$target_prin_graph_dim_1 <- tmp[1, ]
  edge_df$target_prin_graph_dim_2 <- tmp[2, ]
  markers_exprs <- NULL
  if (is.null(markers) == FALSE) {
    markers_fData <- subset(fData(cds), gene_short_name %in% 
                              markers)
    if (nrow(markers_fData) >= 1) {
      markers_exprs <- reshape2::melt(as.matrix(Biobase::exprs(cds[row.names(markers_fData), 
                                                          ])))
      colnames(markers_exprs)[1:2] <- c("feature_id", "cell_id")
      #Set limits for  values
      markers_exprs$value[markers_exprs$value < limits[1]] <- limits[1]
      markers_exprs$value[markers_exprs$value > limits[2]] <- limits[2]
      ##
      markers_exprs <- merge(markers_exprs, markers_fData, 
                             by.x = "feature_id", by.y = "row.names")
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 
      0) {
    data_df <- merge(data_df, markers_exprs, by.x = "sample_name", 
                     by.y = "cell_id")
    g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2 
                                    )) + facet_wrap(~feature_label, nrow = nrow)
  } #remove size adjustment
  else {
    g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
  }
  if (show_tree) {
    g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                     y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                     yend = "target_prin_graph_dim_2"), size = cell_link_size, 
                          linetype = "solid", na.rm = TRUE, data = edge_df)
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 
      0) {
    g <- g + geom_jitter(aes(color = value), alpha = alpha, na.rm = TRUE,size = I(cell_size), 
                         height = 5) + scale_colour_gradient(low = low_col, high = high_col, name = "value")
  } ##set color for expression value
  else {
    g <- g + geom_jitter(aes(color = value), alpha = alpha, size = I(cell_size), 
                         na.rm = TRUE, height = 5) + scale_colour_gradient(low = low_col, high = high_col, name = "value")
  }  ##set color for expression value
  if (show_branch_points && cds@dim_reduce_type == "DDRTree") {
    mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
    branch_point_df <- subset(edge_df, sample_name %in% mst_branch_nodes)[, 
                                                                          c("sample_name", "source_prin_graph_dim_1", "source_prin_graph_dim_2")]
    branch_point_df$branch_point_idx <- match(branch_point_df$sample_name, 
                                              mst_branch_nodes)
    branch_point_df <- branch_point_df[!duplicated(branch_point_df$branch_point_idx), 
                                       ]
    g <- g + geom_point(aes_string(x = "source_prin_graph_dim_1", 
                                   y = "source_prin_graph_dim_2"), size = branch_point_size*cell_size, 
                        na.rm = TRUE, data = branch_point_df) + geom_text(aes_string(x = "source_prin_graph_dim_1", 
                                                                                     y = "source_prin_graph_dim_2", label = "branch_point_idx"), 
                                                                          size = branch_point_size*(0.7)*cell_size, color = "white", na.rm = TRUE, 
                                                                          data = branch_point_df)
  }
  if (show_cell_names) {
    g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
  }
  g <- g + monocle:::monocle_theme_opts() + xlab(paste("Component", x)) + 
    ylab(paste("Component", y)) + theme(legend.position = legend.position, 
                                        legend.key.height = grid::unit(0.35, "in"), legend.text = element_text(size = legend_text_size)) + theme(legend.key = element_blank()) + 
    theme(panel.background = element_rect(fill = "white"), strip.text.x = element_text(face = "italic"))
  g + theme(axis.line.x = element_line(colour = 'black', size = axis_line_size),
              axis.line.y = element_line(colour = 'black', size = axis_line_size),
              axis.ticks = element_line(colour = "black", size = axis_tick_size))
}

###






#plot_complex_cell_trajectory2

plot_complex_cell_trajectory2 <- function (cds, x = 1, y = 2, root_states = NULL, color_by = "State", 
          show_tree = TRUE, show_backbone = TRUE, show_legend = TRUE, backbone_color = "black", 
          markers = NULL, show_cell_names = FALSE, cell_size = 1.5, 
          cell_link_size = 0.75, cell_name_size = 2, legend_text_size = 12, show_branch_points = TRUE, branch_point_size = 2.0,
          theta = 0, axis_line_size = 0.5, axis_tick_size = 0.5) 
{
  gene_short_name <- NA
  sample_name <- NA
  data_dim_1 <- NA
  data_dim_2 <- NA
  lib_info_with_pseudo <- pData(cds)
  if (is.null(cds@dim_reduce_type)) {
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
  }
  if (cds@dim_reduce_type == "ICA") {
    reduced_dim_coords <- reducedDimS(cds)
  }
  else if (cds@dim_reduce_type %in% c("SimplePPT", "DDRTree", 
                                      "SGL-tree")) {
    reduced_dim_coords <- reducedDimK(cds)
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  }
  else {
    stop("Error: unrecognized dimensionality reduction method.")
  }
  if (is.null(reduced_dim_coords)) {
    stop("You must first call reduceDimension() before using this function")
  }
  dp_mst <- minSpanningTree(cds)
  if (is.null(root_states)) {
    if (is.null(lib_info_with_pseudo$Pseudotime)) {
      root_cell <- row.names(lib_info_with_pseudo)[degree(dp_mst) == 
                                                     1][1]
    }
    else root_cell <- row.names(subset(lib_info_with_pseudo, 
                                       Pseudotime == 0))
    if (cds@dim_reduce_type != "ICA") 
      root_cell <- igraph::V(dp_mst)$name[cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[root_cell, 
                                                                                                ]]
  }
  else {
    candidate_root_cells <- row.names(subset(pData(cds), 
                                             State %in% root_states))
    if (cds@dim_reduce_type == "ICA") {
      root_cell <- candidate_root_cells[which(degree(dp_mst, 
                                                     candidate_root_cells) == 1)]
    }
    else {
      Y_candidate_root_cells <- igraph::V(dp_mst)$name[cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[candidate_root_cells, 
                                                                                                             ]]
      root_cell <- Y_candidate_root_cells[which(degree(dp_mst, 
                                                       Y_candidate_root_cells) == 1)]
    }
  }
  tree_coords <- igraph::layout_as_tree(dp_mst, root = root_cell)
  ica_space_df <- data.frame(tree_coords)
  row.names(ica_space_df) <- colnames(reduced_dim_coords)
  colnames(ica_space_df) <- c("prin_graph_dim_1", "prin_graph_dim_2")
  ica_space_df$sample_name <- row.names(ica_space_df)
  if (is.null(dp_mst)) {
    stop("You must first call orderCells() before using this function")
  }
  edge_list <- as.data.frame(igraph::get.edgelist(dp_mst))
  colnames(edge_list) <- c("source", "target")
  edge_df <- merge(ica_space_df, edge_list, by.x = "sample_name", 
                   by.y = "source", all = TRUE)
  edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = "source_prin_graph_dim_1", 
                                     prin_graph_dim_2 = "source_prin_graph_dim_2"))
  edge_df <- merge(edge_df, ica_space_df[, c("sample_name", 
                                             "prin_graph_dim_1", "prin_graph_dim_2")], by.x = "target", 
                   by.y = "sample_name", all = TRUE)
  edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = "target_prin_graph_dim_1", 
                                     prin_graph_dim_2 = "target_prin_graph_dim_2"))
  if (cds@dim_reduce_type == "ICA") {
    S_matrix <- tree_coords[, ]
  }
  else if (cds@dim_reduce_type %in% c("DDRTree", "SimplePPT", 
                                      "SGL-tree")) {
    S_matrix <- tree_coords[closest_vertex, ]
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  }
  data_df <- data.frame(S_matrix)
  row.names(data_df) <- colnames(reducedDimS(cds))
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- row.names(data_df)
  data_df <- merge(data_df, lib_info_with_pseudo, by.x = "sample_name", 
                   by.y = "row.names")
  return_rotation_mat <- function(theta) {
    theta <- theta/180 * pi
    matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 
           nrow = 2)
  }
  tmp <- return_rotation_mat(theta) %*% t(as.matrix(data_df[, 
                                                            c(2, 3)]))
  data_df$data_dim_1 <- tmp[1, ]
  data_df$data_dim_2 <- tmp[2, ]
  tmp <- return_rotation_mat(theta = theta) %*% t(as.matrix(edge_df[, 
                                                                    c("source_prin_graph_dim_1", "source_prin_graph_dim_2")]))
  edge_df$source_prin_graph_dim_1 <- tmp[1, ]
  edge_df$source_prin_graph_dim_2 <- tmp[2, ]
  tmp <- return_rotation_mat(theta) %*% t(as.matrix(edge_df[, 
                                                            c("target_prin_graph_dim_1", "target_prin_graph_dim_2")]))
  edge_df$target_prin_graph_dim_1 <- tmp[1, ]
  edge_df$target_prin_graph_dim_2 <- tmp[2, ]
  markers_exprs <- NULL
  if (is.null(markers) == FALSE) {
    markers_fData <- subset(fData(cds), gene_short_name %in% 
                              markers)
    if (nrow(markers_fData) >= 1) {
      markers_exprs <- reshape2::melt(as.matrix(Biobase::exprs(cds[row.names(markers_fData), 
                                                          ])))
      colnames(markers_exprs)[1:2] <- c("feature_id", "cell_id")
      markers_exprs <- merge(markers_exprs, markers_fData, 
                             by.x = "feature_id", by.y = "row.names")
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 
      0) {
    data_df <- merge(data_df, markers_exprs, by.x = "sample_name", 
                     by.y = "cell_id")
    g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2, 
                                    size = log10(value + 0.1))) + facet_wrap(~feature_label)
  }
  else {
    g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
  }
  if (show_tree) {
    g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                     y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                     yend = "target_prin_graph_dim_2"), size = cell_link_size, 
                          linetype = "solid", na.rm = TRUE, data = edge_df)
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 
      0) {
    g <- g + geom_jitter(aes_string(color = color_by), na.rm = TRUE, 
                         height = 5)
  }
  else {
    g <- g + geom_jitter(aes_string(color = color_by), size = I(cell_size), 
                         na.rm = TRUE, height = 5)
  }
  if (show_branch_points && cds@dim_reduce_type == "DDRTree") {
    mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
    branch_point_df <- subset(edge_df, sample_name %in% mst_branch_nodes)[, 
                                                                          c("sample_name", "source_prin_graph_dim_1", "source_prin_graph_dim_2")]
    branch_point_df$branch_point_idx <- match(branch_point_df$sample_name, 
                                              mst_branch_nodes)
    branch_point_df <- branch_point_df[!duplicated(branch_point_df$branch_point_idx), 
                                       ]
    g <- g + geom_point(aes_string(x = "source_prin_graph_dim_1", 
                                   y = "source_prin_graph_dim_2"), size = branch_point_size * cell_size, 
                        na.rm = TRUE, data = branch_point_df) + geom_text(aes_string(x = "source_prin_graph_dim_1", 
                                                                                     y = "source_prin_graph_dim_2", label = "branch_point_idx"), 
                                                                          size = branch_point_size * cell_size*0.7, color = "white", na.rm = TRUE, 
                                                                          data = branch_point_df)
  }
  if (show_cell_names) {
    g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
  }
  g <- g + monocle:::monocle_theme_opts() + xlab(paste("Component", x)) + ylab(paste("Component", y)) 
  if(show_legend) {
  g <- g + theme(legend.position = "top", legend.key.height = grid::unit(0.35, "in"), legend.text = element_text(size=legend_text_size), legend.key.size = unit(0.5, "cm")) + theme(legend.key = element_blank()) + theme(panel.background = element_rect(fill = "white"))
  }
  
  g + theme(legend.position="none") + theme(axis.line.x = element_line(colour = 'black', size = axis_line_size),
            axis.line.y = element_line(colour = 'black', size = axis_line_size),
            axis.ticks = element_line(colour = "black", size = axis_tick_size))
}




















###



#plot_genes_branched_heatmap2

plot_genes_branched_heatmap2 <- function (cds_subset, branch_point = 1, branch_states = NULL, 
          branch_labels = c("Cell fate 1", "Cell fate 2"), cluster_rows = TRUE, 
          hclust_method = "ward.D2", num_clusters = 6, hmcols = NULL, 
          branch_colors = c("#979797", "#F05662", "#7990C8"), add_annotation_row = NULL, 
          add_annotation_col = NULL, show_rownames = FALSE, use_gene_short_name = TRUE, 
          scale_max = 3, scale_min = -3, norm_method = c("log", "vstExprs"), 
          trend_formula = "~sm.ns(Pseudotime, df=3) * Branch", return_heatmap = FALSE, 
          show_legend = TRUE, show_annotation_legend = TRUE, show_annotation_names_row = TRUE,
          show_annotation_names_col = TRUE, font_size = 12, cores = 1, ...) 
{
  cds <- NA
  new_cds <- buildBranchCellDataSet(cds_subset, branch_states = branch_states, 
                                    branch_point = branch_point, progenitor_method = "duplicate", 
                                    ...)
  new_cds@dispFitInfo <- cds_subset@dispFitInfo
  if (is.null(branch_states)) {
    progenitor_state <- subset(pData(cds_subset), Pseudotime == 
                                 0)[, "State"]
    branch_states <- setdiff(pData(cds_subset)$State, progenitor_state)
  }
  col_gap_ind <- 101
  newdataA <- data.frame(Pseudotime = seq(0, 100, length.out = 100), 
                         Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[1]))
  newdataB <- data.frame(Pseudotime = seq(0, 100, length.out = 100), 
                         Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[2]))
  BranchAB_exprs <- genSmoothCurves(new_cds[, ], cores = cores, 
                                    trend_formula = trend_formula, relative_expr = T, new_data = rbind(newdataA, 
                                                                                                       newdataB))
  BranchA_exprs <- BranchAB_exprs[, 1:100]
  BranchB_exprs <- BranchAB_exprs[, 101:200]
  common_ancestor_cells <- row.names(pData(new_cds)[pData(new_cds)$State == 
                                                      setdiff(pData(new_cds)$State, branch_states), ])
  BranchP_num <- (100 - floor(max(pData(new_cds)[common_ancestor_cells, 
                                                 "Pseudotime"])))
  BranchA_num <- floor(max(pData(new_cds)[common_ancestor_cells, 
                                          "Pseudotime"]))
  BranchB_num <- BranchA_num
  norm_method <- match.arg(norm_method)
  if (norm_method == "vstExprs") {
    BranchA_exprs <- vstExprs(new_cds, expr_matrix = BranchA_exprs)
    BranchB_exprs <- vstExprs(new_cds, expr_matrix = BranchB_exprs)
  }
  else if (norm_method == "log") {
    BranchA_exprs <- log10(BranchA_exprs + 1)
    BranchB_exprs <- log10(BranchB_exprs + 1)
  }
  heatmap_matrix <- cBind(BranchA_exprs[, (col_gap_ind - 1):1], 
                          BranchB_exprs)
  heatmap_matrix = heatmap_matrix[!apply(heatmap_matrix, 1, 
                                         sd) == 0, ]
  heatmap_matrix = Matrix::t(scale(Matrix::t(heatmap_matrix), 
                                   center = TRUE))
  heatmap_matrix = heatmap_matrix[is.na(row.names(heatmap_matrix)) == 
                                    FALSE, ]
  heatmap_matrix[is.nan(heatmap_matrix)] = 0
  heatmap_matrix[heatmap_matrix > scale_max] = scale_max
  heatmap_matrix[heatmap_matrix < scale_min] = scale_min
  heatmap_matrix_ori <- heatmap_matrix
  heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[, 
                                                            1]) & is.finite(heatmap_matrix[, col_gap_ind]), ]
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1
  exp_rng <- range(heatmap_matrix)
  bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by = 0.1)
  if (is.null(hmcols)) {
    hmcols <- monocle:::blue2green2red(length(bks) - 1)
  }
  ph <- pheatmap(heatmap_matrix, useRaster = T, cluster_cols = FALSE, 
                 cluster_rows = TRUE, show_rownames = F, show_colnames = F, 
                 clustering_distance_rows = row_dist, clustering_method = hclust_method, 
                 cutree_rows = num_clusters, silent = TRUE, filename = NA, 
                 breaks = bks, color = hmcols, legend = show_legend, annotation_legend = show_annotation_legend,
                 annotation_names_row = show_annotation_names_row, annotation_names_col = show_annotation_names_col,
                 fontsize = font_size)
  annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 
                                                       num_clusters)))
  if (!is.null(add_annotation_row)) {
    annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), 
                                                               ])
  }
  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
  annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)), 
                               `Cell Type` = c(rep(branch_labels[1], BranchA_num), rep("Pre-branch", 
                                                                                       2 * BranchP_num), rep(branch_labels[2], BranchB_num)))
  colnames(annotation_col) <- "Cell Type"
  if (!is.null(add_annotation_col)) {
    annotation_col <- cbind(annotation_col, add_annotation_col[fData(cds[row.names(annotation_col), 
                                                                         ])$gene_short_name, 1])
  }
  names(branch_colors) <- c("Pre-branch", branch_labels[1], 
                            branch_labels[2])
  annotation_colors = list(`Cell Type` = branch_colors)
  names(annotation_colors$`Cell Type`) = c("Pre-branch", branch_labels)
  if (use_gene_short_name == TRUE) {
    if (is.null(fData(cds_subset)$gene_short_name) == FALSE) {
      feature_label <- as.character(fData(cds_subset)[row.names(heatmap_matrix), 
                                                      "gene_short_name"])
      feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
      row_ann_labels <- as.character(fData(cds_subset)[row.names(annotation_row), 
                                                       "gene_short_name"])
      row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
    }
    else {
      feature_label <- row.names(heatmap_matrix)
      row_ann_labels <- row.names(annotation_row)
    }
  }
  else {
    feature_label <- row.names(heatmap_matrix)
    row_ann_labels <- row.names(annotation_row)
  }
  row.names(heatmap_matrix) <- feature_label
  row.names(annotation_row) <- row_ann_labels
  ph_res <- pheatmap(heatmap_matrix[, ], useRaster = T, cluster_cols = FALSE, 
                     cluster_rows = TRUE, show_rownames = show_rownames, show_colnames = F, 
                     clustering_distance_rows = row_dist, clustering_method = hclust_method, 
                     cutree_rows = num_clusters, annotation_row = annotation_row, 
                     annotation_col = annotation_col, annotation_colors = annotation_colors, 
                     gaps_col = col_gap_ind, treeheight_row = 20, breaks = bks, 
                     fontsize = font_size, color = hmcols, silent = TRUE, legend = show_legend,
                     annotation_legend = show_annotation_legend ,  annotation_names_row = show_annotation_names_row, 
                     annotation_names_col = show_annotation_names_col)
  grid::grid.rect(gp = grid::gpar("fill", col = NA))
  grid::grid.draw(ph_res$gtable)
  if (return_heatmap) {
    return(list(BranchA_exprs = BranchA_exprs, BranchB_exprs = BranchB_exprs, 
                heatmap_matrix = heatmap_matrix, heatmap_matrix_ori = heatmap_matrix_ori, 
                ph = ph, col_gap_ind = col_gap_ind, row_dist = row_dist, 
                hmcols = hmcols, annotation_colors = annotation_colors, 
                annotation_row = annotation_row, annotation_col = annotation_col, 
                ph_res = ph_res))
  }
}

###
