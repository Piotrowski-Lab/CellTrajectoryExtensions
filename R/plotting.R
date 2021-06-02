#' Plotting a line graph expression dynamic grouped by Modules
#'
#' This function loads a CDS obj from Monocle3 and builds a df
#' that contains cells ordered by pseudotime and calculates the average expression
#' of genes in each module along pseudotime, returns a line graph
#'
#'
#' @param cds CDS obj generated from Monocle3
#' @param module_tbl_complete meta table of genes in each module
#' @return a line graph of class ggplot2
#' @export

module_line_dynamic_plot <- function(cds, module_tbl_complete){
  meta <- as.data.frame(monocle3::pData(cds))
  meta$Cell <- rownames(meta)

  cds_sub <- cds[rownames(cds) %in% module_tbl_complete$Gene.name.uniq,]
  cds_exprs <- as.matrix(SingleCellExperiment::counts(cds_sub))

  #perserve gene order from module table
  count_mtx_sub <- cds_exprs[
    match(module_tbl_complete$Gene.name.uniq, rownames(cds_exprs)),]
  #perserve cell order based off pseudomtime
  count_mtx_sub <- count_mtx_sub[,match(
    names(pseudotime(cds)), colnames(count_mtx_sub))]

  #remove genes with zero expression across all cells
  #will cause NaN during scaling it unomitted
  count_mtx_sub=count_mtx_sub[!apply(count_mtx_sub,1,sum)==0,]
  count_mtx_sub <- t(as.matrix(count_mtx_sub))
  count_mtx_sub <- scale(log1p(count_mtx_sub)) #natural log with pseudocount of 1
  count_mtx_sub <- as.data.frame(Seurat::MinMax(count_mtx_sub , min = -2.5, max = 2.5))
  mod1_df <- as.data.frame(count_mtx_sub)
  mod1_df$Cell <- rownames(mod1_df)

  mod1_df <- reshape2::melt(mod1_df)
  colnames(mod1_df)[2:3] <- c("Gene.name.uniq","expression")

  mod1_df <- merge(mod1_df, module_tbl_complete[,1:2], by = "Gene.name.uniq", all.x =TRUE, all.y = FALSE, sort = FALSE)


  #modified module_tbl_complete xslx
  genes_to_keep <- unique(mod1_df$Gene.name.uniq)
  rownames(module_tbl_complete) <- module_tbl_complete$Gene.name.uniq
  module_tbl_complete <- module_tbl_complete[module_tbl_complete$Gene.name.uniq %in% genes_to_keep,]

  ptime_df <- data.frame(pseudotime = t(pseudotime(cds))[1,],
                         Cell = names(pseudotime(cds)),
                         cell_group = colData(cds)$cell.type.and.trt)

  #order from start of pseudotime to end
  ptime_df <- ptime_df[order(ptime_df$pseudotime),]

  #preserve ordered pseudotime
  ptime_df <- dplyr::inner_join(ptime_df, meta[,c("Cell", "trajectory")],
                                by = "Cell")

  #order cells in correspondence to ptime
  mod1_df$Cell <- factor(mod1_df$Cell, levels = ptime_df$Cell)

  plot_dt <- dplyr::inner_join(mod1_df, ptime_df)

  #count how many genes map to each module
  module_count <- module_tbl_complete %>% dplyr::count(
    module, sort =TRUE, name = "n_module")

  plot_dt <- data.table::setDT(inner_join(plot_dt, module_count))

  plot_dt$mod_label <- factor(paste0(
    "Module ", plot_dt$module, " (", plot_dt$n_module, " Genes)"))

  plot_dt$mod_label <- forcats::fct_reorder(plot_dt$mod_label, as.numeric(plot_dt$module), min)

  #means_dt  <-   plot_dt[, list(value=mean(expression)), by=list(pseudotime, mod_label)]

  means_dt <- plot_dt %>%
    dplyr::group_by(mod_label, pseudotime, Cell, cell_group, trajectory,module) %>%
    dplyr::summarize(mean = mean(expression)) #gives same output as line above

  lg <- ggplot(means_dt) +
    geom_smooth(
      colour='black',
      span=0.2,
      method='loess',
      aes( x=pseudotime, y=mean))

  #rug_dt    <- unique(plot_dt[, list(pseudotime, Cell, cell_group)])

  rug_dt <- plot_dt %>% dplyr::distinct(pseudotime, Cell, cell_group) #same result as above

  lg <- lg + geom_rug(data=rug_dt, sides='b', alpha=.10, aes(x=pseudotime,
                                                             color = cell_group) )
  cell_type_trt <- levels(cds@colData$cell.type.and.trt)
  type_trt_cols <- gg_color_hue(length(cell_type_trt))
  # change z axis tick order, correspond with umap colors
  lg <- lg + scale_color_manual(values = type_trt_cols)

  lg <- lg + facet_grid( mod_label ~ ., scales='free_y' ) +
    theme_bw() +
    theme(
    ) +
    labs(
      x     = 'pseudotime'
      ,y    = 'z-scored gene expression'
    )

  # change legend layout, override alpha gradient from .1 to 1.0
  lg <- lg +theme(legend.position="bottom", legend.title=element_blank()) +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))

  return(lg)
}

#build ggplot2 plot
make_plot_df <- function(cds_sub, gene){
  meta <- as.data.frame(monocle3::pData(cds_sub))
  meta$Cell <- rownames(meta)

  #remove genes that are not in the expr mtx
  gene <- gene[gene %in% rownames(cds_sub@assays@data@listData$counts)]

  #check for duplicates
  gene <- gene[!duplicated(gene)]

  cds_exprs_mtx <- cds_sub@assays@data@listData$counts
  cds_exprs_mtx <- cds_exprs_mtx[rownames(cds_exprs_mtx) %in% gene,]


  ptime_df <- data.frame(pseudotime = t(pseudotime(cds_sub))[1,],
                         Cell = names(pseudotime(cds_sub)),
                         cell_group = colData(cds_sub)$cell.type.and.trt)

  #order from start of pseudotime to end
  ptime_df <- ptime_df[order(ptime_df$pseudotime),]

  #preserve ordered pseudotime
  ptime_df <- dplyr::inner_join(ptime_df, meta[,c("Cell", "trajectory")], by = "Cell")

  #perserve gene order from module table
  count_mtx_sub <- cds_exprs_mtx[
    match(gene, rownames(cds_exprs_mtx)),]
  #perserve cell order based off pseudomtime
  count_mtx_sub <- count_mtx_sub[,match(
    ptime_df$Cell, colnames(count_mtx_sub))]


  #remove genes with zero expression across all cells
  #will cause NaN during scaling it unomitted
  count_mtx_sub <- count_mtx_sub[!apply(count_mtx_sub,1,sum)==0,]
  count_mtx_sub <- t(as.matrix(count_mtx_sub))
  count_mtx_sub <- base::scale(log1p(count_mtx_sub))
  count_mtx_sub <- as.data.frame(Seurat::MinMax(count_mtx_sub , min = -2.5, max = 2.5))
  mod1_df <- as.data.frame(count_mtx_sub)
  mod1_df$Cell <- rownames(mod1_df)


  mod1_df <- reshape2::melt(mod1_df) #428736
  colnames(mod1_df)[2:3] <- c("Gene.name.uniq","expression")

  plot_dt <- dplyr::inner_join(mod1_df,ptime_df)


  return(plot_dt)
}

restore_traj_line <-function(cds, plot_obj){
  #  import principle graph df, root node, branch df,
  # leaf nodes df
  x <- 1
  y <- 2
  ica_space_df <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst) %>%
    as.data.frame() %>%
    dplyr::select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>%
    dplyr::mutate(sample_name = rownames(.),
                  sample_state = rownames(.))

  dp_mst <- cds@principal_graph[["UMAP"]]

  edge_df <- dp_mst %>%
    igraph::as_data_frame() %>%
    dplyr::select_(source = "from", target = "to") %>%
    dplyr::left_join(ica_space_df %>%
                       dplyr::select_(
                         source="sample_name",
                         source_prin_graph_dim_1="prin_graph_dim_1",
                         source_prin_graph_dim_2="prin_graph_dim_2"),
                     by = "source") %>%
    dplyr::left_join(ica_space_df %>%
                       dplyr::select_(
                         target="sample_name",
                         target_prin_graph_dim_1="prin_graph_dim_1",
                         target_prin_graph_dim_2="prin_graph_dim_2"),
                     by = "target")

  mst_branch_nodes <- CellTrajectoryExtensions::branch_nodes(cds, reduction_method ="UMAP")
  branch_point_df <- ica_space_df %>%
    dplyr::slice(match(names(mst_branch_nodes), sample_name)) %>%
    dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))

  mst_leaf_nodes <- CellTrajectoryExtensions::leaf_nodes(cds, reduction_method = "UMAP")
  leaf_df <- ica_space_df %>%
    dplyr::slice(match(names(mst_leaf_nodes), sample_name)) %>%
    dplyr::mutate(leaf_idx = seq_len(dplyr::n()))

  mst_root_nodes <- CellTrajectoryExtensions::root_nodes(cds, reduction_method = "UMAP")
  root_df <- ica_space_df %>%
    dplyr::slice(match(names(mst_root_nodes), sample_name)) %>%
    dplyr::mutate(root_idx = seq_len(dplyr::n()))

  #import, plot traj line
  plot_obj <- plot_obj +
    geom_segment(
      aes_string(x="source_prin_graph_dim_1",
                 y="source_prin_graph_dim_2",
                 xend="target_prin_graph_dim_1",
                 yend="target_prin_graph_dim_2"),
      linetype="solid",
      na.rm=TRUE,
      data=edge_df)

  #plot branching
  plot_obj <- plot_obj +
    geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
               shape = 21, stroke=I(0.75),
               color="white",
               fill="black",
               size=I(as.numeric(2) * 1.5),
               na.rm=TRUE, branch_point_df) +
    geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                         label="branch_point_idx"),
              size=I(as.numeric(5)),
              color="white",
              na.rm=TRUE,
              branch_point_df)

  #plot leaves
  plot_obj <- plot_obj +
    geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
               shape = 21, stroke=I(0.75),
               color="black",
               fill="lightgray",
               size=I(as.numeric(2) * 1.5),
               na.rm=TRUE,
               leaf_df) +
    geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                         label="leaf_idx"),
              size=I(as.numeric(2)),
              color="black",
              na.rm=TRUE,
              leaf_df)

  plot_obj <- plot_obj + #plot root node
    geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
               shape = 21, stroke=I(0.75),
               color="black",
               fill="white",
               size=I(as.numeric(2) * 1.5),
               na.rm=TRUE,
               root_df) +
    geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                         label="root_idx"),
              size=I(as.numeric(2)),
              color="black",
              na.rm=TRUE,
              root_df)

  return(plot_obj)
}
