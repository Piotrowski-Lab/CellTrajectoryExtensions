branch_nodes <- function(cds,reduction_method="UMAP"){
  g = principal_graph(cds)[[reduction_method]]
  branch_points <- which(igraph::degree(g) > 2)
  branch_points = branch_points[branch_points %in% root_nodes(cds, reduction_method) == FALSE]
  return(branch_points)
}

leaf_nodes <- function(cds,reduction_method="UMAP"){
  g = principal_graph(cds)[[reduction_method]]
  leaves <- which(igraph::degree(g) == 1)
  leaves = leaves[leaves %in% root_nodes(cds, reduction_method) == FALSE]
  return(leaves)
}

root_nodes <- function(cds, reduction_method="UMAP"){
  g = principal_graph(cds)[[reduction_method]]
  root_pr_nodes <- which(names(igraph::V(g)) %in%
                           cds@principal_graph_aux[[reduction_method]]$root_pr_nodes)
  names(root_pr_nodes) <-
    cds@principal_graph_aux[[reduction_method]]$root_pr_nodes
  return(root_pr_nodes)
}

# returns the pseudotime value of the branching point
get_branching_point <- function(seurat_obj, cds){
  x <- 1
  y <- 2
  ica_space_df <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst) %>%
    as.data.frame() %>%
    dplyr::select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>%
    dplyr::mutate(sample_name = rownames(.),
                  sample_state = rownames(.))

  #get embedding info
  embeddings <- as.data.frame(seurat_obj[["umap"]]@cell.embeddings)
  embeddings$Cell <- rownames(embeddings)
  #round
  #embeddings[,1:2] <- round(embeddings[,1:2],digits = 2)
  embeddings[,1:2] <- signif(embeddings[,1:2],digits = 3)

  #get branching info from monocle princle graph
  mst_branch_nodes <- branch_nodes(cds, reduction_method ="UMAP")
  branch_point_df <- ica_space_df %>%
    dplyr::slice(match(names(mst_branch_nodes), sample_name)) %>%
    dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))

  #identify branching point bt central + HC Lineage, #1
  branch_point_df <- branch_point_df %>% filter(branch_point_idx == "1")
  #round
  #branch_point_df[,1:2] <- round(branch_point_df[,1:2], digits = 2)
  branch_point_df[,1:2] <- signif(branch_point_df[,1:2], digits = 3)

  #retrieve pseudotime, cell, cell group info
  ptime_df <- data.frame(pseudotime = t(pseudotime(cds))[1,],
                         Cell = names(pseudotime(cds)),
                         cell_group = colData(cds)$cell.type.and.trt)
  #join embedding info with pseudotime info
  ptime_df <- inner_join(ptime_df, embeddings)

  branching_embedding <- ptime_df %>% filter(abs(UMAP_1 - as.numeric(branch_point_df$prin_graph_dim_1)) == min(abs(UMAP_1 - as.numeric(branch_point_df$prin_graph_dim_1)))) %>%
    filter(abs(UMAP_2 - as.numeric(branch_point_df$prin_graph_dim_2)) == min(abs(UMAP_2 - as.numeric(branch_point_df$prin_graph_dim_2))))

  #extract first option, retrieve only the pseudotime val
  branching_ptime <- branching_embedding[1,]

  return(branching_ptime)
}
