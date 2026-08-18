# STRING protein-protein interaction networks
#
# Part of the 'proteoSE' package. Grouped by theme during the 2026 refactor
# (see AUDIT.md). Function bodies are preserved verbatim from the original
# main.R / mlasse.txt; renames and cleanup follow in later phases.



#' Strind-db network for goi
#' 
#' Creates a network with known protein-protein interactions and maps expression values to nodes.
#'
#' @param genes List or vector of gene symbols
#' @param species ID of species (human = 9606)
#' @param expression_data dataframe containing a column "gene_names" and additional columns containing numeric values
#' @param expand if TRUE, the network is expanded with additional known interactors
#' @param add_nodes Number of additional interactor nodes to add when expanding
#'   the network (passed to the STRING query); \code{NA} to add none.
#' @param network_type either "functional" or "physical". If "physical", only genes with a physical interaction are connected (e.g. complexes)
#' @param score minimal score for a drawn interaction.
#' @param common_legend If TRUE, minimal and maximal values for the legend are calculated across all columns
#' @param node_deg_above Number specifiying which genes should be deleted from the network (e.g. if 0, all genes with zero interaction partners are removed)
#'
#' @importFrom dplyr mutate select distinct rename 
#' @importFrom rbioapi rba_string_map_ids rba_string_interaction_partners rba_string_interactions_network
#' @importFrom igraph graph_from_data_frame delete_vertices degree
#' @importFrom cowplot plot_grid
#'
#' @return string network with integrated numeric values (fold-changes, lfq-values, ...)
#' @export
#'
get_network <- function(genes, species = 9606, expression_data = NA, expand = FALSE, add_nodes = NA, network_type = "functional", score = 900, common_legend = FALSE, node_deg_above = NA){
  .require_pkg("ggraph", "draw STRING networks")

  # calculate min/max
  if(common_legend){
    min = min(expression_data[, -!is.numeric(expression_data)], na.rm = TRUE)
    max = max(expression_data[, -!is.numeric(expression_data)], na.rm = TRUE)
  } 
  
  # map gene names to stringIDs 
  prots <- genes %>% rbioapi::rba_string_map_ids(species=species)
  
  # If expand = TRUE, all interactions between input proteins and every other STRING protein are returned.
  # If expand = FALSE, interactions among the input set are retrieved. Additional proteins can be added via the add_nodes argument
  if(expand){
    int_net <- rbioapi::rba_string_interaction_partners(prots$stringId, 
                                                        species = 9606, 
                                                        required_score = score,
                                                        network_type = network_type)
  }else{
    int_net <- rbioapi::rba_string_interactions_network(prots$stringId, 
                                                        species = 9606, 
                                                        required_score = score, 
                                                        network_type = network_type, 
                                                        add_nodes = ifelse(is.na(add_nodes),
                                                                           0,
                                                                           add_nodes))
  } 
  
  if(ncol(int_net > 0)){
  
  # get stringIDs for all proteins of the expression table
  all_prots <- expression_data$gene_names %>% unlist() %>%
    rbioapi::rba_string_map_ids(species=species) %>% 
    merge(expression_data, 
          by.x = "queryItem", 
          by.y = "gene_names")
  
  # Create node table
  node_tbl <- int_net[,c("stringId_A", "preferredName_A")] %>% 
    dplyr::rename("stringId_B" = "stringId_A", "preferredName_B" = "preferredName_A") %>%
    rbind(int_net[,c("stringId_B", "preferredName_B")]) %>% 
    dplyr::distinct() %>% 
    dplyr::rename("clean_name" = "preferredName_B", "name" = "stringId_B") %>%
    merge(all_prots, by.x = "name", by.y = "stringId", all.x = TRUE) %>%
    dplyr::rename("gene_names" = "queryItem") %>%
    dplyr::select(name, clean_name, colnames(expression_data))
  
  # Create graph object
  g_obj <- igraph::graph_from_data_frame(int_net, node_tbl, directed = TRUE)
  
  if(!is.na(node_deg_above)){
    g_obj <- igraph::delete_vertices(g_obj , which(igraph::degree(g_obj)<=node_deg_above))
  }
  
  #Plot networks with expression data
  ps <- list()
  
  for(i in  colnames(expression_data[-grep("gene_names", colnames(expression_data))])){
    
    if(!common_legend){
      min = min(expression_data[,i], na.rm = TRUE)
      max = max(expression_data[,i], na.rm = TRUE)
    }
    
    ps[[i]] <- ggraph::ggraph(g_obj, layout = "stress")+
      ggraph::geom_edge_link0(aes(edge_width = score),
                      edge_colour="grey",
                      alpha=0.7) +
      ggraph::geom_node_point(aes_string(fill = i),
                      shape = 21, size=12) +
      ggraph::geom_node_text(aes(label = clean_name))+
      scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                           limits = c(min,max))+
      ggraph::scale_edge_width_continuous(range = c(0,1))+
      theme_void()+
      labs(title = i)
  }
  
  p <- cowplot::plot_grid(plotlist=ps)
  
  result = list("node_tbl" = node_tbl,
                "graph_obj" = g_obj,
                "plotlist" = ps, 
                "plot" = p,
                "prots" = prots,
                "params" = list("species" = species, 
                                "network_type" = network_type, 
                                "score" = score,
                                "common_legend" = common_legend, 
                                "node_deg_above" = node_deg_above,
                                "add_nodes" = add_nodes))
  
  return(result)
  }else{
    return(NULL)
  }
}

