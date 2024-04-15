#' Convert a pairwise metabolite exchange matrix into a network (igraph) object
#' @title Convert a pairwise metabolite exchange matrix into a network (igraph) object
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Convert cross-feeding database (cfdb) to network (igraph) object
#' @param pair A pairwise metabolite exchange matrix or a list of pairwise metabolite exchange matrices
#' @param mode Whether to generate a directed ("forward" or "reverse") or undirected network. Default is mode="undirected".
#' @import tidyverse igraph
#' @examples
#' pair2igraph(allgenomes_pair_total)
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

pair2igraph <- function(pair, mode="undirected") {

  # Declare function
  pair2igraph_func <- function(mat, mode) {
    if(colnames(mat)[1] == "genomes"){
      mat <- mat %>% column_to_rownames(var="genomes")
    }
    igraph <- mat %>%
      as.matrix() %>%
      graph_from_adjacency_matrix(weighted = TRUE, mode = "upper", diag = FALSE)
    # Set directionality
    if (mode == "forward") {
      igraph <- as.directed(igraph, mode = "arbitrary")
    }
    if (mode == "reverse") {
      igraph <- as.directed(igraph, mode = "arbitrary")
      igraph <- reverse_edges(igraph, eids = E(igraph))
    }
    return(igraph)
  }


  if (!inherits(pair, "list")) {
    # If input is a single matrix
    igraph <- pair2igraph_func(pair,mode=mode)
  } else {
    # If input is a list of matrices
    igraph <- map(pair, pair2igraph_func, mode = mode)
  }

  #Output igraph network(s)
  return(igraph)
}

