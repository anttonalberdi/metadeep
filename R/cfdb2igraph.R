#' Convert cross-feeding database to network (igraph) object
#' @title Convert cross-feeding database (cfdb) to network (igraph) object
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Convert cross-feeding database (cfdb) to network (igraph) object
#' @param cfdb A cross-feeding database (cfdb) generated using mdb2cfdb()
#' @param mode Whether to calculate forward (first to second), reverse (second to first) or total metabolite exchanges. Default is mode="total".
#' @import tidyverse igraph
#' @examples
#' cfdb2igraph(allgenomes_cfdb)
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

cfdb2igraph <- function(cfdb, mode="total") {
  # Input check
  if (inherits(cfdb, "cfdb")) {
    cfdb <- cfdb
  } else {
    stop("Input is not a valid cfdb object created by mdb2cfdb().")
  }

  #Generate pairwise matrix
  igraph <- cfdb2pair(cfdb, mode=mode) %>%
    select(-1) %>%
    as.matrix() %>%
    graph_from_adjacency_matrix(., weighted = TRUE, mode = "upper", diag = FALSE)

  #Output cross-feeding matrix
  return(igraph)
}

