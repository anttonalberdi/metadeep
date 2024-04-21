#' Convert a reaction database into a network (igraph) object
#' @title Convert a reaction database into a network (igraph) object
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Convert cross-feeding database (cfdb) to network (igraph) object
#' @param rdb A reaction database (rdb object) produced by sbml2rdb().
#' @import tidyverse igraph
#' @examples
#' rdb2igraph(allgenomes_rdb)
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

rdb2igraph <- function(rdb) {

  # Declare function
  rdb2igraph_func <- function(rdb) {
    #Generate metabolite-to-metabolite relation matrix
    relations <- full_join(rdb %>% unnest(cols = reactants) %>% select(-products),
                           rdb %>% unnest(cols = products) %>% select(-reactants),
                           by="reaction",
                           relationship = "many-to-many") %>%
                           select(-reaction)

    #Calculate metabolite types
    gedb <- rdb2gedb(rdb)

    #Generate vector of vertices
    vertices <- tibble(name=c(relations$reactants,relations$products) %>% sort() %>% unique()) %>%
    mutate(type = case_when(
      name %in% gedb$sources[[1]] ~ "source",
      name %in% gedb$transits[[1]] ~ "transit",
      name %in% gedb$sinks[[1]] ~ "sink",
      TRUE ~ NA_character_  # if name does not match any condition, assign NA
    )) %>%
    mutate(color = case_when(
        type == "source" ~ "#D81B61",
        type == "transit" ~ "#3D82C4",
        type == "sink" ~ "#FEC111",
        TRUE ~ NA_character_  # for other cases, assign NA
      ))

    #Generate igraph network
    igraph <- graph_from_data_frame(relations, directed=TRUE, vertices=vertices)

    return(igraph)
  }

  #Apply function
  if (!inherits(rdb, "list")) {
    # If input is a single matrix
    igraph <- rdb2igraph_func(rdb)
  } else {
    # If input is a list of matrices
    igraph <- map(rdb, rdb2igraph_func)
  }

  return(igraph)
}
