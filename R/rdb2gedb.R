#' Generate genome database
#' @title Generate a genome database from a reaction database (rdb)
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product rdb
#' @description Classification of metabolites in a single genome into source, transit and sink metabolites
#' @param rdb An rdb object produced by sbml2rdb().
#' @import tidyverse
#' @examples
#' rdb2gedb(genome1_rdb)
#' sbml2rdb("data/genome1.sbml") %>% rdb2gedb()
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

rdb2gedb <- function(rdb) {

  #Declare function
  rdb2gedb_func <- function(rdb){
    # Extract metabolite types
    sources <- unique(unlist(rdb$reactants))[!unique(unlist(rdb$reactants)) %in% unique(unlist(rdb$products))]
    transits <- unique(intersect(unlist(rdb$reactants), unlist(rdb$products)))
    sinks <- unique(unlist(rdb$products))[!unique(unlist(rdb$products)) %in% unique(unlist(rdb$reactants))]
    rdb <- nrow(rdb)
    metabolites <- length(sources) + length(transits) + length(sinks)

    # Compile in a tibble
    gedb <- tibble(genome=deparse(substitute(rdb)),
                  sources=list(sources),
                  transits=list(transits),
                  sinks=list(sinks),
                  reactions=rdb,
                  metabolites=metabolites)

    return(gedb)
  }

  # Input check
  if (!inherits(rdb, "list")) {
    #If input is rdb
    gedb <- rdb2gedb_func(rdb)
    #Add class
    class(gedb) <- c("gedb", class(gedb))
  } else {
    #If input is rdbs
    gedb <- map_df(rdb, rdb2gedb_func) %>%
      mutate(genome=names(rdb))
  }

  #Add class
  class(gedb) <- c("gedb", class(gedb))

  # Output metabolite database
  return(gedb)
}
