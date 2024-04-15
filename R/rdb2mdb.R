#' Classification of metabolites
#' @title Classification of metabolites of a single genome
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product rdb
#' @description Classification of metabolites in a single genome into source, transit and sink metabolites
#' @param rdb An rdb or rdbs object produced by sbml2rdb().
#' @import tidyverse
#' @examples
#' rdb2mdb(genome1_rdb)
#' sbml2rdb("data/genome1.sbml") %>% rdb2mdb()
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

rdb2mdb <- function(rdb) {

  #Declare function
  rdb2mdb_func <- function(rdb){
    # Extract metabolite types
    sources <- unique(unlist(rdb$reactants))[!unique(unlist(rdb$reactants)) %in% unique(unlist(rdb$products))]
    transits <- unique(intersect(unlist(rdb$reactants), unlist(rdb$products)))
    sinks <- unique(unlist(rdb$products))[!unique(unlist(rdb$products)) %in% unique(unlist(rdb$reactants))]
    rdb <- nrow(rdb)
    metabolites <- length(sources) + length(transits) + length(sinks)

    # Compile in a tibble
    mdb <- tibble(genome=deparse(substitute(rdb)),
                  sources=list(sources),
                  transits=list(transits),
                  sinks=list(sinks),
                  reactions=rdb,
                  metabolites=metabolites)

    return(mdb)
  }

  # Input check
  if (!inherits(rdb, "list")) {
    #If input is rdb
    mdb <- rdb2mdb_func(rdb)
    #Add class
    class(mdb) <- c("mdb", class(mdb))
  } else {
    #If input is rdbs
    mdb <- map_df(rdb, rdb2mdb_func) %>%
      mutate(genome=names(rdb))
  }

  #Add class
  class(mdb) <- c("mdb", class(mdb))

  # Output metabolite database
  return(mdb)
}
