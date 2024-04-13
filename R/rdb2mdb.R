#' Classification of metabolites
#' @title Classification of metabolites of a single genome
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product rdb
#' @description Classification of metabolites in a single genome into source, transit and sink metabolites
#' @param rdb An rdb object produced by sbml2rdb().
#' @import tidyverse
#' @examples
#' rdb2mdb(genome1_rdb)
#' sbml2rdb("data/genome1.sbml") %>% rdb2mdb()
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

rdb2mdb <- function(rdb) {
  # Input check
  if (inherits(rdb, "rdb")) {
    rdb <- rdb
  } else if (inherits(rdb, "rdbs")) {
    stop("Input is a multi-genome object, use rdbs2mdbs() instead.")
  } else {
    stop("Input is not a valid rdb object created by sbml2rdb().")
  }

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

  #Add class
  class(mdb) <- c("mdb", class(mdb))

  # Output reaction database
  return(mdb)
}
