#' Classification of metabolites
#' @title Classification of metabolites of multiple genomes
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product rdb
#' @description Classification of metabolites in multiple genomes into source, transit and sink metabolites
#' @param rdb An rdb object produced by sbml2rdb().
#' @import tidyverse SBMLR
#' @examples
#' rdbs2mdbs(allgenomes_rdbs)
#' list.files(path = "data", pattern = "\\.sbml$", full.names = TRUE) %>% sbml2rdb() %>% rdbs2mdbs()
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export


rdbs2mdbs <- function(rdbs) {
  # Input check
  if (inherits(rdbs, "rdbs")) {
    rdbs <- rdbs
  } else if (inherits(rdbs, "rdb")) {
    stop("Input is a single-genome object, use rdb2mdb() instead.")
  } else {
    stop("Input is not a valid rdbs object created by sbmls2rdbs().")
  }

  #Create reaction database
  mdb <- map_df(rdbs, rdb2mdb) %>%
    mutate(genome=names(rdbs))

  # Output reaction database
  return(mdb)
}
