#' Convert SBML to database
#' @title Convert SBML to reaction database
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Convert a SBML file into a tibble contain details of all reactions in a genome
#' @param sbml An object of class SBMLR or an sbml file.
#' @import tidyverse SBMLR
#' @examples
#' sbml2rdb("data/genome1.sbml")
#' readSBML("data/genome1.sbml") %>% sbml2rdb()
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

sbml2rdb <- function(sbml) {

  sbml2rdb_func <- function(sbml){
    reactions <- sbml$reactions
    rdb <- tibble(
      reaction = sapply(reactions, function(x) x$id),
      reactants = lapply(reactions, function(x) x$reactants),
      products = lapply(reactions, function(x) x$products)
    )
    return(rdb)
  }

  # Input check
  if (class(sbml) == "SBMLR") {
    rdb <- sbml2rdb_func(sbml)
  } else if (is.vector(sbml) && all(is.character(sbml)) && all(file.exists(sbml))) {
    sbml2 <- lapply(sbml, readSBML)
    rdb <- lapply(sbml2, sbml2rdb_func)
    names(rdb) <- gsub(".*/([^/]+)\\.sbml", "\\1", sbml)
  } else if (is.character(sbml) && file.exists(sbml)) {
    sbml <- readSBML(sbml)
    rdb <- sbml2rdb_func(sbml)
  } else {
    stop("Input is neither an SBMLR object nor, a vector of SBML files, or a validSBML file.")
  }

  # Output reaction database
  return(rdb)
}
