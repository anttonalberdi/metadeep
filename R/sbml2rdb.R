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
  # Input check
  if (class(sbml) == "SBMLR") {
    reactions <- sbml$reactions
  } else if (is.character(sbml) && file.exists(sbml)) {
    reactions <- readSBML(sbml)$reactions
  } else {
    stop("Input is neither an SBMLR object nor a valid file.")
  }

  #Convert sbml to reaction database
  rdb <- tibble(
    reaction = sapply(reactions, function(x) x$id),
    reactants = lapply(reactions, function(x) x$reactants),
    products = lapply(reactions, function(x) x$products)
  )

  #Add class
  class(rdb) <- c("rdb", class(rdb))

  # Output reaction database
  return(rdb)
}
