#' Convert multiple SBMLs to a reaction database
#' @title Convert multiple SBMLs into a multi-genome reaction database (rdbs)
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Convert multiple SBML files into a list of tibbles containing details of all reactions in all genomes
#' @param sbml A list of sbml files.
#' @import tidyverse SBMLR
#' @examples
#' sbmls2rdbs(c("data/genome1.sbml","data/genome2.sbml"))
#' list.files(path = "data", pattern = "\\.sbml$", full.names = TRUE) %>% sbmls2rdb()
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

sbmls2rdbs <- function(sbmls) {
  # Input check
  if (is.vector(sbmls) && all(is.character(sbmls)) && all(file.exists(sbmls))) {
    sbmls <- sbmls
  } else {
    stop("At least some of the input files do not exist.")
  }

  rdbs <- lapply(sbmls, sbml2rdb)
  names(rdbs) <- gsub(".*/([^/]+)\\.sbml", "\\1", sbmls)

  #Add class
  class(rdbs) <- c("rdbs", class(rdbs))

  # Output multi-genome reaction database
  return(rdbs)

}
