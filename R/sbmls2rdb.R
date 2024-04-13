#' Convert multiple SBMLs to a reaction database
#' @title Convert multiple SBMLs to a reaction database
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Convert multiple SBML files into a list of tibbles containing details of all reactions in each genome
#' @param sbml A list of sbml files.
#' @import tidyverse SBMLR
#' @examples
#' sbmls2rdb(c("data/genome1.sbml","data/genome2.sbml"))
#' list.files(path = "data", pattern = "\\.sbml$", full.names = TRUE) %>% sbmls2rdb()
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

sbmls2rdb <- function(sbmls) {
  # Input check
  if (is.vector(sbmls) && all(is.character(sbmls)) && all(file.exists(sbmls))) {
    sbmls <- sbmls
  } else {
    stop("At least some of the input files do not exist.")
  }

  rdb <- lapply(sbmls, sbml2rdb)
  names(rdb) <- gsub(".*/([^/]+)\\.sbml", "\\1", sbmls)

  # Output multi-genome reaction database
  return(rdb)

}
