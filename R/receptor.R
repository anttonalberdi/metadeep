#' Calculation of receptor potential
#' @title Calculation of receptor potential
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Calculation of receptor potential of each genome
#' @param cfdb A cross-feeding database (cfdb) generated using mdb2cfdb()
#' @import tidyverse
#' @examples
#' receptor(allgenomes_cfdb)
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

receptor <- function(cfdb, mode="total") {
  # Input check
  if (inherits(cfdb, "cfdb")) {
    cfdb <- cfdb
  } else {
    stop("Input is not a valid cfdb object created by mdb2cfdb().")
  }

  forward <- cfdb %>%
    filter(length(reverse) > 0) %>%
    group_by(first) %>%
    summarize(metabolites = list(unique(unlist(reverse))), donors_n=n()) %>%
    mutate(metabolites_n = map_int(metabolites, length)) %>%
    rename(genome=1)

  reverse <- cfdb %>%
    filter(length(forward) > 0) %>%
    group_by(second) %>%
    summarize(metabolites = list(unique(unlist(forward))), donors_n=n()) %>%
    mutate(metabolites_n = map_int(metabolites, length)) %>%
    rename(genome=1)

  receptor_potential <- bind_rows(forward,reverse) %>%
    rowwise() %>%
    group_by(genome) %>%
    summarize(metabolites = list(unique(unlist(metabolites))), donors_n=sum(donors_n)) %>%
    mutate(metabolites_n = map_int(metabolites, length))

  #Output cross-feeding matrix
  return(receptor_potential)
}
