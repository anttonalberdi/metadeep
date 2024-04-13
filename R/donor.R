#' Calculation of donor potential
#' @title Calculation of donor potential
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Calculation of donor potential of each genome
#' @param cfdb A cross-feeding database (cfdb) generated using mdb2cfdb()
#' @import tidyverse
#' @examples
#' donor(allgenomes_cfdb)
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

donor <- function(cfdb, mode="total") {
  # Input check
  if (inherits(cfdb, "cfdb")) {
    cfdb <- cfdb
  } else {
    stop("Input is not a valid cfdb object created by mdb2cfdb().")
  }

  forward <- cfdb %>%
    filter(length(forward) > 0) %>%
    group_by(first) %>%
    summarize(metabolites = list(unique(unlist(forward))), receptors_n=n()) %>%
    mutate(metabolites_n = map_int(metabolites, length)) %>%
    rename(genome=1)

  reverse <- cfdb %>%
    filter(length(reverse) > 0) %>%
    group_by(second) %>%
    summarize(metabolites = list(unique(unlist(reverse))), receptors_n=n()) %>%
    mutate(metabolites_n = map_int(metabolites, length)) %>%
    rename(genome=1)

  donor_potential <- bind_rows(forward,reverse) %>%
    rowwise() %>%
    group_by(genome) %>%
    summarize(metabolites = list(unique(unlist(metabolites))), receptors_n=sum(receptors_n)) %>%
    mutate(metabolites_n = map_int(metabolites, length))

  #Output cross-feeding matrix
  return(donor_potential)
}
