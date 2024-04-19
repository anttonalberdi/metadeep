#' Generate exchange database
#' @title Generate exchange database from a multi-genome metabolite database (gedb)
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Generate exchange database from a multi-genome metabolite database (gedb)
#' @param gedb A metabolite database (gedb) produced by rdb2gedb
#' @param mode Whether to calculate strict or loose metabolite exchange capacities. Default is mode="strict".
#' @import tidyverse
#' @examples
#' gedb2exdb(allgenomes_gedb)
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

gedb2exdb <- function(gedb, mode="strict") {
  # Input check
  if (inherits(gedb, "gedb")) {
    gedb <- gedb
  } else if (inherits(gedb, "rdb")) {
    stop("The input file is a rdb object, run rdb2gedb() first to generate a gedb obnject")
  } else if (inherits(gedb, "rdbs")) {
    stop("The input file is a rdbs object, run rdbs2gedb() first to generate a gedb obnject")
  } else {
    stop("Input is not a valid rdbs object created by sbmls2rdbs().")
  }

  #Calculate strict cross-feeding
  if(mode == "strict"){
    exdb <- combn(gedb$genome, 2, simplify = TRUE) %>%
      t() %>%
      as.data.frame() %>%
      tibble %>%
      rename(first=1,second=2) %>%
      rowwise() %>%
      mutate(forward=list(intersect(gedb %>% filter(genome == first) %>% select(sinks) %>% unlist(), gedb %>% filter(genome == second) %>% select(sources) %>% unlist()))) %>%
      mutate(reverse=list(intersect(gedb %>% filter(genome == second) %>% select(sinks) %>% unlist(), gedb %>% filter(genome == first) %>% select(sources) %>% unlist()))) %>%
      mutate(total = list(unique(c(forward, reverse))))
  }

  #Calculate loose cross-feeding
  if(mode == "loose"){
    exdb <- combn(gedb$genome, 2, simplify = TRUE) %>%
      t() %>%
      as.data.frame() %>%
      tibble %>%
      rename(first=1,second=2) %>%
      rowwise() %>%
      mutate(forward=list(intersect(gedb %>% filter(genome == first) %>% select(transits,sinks) %>% unlist(), gedb %>% filter(genome == second) %>% select(sources) %>% unlist()))) %>%
      mutate(reverse=list(intersect(gedb %>% filter(genome == second) %>% select(sinks) %>% unlist(), gedb %>% filter(genome == first) %>% select(transits,sources) %>% unlist()))) %>%
      mutate(total = list(unique(c(forward, reverse))))
  }

  #Add class
  class(exdb) <- c("exdb", class(exdb))

  #Output cross-feeding database
  return(exdb)
}
