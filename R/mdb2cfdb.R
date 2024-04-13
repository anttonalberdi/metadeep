#' Generate crossfeeding database
#' @title Generate crossfeeding database from a multi-genome metabolite database (mdb)
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Generate crossfeeding database from a multi-genome metabolite database (mdb)
#' @param mdb A metabolite database (mdb) produced by rdb2mdb or rdbs2mdb
#' @param mode Whether to calculate strict or loose metabolite exchange capacities. Default is mode="strict".
#' @import tidyverse
#' @examples
#' mdb2cfdb(allgenomes_mdb)
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

mdb2cfdb <- function(mdb, mode="strict") {
  # Input check
  if (inherits(mdb, "mdb")) {
    mdb <- mdb
  } else if (inherits(mdb, "rdb")) {
    stop("The input file is a rdb object, run rdb2mdb() first to generate a mdb obnject")
  } else if (inherits(mdb, "rdbs")) {
    stop("The input file is a rdbs object, run rdbs2mdb() first to generate a mdb obnject")
  } else {
    stop("Input is not a valid rdbs object created by sbmls2rdbs().")
  }

  #Calculate strict cross-feeding
  if(mode == "strict"){
    cfdb <- combn(mdb$genome, 2, simplify = TRUE) %>%
      t() %>%
      as.data.frame() %>%
      tibble %>%
      rename(first=1,second=2) %>%
      rowwise() %>%
      mutate(forward=list(intersect(mdb %>% filter(genome == first) %>% select(sinks) %>% unlist(), mdb %>% filter(genome == second) %>% select(sources) %>% unlist()))) %>%
      mutate(reverse=list(intersect(mdb %>% filter(genome == second) %>% select(sinks) %>% unlist(), mdb %>% filter(genome == first) %>% select(sources) %>% unlist()))) %>%
      mutate(total = list(unique(c(forward, reverse))))
  }

  #Calculate loose cross-feeding
  if(mode == "loose"){
    cfdb <- combn(mdb$genome, 2, simplify = TRUE) %>%
      t() %>%
      as.data.frame() %>%
      tibble %>%
      rename(first=1,second=2) %>%
      rowwise() %>%
      mutate(forward=list(intersect(mdb %>% filter(genome == first) %>% select(transits,sinks) %>% unlist(), mdb %>% filter(genome == second) %>% select(sources) %>% unlist()))) %>%
      mutate(reverse=list(intersect(mdb %>% filter(genome == second) %>% select(sinks) %>% unlist(), mdb %>% filter(genome == first) %>% select(transits,sources) %>% unlist()))) %>%
      mutate(total = list(unique(c(forward, reverse))))
  }

  #Add class
  class(cfdb) <- c("cfdb", class(cfdb))

  #Output cross-feeding database
  return(cfdb)
}
