#' Convert cross-feeding database to pairwise matrix
#' @title Convert cross-feeding database (cfdb) to pairwise matrix
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Convert cross-feeding database (cfdb) to pairwise matrix
#' @param cfdb A cross-feeding database (cfdb) generated using mdb2cfdb()
#' @param mode Whether to calculate forward (columns to rows), reverse (rows to columns) or total metabolite exchanges. Default is mode="total".
#' @import tidyverse
#' @examples
#' cfdb2pair(allgenomes_cfdb)
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

cfdb2pair <- function(cfdb, mode="total") {
  # Input check
  if (inherits(cfdb, "cfdb")) {
    cfdb <- cfdb
  } else {
    stop("Input is not a valid cfdb object created by mdb2cfdb().")
  }

  #Identify first and last genomes
  firstgenome <- allgenomes_cfdb %>% arrange(first) %>% head(1) %>% pull(first)
  lastgenome <- allgenomes_cfdb %>% arrange(second) %>% tail(1) %>% pull(second)

  #Generate exchange matrix
  pair <- cfdb %>%
    select(first,second,{{ mode }}) %>%
    rename(exchange=3) %>%
    mutate(exchange = length(exchange)) %>%
    arrange(first,second) %>%
    pivot_wider(names_from = second, values_from = exchange, values_fill = NA) %>%
    mutate(firstcol=NA, .before = 2) %>%
    rename(genomes=1) %>%
    add_row(genomes = NA) %>%
    rename(!! firstgenome := 2) %>%
    mutate(genomes = if_else(row_number() == n(), {{ lastgenome }}, genomes))

  #Output cross-feeding matrix
  return(pair)
}
