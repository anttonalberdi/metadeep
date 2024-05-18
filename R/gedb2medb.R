#' Generate metabolite database
#' @title Generate a genome database from a reaction database (rdb)
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product rdb
#' @description Classification of metabolites in a single genome into source, transit and sink metabolites
#' @param rdb An rdb object produced by sbml2rdb().
#' @import tidyverse
#' @examples
#' gedb2medb(allgenomes_gedb)
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

gedb2medb <- function(gedb) {

  medb <- gedb %>%
    select(genome,sources,transits,sinks) %>%
    pivot_longer(cols = c(sources, transits, sinks), names_to = "type", values_to = "metabolite") %>%
    unnest(cols = metabolite) %>%
    group_by(metabolite, type) %>%
    summarise(genomes=list(unique(genome)), .groups = "drop") %>%
    pivot_wider(names_from = type, values_from = genomes) %>%
    mutate(exchangeable = ifelse(lengths(sinks) > 0 & lengths(sources) > 0, "strict",
                                ifelse((lengths(transits) > 0 | lengths(sinks) > 0) & lengths(sources) > 0, "loose", "no"))) %>%
    mutate(competition=ifelse(lengths(sources) > 1,"yes","no"))

  return(medb)

}
