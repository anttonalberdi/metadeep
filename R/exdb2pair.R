#' Convert exchange database to pairwise matrix
#' @title Convert exchange database (exdb) to pairwise matrix
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Convert exchange database (exdb) to pairwise matrix
#' @param exdb An exchange database (exdb) generated using medb2exdb()
#' @param metabolites An optional vector of metabolites for which to run the calculations. Default is all metabolites.
#' @param genomes An optional vector of genomes for which to run the calculations. Default is all metabolites.
#' @param samples An optional vector of samples for which to run the calculations. Default is all samples
#' @import tidyverse
#' @examples
#' exdb2pair(allgenomes_exdb)
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

exdb2pair <- function(exdb, metabolites, genomes, samples) {

  #If focal genomes are not provided, then use all genomes from medb
  if(missing(genomes)){genomes=unique(unlist(lapply(exdb, function(x) unique(c(x$donor,x$receptor))))) %>% sort()}

  #If focal metabolites are not provided, then use all metabolites from medb
  if(missing(metabolites)){metabolites=names(exdb)}

  #If focal samples are not provided, then use all samples from medb
  if(missing(samples)){samples=exdb[[1]] %>% names() %>% .[-c(1, 2)] }

  #Unlist, filter, group, split and convert to pairwise
  pair <- exdb %>%
    bind_rows(.id = "metabolite") %>%
    filter(metabolite %in% metabolites) %>%
    filter(donor %in% genomes) %>%
    filter(receptor %in% genomes) %>%
    select(all_of(c("donor","receptor",samples))) %>%
    pivot_longer(cols = starts_with("sample"), names_to = "sample", values_to="exchange") %>%
    group_by(sample,donor,receptor) %>%
    summarise(exchange=sum(exchange), .groups="drop") %>% #summarise per donor and receptor
    group_by(sample) %>%
    group_split() %>%
    set_names(samples) %>%
    map(~select(.x, -sample)) %>%
    map(function(tibble) {
      tibble %>%
        complete(donor, receptor, fill = list(exchange = 0)) %>%
        pivot_wider(names_from = receptor, values_from = exchange, values_fill = NA)
    })

  if(length(pair)==1){pair=pair[[1]]}

  #Output pairwise exchange matrix or matrices
  return(pair)

}
