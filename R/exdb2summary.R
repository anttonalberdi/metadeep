#' Summarise pairwise matrix or list of pairwise matrices
#' @title Summarise pairwise matrix or list of pairwise matrices
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Summarise pairwise matrix or list of pairwise matrices
#' @param exdb A exchange database (exdb) generated using medb2exdb()
#' @param summarise Whether to summarise genomes, metabolites or samples. Default is all three.
#' @param exchange Whether to calculate genome summaries as average, donor, receptor metabolite exchange capacities. Default is average.
#' @param metabolites An optional vector of metabolites for which to run the calculations. Default is all metabolites.
#' @param genomes An optional vector of genomes for which to run the calculations. Default is all genomes
#' @param samples An optional vector of samples for which to run the calculations. Default is all samples
#' @import tidyverse
#' @examples
#' exdb2summary(allgenomes_exchange_total)
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

exdb2summary <- function(exdb,summarise=c("genomes","metabolites","samples"), exchange="average", metabolites, genomes, samples) {

  #If focal genomes are not provided, then use all genomes from medb
  if(missing(genomes)){genomes=unique(unlist(lapply(exdb, function(x) unique(c(x$donor,x$receptor))))) %>% sort()}

  #If focal metabolites are not provided, then use all metabolites from medb
  if(missing(metabolites)){metabolites=names(exdb)}

  #If focal samples are not provided, then use all samples from medb
  if(missing(samples)){samples=exdb[[1]] %>% names() %>% .[-c(1, 2)] }

  #Declare output list
  summary <- list()

  #Unlist and filter
  exdb_filt <- exdb %>%
    bind_rows(.id = "metabolite") %>%
    filter(metabolite %in% metabolites) %>%
    filter(donor %in% genomes) %>%
    filter(receptor %in% genomes) %>%
    select(all_of(c("metabolite","donor","receptor",samples)))

  #If genomes need to be summarised
  if("genomes" %in% summarise){

    #If donor statistics are needed
    if("donor" %in% exchange){
      genomes_donor <- exdb_filt %>%
        group_by(donor) %>%
        summarise(across(where(is.double), \(x) sum(x, na.rm = TRUE)))
      if("receptor" %in% exchange | "average" %in% exchange){ #append to output
        summary$genomes$donor <- genomes_donor
      }else{
        summary$genomes <- genomes_donor
        }
    }

    #If receptor statistics are needed
    if("receptor" %in% exchange){
      genomes_receptor <- exdb_filt %>%
        group_by(receptor) %>%
        summarise(across(where(is.double), \(x) sum(x, na.rm = TRUE)))
      if("donor" %in% exchange | "average" %in% exchange){ #append to output
        summary$genomes$receptor <- genomes_receptor
      }else{
        summary$genomes <- genomes_receptor
      }
    }

    #If average exchange statistics are needed
    if("average" %in% exchange){
      genomes_average <- exdb_filt %>%
        group_by(donor,receptor) %>%
        summarise(across(where(is.double), \(x) sum(x, na.rm = TRUE)), .groups = "drop") %>%
        pivot_longer(cols=c(donor,receptor), names_to = "role",values_to = "genome") %>%
        group_by(genome) %>%
        summarise(across(where(is.double), \(x) sum(x, na.rm = TRUE) / 2), .groups = "drop")
      if("donor" %in% exchange | "receptor" %in% exchange){ #append to output
        summary$genomes$average <- genomes_average
      }else{
        summary$genomes <- genomes_average
      }
    }
  }

  #If samples need to be summarised
  if("metabolites" %in% summarise | "samples" %in% summarise){

    metabolites_total <- exdb_filt %>%
      group_by(metabolite)  %>%
      summarise(across(where(is.double), \(x) sum(x, na.rm = TRUE)))
    if("metabolites" %in% summarise){
      summary$metabolites <- metabolites_total
    }

  }

  #If samples need to be summarised
  if("samples" %in% summarise){

   samples_total <- metabolites_total %>%
      summarise(across(where(is.double), \(x) sum(x, na.rm = TRUE))) %>%
     pivot_longer(cols = everything(), names_to = "sample", values_to = "effective_exchange") %>%
     mutate(maximum_exchange=exdb_filt %>% summarise(across(where(is.double), ~ sum(. != 0))) %>% c() %>% unlist()) %>%
     mutate(exchange_maximisation=effective_exchange/maximum_exchange) %>%
     mutate(n_genomes=genomes_average %>% summarise(across(where(is.double), ~ sum(. != 0))) %>% c() %>% unlist()) %>%
     mutate(int_per_genome=maximum_exchange/n_genomes)
    summary$samples <- samples_total

  }

  if(length(summary)==1){summary=summary[[1]]}

  return(summary)

}
