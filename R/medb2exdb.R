#' Generate exchange database
#' @title Generate exchange database from a metabolite database (medb)
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Generate exchange database from a multi-genome metabolite database (gedb)
#' @param medb A metabolite database (medb) produced by gedb2medb
#' @param mode Whether to calculate strict or loose metabolite exchange capacities. Default is mode="strict".
#' @param abundance An optional dataframe containing relative abundance data of bacteria
#' @param metabolites An optional vector of metabolites for which to run the calculations. Default is all metabolites.
#' @param genomes An optional vector of genomes for which to run the calculations. Default is all metabolites.
#' @param samples An optional vector of samples for which to run the calculations. Default is all samples
#' @param verbosity Whether to print the progress of the computation when using abundance data. Default=TRUE
#' @import tidyverse
#' @examples
#' medb2exdb(allgenomes_gedb)
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

medb2exdb <- function(medb, mode="strict", abundance, metabolites, genomes, samples, verbosity=TRUE){

  #If focal genomes are not provided, then use all genomes from medb
  if(missing(genomes)){genomes=c(unlist(medb$sinks),unlist(medb$transits),unlist(medb$sources)) %>%
      unique() %>% sort()}

  # If no abundance data is provided
  if(missing(abundance)){
    abundance <- tibble(genome=genomes,exchange=rep(1/length(genomes),length(genomes)))
  }

  # If no sample information is provided
  if(missing(samples)){
  samples <- names(abundance)[2:length(names(abundance))]
  }

  #Select strict or loose mode
  if(mode == "strict"){medb <- medb %>% filter(exchangeable == "strict")}
  if(mode == "loose"){medb <- medb %>% filter(exchangeable != "no")}

  #If focal metabolites are not provided, then use all metabolites from medb
  if(missing(metabolites)){metabolites=c(medb$metabolite)}

  filter_genomes <- function(medb, genomes) {
        medb %>%
        mutate(donor = keep(donor, ~ .x %in% genomes),
               receptor = purrr::keep(receptor, ~ .x %in% genomes))
    }

  #Create metabolite exchange list
  exchange_list <- medb %>%
    filter(metabolite %in% metabolites) %>%
    rowwise() %>%
    mutate(donor=list(unique(unlist(list(sinks,transits))))) %>%
    select(-c(exchangeable,sinks,transits)) %>% # not necessary anymore
    rename(receptor=sources) %>%
    select(metabolite,donor,receptor) %>%
    group_by(metabolite) %>%
    group_split() %>%
    set_names(medb %>% pull(metabolite)) %>%
    map(~select(.x, -metabolite))

  #Declare function
  calculate_exchange <- function(exchange, sample){
    exchange %>%
      unnest(donor) %>%
      unnest(receptor) %>%
      arrange(donor,receptor) %>%
      #get relative abundances of donors
      left_join(abundance %>% select(genome,{{ sample }}) %>% rename(abun_donor=2),
                by=join_by(donor==genome)) %>%
      #get relative abundances of receptors
      left_join(abundance %>% select(genome,{{ sample }}) %>% rename(abun_receptor=2),
                by=join_by(receptor==genome)) %>%
      #calculate number of donors
      mutate(donor_n=length(unique(donor))) %>%
      #calculate cumulative relative abundance of donors
      mutate(abun_total_donor=
               abundance %>% select(genome,{{ sample }}) %>% rename(abun=2) %>% filter(genome %in% unique(donor)) %>% summarise(sum(abun)) %>% pull()) %>%
      #calculate cumulative relative abundance of receptors
      mutate(abun_total_receptor=
               abundance %>% select(genome,{{ sample }}) %>% rename(abun=2) %>% filter(genome %in% unique(receptor)) %>% summarise(sum(abun)) %>% pull()) %>%
      #calculate exchange ratio to penalise flow estimates
      mutate(proportion=ifelse(abun_total_donor/abun_total_receptor>1,abun_total_donor/abun_total_receptor,abun_total_receptor/abun_total_donor)) %>%
      rowwise() %>%
      # NUMBER_OF_DONORS * RELATIVE_PROPORTION_OF_DONOR * RELATIVE_PROPORTION_OF_RECEPTOR / |DONOR-RECEPTOR|_PROPORTION
      mutate(!!sample := (donor_n * (abun_donor/abun_total_donor) * (abun_receptor/abun_total_receptor)) / proportion) %>%
      select(donor,receptor,{{sample}}) %>%
      ungroup()
  }

  #Calculate exchanges
  n=0
  for(sample in samples){
    n=n+1
    if(verbosity==TRUE){message(paste0("Processing ",sample," (",n,"/",length(samples),")"))}
    if(n==1){
      exchange_list=map(exchange_list, calculate_exchange, sample = sample)
    }else{
      exchange_list <- map2(exchange_list, map(exchange_list,calculate_exchange, sample = sample), ~ left_join(.x, .y, by = join_by(donor==donor,receptor==receptor)))
    }
  }

  return(exchange_list)
}
