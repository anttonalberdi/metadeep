#' Generate competition database
#' @title Generate competition database from a metabolite database (medb)
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Generate competition database from a multi-genome metabolite database (gedb)
#' @param medb A metabolite database (medb) produced by gedb2medb
#' @param abundance An optional dataframe containing relative abundance data of bacteria
#' @param weighted Whether to calculate unweighted or weighted exchanges. Default is TRUE if abundances are provided, and FALSE otherwise.
#' @param metabolites An optional vector of metabolites for which to run the calculations. Default is all metabolites.
#' @param samples An optional vector of samples for which to run the calculations. Default is all samples
#' @param verbosity Whether to print the progress of the computation when using abundance data. Default=TRUE
#' @import tidyverse
#' @examples
#' medb2codb(allgenomes_gedb)
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

medb2codb <- function(medb, abundance, weighted, metabolites, samples, verbosity=TRUE){

  #Declare weighted if not declared
  if(missing(weighted)){
    if(!missing(abundance)){
      weighted=TRUE
    }else{
      weighted=FALSE
    }
  }

  #Force relative abundance transformation
  if(!missing(abundance)){
    abundance <- abundance %>%
      mutate(across(where(is.double), ~ ./sum(.)))
  }

  #Convert relative abundances into presence/absence
  if(!missing(abundance) & weighted==FALSE){
    occurrence <- abundance %>%
      mutate_if(is.double, ~ if_else(. > 0, 1, 0))
  }

  # If no sample information is provided
  if(missing(samples) & !missing(abundance)){
  samples <- names(abundance)[2:length(names(abundance))]
  }

  #If focal metabolites are not provided, then use all metabolites from medb
  if(missing(metabolites)){metabolites=c(medb$metabolite)}

  #Create metabolite exchange list
  competition_list <- medb %>%
    filter(metabolite %in% metabolites) %>%
    rename(competitor=sources) %>%
    select(metabolite,competitor) %>%
    filter(lengths(competitor) > 1) %>% #keep only metabolites with >1 competing genomes
    group_by(metabolite) %>%
    {unique_metabolites <- group_keys(.)$metabolite  # Extract unique metabolite names
      group_split(.) %>%                             # Split the data
        set_names(unique_metabolites)} %>%           # Assign metabolite names to list
    map(~select(.x, -metabolite))                    # Remove redundant column

  #Declare genome filtering function
  filter_genomes <- function(competition, genomes) {
      competition %>%
      mutate(competitor = keep(competitor, ~ .x %in% genomes))
  }

  #Declare unweighted function without occurrence data (single calculation)
  Omc <- function(competition){
    competition %>%
      unnest(competitor) %>%
      #Cm: number of competitors
      mutate(Cm=competitor %>% unique() %>% length()) %>%
      #Effective competition parameters
      mutate(Pmc=1/Cm) %>%
      #Effective competition
      mutate(competition=Cm*Pmc) %>%
      select(competitor,competition) %>%
      ungroup()
  }

  #Declare unweighted function with occurrence data (multiple calculations)
  Omc_occurrence <- function(competition, sample){
    suppressWarnings({
      competition %>%
      unnest(competitor) %>%
      #Occurrence of competitors
      left_join(occurrence %>% select(genome,{{ sample }}) %>% rename(Oc=2),
                by=join_by(competitor==genome)) %>%
      #Remove genomes with no occurrence
      filter(Oc > 0) %>%
      #Cm: number of competitors
      mutate(Cm=competitor %>% unique() %>% length()) %>%
      #Effective competition parameters
      mutate(Pmc=1/Cm) %>%
      #Effective exchange
      mutate(!!sample :=Cm*Pmc) %>%
      select(competitor,{{ sample }}) %>%
      ungroup()
    })
  }

  #Declare weighted function
  wOmc <- function(competition, sample){
    suppressWarnings({
      competition %>%
        unnest(competitor) %>%
        #Occurrence of competitors
        left_join(abundance %>% select(genome,{{ sample }}) %>% rename(Ac=2),
                  by=join_by(competitor==genome)) %>%
        #Remove genomes with no occurrence
        filter(Ac > 0) %>%
        #Cm: number of competitors
        mutate(Cm=competitor %>% unique() %>% length()) %>%
        #cm: cumulative relative abundances of competitors
        mutate(cm=
               abundance %>% select(genome,{{ sample }}) %>% rename(abun=2) %>% filter(genome %in% unique(competitor)) %>% summarise(sum(abun)) %>% pull()) %>%
        #Effective competition parameters
        mutate(wCm=Cm*cm) %>%
        mutate(wPmc=1/(Ac/cm) / sum(1/(Ac/cm))) %>%
        #Effective competition
        mutate(!!sample :=ifelse(Cm>1,wCm*wPmc,0)) %>% # set to zero if no competitors remain after filtering
        select(competitor,{{sample}}) %>%
        ungroup()
    })
  }

  #Calculate exchanges
  if(missing(abundance)){
    competition_list=map(competition_list, Omc)
  }else{
    n=0
    for(sample in samples){
      n=n+1
      if(verbosity==TRUE){message(paste0("Processing ",sample," (",n,"/",length(samples),")"))}
      if(n==1){
        if(weighted==FALSE){
          competition_list=map(competition_list, Omc_occurrence, sample = sample)
        }else{
          competition_list=map(competition_list, wOmc, sample = sample)
        }
      }else{
        if(weighted==FALSE){
          competition_list <- map2(competition_list, map(competition_list, Omc_occurrence, sample = sample), ~ left_join(.x, .y, by = "competitor"))
        }else{
          competition_list <- map2(competition_list, map(competition_list, wOmc, sample = sample), ~ left_join(.x, .y, by = "competitor"))
        }
      }
    }
  }

  #Remove metabolites with no exchanges
  competition_list <- competition_list %>%
    discard(~ nrow(.) == 0)

  #Convert NAs to zero
  competition_list <- competition_list %>%
    map(~ .x %>% replace(is.na(.), 0))

  return(competition_list)
}
