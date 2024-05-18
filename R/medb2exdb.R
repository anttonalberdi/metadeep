#' Generate exchange database
#' @title Generate exchange database from a metabolite database (medb)
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Generate exchange database from a multi-genome metabolite database (gedb)
#' @param medb A metabolite database (medb) produced by gedb2medb
#' @param mode Whether to calculate strict or loose metabolite exchange capacities. Default is mode="strict".
#' @param abundance An optional dataframe containing relative abundance data of bacteria
#' @param weighted Whether to calculate unweighted or weighted exchanges. Default is TRUE if abundances are provided, and FALSE otherwise.
#' @param metabolites An optional vector of metabolites for which to run the calculations. Default is all metabolites.
#' @param samples An optional vector of samples for which to run the calculations. Default is all samples
#' @param verbosity Whether to print the progress of the computation when using abundance data. Default=TRUE
#' @import tidyverse
#' @examples
#' medb2exdb(allgenomes_gedb)
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

medb2exdb <- function(medb, mode="strict", abundance, weighted, metabolites, samples, verbosity=TRUE){

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
  exchange_list <- medb %>%
    filter(metabolite %in% metabolites) %>%
    rowwise() %>%
    #select donors depending on exchange mode
    mutate(donor=case_when(
        mode=="strict" ~ list(unique(unlist(list(sinks)))),
        mode=="loose" ~ list(unique(unlist(list(sinks,transits)))),
        TRUE ~ NA)) %>%
    rename(receptor=sources) %>%
    select(metabolite,donor,receptor) %>%
    filter(length(donor) > 0 & length(receptor) > 0) %>% #keep only exchangable metabolites
    group_by(metabolite) %>%
    {unique_metabolites <- group_keys(.)$metabolite  # Extract unique metabolite names
      group_split(.) %>%                             # Split the data
        set_names(unique_metabolites)} %>%           # Assign metabolite names to list
    map(~select(.x, -metabolite))                    # Remove redundant column

  #Declare genome filtering function
  filter_genomes <- function(exchange, genomes) {
      exchange %>%
      mutate(donor = keep(donor, ~ .x %in% genomes),
             receptor = purrr::keep(receptor, ~ .x %in% genomes))
  }

  #Declare unweighted function without occurrence data (single calculation)
  Xmdr <- function(exchange){
    exchange %>%
      unnest(donor) %>%
      unnest(receptor) %>%
      arrange(donor,receptor) %>%
      #Dm: number of donors
      mutate(Dm=donor %>% unique() %>% length()) %>%
      #Rm: number of receptors
      mutate(Rm=receptor %>% unique() %>% length()) %>%
      #Effective exchange parameters
      mutate(Im=Dm*Rm) %>%
      mutate(Em=min(Dm/Rm,Rm/Dm)) %>%
      mutate(Pmdr=1/(Dm*Rm)) %>%
      #Effective exchange
      mutate(exchange=Im*Em*Pmdr) %>%
      select(donor,receptor,exchange) %>%
      ungroup()
  }

  #Declare unweighted function with occurrence data (multiple calculations)
  Xmdr_occurrence <- function(exchange, sample){
    suppressWarnings({
    exchange %>%
      unnest(donor) %>%
      unnest(receptor) %>%
      arrange(donor,receptor) %>%
      #Occurrence of donors
      left_join(occurrence %>% select(genome,{{ sample }}) %>% rename(Od=2),
                by=join_by(donor==genome)) %>%
      #Occurrence of receptors
      left_join(occurrence %>% select(genome,{{ sample }}) %>% rename(Or=2),
                by=join_by(receptor==genome)) %>%
      #Remove genomes with no occurrence
      filter(Od > 0 & Or > 0) %>%
      #Dm: number of donors
      mutate(Dm=donor %>% unique() %>% length()) %>%
      #Rm: number of receptors
      mutate(Rm=receptor %>% unique() %>% length()) %>%
      #Effective exchange parameters
      mutate(Im=Dm*Rm) %>%
      mutate(Em=min(Dm/Rm,Rm/Dm)) %>%
      mutate(Pmdr=1/(Dm*Rm)) %>%
      #Effective exchange
      mutate(!!sample :=Im*Em*Pmdr) %>%
      select(donor,receptor,{{ sample }}) %>%
      ungroup()
    })
  }

  #Declare weighted function
  wXmdr <- function(exchange, sample){
    suppressWarnings({
    exchange %>%
      unnest(donor) %>%
      unnest(receptor) %>%
      arrange(donor,receptor) %>%
      #Ad: relative abundances of donors
      left_join(abundance %>% select(genome,{{ sample }}) %>% rename(Ad=2),
                by=join_by(donor==genome)) %>%
      #Ad: relative abundances of receptors
      left_join(abundance %>% select(genome,{{ sample }}) %>% rename(Ar=2),
                by=join_by(receptor==genome)) %>%
      #Remove genomes with no occurrence
      filter(Ad > 0 & Ar > 0) %>%
      #Dm: number of donors
      mutate(Dm=donor %>% unique() %>% length()) %>%
      #Rm: number of receptors
      mutate(Rm=receptor %>% unique() %>% length()) %>%
      #dm: cumulative relative abundances of donors
      mutate(dm=
               abundance %>% select(genome,{{ sample }}) %>% rename(abun=2) %>% filter(genome %in% unique(donor)) %>% summarise(sum(abun)) %>% pull()) %>%
      #rm: cumulative relative abundances of receptors
      mutate(rm=
               abundance %>% select(genome,{{ sample }}) %>% rename(abun=2) %>% filter(genome %in% unique(receptor)) %>% summarise(sum(abun)) %>% pull()) %>%
      #Effective exchange parameters
      mutate(wIm=Dm*Rm*(dm+rm)) %>%
      mutate(wEm=min((Dm*dm)/(Rm*rm),(Rm*rm)/(Dm*dm))) %>%
      mutate(wPmdr=(Ad/dm) * (Ar/rm)) %>%
      #Effective exchange
      mutate(!!sample :=wIm*wEm*wPmdr) %>%
      select(donor,receptor,{{sample}}) %>%
      ungroup()
    })
  }

  #Calculate exchanges
  if(missing(abundance)){
    exchange_list=map(exchange_list, Xmdr)
  }else{
    n=0
    for(sample in samples){
      n=n+1
      if(verbosity==TRUE){message(paste0("Processing ",sample," (",n,"/",length(samples),")"))}
      if(n==1){
        if(weighted==FALSE){
          exchange_list=map(exchange_list, Xmdr_occurrence, sample = sample)
        }else{
          exchange_list=map(exchange_list, wXmdr, sample = sample)
        }
      }else{
        if(weighted==FALSE){
          exchange_list <- map2(exchange_list, map(exchange_list, Xmdr_occurrence, sample = sample), ~ left_join(.x, .y, by = join_by(donor==donor,receptor==receptor)))
        }else{
          exchange_list <- map2(exchange_list, map(exchange_list, wXmdr, sample = sample), ~ left_join(.x, .y, by = join_by(donor==donor,receptor==receptor)))
        }
      }
    }
  }

  #Remove metabolites with no exchanges
  exchange_list <- exchange_list %>%
    discard(~ nrow(.) == 0)

  #Convert NAs to zero
  exchange_list <- exchange_list %>%
    map(~ .x %>% replace(is.na(.), 0))

  return(exchange_list)
}
