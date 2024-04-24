#' Calculation of metabolite exchange potential of each genome
#' @title Calculation of metabolite exchange potential of each genome
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Calculation of metabolite exchange potential of each genome
#' @param exdb A exchange database (exdb) generated using mdb2exdb()
#' @param exchange Whether to calculate donor, receptor or total metabolite exchange capacities. Default is mode="total".
#' @param abundance An optional data frame containing relative abundance data of bacteria
#' @param focal An optional focal genome name or vector of genome names to calculate donor potential
#' @param verbosity Whether to print the progress of the computation when using abundance data. Default=TRUE
#' @import tidyverse
#' @examples
#' exdb2focal(allgenomes_exdb)
#' exdb2focal(allgenomes_exdb, abundance=genome_abundances)
#' exdb2focal(allgenomes_exdb, abundance=genome_abundances, focal="genome1")
#' exdb2focal(allgenomes_exdb, focal=c("genome1","genome2"))
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

exdb2focal <- function(exdb, exchange="total", abundance, focal, verbosity=TRUE) {

  # Input check
  if (inherits(exdb, "exdb")) {
    exdb <- exdb
  } else {
    stop("Input is not a valid exdb object created by mdb2exdb().")
  }

  #If focal genomes are not provided, then select all genomes in exdb
  if(missing(focal)){focal=c(exdb$first,exdb$second) %>% unique() %>% sort()}

  # If no abundance data is provided
  if(missing(abundance)){

    #Calculate donor potential
    if(exchange %in% c("donor","total")){
      #Metabolites each genome in column first can provide to the rest of genomes in column second
      first2seconds <- exdb %>%
        filter(length(forward) > 0) %>% # remove pairings with no exchange of metabolites
        filter(first %in% focal) %>%
        group_by(first) %>%
        summarize(metabolites = list(unique(unlist(forward))), receptors=list(second)) %>%
        rename(genome=1)

      #Metabolites each genome in column second can provide to the rest of genomes in column first
      second2firsts <- exdb %>%
        filter(length(reverse) > 0) %>% # remove pairings with no exchange of metabolites
        filter(second %in% focal) %>%
        group_by(second) %>%
        summarize(metabolites = list(unique(unlist(reverse))), receptors=list(first)) %>%
        rename(genome=1)

      donor_potential <- bind_rows(first2seconds,second2firsts) %>%
        rowwise() %>%
        group_by(genome) %>%
        summarize(metabolites = list(unique(unlist(metabolites))), receptors=list(unique(unlist(receptors))))
    }

    #Calculate receptor potential
    if(exchange %in% c("receptor","total")){
      #Metabolites each genome in column first can receive from the rest of genomes in column second
      seconds2first <- exdb %>%
        filter(length(reverse) > 0) %>% # remove pairings with no exchange of metabolites
        filter(first %in% focal) %>%
        group_by(first) %>%
        summarize(metabolites = list(unique(unlist(reverse))), donors=list(second)) %>%
        rename(genome=1)

      #Metabolites each genome in column second can receive from the rest of genomes in column first
      firsts2second <- exdb %>%
        filter(length(forward) > 0) %>% # remove pairings with no exchange of metabolites
        filter(second %in% focal) %>%
        group_by(second) %>%
        summarize(metabolites = list(unique(unlist(forward))), donors=list(first)) %>%
        rename(genome=1)

      receptor_potential <- bind_rows(seconds2first,firsts2second) %>%
        rowwise() %>%
        group_by(genome) %>%
        summarize(metabolites = list(unique(unlist(metabolites))), donors=list(unique(unlist(donors))))
    }

    if(exchange == "total"){
      exchange_potential <- donor_potential %>%
        left_join(receptor_potential,by=join_by(genome==genome)) %>%
        rename(metabolites_donated=metabolites.x,metabolites_received=metabolites.y)  %>%
        rowwise() %>%
        mutate(metabolites_exchanged=list(unique(unlist(list(metabolites_donated,metabolites_received))))) %>%
        mutate(exchangers=list(unique(unlist(receptors),unlist(donors)))) %>%
        ungroup()
    }

    #Generate output
    if(exchange == "donor"){focal <- donor_potential}
    if(exchange == "receptor"){focal <- receptor_potential}
    if(exchange == "total"){focal <- exchange_potential}

  }

  # If abundance data is provided
  if(!missing(abundance)){

    #Force relative abundance transformation
    abundance <- abundance %>%
      mutate(across(where(is.double), ~ ./sum(.)))

    #Identify samples
    samples <- names(abundance)[2:length(names(abundance))]
    #Per-sample loop is used, because in real datasets multi-sample calculations become too memory-exhaustive

    #Declare tibbles
    donor_potential <- tibble(genome=focal)
    receptor_potential <- tibble(genome=focal)

    #Iteracte across samples
    n=0
    for (sample in samples){
      n=n+1
      if(verbosity==TRUE){message(paste0("Processing ",sample," (",n,"/",length(samples),")"))}

      #Filter by focal taxa and append abundances to exdb
      exdb_abun <- exdb %>%
        left_join(abundance %>% select(genome,{{ sample }}) %>% rename(abun_first=2),
                  by=join_by(first==genome)) %>%
        left_join(abundance %>% select(genome,{{ sample }}) %>% rename(abun_second=2),
                  by=join_by(second==genome))

      #Calculate donor potential
      if(exchange %in% c("donor","total")){
        #Effective number of metabolites each genome in column first can provide to the rest of genomes in column second
        first2seconds <- exdb_abun %>%
          filter(first %in% focal) %>%
          select(-reverse,-total) %>% 		# Drop reverse and total columns
          rename(abun_donor=abun_first,abun_receptor=abun_second) %>% #Rename abundances of first/second to donor/receptor
          unnest(cols = forward) %>% 			# Expand metabolites
          pivot_longer(cols = c("abun_donor", "abun_receptor"), 	# Pivot longer to perform calculations
                       names_to = c(".value", "genome"),
                       names_pattern = "(.+)_(donor|receptor)",
                       values_drop_na = TRUE) %>%
          rename(donor=first,receptor=second,metabolite=3)

        #Effective number of metabolites each genome in column second can provide to the rest of genomes in column first
        second2firsts <- exdb_abun %>%
          filter(second %in% focal) %>%
          select(-forward,-total) %>% 			# Drop reverse and total columns
          rename(abun_receptor=abun_first,abun_donor=abun_second) %>% #Rename abundances of first/second to receptor/donor
          unnest(cols = reverse) %>% 			# Expand metabolites
          pivot_longer(cols = c("abun_donor", "abun_receptor"), 	# Pivot longer to perform calculations
                       names_to = c(".value", "genome"),
                       names_pattern = "(.+)_(donor|receptor)",
                       values_drop_na = TRUE) %>%
          rename(donor=second,receptor=first,metabolite=3)

        #Merge tables and calculate metabolite exchanges
        donor_potential_sample <- bind_rows(first2seconds,second2firsts) %>%
          group_by(donor,metabolite) %>%
          summarise(potential = pmin(1, max(abun[genome == "donor"]) / sum(abun[genome == "receptor"])), .groups = "drop") %>%
          mutate(potential = if_else(is.na(potential), 0, potential)) %>% #convert NaNs derived from n/0 to 0
          rename(genome=1) %>%
          rowwise() %>%
          group_by(genome) %>%
          summarise(exchange=sum(potential), .groups = "drop")

        #Append sample to table
        donor_potential <- donor_potential %>% mutate(!!sample :=donor_potential_sample$exchange)
      }

      if(exchange %in% c("receptor","total")){
        #Effective number of metabolites each genome in column first can receive from the rest of genomes in column second
        seconds2first <- exdb_abun %>%
          filter(first %in% focal) %>%
          select(-forward,-total) %>% 			# Drop reverse and total columns
          rename(abun_donor=abun_second,abun_receptor=abun_first) %>% #Rename abundances of first/second to donor/receptor
          unnest(cols = reverse) %>% 			# Expand metabolites
          pivot_longer(cols = contains(c("abun_donor", "abun_receptor")), 	# Pivot longer to create one column per sample
                       names_to = c(".value", "genome"),
                       names_pattern = "(.+)_(donor|receptor)",
                       values_drop_na = TRUE) %>%
          rename(donor=second,receptor=first,metabolite=3)

        #Effective number of metabolites each genome in column second can receive from the rest of genomes in column first
        firsts2second <- exdb_abun %>%
          filter(second %in% focal) %>%
          select(-reverse,-total) %>% 			# Drop reverse and total columns
          rename(abun_receptor=abun_second,abun_donor=abun_first) %>% #Rename abundances of first/second to receptor/donor
          unnest(cols = forward) %>% 			# Expand metabolites
          pivot_longer(cols = contains(c("abun_donor", "abun_receptor")), 	# Pivot longer to create one column per sample
                       names_to = c(".value", "genome"),
                       names_pattern = "(.+)_(donor|receptor)",
                       values_drop_na = TRUE) %>%
          rename(donor=first,receptor=second,metabolite=3)

        #Merge tables and calculate metabolite exchanges
        receptor_potential_sample <- bind_rows(seconds2first,firsts2second) %>%
          group_by(receptor,metabolite) %>%
          summarise(potential = pmin(1, sum(abun[genome == "donor"]) / max(abun[genome == "receptor"]))) %>%
          mutate(potential = if_else(is.na(potential), 0, potential)) %>% #convert NaNs derived from n/0 to 0
          rename(genome=1) %>%
          rowwise() %>%
          group_by(genome) %>%
          summarise(exchange=sum(potential), .groups = "drop")

        #Append sample to table
        receptor_potential <- receptor_potential %>% mutate(!!sample := receptor_potential_sample$exchange)
      }

      if(exchange == "total"){
        exchange_potential <- donor_potential %>%
          mutate(!!sample := donor_potential %>% pull({{ sample }}) + receptor_potential %>% pull({{ sample }}))

      }

    }

    #Generate output
    if(exchange == "donor"){focal <- donor_potential}
    if(exchange == "receptor"){focal <- receptor_potential}
    if(exchange == "total"){focal <- exchange_potential}

  }

  return(focal)

}

