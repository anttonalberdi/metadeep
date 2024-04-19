#' Calculation of receptor potential
#' @title Calculation of receptor potential
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Calculation of receptor potential of each genome
#' @param exdb A exchange database (exdb) generated using mdb2exdb()
#' @param abundance An optional data frame containing relative abundance data of bacteria
#' @param focal An optional focal genome name or vector of genome names to calculate receptor potential
#' @param verbosity Whether to print the progress of the computation when using abundance data. Default=TRUE
#' @import tidyverse
#' @examples
#' receptor(allgenomes_exdb)
#' receptor(allgenomes_exdb, abundance=genome_abundances)
#' receptor(allgenomes_exdb, abundance=genome_abundances, focal="genome1")
#' receptor(allgenomes_exdb, focal=c("genome1","genome2"))
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

receptor <- function(exdb, abundance, focal, verbosity=TRUE) {

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

    #Metabolites each genome in column first can receive from the rest of genomes in column second
      seconds2first <- exdb %>%
        filter(length(reverse) > 0) %>% # remove pairings with no exchange of metabolites
        filter(first %in% focal) %>%
        group_by(first) %>%
        summarize(metabolites = list(unique(unlist(reverse))), donors_n=n()) %>%
        mutate(metabolites_n = map_int(metabolites, length)) %>%
        rename(genome=1)

      #Metabolites each genome in column second can receive from the rest of genomes in column first
      firsts2second <- exdb %>%
        filter(length(forward) > 0) %>% # remove pairings with no exchange of metabolites
        filter(second %in% focal) %>%
        group_by(second) %>%
        summarize(metabolites = list(unique(unlist(forward))), donors_n=n()) %>%
        mutate(metabolites_n = map_int(metabolites, length)) %>%
        rename(genome=1)

      receptor_potential <- bind_rows(seconds2first,firsts2second) %>%
        rowwise() %>%
        group_by(genome) %>%
        summarize(metabolites = list(unique(unlist(metabolites))), donors_n=sum(donors_n)) %>%
        mutate(metabolites_n = map_int(metabolites, length))
  }

  # If abundance data is provided
  if(!missing(abundance)){

    #Force relative abundance transformation
    abundance <- abundance %>%
      mutate(across(where(is.double), ~ ./sum(.)))

    #Identify samples
    samples <- names(abundance)[2:length(names(abundance))]

    #Per-sample loop is used, because in real datasets multi-sample calculations become too memory-exhaustive
    receptor_potential <- tibble(genome=focal)
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
      receptor_potential <- receptor_potential %>% mutate(!!sample :=receptor_potential_sample$exchange)

      }
  }

  #Output receptor potential table
  return(receptor_potential)
}
