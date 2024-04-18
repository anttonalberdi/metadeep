#' Calculation of donor potential
#' @title Calculation of donor potential
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Calculation of donor potential of each genome
#' @param exdb A exchange database (exdb) generated using mdb2exdb()
#' @param abundance An optional data frame containing relative abundance data of bacteria
#' @param focal An optional focal genome name or vector of genome names to calculate donor potential
#' @param verbosity Whether to print the progress of the computation when using abundance data. Default=TRUE
#' @import tidyverse
#' @examples
#' donor(allgenomes_exdb)
#' donor(allgenomes_exdb, abundance=genome_abundances)
#' donor(allgenomes_exdb, abundance=genome_abundances, focal="genome1")
#' donor(allgenomes_exdb, focal=c("genome1","genome2"))
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

donor <- function(exdb, abundance, focal, verbosity=TRUE) {

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

      #Metabolites each genome in column first can provide to the rest of genomes in column second
      first2seconds <- exdb %>%
        filter(length(forward) > 0) %>% # remove pairings with no exchange of metabolites
        filter(first %in% focal) %>%
        group_by(first) %>%
        summarize(metabolites = list(unique(unlist(forward))), receptors_n=n()) %>%
        mutate(metabolites_n = map_int(metabolites, length)) %>%
        rename(genome=1)

      #Metabolites each genome in column second can provide to the rest of genomes in column first
      second2firsts <- exdb %>%
        filter(length(reverse) > 0) %>% # remove pairings with no exchange of metabolites
        filter(second %in% focal) %>%
        group_by(second) %>%
        summarize(metabolites = list(unique(unlist(reverse))), receptors_n=n()) %>%
        mutate(metabolites_n = map_int(metabolites, length)) %>%
        rename(genome=1)

    donor_potential <- bind_rows(first2seconds,second2firsts) %>%
      rowwise() %>%
      group_by(genome) %>%
      summarize(metabolites = list(unique(unlist(metabolites))), receptors_n=sum(receptors_n)) %>%
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
    donor_potential <- tibble(genome=focal)
    n=0
    for (sample in samples){
      n=n+1
      if(verbosity==TRUE){message(paste0("Processing ",sample," (",n,"/",length(samples),")"))}

      #Filter by focal taxa and append abundances to exdb
      exdb_abun <- exdb %>%
        filter(first %in% focal, second %in% focal) %>%
        left_join(abundance %>% select(genome,{{ sample }}) %>% rename(abun_first=2),
                  by=join_by(first==genome)) %>%
        left_join(abundance %>% select(genome,{{ sample }}) %>% rename(abun_second=2),
                  by=join_by(second==genome))

      #Effective number of metabolites each genome in column first can provide to the rest of genomes in column second
      first2seconds <- exdb_abun %>%
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
        summarise(potential = pmin(1, max(abun[genome == "donor"]) / sum(abun[genome == "receptor"]))) %>%
        mutate(potential = if_else(is.na(potential), 0, potential)) %>% #convert NaNs derived from n/0 to 0
        rename(genome=1) %>%
        rowwise() %>%
        group_by(genome) %>%
        summarise(exchange=sum(potential), .groups = "drop")

      #Append sample to table
      donor_potential <- donor_potential %>% mutate(!!sample :=donor_potential_sample$exchange)
    }
  }

  #Output cross-feeding matrix
  return(donor_potential)
}
