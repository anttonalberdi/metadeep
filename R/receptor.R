#' Calculation of receptor potential
#' @title Calculation of receptor potential
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Calculation of receptor potential of each genome
#' @param exdb A exchange database (exdb) generated using mdb2exdb()
#' @param abundance An optional data frame containing relative abundance data of bacteria
#' @param focal An optional focal genome name or vector of genome names to calculate receptor potential
#' @import tidyverse
#' @examples
#' receptor(allgenomes_exdb)
#' receptor(allgenomes_exdb, abundance=genome_abundances)
#' receptor(allgenomes_exdb, abundance=genome_abundances, focal="genome1")
#' receptor(allgenomes_exdb, focal=c("genome1","genome2"))
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

receptor <- function(exdb, abundance, focal) {
  # Input check
  if (inherits(exdb, "exdb")) {
    exdb <- exdb
  } else {
    stop("Input is not a valid exdb object created by mdb2exdb().")
  }

  # If no abundance data is provided
  if(missing(abundance)){
    if(missing(focal)){
      forward <- exdb %>%
        filter(length(reverse) > 0) %>%
        group_by(first) %>%
        summarize(metabolites = list(unique(unlist(reverse))), donors_n=n()) %>%
        mutate(metabolites_n = map_int(metabolites, length)) %>%
        rename(genome=1)

      reverse <- exdb %>%
        filter(length(forward) > 0) %>%
        group_by(second) %>%
        summarize(metabolites = list(unique(unlist(forward))), donors_n=n()) %>%
        mutate(metabolites_n = map_int(metabolites, length)) %>%
        rename(genome=1)
    } else {
      forward <- exdb %>%
        filter(length(reverse) > 0) %>%
        filter(first %in% focal) %>%
        group_by(first) %>%
        summarize(metabolites = list(unique(unlist(reverse))), donors_n=n()) %>%
        mutate(metabolites_n = map_int(metabolites, length)) %>%
        rename(genome=1)

      reverse <- exdb %>%
        filter(length(forward) > 0) %>%
        filter(second %in% focal) %>%
        group_by(second) %>%
        summarize(metabolites = list(unique(unlist(forward))), donors_n=n()) %>%
        mutate(metabolites_n = map_int(metabolites, length)) %>%
        rename(genome=1)
    }

    receptor_potential <- bind_rows(forward,reverse) %>%
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

    if(missing(focal)){
      suppressWarnings({
        forward <- exdb %>%
          # Append abundance values of first and second genomes
          left_join(abundance %>%
                      rename_with(~ paste0(., "_firstX"), -1),
                    by=join_by(first==genome)) %>%
          left_join(abundance %>%
                      rename_with(~ paste0(., "_secondX"), -1),
                    by=join_by(second==genome)) %>%
          select(-forward,-total) %>% 			# Drop reverse and total columns
          unnest(cols = reverse) %>% 			# Expand metabolites
          pivot_longer(cols = contains(c("firstX", "secondX")), 	# Pivot longer to create one column per sample
                       names_to = c(".value", "sample"),
                       names_pattern = "(.+)_(firstX|secondX)",
                       values_drop_na = TRUE) %>%
          rename(metabolite=3)

        reverse <- exdb %>%
          # Append abundance values of first and second genomes
          left_join(abundance %>%
                      rename_with(~ paste0(., "_secondX"), -1),
                    by=join_by(first==genome)) %>%
          left_join(abundance %>%
                      rename_with(~ paste0(., "_firstX"), -1),
                    by=join_by(second==genome)) %>%
          select(-reverse,-total) %>% 			# Drop reverse and total columns
          unnest(cols = forward) %>% 			# Expand metabolites
          pivot_longer(cols = contains(c("firstX", "secondX")), 	# Pivot longer to create one column per sample
                       names_to = c(".value", "sample"),
                       names_pattern = "(.+)_(firstX|secondX)",
                       values_drop_na = TRUE) %>%
          rename(first=second,second=first,metabolite=3)
      })
    } else {
      suppressWarnings({
      forward <- exdb %>%
        filter(first %in% focal) %>%
        # Append abundance values of first and second genomes
        left_join(abundance %>%
                    rename_with(~ paste0(., "_firstX"), -1),
                  by=join_by(first==genome)) %>%
        left_join(abundance %>%
                    rename_with(~ paste0(., "_secondX"), -1),
                  by=join_by(second==genome)) %>%
        select(-forward,-total) %>% 			# Drop reverse and total columns
        unnest(cols = reverse) %>% 			# Expand metabolites
        pivot_longer(cols = contains(c("firstX", "secondX")), 	# Pivot longer to create one column per sample
                     names_to = c(".value", "sample"),
                     names_pattern = "(.+)_(firstX|secondX)",
                     values_drop_na = TRUE) %>%
        rename(metabolite=3)

      reverse <- exdb %>%
        filter(second %in% focal) %>%
        # Append abundance values of first and second genomes
        left_join(abundance %>%
                    rename_with(~ paste0(., "_secondX"), -1),
                  by=join_by(first==genome)) %>%
        left_join(abundance %>%
                    rename_with(~ paste0(., "_firstX"), -1),
                  by=join_by(second==genome)) %>%
        select(-reverse,-total) %>% 			# Drop reverse and total columns
        unnest(cols = forward) %>% 			# Expand metabolites
        pivot_longer(cols = contains(c("firstX", "secondX")), 	# Pivot longer to create one column per sample
                     names_to = c(".value", "sample"),
                     names_pattern = "(.+)_(firstX|secondX)",
                     values_drop_na = TRUE) %>%
        rename(first=second,second=first,metabolite=3)
      })
    }

    #Merge tables and calculate metabolite exchanges
    receptor_potential <- bind_rows(forward,reverse) %>%
      group_by(first,metabolite) %>%
      #summarise(across(where(is.double), ~sum(.x[sample == "secondX"]) / max(.x[sample == "firstX"]))) %>% #without topping
      summarise(across(where(is.double), ~ pmin(1, sum(.x[sample == "secondX"]) / max(.x[sample == "firstX"]))), .groups = "drop") %>%
      mutate_if(is.double, ~ifelse(is.nan(.), 0, .)) %>% #convert NaNs derived from n/0 to 0
      rename(genome=1) %>%
      rowwise() %>%
      group_by(genome) %>%
      summarise(across(where(is.double), sum))
  }

  #Output receptor potential table
  return(receptor_potential)
}
