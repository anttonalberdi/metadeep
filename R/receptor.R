#' Calculation of receptor potential
#' @title Calculation of receptor potential
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Calculation of receptor potential of each genome
#' @param cfdb A cross-feeding database (cfdb) generated using mdb2cfdb()
#' @param abundance An optional data frame containing relative abundance data of bacteria
#' @param focal An optional focal genome name or vector of genome names to calculate receptor potential
#' @import tidyverse
#' @examples
#' receptor(allgenomes_cfdb)
#' receptor(allgenomes_cfdb, abundance=genome_abundances)
#' receptor(allgenomes_cfdb, abundance=genome_abundances, focal="genome1")
#' receptor(allgenomes_cfdb, focal=c("genome1","genome2"))
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

receptor <- function(cfdb, abundance, focal) {
  # Input check
  if (inherits(cfdb, "cfdb")) {
    cfdb <- cfdb
  } else {
    stop("Input is not a valid cfdb object created by mdb2cfdb().")
  }

  # If no abundance data is provided
  if(missing(abundance)){
    if(missing(focal)){
      forward <- cfdb %>%
        filter(length(reverse) > 0) %>%
        group_by(first) %>%
        summarize(metabolites = list(unique(unlist(reverse))), donors_n=n()) %>%
        mutate(metabolites_n = map_int(metabolites, length)) %>%
        rename(genome=1)

      reverse <- cfdb %>%
        filter(length(forward) > 0) %>%
        group_by(second) %>%
        summarize(metabolites = list(unique(unlist(forward))), donors_n=n()) %>%
        mutate(metabolites_n = map_int(metabolites, length)) %>%
        rename(genome=1)
    } else {
      forward <- cfdb %>%
        filter(length(reverse) > 0) %>%
        filter(first %in% focal) %>%
        group_by(first) %>%
        summarize(metabolites = list(unique(unlist(reverse))), donors_n=n()) %>%
        mutate(metabolites_n = map_int(metabolites, length)) %>%
        rename(genome=1)

      reverse <- cfdb %>%
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
      forward <- cfdb %>%
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
        # Calculate the ratio between donor and pool of receptors for each metabolite
        group_by(first,reverse) %>%
        #summarise(across(where(is.double), ~sum(.x[sample == "secondX"]) / max(.x[sample == "firstX"]))) %>% #without topping
        summarise(across(where(is.double), ~ pmin(1, sum(.x[sample == "secondX"]) / max(.x[sample == "firstX"]))), .groups = "drop") %>%

        # Calculate the donor potential per genome
        group_by(first) %>%
        summarise(across(where(is.double), sum)) %>%
        rename(genome=1)

      reverse <- cfdb %>%
        # Append abundance values of first and second genomes
        left_join(abundance %>%
                    rename_with(~ paste0(., "_firstX"), -1),
                  by=join_by(first==genome)) %>%
        left_join(abundance %>%
                    rename_with(~ paste0(., "_secondX"), -1),
                  by=join_by(second==genome)) %>%
        select(-reverse,-total) %>% 			# Drop reverse and total columns
        unnest(cols = forward) %>% 			# Expand metabolites
        pivot_longer(cols = contains(c("firstX", "secondX")), 	# Pivot longer to create one column per sample
                     names_to = c(".value", "sample"),
                     names_pattern = "(.+)_(firstX|secondX)",
                     values_drop_na = TRUE) %>%
        # Calculate the ratio between donor and pool of receptors for each metabolite
        group_by(second, forward) %>%
        #summarise(across(where(is.double), ~sum(.x[sample == "firstX"]) / max(.x[sample == "secondX"]))) %>% #without topping
        summarise(across(where(is.double), ~ pmin(1, sum(.x[sample == "firstX"]) / max(.x[sample == "secondX"]))), .groups = "drop") %>%
        # Calculate the donor potential per genome
        group_by(second) %>%
        summarise(across(where(is.double), sum)) %>%
        rename(genome=1)
    } else {
      suppressWarnings({
      forward <- cfdb %>%
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
        # Calculate the ratio between donor and pool of receptors for each metabolite
        group_by(first,reverse) %>%
        #summarise(across(where(is.double), ~sum(.x[sample == "secondX"]) / max(.x[sample == "firstX"]))) %>% #without topping
        summarise(across(where(is.double), ~ pmin(1, sum(.x[sample == "secondX"]) / max(.x[sample == "firstX"]))), .groups = "drop") %>%
        # Calculate the donor potential per genome
        group_by(first) %>%
        summarise(across(where(is.double), sum)) %>%
        rename(genome=1)

      reverse <- cfdb %>%
        filter(second %in% focal) %>%
        # Append abundance values of first and second genomes
        left_join(abundance %>%
                    rename_with(~ paste0(., "_firstX"), -1),
                  by=join_by(first==genome)) %>%
        left_join(abundance %>%
                    rename_with(~ paste0(., "_secondX"), -1),
                  by=join_by(second==genome)) %>%
        select(-reverse,-total) %>% 			# Drop reverse and total columns
        unnest(cols = forward) %>% 			# Expand metabolites
        pivot_longer(cols = contains(c("firstX", "secondX")), 	# Pivot longer to create one column per sample
                     names_to = c(".value", "sample"),
                     names_pattern = "(.+)_(firstX|secondX)",
                     values_drop_na = TRUE) %>%
        # Calculate the ratio between donor and pool of receptors for each metabolite
        group_by(second, forward) %>%
       # summarise(across(where(is.double), ~sum(.x[sample == "firstX"]) / max(.x[sample == "secondX"]))) %>% #without topping
        summarise(across(where(is.double), ~ pmin(1, sum(.x[sample == "firstX"]) / max(.x[sample == "secondX"]))), .groups = "drop") %>%
        # Calculate the donor potential per genome
        group_by(second) %>%
        summarise(across(where(is.double), sum)) %>%
        rename(genome=1)
      })
    }

    receptor_potential <- bind_rows(forward,reverse) %>%
      rowwise() %>%
      group_by(genome) %>%
      summarise(across(where(is.double), sum))

  }

  #Output receptor potential table
  return(receptor_potential)
}
