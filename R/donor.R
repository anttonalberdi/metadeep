#' Calculation of donor potential
#' @title Calculation of donor potential
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Calculation of donor potential of each genome
#' @param cfdb A cross-feeding database (cfdb) generated using mdb2cfdb()
#' @param abundance A data frame containing relative abundance data of bacteria
#' @import tidyverse
#' @examples
#' donor(allgenomes_cfdb)
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

donor <- function(cfdb, abundance) {
  # Input check
  if (inherits(cfdb, "cfdb")) {
    cfdb <- cfdb
  } else {
    stop("Input is not a valid cfdb object created by mdb2cfdb().")
  }

  # If no abundance data is provided
  if(missing(abundance)){
    forward <- cfdb %>%
      filter(length(forward) > 0) %>%
      group_by(first) %>%
      summarize(metabolites = list(unique(unlist(forward))), receptors_n=n()) %>%
      mutate(metabolites_n = map_int(metabolites, length)) %>%
      rename(genome=1)

    reverse <- cfdb %>%
      filter(length(reverse) > 0) %>%
      group_by(second) %>%
      summarize(metabolites = list(unique(unlist(reverse))), receptors_n=n()) %>%
      mutate(metabolites_n = map_int(metabolites, length)) %>%
      rename(genome=1)

    donor_potential <- bind_rows(forward,reverse) %>%
      rowwise() %>%
      group_by(genome) %>%
      summarize(metabolites = list(unique(unlist(metabolites))), receptors_n=sum(receptors_n)) %>%
      mutate(metabolites_n = map_int(metabolites, length))
  }

  # If abundance data is provided
  if(!missing(abundance)){
    forward <- cfdb %>%
      # Append abundance values of first and second genomes
      left_join(genome_abundances %>%
                  rename_with(~ paste0(., "_first"), -1),
                by=join_by(first==genome)) %>%
      left_join(genome_abundances %>%
                  rename_with(~ paste0(., "_second"), -1),
                by=join_by(second==genome)) %>%
      select(-reverse,-total) %>% 			# Drop reverse and total columns
      unnest(cols = forward) %>% 			# Expand metabolites
      pivot_longer(cols = starts_with("sample"), 	# Pivot longer to create one column per sample
                   names_to = c(".value", "sample"),
                   names_sep = "_") %>%
      # Calculate the ratio between donor and pool of receptors for each metabolite
      group_by(first,forward) %>%
      summarise(across(where(is.double), ~max(.x[sample == "first"]) / sum(.x[sample == "second"]))) %>%
      # Calculate the donor potential per genome
      group_by(first) %>%
      summarise(across(where(is.double), sum)) %>%
      rename(genome=1)

    reverse <- cfdb %>%
      # Append abundance values of first and second genomes
      left_join(genome_abundances %>%
                  rename_with(~ paste0(., "_first"), -1),
                by=join_by(first==genome)) %>%
      left_join(genome_abundances %>%
                  rename_with(~ paste0(., "_second"), -1),
                by=join_by(second==genome)) %>%
      select(-forward,-total) %>% 			# Drop reverse and total columns
      unnest(cols = reverse) %>% 			# Expand metabolites
      pivot_longer(cols = starts_with("sample"), 	# Pivot longer to create one column per sample
                   names_to = c(".value", "sample"),
                   names_sep = "_") %>%
      # Calculate the ratio between donor and pool of receptors for each metabolite
      group_by(second, reverse) %>%
      summarise(across(where(is.double), ~max(.x[sample == "second"]) / sum(.x[sample == "first"]))) %>%
      # Calculate the donor potential per genome
      group_by(second) %>%
      summarise(across(where(is.double), sum)) %>%
      rename(genome=1)


    donor_potential <- bind_rows(forward,reverse) %>%
      rowwise() %>%
      group_by(genome) %>%
      summarise(across(where(is.double), sum))
  }

  #Output cross-feeding matrix
  return(donor_potential)
}
