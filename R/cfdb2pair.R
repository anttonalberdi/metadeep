#' Convert cross-feeding database to pairwise matrix
#' @title Convert cross-feeding database (cfdb) to pairwise matrix
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Convert cross-feeding database (cfdb) to pairwise matrix
#' @param cfdb A cross-feeding database (cfdb) generated using mdb2cfdb()
#' @param mode Whether to calculate forward (columns to rows), reverse (rows to columns) or total metabolite exchanges. Default is mode="total".
#' @param abundance An optional data frame containing relative abundance data of bacteria
#' @import tidyverse
#' @examples
#' cfdb2pair(allgenomes_cfdb)
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

cfdb2pair <- function(cfdb, mode="total", abundance) {
  # Input check
  if (inherits(cfdb, "cfdb")) {
    cfdb <- cfdb
  } else {
    stop("Input is not a valid cfdb object created by mdb2cfdb().")
  }

  #Identify first and last genomes
  firstgenome <- cfdb %>% arrange(first) %>% head(1) %>% pull(first)
  lastgenome <- cfdb %>% arrange(second) %>% tail(1) %>% pull(second)

  if(missing(abundance)){
    #Generate exchange matrix
    pair <- cfdb %>%
      select(first,second,{{ mode }}) %>%
      rename(exchange=3) %>%
      mutate(exchange = length(exchange)) %>%
      arrange(first,second) %>%
      pivot_wider(names_from = second, values_from = exchange, values_fill = NA) %>%
      mutate(firstcol=NA, .before = 2) %>%
      rename(genomes=1) %>%
      add_row(genomes = NA) %>%
      rename(!! firstgenome := 2) %>%
      mutate(genomes = if_else(row_number() == n(), {{ lastgenome }}, genomes))
  } else {
    #Force relative abundance transformation
    abundance <- abundance %>%
        mutate(across(where(is.double), ~ ./sum(.)))

    #Identify samples
    samples <- names(abundance)[2:length(names(abundance))]

    #Generate forward list
    forward <- cfdb %>%
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
      rename(metabolite=3) %>%
      group_by(first,second,metabolite) %>%
      summarise(across(where(is.double), ~ pmin(1, max(.x[sample == "firstX"]) / sum(.x[sample == "secondX"]))), .groups = "drop") %>%
      mutate_if(is.double, ~ifelse(is.nan(.), 0, .)) %>% #convert NaNs derived from n/0 to 0
      rowwise() %>%
      group_by(first,second) %>%
      summarise(across(where(is.double), sum)) %>%
      ungroup()

    create_pair_forward <- function(numeric_col) {
      forward %>%
        select(first, second, {{ numeric_col }}) %>%
        rename(sample = {{ numeric_col }}) %>%
        arrange(second) %>%
        pivot_wider(names_from = second, values_from = sample, values_fill = NA) %>%
        arrange(first) %>%
        mutate(firstcol=NA, .before = 2) %>%
        rename(genomes=1) %>%
        add_row(genomes = NA) %>%
        rename(!! firstgenome := 2) %>%
        mutate(genomes = if_else(row_number() == n(), {{ lastgenome }}, genomes)) %>%
        column_to_rownames(var="genomes")
    }

    forward_list <- map(samples, create_pair_forward)
    names(forward_list) <- samples

    #Generate reverse list
    reverse <- cfdb %>%
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
      rename(metabolite=3) %>%
      group_by(second,first,metabolite) %>%
      summarise(across(where(is.double), ~ pmin(1, max(.x[sample == "secondX"]) / sum(.x[sample == "firstX"]))), .groups = "drop") %>%
      mutate_if(is.double, ~ifelse(is.nan(.), 0, .)) %>% #convert NaNs derived from n/0 to 0
      rowwise() %>%
      group_by(first,second) %>%
      summarise(across(where(is.double), sum)) %>%
      ungroup()


    create_pair_reverse <- function(numeric_col) {
      reverse %>%
        select(first, second, {{ numeric_col }}) %>%
        rename(sample = {{ numeric_col }}) %>%
        arrange(second) %>%
        pivot_wider(names_from = second, values_from = sample, values_fill = NA) %>%
        arrange(first) %>%
        mutate(firstcol=NA, .before = 2) %>%
        rename(genomes=1) %>%
        add_row(genomes = NA) %>%
        rename(!! firstgenome := 2) %>%
        mutate(genomes = if_else(row_number() == n(), {{ lastgenome }}, genomes)) %>%
        column_to_rownames(var="genomes")
    }

    reverse_list <- map(samples, create_pair_reverse)
    names(reverse_list) <- samples

    #Generate total list
    total_list <- Map(`+`, forward_list, reverse_list)

    #Select result
    if(mode=="total"){
      pair <- total_list
    }else if(mode=="forward"){
      pair <- forward_list
    }else{
      pair <- reverse_list
    }
  }

  #Output cross-feeding matrix
  return(pair)
}
