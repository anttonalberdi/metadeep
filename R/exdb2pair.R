#' Convert exchange database to pairwise matrix
#' @title Convert exchange database (exdb) to pairwise matrix
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Convert exchange database (exdb) to pairwise matrix
#' @param exdb An exchange database (exdb) generated using mdb2exdb()
#' @param exchange Whether to calculate forward (columns to rows), reverse (rows to columns) or total metabolite exchanges. Default is mode="total".
#' @param abundance An optional data frame containing relative abundance data of bacteria
#' @param verbosity Whether to print the progress of the computation when using abundance data. Default=TRUE
#' @import tidyverse
#' @examples
#' exdb2pair(allgenomes_exdb)
#' exdb2pair(allgenomes_exdb, exchange="total")
#' exdb2pair(allgenomes_exdb, exchange=c("forward","reverse"))
#' exdb2pair(allgenomes_exdb, exchange="total", abundance=genome_abundances)
#' exdb2pair(allgenomes_exdb, exchange="total", abundance=genome_abundances, verbosity=FALSE)
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

exdb2pair <- function(exdb, exchange=c("forward","reverse","total"), abundance, verbosity=TRUE) {
  # Input check
  if (inherits(exdb, "exdb")) {
    exdb <- exdb
  } else {
    stop("Input is not a valid exdb object created by mdb2exdb().")
  }

  # General function to turn the exchange tibble into matrix
  tibble2pair <- function(tibble) {
    #Identify first and last genomes
    firstgenome <- tibble %>% arrange(first) %>% head(1) %>% pull(first)
    lastgenome <- tibble %>% arrange(second) %>% tail(1) %>% pull(second)

    #Perform conversion
    tibble  %>%
        arrange(second) %>%
        pivot_wider(names_from = second, values_from = exchange, values_fill = NA) %>%
        arrange(first) %>%
        mutate(firstcol=NA, .before = 2) %>%
        rename(genomes=1) %>%
        add_row(genomes = NA) %>%
        rename(!! firstgenome := 2) %>%
        mutate(genomes = if_else(row_number() == n(), {{ lastgenome }}, genomes)) %>%
        column_to_rownames(var="genomes")
    }

  if(missing(abundance)){
    #Generate exchange matrix
    first2second <- exdb %>%
      select(first,second,forward) %>%
      rename(exchange=3) %>%
      mutate(exchange = length(exchange)) %>%
      tibble2pair()

    second2first <- exdb %>%
      select(first,second,reverse) %>%
      rename(exchange=3) %>%
      mutate(exchange = length(exchange)) %>%
      tibble2pair()

    total <- exdb %>%
      select(first,second,total) %>%
      rename(exchange=3) %>%
      mutate(exchange = length(exchange)) %>%
      tibble2pair()

    #Print matrices into the output list
    pair <- list()
    if("forward" %in% exchange){pair[["forward"]] <- first2second}
    if("reverse" %in% exchange){pair[["reverse"]] <- second2first}
    if("total" %in% exchange){pair[["total"]] <- total}
    if(length(pair) == 1){pair <- pair[[1]]} #unlist if a single exchange type

  } else {
    #Force relative abundance transformation
    abundance <- abundance %>%
        mutate(across(where(is.double), ~ ./sum(.)))

    #Identify samples
    samples <- names(abundance)[2:length(names(abundance))]

    #Per-sample loop is used, because in real datasets multi-sample calculations become too memory-exhaustive
    pair <- list()
    n=0
    for (sample in samples){
      n=n+1
      if(verbosity==TRUE){message(paste0("Processing ",sample," (",n,"/",length(samples),")"))}
      #Append abundances to exdb
      exdb_abun <- exdb %>%
        left_join(abundance %>% select(genome,{{ sample }}) %>% rename(abun_first=2),
                  by=join_by(first==genome)) %>%
        left_join(abundance %>% select(genome,{{ sample }}) %>% rename(abun_second=2),
                  by=join_by(second==genome))

      #Compute first2second exchange
      first2second <- exdb_abun %>%
        select(-reverse,-total) %>% 		# Drop non-required columns
        unnest(cols = forward) %>% 			# Expand metabolites
        pivot_longer(cols = c("abun_first", "abun_second"), 	# Pivot longer to perform calculations
                     names_to = c(".value", "genome"),
                     names_pattern = "(.+)_(first|second)",
                     values_drop_na = TRUE) %>%
        rename(metabolite=3) %>%
        group_by(first,second,metabolite) %>%
        summarise(across(where(is.double), ~ pmin(1, max(.x[genome == "first"]) / sum(.x[genome == "second"]))), .groups = "drop") %>%
        mutate_if(is.double, ~ifelse(is.nan(.), 0, .)) %>% #convert NaNs derived from n/0 to 0
        rowwise() %>%
        group_by(first,second) %>%
        summarise(exchange=sum(abun), .groups = "drop") %>%
        ungroup()

      #Compute second2first exchange
      second2first <- exdb_abun %>%
        select(-forward,-total) %>% 		# Drop non-required columns
        unnest(cols = reverse) %>% 			# Expand metabolites
        pivot_longer(cols = c("abun_first", "abun_second"), 	# Pivot longer to perform calculations
                       names_to = c(".value", "genome"),
                       names_pattern = "(.+)_(first|second)",
                       values_drop_na = TRUE) %>%
        rename(metabolite=3) %>%
        group_by(second,first,metabolite) %>%
        summarise(across(where(is.double), ~ pmin(1, max(.x[genome == "second"]) / sum(.x[genome == "first"]))), .groups = "drop") %>%
        mutate_if(is.double, ~ifelse(is.nan(.), 0, .)) %>% #convert NaNs derived from n/0 to 0
        rowwise() %>%
        group_by(first,second) %>%
        summarise(exchange=sum(abun), .groups = "drop") %>%
        ungroup()

      #Compute total exchange (first2second + second2first)
      total <- first2second %>%
          mutate(exchange=exchange+second2first$exchange)

      #Print matrices into the output list
      pair_sample <- list()
      if("forward" %in% exchange){pair_sample[["forward"]] <- tibble2pair(first2second)}
      if("reverse" %in% exchange){pair_sample[["reverse"]] <- tibble2pair(second2first)}
      if("total" %in% exchange){pair_sample[["total"]] <- tibble2pair(total)}
      if(length(pair_sample) == 1){pair_sample <- pair_sample[[1]]} #unlist if a single exchange type

      #Append output to sample list
      pair[[sample]] <- pair_sample
    }
  }

  #Output pairwise exchange matrix or matrices
  return(pair)
}
