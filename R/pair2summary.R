#' Summarise pairwise matrix or list of pairwise matrices
#' @title Summarise pairwise matrix or list of pairwise matrices
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Summarise pairwise matrix or list of pairwise matrices
#' @param pair A cross-feeding database (cfdb) generated using mdb2cfdb()
#' @import tidyverse
#' @examples
#' pair2summary(allgenomes_exchange_total)
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

pair2summary <- function(pair) {

  if (!inherits(pair, "list")) {

    #Remove genome names column
    mat <- pair %>% select(-genomes)

    #If pair object is a single paiwrise matrix
      sum_val <- sum(as.matrix(mat), na.rm = TRUE)
      median_val <- median(as.matrix(mat), na.rm = TRUE)
      mean_val <- mean(as.matrix(mat), na.rm = TRUE)
      sd_val <- sd(as.matrix(mat), na.rm = TRUE)

      summary <- tibble(
        sum = sum_val,
        median= median_val,
        mean = mean_val,
        sd = sd_val)

    } else{
    #If pair object is a list of paiwrise matrix
      summary <- map_dfr(names(pair), ~ {
        # Extract matrix by name
        mat <- pair[[.x]]

        # Calculate mean and sd for each matrix
        sum_val <- sum(as.matrix(mat), na.rm = TRUE)
        median_val <- median(as.matrix(mat), na.rm = TRUE)
        mean_val <- mean(as.matrix(mat), na.rm = TRUE)
        sd_val <- sd(as.matrix(mat), na.rm = TRUE)

        # Return a data frame with list name, mean, and sd
        tibble(
          sample = .x,
          sum = sum_val,
          median= median_val,
          mean = mean_val,
          sd = sd_val)
      })
    }

  return(summary)
  }
