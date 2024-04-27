#' Quantify metabolite exchanges in microbiomes
#' @title Quantify metabolite exchanges in microbiomes
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords SBML tibble reaction reactant product
#' @description Quantify metabolite exchanges in microbiome. Note this is a wrapper of many functions in metaDEEP.
#'  Run the pipeline step by step to achieve more customisation capacity
#' @param sbml An object of class SBMLR, an sbml file or, most commonly, a vector of smbl files.
#' @param summarise Whether to summarise genomes, metabolites or samples. Default is all three.
#' @param exchange Whether to calculate genome summaries as average, donor, receptor metabolite exchange capacities. Default is average.
#' @param metabolites An optional vector of metabolites for which to run the calculations. Default is all metabolites.
#' @param genomes An optional vector of genomes for which to run the calculations. Default is all genomes
#' @param samples An optional vector of samples for which to run the calculations. Default is all samples
#' @import tidyverse
#' @examples
#' metadeep(sbml_files)
#' @references
#' Keating, S.M. et al. (2020). SBML Level 3: an extensible format for the exchange and reuse of biological models. Molecular Systems Biology 16: e9110
#' @export

metadeep <- function(sbml, mode="strict", abundance, summarise=c("genomes","metabolites","samples"), exchange="average", metabolites, genomes, samples, verbosity=TRUE){

  # Load SBML files
  message("Loading SBML files. This might take some time...")
  rdb <- sbml2rdb(sbml)

  # Generate genome database
  message("Generating genome database...")
  gedb <- rdb2gedb(rdb)

  # Generate metabolite database
  message("Generating metabolite database...")
  medb <- gedb2medb(gedb)

  # Generate exchange database
  if(missing(genomes)){genomes=c(unlist(medb$sinks),unlist(medb$transits),unlist(medb$sources)) %>%
    unique() %>% sort()}

  if(missing(metabolites)){metabolites=c(medb$metabolite)}

  if(missing(abundance)){
    abundance <- tibble(genome=genomes,exchange=rep(1/length(genomes),length(genomes)))
  }else{
    #Force relative abundance transformation
    abundance <- abundance %>%
      mutate(across(where(is.double), ~ ./sum(.)))
  }

  exdb <- medb2exdb(medb, mode=mode, abundance, metabolites=metabolites, genomes=genomes, samples=samples, verbosity=verbosity)
  message("Generating metabolite exchange database. This might take a while...")

  # Generate summary
  message("Generating exchange summary...")
  summary <- exdb2summary(exdb, summarise=summarise, exchange=exchange, metabolites=metabolites, genomes=genomes, samples=samples)

  return(summary)

}


