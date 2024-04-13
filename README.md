# MetaDEEP

etaDEEP is an R package to analyse metabolic dependence and exchange metrics from genome-scale metabolic networks (GSMN).

## Installation

MetaDEEP can be installed from this Github repository using devtools.

```r
install.packages("devtools")
library(devtools)
install_github("anttonalberdi/metadeep")
library(metadeep)
```

## Usage
Basic usage of MetaDEEP package

### Load a single SBML
Load and convert a SBML file into MetaDEEP reaction database (rdb).

```r
# Read SBML directly into MetaDEEP reaction database (rdb)
genome1_rdb <- sbml2rdb("data/genome1.sbml")

# Read SBML in two steps into MetaDEEP reaction database (rdb)
genome1_rdb <- readSBML("data/genome1.sbml") %>% sbml2rdb()
```

### Load multiple SBML
Convert multiple SBML files into a list of MetaDEEP reaction databases (rdb).

```r
sbml_files <- list.files(path = "data", pattern = "\\.sbml$", full.names = TRUE)
sbmls2rdb(sbml_files)
```

### Classify metabolite types
Classify metabolites in a reaction database (rdb) into source, transit and sink metabolites stored in a metabolite database (mdb).

```r
genome1_mdb <- rdb2mdb(genome1_rdb)
```

```r
allgenomes_mdb <- map_df(allgenomes_rdb, rdb2mdb)
```
