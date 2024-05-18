# MetaDEEP

MetaDEEP, the **Meta**bolic **D**ependence and **E**xchange **E**valuation **P**ackage, is an R package to analyse metabolic dependence and exchange metrics of microbiomes based on genome-scale metabolic networks (GSMN) and genome abundance data. MetaDEEP identifies metabolites that can be exchanged among bacteria in a community based on GSMNs, and calculates effective metabolite exchange values between genomes, accounting for the number of donors and receptors for each metabolite and their respective relative abundances. The summarisation of pairwise exchanges yields genome-, metabolite- and sample-level metabolite exchange metrics that allow characterising microbiome structures in terms of effective metabolite exchange capacities.  

## Installation

MetaDEEP can be installed from this Github repository using devtools.

```r
install.packages("devtools")
library(devtools)
install_github("anttonalberdi/metadeep")
library(metadeep)
```

### Dependencies

MetaDEEP only has three strict dependencies:

- [tidyverse](https://www.tidyverse.org)
- [SBMLR](https://www.bioconductor.org/packages/release/bioc/html/SBMLR.html)
- [igraph](https://igraph.org)

## Worflow

MetaDEEP can simply be used by running the master function **metadeep()**, which pipes all the intermediate steps and yields summary metabolite exchange metrics for microbiomes.

```r
sbml_files <- list.files(path = "data", pattern = "\\.sbml$", full.names = TRUE)
metadeep(sbml_files, abundance=genome_abundances, summarise="samples")
```

The function**metadeep()** runs the following workflow:

1. **sbml2rdb()** to load SBML files into a reaction database.
2. **rdb2igraph()** to generate genome-specific igraph network objects.
3. **rdb2gedb()** to generate a genome-metabolite database.
4. **gedb2medb()** to generate a metabolite-genome database.
5. **medb2exdb()** to generate a metabolite exchange database.
6. **exdb2summary()** to summarise metabolite exchange metrics in diverse ways.
7. **exdb2pair()** to display exchange metrics into pairwise matrices.

If working with large data sets, or exploring different metabolite exchange characteristics, it is adviced to run the workflow step-by-step, to better customise the analyses.

## Step-by-step usage

Step-by-step usage and interpretation of the MetaDEEP package and its contents.

### Load a single SBML (sbml2rdb)
Load and convert a SBML file into MetaDEEP reaction database (rdb). The resulting object is a tibble containing character lists of reactants and products of each reaction. 

```r
# Read SBML directly into MetaDEEP reaction database (rdb)
genome1_rdb <- sbml2rdb("data/genome1.sbml")

# Read SBML in two steps into MetaDEEP reaction database (rdb)
genome1_rdb <- readSBML("data/genome1.sbml") %>% sbml2rdb()
```

| reaction                          | reactants    | products     |
|-----------------------------------|--------------|--------------|
| R_RXN__45__18707                  | <chr [2]>    | <chr [2]>    |
| R_RXN__45__14267                  | <chr [2]>    | <chr [2]>    |
| R_3__46__1__46__26__46__4__45__RXN| <chr [2]>    | <chr [2]>    |
| R_CARDIOLIPSYN__45__RXN           | <chr [1]>    | <chr [2]>    |

#### Explore a specific reaction

It is possible to visualise the metabolites involved in each reaction.

```r
genome1_rdb %>% 
    filter(reaction == "R_RXN__45__18707") %>% 
    unnest(cols = c(reactants, products))
```

| reaction         | reactants                                           | products                                      |
|------------------|-----------------------------------------------------|-----------------------------------------------|
| R_RXN__45__18707 | M_L__45__Cysteine__45__Desulfurase__45__persulfide_c | M_Cysteine__45__Desulfurase__45__L__45__cystei… |
| R_RXN__45__18707 | M_TusA__45__L__45__cysteine_c                       | M_TusA__45__Persulfides_c                    |

### Load multiple SBML (sbml2rdb)
Convert multiple SBML files into a list of MetaDEEP reaction databases (rdb). The resulting object is a list of tibbles containing character lists of reactants and products of each reaction. This step might take a few minutes (if working with hundreds of genomes) or even a few hours (if working with thousands of genomes) to accomplish. The resulting rdbs object can be considerably large (several GBs if working with thousands of genomes), which might require large memory allocation.

```r
sbml_files <- list.files(path = "data", pattern = "\\.sbml$", full.names = TRUE)
allgenomes_rdb <- sbml2rdb(sbml_files)
```

$genome1
| reaction                          | reactants    | products     |
|-----------------------------------|--------------|--------------|
| R_RXN__45__18707                  | <chr [2]>    | <chr [2]>    |
| R_RXN__45__14267                  | <chr [2]>    | <chr [2]>    |
| R_3__46__1__46__26__46__4__45__RXN| <chr [2]>    | <chr [2]>    |
| R_CARDIOLIPSYN__45__RXN           | <chr [1]>    | <chr [2]>    |

$genome2
| reaction                                 | reactants    | products     |
|------------------------------------------|--------------|--------------|
| R_THREOSPON__45__RXN                     | <chr [2]>    | <chr [2]>    |
| R_UDPNACETYLGLUCOSAMENOLPYRTRANS__45__RXN| <chr [2]>    | <chr [2]>    |
| R_RXN__45__22610                         | <chr [2]>    | <chr [2]>    |
| R_RXN__45__15920                         | <chr [2]>    | <chr [3]>    |

### Convert a reaction database into an igraph network (rdb2igraph)

```r
genome1_igraph <- rdb2igraph(genome1_rdb)
allgenomes_igraph <- rdb2igraph(allgenomes_rdb)
```

#### Visualise the network

```r
genome1_igraph %>%
  as.undirected() %>%
  plot(., 
        layout = layout_with_fr(genome1_igraph), 
        #vertex attributes
        vertex.color=V(genome1_igraph)$color,
        vertex.frame.color=NA,
        vertex.label=NA,
        vertex.size=1.5,
        #edge attributes
        edge.label = NA, 
        edge.curved = 0.1,
        edge.color="#cccccc70",
        edge.label.color="#8E8E1E"
        )
```

![Metabolic networks](figures/metabolite_networks.png)

### Classify metabolite types into a genome database (rdb2gedb)
Classify metabolites in a reaction database (rdb) into source, transit and sink metabolites stored in a genome database (gedb).

- **Source metabolites:** those that the bacterium is able to use but not to produce.
- **Transit metabolites:** those that the bacterium is able to use and produce.
- **Sink metabolites:** those that the bacterium is able to produce but not to use.

#### In a single genome

```r
genome1_gedb <- rdb2gedb(genome1_rdb)
```

| genome  | sources    | transits   | sinks      | reactions | metabolites |
|---------|------------|------------|------------|-----------|-------------|
| genome1 | <chr [174]>| <chr [95]> | <chr [179]>| 267       | 448         |

#### In multiple genomes

```r
allgenomes_gedb <- rdb2gedb(allgenomes_rdb)
```

| genome  | sources    | transits   | sinks      | reactions | metabolites |
|---------|------------|------------|------------|-----------|-------------|
| genome1 | <chr [174]>| <chr [95]> | <chr [179]>| 267       | 448         |
| genome2 | <chr [442]>| <chr [407]>| <chr [441]>| 1011      | 1290        |
| genome3 | <chr [233]>| <chr [159]>| <chr [245]>| 395       | 637         |
| genome4 | <chr [262]>| <chr [196]>| <chr [282]>| 508       | 740         |

### Convert a genome database into a metabolite database (gedb2medb)
Transform the genome database (gedb) into a metabolite database that displays the role (sink, transit or source) each metabolite can play across genomes. The last column indicates whether that metabolite can be exchanged and the mode in which it can be exchanged (strict or loose).

- **Strict mode**: only sink (donor) and source (receptor) metabolites are considered exchangeable.
- **Loose mode**: sink (donor), transit (donor) and source (receptor) metabolites are considered exchangeable.

```r
allgenomes_medb <- gedb2medb(allgenomes_gedb)
```

| metabolites                                                       | sinks       | transits | sources   | exchangeable |
|-------------------------------------------------------------------|-------------|----------|-----------|-------------|
| M_10__45__FORMYL__45__DIHYDROFOLATE__45__GLU__45__N_c             | `<chr>`     | `<NULL>` | `<NULL>`  | no          |
| M_11Z__45__3__45__oxo__45__icos__45__11__45__enoyl__45__ACPs_c    | `<NULL>`    | `<chr>`  | `<NULL>`  | no          |
| M_11Z__45__icos__45__11__45__enoyl__45__ACPs_c                    | `<chr>`     | `<NULL>` | `<NULL>`  | no          |
| M_ACETYL__45__D__45__GLUCOSAMINYLDIPHOSPHO__45__UNDECAPRE_c       | `<chr [2]>` | `<NULL>` | `<chr [1]>` | strict      |
| M_AGMATHINE_c                                                      | `<NULL>`    | `<chr [1]>`| `<chr [1]>` | loose       |
| M_13__45__HYDROXY__45__MAGNESIUM__45__PROTOPORP_c                 | `<NULL>`    | `<chr>`  | `<NULL>`  | no          |
| M_16S__45__rRNA__45__2__45__O__45__methylcytidine1402_c           | `<chr>`     | `<NULL>` | `<NULL>`  | no          |
| M_CL__45___c                                                       | `<chr [1]>` | `<chr [1]>`| `<chr [1]>` | strict      |
| M_16S__45__rRNA__45__5__45__O__45__methylcytosine967_c            | `<chr>`     | `<NULL>` | `<NULL>`  | no          |
| M_CDPDIACYLGLYCEROL_c                                              | `<NULL>`    | `<chr [3]>`| `<chr [1]>` | loose       |

### Calculate metabolite exchange potential (medb2exdb)

Calculate the exchangeability of each metabolite across genomes. The resulting object is a list of tibbles each indicating pairwise exchange capacities for a specific metabolite across genomes. 

$M_3__45__5__45__ADP_c
| donor   | receptor | exchange |
|---------|----------|---------|
| genome1 | genome3  | 0.333   |
| genome2 | genome3  | 0.333   |
| genome4 | genome3  | 0.333   |

$M_3__45__OXOADIPATE__45__ENOL__45__LACTONE_c
| donor   | receptor | exchange |
|---------|----------|---------|
| genome2 | genome4  | 1       |

$M_C6_c
| donor   | receptor | exchange |
|---------|----------|---------|
| genome2 | genome4  | 0.5     |
| genome3 | genome4  | 0.5     |

#### Relative abundance-weighed exchanges

If abundance information is provided, exchangeability is calculated considering the relative abundances of the genomes in each sample.

```r
genome_abundances <- data.frame(genome=c("genome1","genome2","genome3","genome4"),
        sample1=c(0.25,0.25,0.25,0.25),
        sample2=c(0.40,0.40,0.10,0.10),
        sample3=c(0.85,0.05,0.05,0.05))

allgenomes_exdb <- medb2exdb(allgenomes_medb, abundance=genome_abundances)
```

$M_3__45__5__45__ADP_c
| donor   | receptor | sample1 | sample2 | sample3 |
|---------|----------|---------|---------|---------|
| genome1 | genome3  | 0.333   | 0.148   | 0.141   |
| genome2 | genome3  | 0.333   | 0.148   | 0.00831 |
| genome4 | genome3  | 0.333   | 0.0370  | 0.00831 |

$M_3__45__OXOADIPATE__45__ENOL__45__LACTONE_c
| donor   | receptor | sample1 | sample2 | sample3 |
|---------|----------|---------|---------|---------|
| genome2 | genome4  | 1       | 0.25    | 1       |

$M_C6_c
| donor   | receptor | sample1 | sample2 | sample3 |
|---------|----------|---------|---------|---------|
| genome2 | genome4  | 0.5     | 0.32    | 0.5     |
| genome3 | genome4  | 0.5     | 0.08    | 0.5     |

### Summarise exchanges (exdb2summary)

The exchange database (exdb) can be summarised into multiple metrics using exdb2summary. If exchange database has been generated using abundance data, summary statistics are generated for all samples. The default option outputs all three summary types:

- **Genomes**: metabolite exchange metrics summarised by genome.
- **Metabolites**: metabolite exchange metrics summarised by metabolite.
- **Samples**: overall metabolite exchange metrics summarised by sample.

```r
allgenomes_summary <- exdb2summary(allgenomes_exdb)
```

$genomes
| genome  | sample1 | sample2 | sample3 |
|---------|---------|---------|---------|
| genome1 | 4.25    | 3.76    | 1.03    |
| genome2 | 8.29    | 5.33    | 5.64    |
| genome3 | 5.75    | 2.21    | 3.59    |
| genome4 | 8.71    | 3.45    | 5.88    |

$metabolites
| metabolite                                                  | sample1 | sample2 | sample3 |
|-------------------------------------------------------------|---------|---------|---------|
| M_3__45__5__45__ADP_c                                       | 1       | 0.333   | 0.158   |
| M_3__45__OXOADIPATE__45__ENOL__45__LACTONE_c                | 1       | 0.25    | 1       |
| M_4__45__AMINO__45__4__45__DEOXYCHORISMATE_c                | 1       | 1       | 0.0588  |
| M_ACETYL__45__D__45__GLUCOSAMINYLDIPHOSPHO__45__UNDECAPRE_c | 1       | 1       | 1       |
| M_ALA__45__tRNAs_c                                          | 2       | 0.5     | 0.222   |
[...]

$samples
| sample  | exchange |
|---------|----------|
| sample1 | 27       |
| sample2 | 14.8     |
| sample3 | 16.1     |

#### Genome summary

The default genome summary yields average exchange values. However, genomes can be summarised by both donor and receptor capacities. This can be specified through the ***exchange*** parameter.

```r
allgenomes_summary <- exdb2summary(allgenomes_exdb, summarise="genomes", exchange=c("donor","receptor","average"))
```

$donor
| donor   | sample1 | sample2 | sample3 |
|---------|---------|---------|---------|
| genome1 | 2.5     | 1.57    | 0.588   |
| genome2 | 9.33    | 7.24    | 6.94    |
| genome3 | 6       | 1.90    | 3.42    |
| genome4 | 9.17    | 4.04    | 5.20    |

$receptor
| receptor | sample1 | sample2 | sample3 |
|----------|---------|---------|---------|
| genome1  | 6       | 5.95    | 1.48    |
| genome2  | 7.25    | 3.41    | 4.33    |
| genome3  | 5.5     | 2.52    | 3.77    |
| genome4  | 8.25    | 2.87    | 6.57    |

$average
| genome  | sample1 | sample2 | sample3 |
|---------|---------|---------|---------|
| genome1 | 4.25    | 3.76    | 1.03    |
| genome2 | 8.29    | 5.33    | 5.64    |
| genome3 | 5.75    | 2.21    | 3.59    |
| genome4 | 8.71    | 3.45    | 5.88    |

#### Focal genomes, metabolites and samples

One can also limit the summaries to certain focal genomes, metabolites and/or samples.

```r
allgenomes_summary <- exdb2summary(allgenomes_exdb, genomes=c("genome1","genome2"), metabolites="M_PROTON_e", samples=c("sample2","sample3"))
```

$genomes
| genome  | sample2 | sample3 |
|---------|---------|---------|
| genome1 | 0.64    | 0.0588  |
| genome2 | 0.64    | 0.0588  |

$metabolites
| metabolite | sample2 | sample3 |
|------------|---------|---------|
| M_PROTON_e | 1.28    | 0.118   |

$samples
| sample  | exchange |
|---------|----------|
| sample2 | 1.28     |
| sample3 | 0.118    |

#### Diversity-weighed exchange

```r
library(hilldiv2)
diversity <- genome_abundances %>%
  column_to_rownames(var="genome") %>%
  hilldiv(.,q=1) %>% t() %>%
  as.data.frame() %>%
  rownames_to_column(var="sample")
  
exdb2summary(allgenomes_exdb, summarise="samples") %>%
  left_join(diversity,by="sample") %>%
  mutate(exchange_w=exchange/q1)
```

### Pairwise exchange matrix (exdb2pair)

```r
allgenomes_pair <- exdb2pair(allgenomes_exdb)
```
