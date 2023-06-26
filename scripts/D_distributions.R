setwd('~/LD/Paper')
suppressPackageStartupMessages(library(AlphaSimR))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(Rcpp))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(broom))


sourceCpp('scripts/cpp_functions.cpp')
source('scripts/R_functions.R')

## Set up parameter space

## Load data and functions
pop_foundation <- readRDS('input/founding_pop.RData')

qtl_per_chr <- sapply(pop_foundation$genMap, length)
# 


## Calculate how many of the QTL positions are on each chromosome
chr_qtl_map <-
  split(1:2500, do.call(c, sapply(1:10, function(chr) {
    rep(chr, qtl_per_chr[chr])
  })))

## Set elite and landrace together as the founding population
founder_pop <- newMapPop(
  genMap = pop_foundation$genMap,
  haplotypes = pop_foundation$haplotypes,
  inbred = TRUE
)

SP <- SimParam$new(founder_pop)

## Set crossover interference to zero
SP$v <- 1

## Divide the founding population into elite and landrace again
pop <-  newPop(founder_pop, simParam = SP)

EL <-
  selectInd(
    pop,
    nInd = 115,
    candidates = 1:115,
    use = 'rand',
    simParam = SP
  )

X_EL <- pullSegSiteGeno(EL, simParam = SP)

LR <-
  selectInd(
    pop,
    nInd = 115,
    candidates = 116:230,
    use = 'rand',
    simParam = SP
  )

X_LR <- pullSegSiteGeno(LR, simParam = SP)


result <- lapply(list(X_EL, X_LR), function(X) {
  D <- D_calc(X)
  

  diagonal <- diag(D)
  
  off_diag <- rep(D[upper.tri(D)], 2)
  

  tmp_within_chr <- lapply(1:length(chr_qtl_map), function(chr){
  D_chr <- D[chr_qtl_map[[chr]],
    chr_qtl_map[[chr]]]
  
  return(rep(D_chr[upper.tri(D_chr)], 2))
  })
  
  within_chr <- do.call(c, tmp_within_chr)
  
  
  chr <- 1
  
  tmp_between_chr <- lapply(1:length(chr_qtl_map), function(chr){
    return(D[chr_qtl_map[[chr]], -chr_qtl_map[[chr]]])
  })
  
  between_chr <- do.call(c, tmp_between_chr)
  
  
  
  result_poly <- list(
    diag = as_tibble(diagonal),
    off_diag = as_tibble(off_diag),
    within_chr = as_tibble(within_chr),
    between_chr = as_tibble(between_chr)
  )
})

names(result) <- c('Elite', 'Landrace')

result_df <- bind_rows(map(result, bind_rows, .id = 'elements'), .id = 'population')

saveRDS(result_df, file = 'output/D_ancestral.RData')


