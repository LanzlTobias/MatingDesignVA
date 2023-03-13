## Calculates the off-diagonal values of the matrix D in the ancestral populations.

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
# ## Sample the QTL positions for each trait out out the pool of 2500 possibilities
# n_qtl = 1000
# set.seed(1)
# selected_qtl <- sort(sample(1:2500, n_qtl))
# 
## Calculate how many of the selected QTL positions are on each chromosome
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
  
  D_poly <- D[diag(D) > 0,
              diag(D) > 0]
  
  diagonal <- diag(D_poly)
  
  off_diag <- D_poly[upper.tri(D_poly)]
  
  qtl_per_chr_poly <-  sapply(chr_qtl_map, function(chr) {
    sum((1:ncol(D))[diag(D) > 0] %in% chr)
  })
  
  which_chr <-
    do.call(c, lapply(1:length(qtl_per_chr_poly), function(chr) {
      return(rep(chr, qtl_per_chr_poly[chr]))
    }))
  
  within_chr <-
    do.call(c, lapply(1:length(qtl_per_chr_poly), function(chr) {
      D_tmp <- D_poly[which_chr == chr,
                      which_chr == chr]
      
      return(D_tmp[upper.tri(D_tmp)])
    }))
  
  between_chr <-
    do.call(c, lapply(1:length(qtl_per_chr_poly), function(chr) {
      D_tmp <- D_poly[which_chr == chr,
                      which_chr != chr]
      
      return(D_tmp[upper.tri(D_tmp)])
    }))
  
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



