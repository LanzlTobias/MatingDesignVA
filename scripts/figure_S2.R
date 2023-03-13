## Simulation for the data used in figure S2

setwd('~/LD/Paper')

## Load packages
suppressPackageStartupMessages(library(AlphaSimR))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(Rcpp))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(future.apply))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(broom))

## Set up workers for parallel running jobs
plan(multicore, workers = 80)

## Load functions
sourceCpp('scripts/cpp_functions.cpp')
source('scripts/R_functions.R')



## Set up parameter space

traits <- 10000
P_vector <-  c(2, 4)
N <- 1000

## Load data 
pop_foundation <- readRDS('input/founding_pop.RData')

## Sample number of potential QTL positions per chromosome
qtl_per_chr <- sapply(pop_foundation$genMap, length)

## Sample the QTL positions for one set out out the pool of 2500 possibilities
n_qtl = 1000
set.seed(1)
selected_qtl <- sort(sample(1:2500, n_qtl))

## Create a list with the index of the position over the chromosomes
chr_qtl_map <-
  split(1:2500, do.call(c, sapply(1:10, function(chr) {
    rep(chr, qtl_per_chr[chr])
  })))

## Calculate how many of the selected QTL positions are on each chromosome
selected_qtl_per_chr <-  sapply(chr_qtl_map, function(chr) {
  sum(selected_qtl %in% chr)
})


## Sample the QTL effects 
alpha <- lapply(1:traits, function(x) {
  set.seed(x)
  
  alpha <-  rnorm(n_qtl, mean = 0, sd = 1)
  
  return(alpha)
})
set.seed(1)

start <-  Sys.time()

message('Start: ', start)
# Preparation ------------------------------------------------------------

## Set elite and landrace together as the founding population
founder_pop <- newMapPop(
  genMap = pop_foundation$genMap,
  haplotypes = pop_foundation$haplotypes,
  inbred = TRUE
)

SP <- SimParam$new(founder_pop)

## Set crossover interference to zero (according to Haldane's mapping function)
SP$v <- 1

## Filter for the Elite population
pop <-  newPop(founder_pop, simParam = SP)

EL <-
  selectInd(
    pop,
    nInd = 115,
    candidates = 1:115,
    use = 'rand',
    simParam = SP
  )

result <- lapply(P_vector, function(P) {
  set.seed(1)
  
  ## Sample the P parental lines
  EL_sampled_lines <-  sample(1:EL@nInd, P)
  EL_parents <-  selectInd(
    pop = EL,
    nInd = P,
    use = 'rand',
    candidates = EL_sampled_lines,
    simParam = SP
  )

  ## Set the crossing plan according to disjoint crosses scheme (DC)
  p1 <- EL_parents@id[1:EL_parents@nInd %% 2 == 0]
  p2 <- EL_parents@id[1:EL_parents@nInd %% 2 != 0]
  
  crossplan <-  cbind(p1, p2)
  
  if (N <= nrow(crossplan)) {
    crossplan = crossplan[sample(1:nrow(crossplan), size = N),]
  } else {
    tmp1 <-  rep(crossplan[, 1], ceiling(N / nrow(crossplan)))
    tmp2 <-  rep(crossplan[, 2], ceiling(N / nrow(crossplan)))
    crossplan <-  cbind(tmp1, tmp2)[sample(1:length(tmp1), size = N), ]
  }  
  
  ## G1
  G1 <- makeCross(
    pop = EL_parents,
    crossplan,
    nProgeny = 1,
    simParam = SP
  )
  
  ## Generate DH-lines
  G1_DH <-  makeDH(
    pop = G1,
    nDH = 1,
    keepParents = TRUE,
    simParam = SP
  )
  
  ## Extract the matrix of the 2500 genotypic scores
  X <- pullSegSiteGeno(G1_DH, simParam = SP)
  
  ## Calculate matrix D based on X
  D <- D_calc(X)
  
  ## Calculate the conditional expectations and variances based on the combination of replication and set
  expectation_df <- conditional_var_exp(D[selected_qtl,
                                        selected_qtl],
                                      selected_qtl_per_chr)
  
  ## Calculate the parameter for each sample of allelic effects
  simulation_df <- bind_rows(future_lapply(alpha, function(a) {
    VA_decomposition(a = a,
                     mtx = D[selected_qtl,
                             selected_qtl],
                     selected_qtl_per_chr = selected_qtl_per_chr)
    
  }), .id = 'trait')
  
  
  return(list(expectation = expectation_df,
              simulation = simulation_df))
})

names(result) <-  paste0('P=', P_vector)

message(Sys.time() - start)

## Summarize the results
expectation_df <- bind_rows(map(result, 'expectation'), .id = 'P') %>% 
  mutate(P = factor(P, levels = paste0('P=', c(1, 2, 4, 16)),
                    labels = paste0('P=', c(1, 2, 4, 16))))

simulation_df <- bind_rows(map(result, 'simulation'), .id = 'P') %>% 
  mutate(P = factor(P, levels = paste0('P=', P_vector),
                    labels = paste0('P=', P_vector)))

## Save the results
saveRDS(expectation_df, file = 'output/figS2_expectation.RData')
saveRDS(simulation_df, file = 'output/figS2_simulation.RData')


