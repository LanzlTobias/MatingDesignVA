## Script for running the main simulation. 
## Needs the source files in R and C++ as well as the provided input files. 
## Can be parallelized by changing the amount of workers.


## Load packages
setwd('~/LD/Paper')
suppressPackageStartupMessages(library(AlphaSimR))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(future.apply))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(Rcpp))

## Set up parallel computing
plan(multicore, workers = 80)

## Set up parameter space
replications <- 500
sets <- 50
P_vector <-  c(2, 4, 8, 16)
N_vector <-  c(50, 250, 1000)
L_vector <- c(50, 250, 1000)

## Load data and functions
pop_foundation <- readRDS('input/founding_pop.RData')
sourceCpp('scripts/cpp_functions.cpp')
source('scripts/R_functions.R')

qtl_per_chr <- sapply(pop_foundation$genMap, length)

## Start time keeping
start <-  Sys.time()

message('Start: ', start)
# Preparation ------------------------------------------------------------

## Set elite and landrace together as the ancestral population
ancestral_pop <- newMapPop(
  genMap = pop_foundation$genMap,
  haplotypes = pop_foundation$haplotypes,
  inbred = TRUE
)

SP <- SimParam$new(ancestral_pop)

## Set crossover interference to zero
SP$v <- 1


## Divide the founding population into elite and landrace again
pop <-  newPop(ancestral_pop, simParam = SP)

EL <-
  selectInd(
    pop,
    nInd = 115,
    candidates = 1:115,
    use = 'rand',
    simParam = SP
  )
LR <-
  selectInd(
    pop,
    nInd = 115,
    candidates = 116:230,
    use = 'rand',
    simParam = SP
  )



# Start Loop --------------------------------------------------------------

L_result <- lapply(L_vector, function(L) {
  
  ## select samples of L QTL positions
  selected_qtl <-  lapply(1:sets, function(x, L) {
    set.seed(x)
    return(sort(sample(1:2500, size = L)))
  }, L = L)
  
  P_result <-  lapply(P_vector, function(P, selected_qtl) {
    N_result <-  lapply(N_vector, function(N, P, selected_qtl) {
      set.seed(P)
      
      message('L', L, '_P', P, '_N', N)
      print(Sys.time() - start)
      rep_result  <-
        future_lapply(1:replications, function(rep, N, P, selected_qtl) {
          ## Sample the P parental lines from both ancestral populations
          EL_sampled_lines <-  sample(1:EL@nInd, P)
          LR_sampled_lines <-  sample(1:LR@nInd, P)
          
          EL_parents <-  selectInd(
            pop = EL,
            nInd = P,
            use = 'rand',
            candidates = EL_sampled_lines,
            simParam = SP
          )
          
          
          LR_parents <-  selectInd(
            pop = LR,
            nInd = P,
            use = 'rand',
            candidates = LR_sampled_lines,
            simParam = SP
          )
          
          ## Simulation and analysis of the populations after the three mating designs
          FC_result <-  list(
            Elite = cross_simulation_FC(
              parents = EL_parents,
              N = N,
              selected_qtl = selected_qtl,
              SP = SP
            ),
            Landrace = cross_simulation_FC(
              parents = LR_parents,
              N = N,
              selected_qtl = selected_qtl,
              SP = SP
            )
          )
          
          FC_result <-
            bind_rows(map(FC_result, 'parameter_df'), .id = 'ancestral')
          
          HC_result <-  list(
            Elite = cross_simulation_HC(
              parents = EL_parents,
              N = N,
              selected_qtl = selected_qtl,
              SP = SP
            ),
            Landrace = cross_simulation_HC(
              parents = LR_parents,
              N = N,
              selected_qtl = selected_qtl,
              SP = SP
            )
          )
          
          HC_result <-
            bind_rows(map(HC_result, 'parameter_df'), .id = 'ancestral')
          
          DC_result <-  list(
            Elite = cross_simulation_DC(
              parents = EL_parents,
              N = N,
              selected_qtl = selected_qtl,
              SP = SP
            ),
            Landrace = cross_simulation_DC(
              parents = LR_parents,
              N = N,
              selected_qtl = selected_qtl,
              SP = SP
            )
          )
          
          DC_result <-
            bind_rows(map(DC_result, 'parameter_df'), .id = 'ancestral')
          
          
          
          return(bind_rows(
            list(
              DC = DC_result,
              FC = FC_result,
              HC = HC_result
            ),
            .id = 'MD'
          ))
          
        }, N = N, P = P, selected_qtl = selected_qtl, future.seed = TRUE)
      
      
      
      return(bind_rows(rep_result, .id = 'rep'))
    }, P = P, selected_qtl = selected_qtl)
    
    names(N_result) <-  N_vector
    return(bind_rows(N_result, .id = 'N'))
    
  }, selected_qtl = selected_qtl)
  
  names(P_result) <-  P_vector
  
  return(bind_rows(P_result, .id = 'P'))
})

names(L_result) <- L_vector

main_sim_result <- bind_rows(L_result, .id = 'L')

message('Starting to save result')
print(Sys.time() - start)
data.table::fwrite(main_sim_result, 'output/main_sim_result.csv')

message('Finished to save result')
print(Sys.time() - start)
