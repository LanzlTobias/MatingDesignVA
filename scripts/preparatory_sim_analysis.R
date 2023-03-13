## Same as ANOVA_Mean.R for the output of preparatory_sim.R

library(lme4)
library(tidyr)
library(dplyr)
library(purrr)
library(broom)
library(future)
library(future.apply)

setwd('~/LD/Paper')

## How many cores should be used?
plan(multicore, workers = 50)

## Read in the data
main_sim_result <-  data.table::fread('output/preparatory_sim_result.csv') %>%  
  mutate(P = factor(P, levels = c(2, 4, 8, 16),
                    labels = paste0('P=', c(2, 4, 8, 16))),
         N = factor(N, levels = c(50, 250, 1000),
                    labels = paste0('N=', c(50, 250, 1000))),
         G = factor(G, levels = 1:4,
                    labels = paste0('G', 1:4, '_DH')),
         L = factor(L, levels = c(50, 250, 1000),
                    labels = paste0('L=', c(50, 250, 1000))),
         MD = factor(MD,
                     levels = c('DC', 'FC', 'HC')))

## Number of replications and sets
replications <- length(unique(main_sim_result$rep))
sets <- length(unique(main_sim_result$set))

## Function for calculating the variance components and means
varComp <- function(df){
  lmer(value ~ (1|rep) + (1|set), data = df) %>% 
    tidy() %>% 
    filter(term != '(Intercept)') %>% 
    mutate(varComp = estimate**2) %>% 
    select(group, varComp) %>% 
    pivot_wider(values_from = 'varComp', names_from = 'group') %>% 
    mutate(Mean = mean(df$value, na.rm = TRUE))
}

## Setting up a nested data.frame for every factor combination
nested_df <- main_sim_result %>% 
  pivot_longer(var_Vg:var_VA, names_to = 'parameter') %>% 
  group_by(MD, ancestral, P, N, G, L, parameter) %>% 
  nest()

variance_components <- nested_df

## Calculating variance components for replication and set for every factor combination
variance_components$model <- future_lapply(nested_df$data, varComp)

## Calculating standard errors
variance_components <- variance_components %>% 
  unnest(model) %>% 
  select(-data) %>% 
  mutate(SE = sqrt(rep/replications + set/sets + Residual/(replications*sets))) %>% 
  ungroup() 

saveRDS(variance_components, file = 'output/preparatory_summary_df.RData')
message('Finished variance components')

