## Script for calculating the decay of LD in the ancestral populations.


setwd('~/LD/Paper')
library(dplyr)
library(magrittr)

# Functions ---------------------------------------------------------------

## Takes as input a individual x marker matrix (geno),
## a matrix of pairwise distances between the marker positions (marker_dist) 
## or the marker positions (marker_pos),
## and (optionally) a threshold until which distance pairwise LD should be computed
## Returns a the pairwise r2 and distance values
pairwise_LD = function(geno, dist_threshold = NULL, marker_dist = NULL, marker_pos = NULL){
  if(any(is.na(geno))){
    stop('Input contains missing marker information')
  }
  
  if(!is.null(marker_dist)){
    dist_vec = as.vector(marker_dist[lower.tri(marker_dist)])
  } else if(!is.null(marker_pos)){
    dist_vec = as.vector(dist(marker_pos)) 
  } else {
    stop('Give either marker positions or a distance matrix for the marker')
  }
  r2_mtx = cor(geno)
  r2_vec = r2_mtx[lower.tri(r2_mtx)]**2
  
  if(is.null(dist_threshold)){
    return(tibble(r2 = r2_vec[which(dist_vec > 0)],
                  distance = dist_vec[(dist_vec > 0)]))
  } else {
    return(tibble(r2 = r2_vec[which(dist_vec < dist_threshold & dist_vec > 0)],
                  distance = dist_vec[(dist_vec < dist_threshold & dist_vec > 0)]))
  }
}

## Takes as input two vectors of pairwise r2 and distance values (r2 & distance), 
## the number of individuals that were used for their calculation (nInd)
## and a threshold for the decay distance (r2_threshold)
## Returns the distance where the non-linear regression crosses the threshold
## and the p, which can be used to plot the curve
LD_decay = function(r2, distance, nInd, r2_threshold){
  p = smooth.fit(distance, r2, n = nInd)
  
  r2decay = t = try(uniroot(f = fitcurve,
                            n = nInd,
                            interval=c(0, 10000000000),
                            p = p,
                            threshold = r2_threshold))
  
  if (inherits(t, "try-error")){
    r2_decay = NA
  } else {
    r2_decay = r2decay[[1]]
  }
  
  return(tibble(p = p, r2_decay = r2_decay))
}

## Functions necessary inside LD_decay
fitcurve <- function(x,p,n, threshold){
  -1*threshold + ((10 + p*x)) / ((2+p*x) * (11 + p*x) ) *
    ( 1 + ( (3+ p*x) * (12 + 12 * p + p^2*x^2)) / ( n*(2+p*x) * (11 + p*x)))	
}


smooth.fit <- function(overallDist, overallr2, n) {
  nonlinearoverall <- nls(overallr2 ~ ((10 + p * overallDist))/((2 +
                                                                   p * overallDist) * (11 + p * overallDist)) * (1 +
                                                                                                                   ((3 + p * overallDist) * (12 + 12 * p + p^2 * overallDist^2))/(n *
                                                                                                                                                                                    (2 + p * overallDist) * (11 + p * overallDist))),
                          start = list(p = 0.001), contro=list(maxiter=1000),
                          na.action = na.exclude)
  p <- coef(nonlinearoverall)
  return(p)	        
} 


# Calculation -------------------------------------------------------------

full_pop <-  readRDS('input/founding_pop.RData')

EL <-  lapply(1:length(full_pop$haplotypes), function(x){
  full_pop$haplotypes[[x]][grepl('Elite_', 
                                  rownames(full_pop$haplotypes[[x]])),]/2
})

LR <-  lapply(1:length(full_pop$haplotypes), function(x){
  full_pop$haplotypes[[x]][grepl('Landrace', 
                                  rownames(full_pop$haplotypes[[x]])),]/2
})

genMap <- full_pop$genMap


LD_EL <- lapply(1:length(full_pop$haplotypes), function(chr){
  pairwise_LD(geno = EL[[chr]],
              marker_pos = genMap[[chr]])
})


LD_LR <- lapply(1:length(full_pop$haplotypes), function(chr){
  pairwise_LD(geno = LR[[chr]],
              marker_pos = genMap[[chr]])
})


# Average for Elite -------------------------------------------------------
tmp <- bind_rows(LD_EL) %$%
LD_decay(r2 = r2,
         distance = distance,
         nInd = 115,
         r2_threshold = 0.1) %>% 
  pull(r2_decay)

message('Average LD decay for Elite:', round(tmp * 100, 2), 'cM')


# Per chr for Elite -------------------------------------------------------

elite <- lapply(LD_EL, function(df){
  LD_decay(r2 = df$r2,
            distance = df$distance,
            nInd = 115,
            r2_threshold = 0.1)
}) %>% 
  bind_rows(.id = 'chr')


# Average for Landrace ----------------------------------------------------------
tmp <- bind_rows(LD_LR) %$%
  LD_decay(r2 = r2,
            distance = distance,
            nInd = 115,
            r2_threshold = 0.1) %>% 
  pull(r2_decay)

message('Average LD decay for Landrace:', round(tmp * 100, 2), 'cM')

# Per chr for Landrace ----------------------------------------------------

landrace <- lapply(LD_LR, function(df){
  LD_decay(r2 = df$r2,
            distance = df$distance,
            nInd = 115,
            r2_threshold = 0.1)
}) %>% 
  bind_rows(.id = 'chr')


saveRDS(bind_rows(list(Landrace = landrace,
                       Elite = elite), .id = 'ancestral') ,
        file = 'output/summary_LD_decay_chr.RData')

saveRDS(bind_rows(list(Landrace = bind_rows(LD_LR, .id = 'chr'),
                       Elite = bind_rows(LD_EL, .id = 'chr')), 
                  .id = 'ancestral'),
        file = 'output/raw_LD_decay.RData')