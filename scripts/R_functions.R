## Source file for the functions in R

suppressPackageStartupMessages(library(AlphaSimR))
suppressPackageStartupMessages(library(dplyr))



## Calculates D based on X 
D_calc <- function(X){
  n <- nrow(X)
  j <- matrix(1, nrow = n)
  p <- crossprod(j, X) * (1/(2 * n))
  Z <- X - 2 * j %*% p
  D <- crossprod(Z, Z) * (1/n)
  return(D)
}


## Calculates the expectations and variances conditional on X (uses D as input matrix) 
conditional_var_exp = function(mtx, selected_qtl_per_chr) {
  n_fixed <- sum(round(diag(mtx), 10) == 0)
  var_Vg <- var_Vg(mtx)
  var_Cw <- var_Cw(mtx, selected_qtl_per_chr)
  var_Cb <- var_Cb(mtx, selected_qtl_per_chr)
  var_C <-  var_Cw + var_Cb
  var_VA <- var_Vg + var_C
  E_Vg <- E_Vg(mtx)
  var_b <- var_C / (var_Vg + E_Vg**2)
  var_bw <- var_Cw / (var_Vg + E_Vg**2)
  var_bb <- var_Cb / (var_Vg + E_Vg**2)
  
  result_df = tibble(var_Vg, var_Cw, var_Cb, var_C,
                     E_Vg, var_b, var_bw, var_bb, var_VA,
                     nPoly_QTL = sum(selected_qtl_per_chr) - n_fixed,
                     )
  return(result_df)
}



## Calculates the additive genetic variances, its components, as well as the relative contribution of GPD
## For a given matrix D and vector a
VA_decomposition = function(a, mtx, selected_qtl_per_chr) {
  Vg <-  calc_Vg(mtx, a)
  Cw <-  calc_Cw(mtx, selected_qtl_per_chr, a)
  Cb <-  calc_Cb(mtx, selected_qtl_per_chr, a)
  C <-  Cw + Cb
  b <-  C / Vg
  bb <- Cb / Vg
  bw <- Cw / Vg
  VA <- Vg + C
  
  result_df = tibble(Vg, VA, C, Cw, Cb, b, bw, bb)
  return(result_df)
}

# FC ---------------------------------------------------------------------
## Function for the simulation after the factorial crossing scheme
cross_simulation_FC = function(parents, N, selected_qtl, SP) {
  
  parameter_list <- list()
  
  ## Split the parents into two sets
  p1 <- parents@id[1:parents@nInd %% 2 == 0]
  p2 <- parents@id[1:parents@nInd %% 2 != 0]
  
  ## Make a factorial crossing scheme
  crossplan <-  as.matrix(expand.grid(p1, p2))
  
  if (N <= nrow(crossplan)) {
    crossplan = crossplan[sample(1:nrow(crossplan), size = N),]
  } else {
    tmp1 = rep(crossplan[, 1], ceiling(N / nrow(crossplan)))
    tmp2 = rep(crossplan[, 2], ceiling(N / nrow(crossplan)))
    crossplan = cbind(tmp1, tmp2)[sample(1:length(tmp1), size = N), ]
  }


  
  ## G1
  G1 <- makeCross(
    pop = parents,
    crossplan,
    nProgeny = 1,
    simParam = SP
  )
  
  G1_DH <-  makeDH(
    pop = G1,
    nDH = 1,
    keepParents = TRUE,
    simParam = SP
  )
  
  
  ## G2
  G2 <- randCross(G1,
                  nCrosses = N,
                  nProgeny = 1,
                  simParam = SP)
  
  G2_DH <- makeDH(
    pop = G2,
    nDH = 1,
    keepParents = TRUE,
    simParam = SP
  )
  
  
  ## G3
  G3 <-  randCross(G2,
                   nCrosses = N,
                   nProgeny = 1,
                   simParam = SP)
  
  G3_DH <- makeDH(
    pop = G3,
    nDH = 1,
    keepParents = TRUE,
    simParam = SP
  )
  
  
  ## G4
  G4 <- randCross(G3,
                  nCrosses = N,
                  nProgeny = 1,
                  simParam = SP)
  
  G4_DH <- makeDH(
    pop = G4,
    nDH = 1,
    keepParents = TRUE,
    simParam = SP
  )
  
  ## Assign numbers to the QTL positions
  chr_qtl_map <-
    split(1:2500, do.call(c, sapply(1:10, function(chr) {
      rep(chr, qtl_per_chr[chr])
    })))
  
  G_list <- list(G1_DH, G2_DH, G3_DH, G4_DH)
  
  ## Analysis
  
  for (i in 1:4) {
    ## Calculate D based on the 2500 QTL positions
    D <- D_calc(pullSegSiteGeno(pop = G_list[[i]],
                            simParam = SP))
    
    
    ## Subset D and calculate the conditional variances and expectations
    parameter_list[[as.character(i)]] <-
      bind_rows(lapply(1:sets, function(x) {
        selected_qtl_per_chr <-  sapply(chr_qtl_map, function(chr) {
          sum(selected_qtl[[x]] %in% chr)
        })
        
        conditional_var_exp(D[selected_qtl[[x]],
                            selected_qtl[[x]]],
                          selected_qtl_per_chr)
      }), .id = 'set')
    
    
  }
  
  ## Combining results
  
  parameter_df <- bind_rows(parameter_list, .id = 'G')
  
  
  return(list(
    parameter_df = parameter_df
  ))
}


# HC ----------------------------------------------------------------------
## Function for the simulation after the half-diallel crossing scheme
cross_simulation_HC = function(parents, N, selected_qtl, SP) {
  
  parameter_list <- list()
  
  ## make a half-diallel crossing scheme
  crossplan <-  as.matrix(expand.grid(parents@id, parents@id))
  
  crossplan <- crossplan[as.numeric(crossplan[, 1]) < as.numeric(crossplan[, 2]), ,drop = FALSE]
  
  if (N <= nrow(crossplan)) {
    crossplan = crossplan[sample(1:nrow(crossplan), size = N),]
  } else {
    tmp1 = rep(crossplan[, 1], ceiling(N / nrow(crossplan)))
    tmp2 = rep(crossplan[, 2], ceiling(N / nrow(crossplan)))
    crossplan = cbind(tmp1, tmp2)[sample(1:length(tmp1), size = N), ]
  }
  
  
  
  ## G1
  G1 <- makeCross(
    pop = parents,
    crossplan,
    nProgeny = 1,
    simParam = SP
  )
  
  G1_DH <-  makeDH(
    pop = G1,
    nDH = 1,
    keepParents = TRUE,
    simParam = SP
  )
  
  
  ## G2
  G2 <- randCross(G1,
                  nCrosses = N,
                  nProgeny = 1,
                  simParam = SP)
  
  G2_DH <- makeDH(
    pop = G2,
    nDH = 1,
    keepParents = TRUE,
    simParam = SP
  )
  
  
  ## G3
  G3 <-  randCross(G2,
                   nCrosses = N,
                   nProgeny = 1,
                   simParam = SP)
  
  G3_DH <- makeDH(
    pop = G3,
    nDH = 1,
    keepParents = TRUE,
    simParam = SP
  )
  
  
  ## G4
  G4 <- randCross(G3,
                  nCrosses = N,
                  nProgeny = 1,
                  simParam = SP)
  
  G4_DH <- makeDH(
    pop = G4,
    nDH = 1,
    keepParents = TRUE,
    simParam = SP
  )
  
  
  ## Assign numbers to the QTL positions
  chr_qtl_map <-
    split(1:2500, do.call(c, sapply(1:10, function(chr) {
      rep(chr, qtl_per_chr[chr])
    })))
  
  G_list <- list(G1_DH, G2_DH, G3_DH, G4_DH)
  
  ## Analysis
  
  for (i in 1:4) {
    ## Calculate D based on the 2500 QTL positions
    D <- D_calc(pullSegSiteGeno(pop = G_list[[i]],
                                simParam = SP))
    
    
    ## Subset D and calculate the conditional variances and expectations
    parameter_list[[as.character(i)]] <-
      bind_rows(lapply(1:sets, function(x) {
        selected_qtl_per_chr <-  sapply(chr_qtl_map, function(chr) {
          sum(selected_qtl[[x]] %in% chr)
        })
        
        conditional_var_exp(D[selected_qtl[[x]],
                              selected_qtl[[x]]],
                            selected_qtl_per_chr)
      }), .id = 'set')
    
    
  }
  
  ## Combining results
  
  parameter_df <- bind_rows(parameter_list, .id = 'G')
  
  
  return(list(
    parameter_df = parameter_df
  ))
}


# DC ----------------------------------------------------------------------
## Function for the simulation after the disjoint crossing scheme
cross_simulation_DC = function(parents, N, selected_qtl, SP) {
  
  parameter_list <- list()
  
  ## Split the parents into two sets
  p1 <- parents@id[1:parents@nInd %% 2 == 0]
  p2 <- parents@id[1:parents@nInd %% 2 != 0]
  
  
  ## Make a disjoint crossing scheme
  crossplan <-  cbind(p1, p2)
  
  if (N <= nrow(crossplan)) {
    crossplan = crossplan[sample(1:nrow(crossplan), size = N),]
  } else {
    tmp1 = rep(crossplan[, 1], ceiling(N / nrow(crossplan)))
    tmp2 = rep(crossplan[, 2], ceiling(N / nrow(crossplan)))
    crossplan = cbind(tmp1, tmp2)[sample(1:length(tmp1), size = N), ]
  }  

  
  ## G1
  G1 <- makeCross(
    pop = parents,
    crossplan,
    nProgeny = 1,
    simParam = SP
  )
  
  G1_DH <-  makeDH(
    pop = G1,
    nDH = 1,
    keepParents = TRUE,
    simParam = SP
  )
  
  
  ## G2
  G2 <- randCross(G1,
                  nCrosses = N,
                  nProgeny = 1,
                  simParam = SP)
  
  G2_DH <- makeDH(
    pop = G2,
    nDH = 1,
    keepParents = TRUE,
    simParam = SP
  )
  
  
  ## G3
  G3 <-  randCross(G2,
                   nCrosses = N,
                   nProgeny = 1,
                   simParam = SP)
  
  G3_DH <- makeDH(
    pop = G3,
    nDH = 1,
    keepParents = TRUE,
    simParam = SP
  )
  
  
  ## G4
  G4 <- randCross(G3,
                  nCrosses = N,
                  nProgeny = 1,
                  simParam = SP)
  
  G4_DH <- makeDH(
    pop = G4,
    nDH = 1,
    keepParents = TRUE,
    simParam = SP
  )
  
  
  ## Assign numbers to the QTL positions
  chr_qtl_map <-
    split(1:2500, do.call(c, sapply(1:10, function(chr) {
      rep(chr, qtl_per_chr[chr])
    })))
  
  G_list <- list(G1_DH, G2_DH, G3_DH, G4_DH)
  
  ## Analysis
  
  for (i in 1:4) {
    ## Calculate D based on the 2500 QTL positions
    D <- D_calc(pullSegSiteGeno(pop = G_list[[i]],
                                simParam = SP))
    
    
    ## Subset D and calculate the conditional variances and expectations
    parameter_list[[as.character(i)]] <-
      bind_rows(lapply(1:sets, function(x) {
        selected_qtl_per_chr <-  sapply(chr_qtl_map, function(chr) {
          sum(selected_qtl[[x]] %in% chr)
        })
        
        conditional_var_exp(D[selected_qtl[[x]],
                              selected_qtl[[x]]],
                            selected_qtl_per_chr)
      }), .id = 'set')
    
    
  }
  
  ## Combining results
  
  parameter_df <- bind_rows(parameter_list, .id = 'G')
  
  
  return(list(
    parameter_df = parameter_df
  ))
}



