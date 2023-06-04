## Tests for the formulas in Appendix A.

set.seed(1)
library(ggplot2)
library(dplyr)
library(tidyr)
library(moments)

# Number of loci
l <- 10
# Probability of the reference allele being positive
p <- 0.5

# Vector with the number of positive reference alleles
k_vec <- 0:l

# P[K=k]
choose(l, k_vec) * p ** k_vec * (1 - p) ** (l - k_vec)
sum(choose(l, k_vec) * p ** k_vec * (1 - p) ** (l - k_vec))

# Vector of number of coupling phase linkages - repulsion phase linkages
y_vec <- -choose(l, 2):choose(l, 2)

# P[Y=y|K=k]
tmp <- sapply(y_vec, function(y) {
  sapply(k_vec, function(k) {
    y == 2 * k ** 2 - 2 * l * k + 0.5 * l * (l - 1)
  })
}) * 1
rownames(tmp) <- paste0('k_', 1:length(k_vec))
colnames(tmp) <- paste0('y_', 1:length(y_vec))
tmp
rowSums(tmp)
colSums(tmp)
# one y value for each k value
# maximum two values of k for each y value

# P[Y=y]
tmp <- sapply(y_vec, function(y) {
  sum(sapply(k_vec, function(k) {
    P_X_k <- choose(l, k) * p ** k * (1 - p) ** (l - k)
    p_Y_y_P_X_k <- y == 2 * k ** 2 - 2 * l * k + 0.5 * l * (l - 1)
    P_X_k * p_Y_y_P_X_k
  }))
})

sum(tmp)

tibble(y_vec, tmp) %>%
  ggplot(aes(y_vec, tmp)) +
  geom_point() +
  geom_vline(xintercept = -floor(l / 2), colour = 'red') +
  theme_bw(base_size = 15) +
  xlab('y') +
  ylab('P[Y=y]') +
  labs(caption = 'Probability density function of P[Y=y] for l=10. The red line shows -[l/2].') +
  scale_y_continuous(limit = c(0, 0.5))

## Data for fig.S6 with 1000 loci (additionally with p=0.75 and p=0.9)

probability_function <- function(p, l) {
  k_vec <- 0:l
  y_vec <- -choose(l, 2):choose(l, 2)
  
  P_Y_y <- vapply(y_vec, function(y) {
    K_y <- y == 2 * k_vec ** 2 - 2 * l * k_vec + 0.5 * l * (l - 1)
    K_y * choose(l, k_vec) * p ** k_vec * (1 - p) ** (l - k_vec)
  }, FUN.VALUE = numeric(length(k_vec)))
  return(colSums(P_Y_y))
}

l <- 1000
p_df <-
  bind_cols(lapply(c(0.5), probability_function, l = l))

names(p_df) <- paste0('p=', c(0.5))

p_df$y <- -choose(l, 2):choose(l, 2)

saveRDS(p_df, file = 'output/Appendix_A.RData')

suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(future.apply))

plan(multicore, workers = 10)


negative_probability_function <- function(l, p) {
  k_vec <- 0:l
  y_vec <- -choose(l, 2):(-1)
  
  P_Y_y <- vapply(y_vec, function(y) {
    K_y <- y == 2 * k_vec ** 2 - 2 * l * k_vec + 0.5 * l * (l - 1)
    sum(K_y * choose(l, k_vec) * p ** k_vec * (1 - p) ** (l - k_vec))
  }, FUN.VALUE = numeric(1))
  return(sum(P_Y_y))
}

l <- c(seq(100, 2500, 100), 5000)
neg_prob <- future_vapply(l, negative_probability_function, p = 0.5, FUN.VALUE = numeric(1))

result <- cbind(l, neg_prob)

saveRDS(result, "Appendix_A_negative_prob.RData")


