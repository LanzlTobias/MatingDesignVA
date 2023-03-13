Description of the scripts and input file:

LD_decay:
Script for calculating the decay of GPD in the ancestral populations.

R_functions.R:
Source file for the functions I wrote in R.

cpp_functions.cpp:
Source file for the functions I wrote in C++.

main_simulation.R:
Script for running the main simulation. Needs the source files in R and C++ as well as the input created by create_pop_foundation_v10.R. Can be parallelized by changing the amount of workers.

ANOVA_Mean.R:
Script for summarizing the output of the main simulation. Performs an ANOVA for each setting and calculates the means and standard errors. Can be parallelized by changing the amount of workers.

preparatory_sim.R
Basically main_simulation.R for 100 replications and 100 sets of QTL positions

preparatory_sim_analysis.R
Basically ANOVA_Mean.R for the output of preparatory_sim.R

Appendix_A.R:
Tests for the formulas in Appendix A.

figure_S2.R:
Simulation for the data used in figure S2

D_distributions.R
Calculates the off-diagonal values of the matrix D in the ancestral populations.

figures.Rmd: 
R Notebook for the figures in the paper.

founding_pop.RData:
Input file for the main simulation. Contains masked information for the coding of the genotypes

