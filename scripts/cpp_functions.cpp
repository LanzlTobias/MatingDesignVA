// Source file for the functions in C++.

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double var_Vg(NumericMatrix mtx) {
  int N = mtx.nrow();
  double sum_D_sq = 0;
  for (int i = 0; i < N; i++){
    sum_D_sq += pow(mtx(i, i), 2);
  }
  return 2 * sum_D_sq;
}

// [[Rcpp::export]]
double E_Vg(NumericMatrix mtx) {
  int N = mtx.nrow();
  double sum_D = 0;
  for (int i = 0; i < N; i++){
    sum_D += mtx(i, i);
  }
  return sum_D;
}

// [[Rcpp::export]]
double var_Cw(NumericMatrix mtx, NumericVector qtl_per_chr) {
  int nrow = mtx.nrow();
  int n_chr = qtl_per_chr.length();
  IntegerVector start(n_chr);
  IntegerVector stop(n_chr);
  double sum_D_sq = 0;
  
  start[0] = 1;
  for (int i = 1; i < n_chr; i++){
    start[i] = start[i-1] + qtl_per_chr[i-1];
  }
  for (int i = 0; i < n_chr - 1; i++){
    stop[i] = start[i+1] - 1;
  }
  stop[n_chr-1] = nrow;
  
  for (int chr = 0; chr < n_chr; chr++){
    if (qtl_per_chr[chr] != 0){
      NumericMatrix mtx_chr = mtx(Range(start[chr]-1, stop[chr]-1),
                                  Range(start[chr]-1, stop[chr]-1));
      int n_loci_chr = stop[chr] - (start[chr]-1);
      
      for (int i = 0; i < n_loci_chr; i++){
        for (int j = 0; j < n_loci_chr; j++){
          if (i < j){
            sum_D_sq += 2 *pow(mtx_chr(i,j), 2);
          }
        }
      }
    }
  }
  return 2 * sum_D_sq;
}

// [[Rcpp::export]]
double var_Cb(NumericMatrix mtx, NumericVector qtl_per_chr) {
  int nrow = mtx.nrow();
  int n_chr = qtl_per_chr.length();
  IntegerVector start(n_chr);
  IntegerVector stop(n_chr);
  double sum_D_sq = 0;
  
  start[0] = 1;
  for (int i = 1; i < n_chr; i++){
    start[i] = start[i-1] + qtl_per_chr[i-1];
  }
  for (int i = 0; i < n_chr - 1; i++){
    stop[i] = start[i+1] - 1;
  }
  stop[n_chr-1] = nrow;
  
  for (int chr1 = 0; chr1 < n_chr; chr1++){
    if (qtl_per_chr[chr1] != 0){
      for (int chr2 = 0; chr2 < n_chr; chr2++){
        if (qtl_per_chr[chr2] != 0){
          NumericMatrix mtx_chr = mtx(Range(start[chr1]-1, stop[chr1]-1),
                                      Range(start[chr2]-1, stop[chr2]-1));
          int mtx_row = mtx_chr.nrow();
          int mtx_col = mtx_chr.ncol();
          
          for (int i = 0; i < mtx_row; i++){
            for (int j = 0; j < mtx_col; j++){
              if (chr1 < chr2){
                sum_D_sq += 2 * pow(mtx_chr(i, j), 2);
              }
            }
          }
        }
      }
    }
  }
  return 2 * sum_D_sq;
}


// [[Rcpp::export]]
double calc_Vg(NumericMatrix mtx, NumericVector alpha) {
  int nrow = mtx.nrow();
  double sum_var = 0;
  for (int i = 0; i < nrow; i++){
    sum_var += mtx(i,i) * alpha(i) * alpha(i);
  }
  return sum_var;
}

// [[Rcpp::export]]   
double calc_Cw(NumericMatrix mtx, NumericVector qtl_per_chr, NumericVector alpha) {
  int nrow = mtx.nrow();
  int n_chr = qtl_per_chr.length();
  IntegerVector start(n_chr);
  IntegerVector stop(n_chr);
  double Cw = 0;
  int n_loci = 0;
  for (int i = 0; i < n_chr; i++){
    n_loci += (pow(qtl_per_chr[i], 2) - qtl_per_chr[i]);
  }
  
  start[0] = 1;
  for (int i = 1; i < n_chr; i++){
    start[i] = start[i-1] + qtl_per_chr[i-1];
  }
  for (int i = 0; i < n_chr - 1; i++){
    stop[i] = start[i+1] - 1;
  }
  stop[n_chr-1] = nrow;
  
  for (int chr = 0; chr < n_chr; chr++){
    NumericMatrix mtx_chr = mtx(Range(start[chr]-1, stop[chr]-1),
                                Range(start[chr]-1, stop[chr]-1));
    NumericVector alpha_chr = alpha[seq(start[chr]-1, stop[chr]-1)];
    int n_loci_chr = stop[chr] - (start[chr]-1);
    
    for (int i = 0; i < n_loci_chr; i++){
      for (int j = 0; j < n_loci_chr; j++){
        if (i == j) Cw = Cw;
        else Cw += mtx_chr(i,j) * alpha_chr(i) * alpha_chr(j);
      }
    }
  }
  return Cw;
}

// [[Rcpp::export]]   
double calc_Cw_new(NumericMatrix mtx, NumericVector qtl_per_chr, NumericVector alpha) {
  int nrow = mtx.nrow();
  int n_chr = qtl_per_chr.length();
  IntegerVector start(n_chr);
  IntegerVector stop(n_chr);
  double Cw = 0;
  int n_loci = 0;
  for (int i = 0; i < n_chr; i++){
    n_loci += (pow(qtl_per_chr[i], 2) - qtl_per_chr[i]);
  }
  
  start[0] = 1;
  for (int i = 1; i < n_chr; i++){
    start[i] = start[i-1] + qtl_per_chr[i-1];
  }
  for (int i = 0; i < n_chr - 1; i++){
    stop[i] = start[i+1] - 1;
  }
  stop[n_chr-1] = nrow;
  
  for (int chr = 0; chr < n_chr; chr++){
    NumericMatrix mtx_chr = mtx(Range(start[chr]-1, stop[chr]-1),
                                Range(start[chr]-1, stop[chr]-1));
    NumericVector alpha_chr = alpha[seq(start[chr]-1, stop[chr]-1)];
    int n_loci_chr = stop[chr] - (start[chr]-1);
    
    for (int i = 0; i < n_loci_chr; i++){
      for (int j = 0; j < n_loci_chr; j++){
        if (i < j) Cw += 2 * mtx_chr(i,j) * alpha_chr(i) * alpha_chr(j);
      }
    }
  }
  return Cw;
}

// [[Rcpp::export]]          
double calc_Cb(NumericMatrix mtx, NumericVector qtl_per_chr, NumericVector alpha) {
  int nrow = mtx.nrow();
  int n_chr = qtl_per_chr.length();
  IntegerVector start(n_chr);
  IntegerVector stop(n_chr);
  double Cb = 0;
  
  start[0] = 1;
  for (int i = 1; i < n_chr; i++){
    start[i] = start[i-1] + qtl_per_chr[i-1];
  }
  for (int i = 0; i < n_chr - 1; i++){
    stop[i] = start[i+1] - 1;
  }
  stop[n_chr-1] = nrow;
  
  for (int chr1 = 0; chr1 < n_chr; chr1++){
    for (int chr2 = 0; chr2 < n_chr; chr2++){
      NumericMatrix mtx_chr = mtx(Range(start[chr1]-1, stop[chr1]-1),
                                  Range(start[chr2]-1, stop[chr2]-1));
      NumericVector alpha_chr1 = alpha[seq(start[chr1]-1, stop[chr1]-1)];
      NumericVector alpha_chr2 = alpha[seq(start[chr2]-1, stop[chr2]-1)];
      int mtx_row = mtx_chr.nrow();
      int mtx_col = mtx_chr.ncol();
      
      for (int i = 0; i < mtx_row; i++){
        for (int j = 0; j < mtx_col; j++){
          if (chr1 < chr2) Cb += 2 * mtx_chr(i,j) * alpha_chr1(i) * alpha_chr2(j);
        }
      }
    }
  }
  return Cb;
}


