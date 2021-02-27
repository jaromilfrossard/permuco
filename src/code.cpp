#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
IntegerVector get_cluster_matrix(NumericMatrix distribution, double threshold){
    // initialize an accumulator variable
    double acc = 0;

    // initialize the result matrix
    IntegerMatrix res(distribution.nrow(),distribution.ncol());

    for(int rowi = 0; rowi < res.nrow(); rowi++){
      if(distribution(rowi,0)>threshold){
        res(rowi,0) = 1;
        acc = 1;
      }else{
        res(rowi,0) = 0;
        acc = 0;
      }

      for(int coli = 1; coli < res.ncol(); coli++){
        if(distribution(rowi,coli)<=threshold){
          res(rowi,coli) = 0;
        }else{
          if(res(rowi,coli-1)==0){
            acc = acc+1;
          }
          res(rowi,coli) = acc;
        }
      }

    }
    return res;


  }
