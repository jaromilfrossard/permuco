#include <Rcpp.h>
using namespace Rcpp;


// //[[Rcpp::export]]
// NumericVector vector_extend(NumericVector x, double threshold){
//   NumericVector res(x.length());
//   int acc = 0;
//   for(int i = 0; i<x.length();i++){
//     if(x(i) <= threshold){
//       res(i) = 0;
//       acc=0;
//     }
//     if(x(i)>threshold){
//       acc = acc+1;
//       for(int j =0; j<acc;j++){
//         res(i-j) = acc;
//       }
//     }
//
//   }
//   return res;
// }

//[[Rcpp::export]]
NumericVector vector_extend(NumericVector x, double threshold){
  NumericVector res(x.length());
  int acc = 0;
  for(int i = 0; i<x.length();i++){
    if(x(i) <= threshold){
      res(i) = 0;
      acc=0;
    }
    if(x(i)>threshold){
      acc = acc+1;
      res(i) = acc;
    }

  }
  acc = res(res.length()-1);
  for(int i = res.length()-2; i>=0;i--){
    if(res(i) >0){
      if(acc>0){
        res(i) = acc;
      }
    }
    acc = res(i);
  }
  return res;
}




//[[Rcpp::export]]
NumericMatrix tfce_distribution(NumericMatrix distribution, double E, double H, double dh, NumericVector dhi){

  NumericMatrix res(distribution.nrow(),distribution.ncol());


  for(int rowi = 0; rowi < res.nrow(); rowi++){
    int rowmax = max(distribution(rowi, _));
    for(int ddi = 0; dhi(ddi)<rowmax;ddi++){
      res(rowi,_) = res(rowi,_) + pow(vector_extend(distribution(rowi, _),dhi(ddi)),E)*pow(dhi[ddi],H)*dh ;

    }

  }
  return res;


}



