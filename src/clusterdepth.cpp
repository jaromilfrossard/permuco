#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
IntegerVector get_cluster_matrix(NumericMatrix distribution, double threshold, String side){
    // initialize an accumulator variable
    double acc = 0;

    // initialize the result matrix
    IntegerMatrix res(distribution.nrow(),distribution.ncol());

    if((side == "all")||(side == "ending")){
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

    }else if(side == "starting"){
      for(int rowi = 0; rowi < res.nrow(); rowi++){
        acc = 0;
        for(int coli = 0; coli < res.ncol(); coli++){
          if((coli==0)|(distribution(rowi,coli)<threshold)){
            res(rowi,coli) = 0;
          }else if((coli>0) & (distribution(rowi,coli)>threshold)){
            if(distribution(rowi,coli-1)<threshold){
              acc ++;
            }
            res(rowi,coli) = acc;

          }

        }

      }

    }

    if(side == "ending"){
      for(int rowi = 0; rowi < res.nrow(); rowi++){
        int coli = res.ncol()-1;
        while((res(rowi,coli)>0)&(coli>0)){
          res(rowi,coli)=0;
          coli = coli-1;
        }
      }
    }



    return res;


  }




//[[Rcpp::export]]
IntegerMatrix get_clusterdepth_head(IntegerMatrix cluster, String border){

  double acc = 0;


  // initialize the result matrix
  IntegerMatrix res(cluster.nrow(),cluster.ncol());

  for(int rowi = 0; rowi < res.nrow(); rowi++){
    acc = 0;
    for(int coli = 0; coli < res.ncol(); coli++){
      if(cluster(rowi,coli)==0){
        acc=0;
      }else if(cluster(rowi,coli)>0){
        acc = acc+1;
      }
      res(rowi,coli) = acc;
    }
  }

  // handle border

  if(border == "rm"){
    for(int rowi = 0; rowi < res.nrow(); rowi++){
      int coli = 0;
      while((res(rowi,coli)>=(coli+1))&&(coli<res.ncol())){
        res(rowi,coli)=0;
        coli = coli+1;
      }
    }
  }


  if(border == "reverse"){
    for(int rowi = 0; rowi < res.nrow(); rowi++){
      if(res(rowi,0)>0){

        //length of border cluster
        int lcli=0;
        while((res(rowi,lcli)>0) && (lcli<res.ncol())){
          lcli = lcli+1;
        }


        for(int coli = 0; coli < lcli; coli++){
          res(rowi,coli) = lcli-coli;
        }
      }
    }
  }




  return res;

}


//[[Rcpp::export]]
IntegerMatrix get_clusterdepth_tail(IntegerMatrix cluster, String border){

  double acc = 0;


  // initialize the result matrix
  IntegerMatrix res(cluster.nrow(),cluster.ncol());

  for(int rowi = 0; rowi < res.nrow(); rowi++){
    acc = 0;
    for(int coli = res.ncol()-1; coli >= 0; coli--){
      if(cluster(rowi,coli)==0){
        acc=0;
      }else if(cluster(rowi,coli)>0){
        acc = acc+1;
      }
      res(rowi,coli) = acc;
    }
  }

  if(border == "rm"){
    for(int rowi = 0; rowi < res.nrow(); rowi++){
      int coli = res.ncol()-1;
      while( (res(rowi,abs(coli))>0)&(coli>=0)){
        res(rowi,coli)=0;
        coli = coli-1;
      }
    }
  }


  if(border == "reverse"){
    for(int rowi = 0; rowi < res.nrow(); rowi++){
      if(res(rowi,res.ncol()-1)>0){

        //length of border cluster
        int lcli=0;
        while( (res(rowi,abs(res.ncol()-1-lcli))>0) && (lcli<res.ncol())){
          lcli = lcli+1;
        }

        for(int iend = 0; iend < lcli; iend++){
          res(rowi,res.ncol()-iend-1) = lcli-iend;
        }
      }
    }
  }




  return res;

}


//[[Rcpp::export]]
NumericMatrix depth_distribution_head(NumericMatrix distribution, IntegerMatrix head){
  int max = 0;
  for(int rowi = 0; rowi < head.nrow(); rowi++){
    for(int coli = 0; coli < head.ncol(); coli++){
      if(max < head(rowi,coli)){
        max = head(rowi,coli);
      }
    }
  }

  NumericMatrix res(head.nrow(),max);


  for(int rowi = 0; rowi < head.nrow(); rowi++){
    for(int coli = 0; coli < head.ncol(); coli++){
      if(head(rowi,coli)>0){
        if(res(rowi,head(rowi,coli)-1) < distribution(rowi,coli)){
          res(rowi,head(rowi,coli)-1) = distribution(rowi,coli);
        }
      }
    }
  }




  return res;

}


//[[Rcpp::export]]
NumericMatrix depth_distribution_tail(NumericMatrix distribution, IntegerMatrix tail){
    int max = 0;
    for(int rowi = 0; rowi < tail.nrow(); rowi++){
      for(int coli = 0; coli < tail.ncol(); coli++){
        if(max < tail(rowi,coli)){
          max = tail(rowi,coli);
        }
      }
    }

    NumericMatrix res(tail.nrow(),max);

    for(int rowi = 0; rowi < tail.nrow(); rowi++){
      for(int coli = 0; coli < tail.ncol(); coli++){
        if(tail(rowi, coli) > 0){
          if(res(rowi, res.ncol()-1 - (tail(rowi, coli) -1) ) < distribution(rowi,coli)){
            res(rowi, res.ncol()-1 - (tail(rowi, coli) -1) ) = distribution(rowi,coli);
          }
        }
      }
    }

    return res;

}



//[[Rcpp::export]]
NumericMatrix depth_distribution_unique(NumericMatrix distribution, IntegerMatrix head, IntegerMatrix tail){
    int max = 0;
    for(int rowi = 0; rowi < head.nrow(); rowi++){
      for(int coli = 0; coli < head.ncol(); coli++){
        if(max < head(rowi,coli)){
          max = head(rowi,coli);
        }
        if(max < tail(rowi,coli)){
          max = tail(rowi,coli);
        }
      }
    }

    NumericMatrix res(head.nrow(),max);


    for(int rowi = 0; rowi < head.nrow(); rowi++){
      for(int coli = 0; coli < tail.ncol(); coli++){
        if(head(rowi,coli)>0){
          if(res(rowi,head(rowi,coli)-1) < distribution(rowi,coli)){
            res(rowi,head(rowi,coli)-1) = distribution(rowi,coli);
          }
        }
        if(tail(rowi, coli) > 0){
          if(res(rowi, tail(rowi, coli) -1) < distribution(rowi,coli)){
            res(rowi, tail(rowi, coli) -1) = distribution(rowi,coli);
          }
        }

      }
    }




    return res;

  }





