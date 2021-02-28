// #include <Rcpp.h>
// using namespace Rcpp;
//
//
// // [[Rcpp::export]]
// NumericVector sort_rcpp(NumericVector x) {
//   std::vector<double> tmp = Rcpp::as< std::vector<double> > (x);
//   std::sort(tmp.begin(), tmp.end());
//   return wrap(tmp);
// }
//
//
//
// // [[Rcpp::export]]
// NumericVector rev_(NumericVector x) {
//   NumericVector res = clone<NumericVector>(x);
//   std::reverse(res.begin(), res.end());
//   return res;
// }
//
// // [[Rcpp::export]]
// IntegerVector rank_(NumericVector x) {
//   return match(x, sort_rcpp(x));
// }
//
// // [[Rcpp::export]]
// IntegerVector rank_rev(NumericVector x) {
//   return match(rev_(x), rev_(sort_rcpp(x)));
// }
//
// // [[Rcpp::export]]
// IntegerVector rank_up(NumericVector x) {
//   IntegerVector res = rank_rev(x);
//   res = res.length()-rev_(res);
//   return res;
// }

