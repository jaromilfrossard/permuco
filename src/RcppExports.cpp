// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// get_cluster_matrix
IntegerVector get_cluster_matrix(NumericMatrix distribution, double threshold);
RcppExport SEXP _permuco_get_cluster_matrix(SEXP distributionSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type distribution(distributionSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(get_cluster_matrix(distribution, threshold));
    return rcpp_result_gen;
END_RCPP
}
// get_clusterdepth_head
IntegerMatrix get_clusterdepth_head(IntegerMatrix cluster, String border);
RcppExport SEXP _permuco_get_clusterdepth_head(SEXP clusterSEXP, SEXP borderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type cluster(clusterSEXP);
    Rcpp::traits::input_parameter< String >::type border(borderSEXP);
    rcpp_result_gen = Rcpp::wrap(get_clusterdepth_head(cluster, border));
    return rcpp_result_gen;
END_RCPP
}
// get_clusterdepth_tail
IntegerMatrix get_clusterdepth_tail(IntegerMatrix cluster, String border);
RcppExport SEXP _permuco_get_clusterdepth_tail(SEXP clusterSEXP, SEXP borderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type cluster(clusterSEXP);
    Rcpp::traits::input_parameter< String >::type border(borderSEXP);
    rcpp_result_gen = Rcpp::wrap(get_clusterdepth_tail(cluster, border));
    return rcpp_result_gen;
END_RCPP
}
// depth_distribution_head
NumericMatrix depth_distribution_head(NumericMatrix distribution, IntegerMatrix head);
RcppExport SEXP _permuco_depth_distribution_head(SEXP distributionSEXP, SEXP headSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type distribution(distributionSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type head(headSEXP);
    rcpp_result_gen = Rcpp::wrap(depth_distribution_head(distribution, head));
    return rcpp_result_gen;
END_RCPP
}
// depth_distribution_tail
NumericMatrix depth_distribution_tail(NumericMatrix distribution, IntegerMatrix tail);
RcppExport SEXP _permuco_depth_distribution_tail(SEXP distributionSEXP, SEXP tailSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type distribution(distributionSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type tail(tailSEXP);
    rcpp_result_gen = Rcpp::wrap(depth_distribution_tail(distribution, tail));
    return rcpp_result_gen;
END_RCPP
}
// depth_distribution_unique
NumericMatrix depth_distribution_unique(NumericMatrix distribution, IntegerMatrix head, IntegerMatrix tail);
RcppExport SEXP _permuco_depth_distribution_unique(SEXP distributionSEXP, SEXP headSEXP, SEXP tailSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type distribution(distributionSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type head(headSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type tail(tailSEXP);
    rcpp_result_gen = Rcpp::wrap(depth_distribution_unique(distribution, head, tail));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_permuco_get_cluster_matrix", (DL_FUNC) &_permuco_get_cluster_matrix, 2},
    {"_permuco_get_clusterdepth_head", (DL_FUNC) &_permuco_get_clusterdepth_head, 2},
    {"_permuco_get_clusterdepth_tail", (DL_FUNC) &_permuco_get_clusterdepth_tail, 2},
    {"_permuco_depth_distribution_head", (DL_FUNC) &_permuco_depth_distribution_head, 2},
    {"_permuco_depth_distribution_tail", (DL_FUNC) &_permuco_depth_distribution_tail, 2},
    {"_permuco_depth_distribution_unique", (DL_FUNC) &_permuco_depth_distribution_unique, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_permuco(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
