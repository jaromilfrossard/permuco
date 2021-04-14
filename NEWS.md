# permuco 1.1.1 (github)

 * add signflip test
 * add replication for rANOVA
 * change Pmat(): argument "type" and "counting"
 * plot.lmperm(): ... argument can be used (eg xlim, ylim)
 * add testthat, pkgdown
 * rcpp optimization: get_cluster(), compute_tfce()
 * add compute_clusterdepth()
 * add alternative to compute_tfce()
 * add warnings to minp and troendle


# permuco 1.1.0 (cran)


 * correction of compute_troendle: pvalue for all tests
 * New multiple comparisons procedure: min-P.
 * new display of the output of clusterlm():all test or pseudo-clusters
 * delete the table of the output of clusterlm. Access the table using summary(clusterlm(...))
 * user access to functions compute_tfce, compute_clustermass, compute_troendle, compute_minP
 * update vignette permuco_tutorial.pdf

# permuco 1.0.2


 * correction benjaminin to benjamini
 * adding parametric (uncorrected) pvalues for signal
 * add argument "nbbaselinepts" and "nbptsperunit" in plot.clusterlm
 * change argument names: bilateral in two.sided, left in less,right in greater
 * update vignette permuco_tutorial.pdf


# permuco 1.0.1


 *  add vignette permuco_tutorial.pdf


# permuco 1.0.0

 * Initial release
