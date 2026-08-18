#include "clrm1/clrm1.hpp"
#include "Rcpp.h"
#include "tatami/tatami.hpp"

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector clrm1_cpp(Rcpp::NumericMatrix mat) {
    const std::size_t NR = mat.nrow(), NC = mat.ncol();
    tatami::DenseMatrix<double, int, tatami::ArrayView<double> > wrapped(NR, NC, tatami::ArrayView<double>(mat.begin(), NR * NC), false);
    Rcpp::NumericVector output(NC);
    clrm1::compute(wrapped, clrm1::Options(), static_cast<double*>(output.begin()));
    return output;
}
