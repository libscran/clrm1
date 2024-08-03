#include "clrm1.hpp"
#include "Rcpp.h"
#include "Rtatami.h"

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector clrm1_cpp(SEXP raw_ptr) {
    Rtatami::BoundNumericPointer ptr(raw_ptr);
    const auto& mat = *(ptr->ptr);
    Rcpp::NumericVector output(mat.ncol());
    clrm1::compute(*(ptr->ptr), clrm1::Options(), static_cast<double*>(output.begin()));
    return output;
}
