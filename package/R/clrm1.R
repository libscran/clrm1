#' @export
#' @importFrom Matrix rowSums colMeans
clrm1 <- function(mat) {
    keep <- rowSums(mat) > 0
    lmat <- log1p(mat[keep,,drop=FALSE])
    expm1(colMeans(lmat))
}

#' @export
#' @importFrom DelayedArray DelayedArray
#' @importFrom Matrix rowSums colMeans
clrm1.delayed <- function(mat) {
    keep <- rowSums(mat) > 0
    da <- DelayedArray(mat)
    lmat <- log1p(da[keep,,drop=FALSE])
    expm1(colMeans(lmat))
}

#' @export
#' @importFrom beachmat initializeCpp
#' @useDynLib clrm1
clrm1.cpp <- function(mat) {
    clrm1_cpp(initializeCpp(mat))
}
