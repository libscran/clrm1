# library(testthat); library(clrm1); source("test-clrm1.R")

set.seed(1000)
y <- abs(Matrix::rsparsematrix(100, 1000, 0.1))

test_that("basic checks work out", {
    ref <- clrm1(y)
    expect_equal(ref, clrm1.delayed(y))
    expect_equal(ref, clrm1.cpp(y))
})

test_that("methods are the same after filtering", {
    ref <- clrm1(y)

    y2 <- rbind(y, matrix(0, 10, ncol(y)))
    ref2 <- clrm1(y2)
    expect_identical(ref2, ref)

    expect_equal(ref2, clrm1.delayed(y2))
    expect_equal(ref2, clrm1.cpp(y2))
})
