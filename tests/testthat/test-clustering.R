
# Test for Rfast::Dist() with "cosine" -----------------------------------------


test_that("Rfast::Dist with method cosine still returns similarities", {
    mat_test <- matrix(c(1, 2, 3, 1, 2, 3, 3, 2, 1), nrow = 3, byrow = TRUE)
    dist_test <- Rfast::Dist(mat_test, "cosine")
    expect_identical(dist_test[2, 1], 1)
})
