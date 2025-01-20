
# Test for Rfast::Dist() with "cosine" -----------------------------------------


test_that("Rfast::Dist() with method cosine still returns similarities", {
    mat_test <- matrix(c(1, 2, 3, 1, 2, 3, 3, 2, 1), nrow = 3, byrow = TRUE)
    dist_test <- Rfast::Dist(mat_test, "cosine")
    expect_identical(dist_test[2, 1], 1)
})


# Tests for clustering ---------------------------------------------------------

test_that("ClusterMuscadet() returns an updated muscadet object", {
    data(muscadet_obj)
    # remove clustering result
    muscadet_obj@clustering <- list()

    obj <- clusterMuscadet(
        muscadet_obj,
        dist_method = "euclidean",
        hclust_method = "ward.D",
        k_range = 2:5
    )
    expect_s4_class(obj, "muscadet")
    expect_true(length(obj$clustering) != 0)
})

test_that("imputeClusters() returns a list of k elements as named vectors", {
    # Create matrices with some cells missing in one or the other
    mat1 <- matrix(runif(100), nrow = 20)
    mat2 <- matrix(runif(100), nrow = 20)
    rownames(mat1) <- paste0("Cell", 1:20)
    rownames(mat2) <- paste0("Cell", c(1:5, 11:25))
    mat_list <- list(ATAC = mat1, RNA = mat2)

    # Create cluster assignments for common cells
    clusters <- list(
        `2` = setNames(
            sample(1:2, 15, replace = TRUE),
            intersect(rownames(mat1), rownames(mat2))
        ),
        `3` = setNames(
            sample(1:3, 15, replace = TRUE),
            intersect(rownames(mat1), rownames(mat2))
        ),
        `4` = setNames(
            sample(1:4, 15, replace = TRUE),
            intersect(rownames(mat1), rownames(mat2))
        )
    )

    # Impute cluster assignments for missing cells
    imputed_clusters <- imputeClusters(mat_list, clusters, knn_imp = 5)

    expect_type(imputed_clusters, "list")
    expect_length(imputed_clusters, 3)

    expect_length(imputed_clusters[[1]], 25)
    expect_length(imputed_clusters[[2]], 25)
    expect_length(imputed_clusters[[3]], 25)

    expect_type(imputed_clusters[[1]], "integer")
    expect_type(imputed_clusters[[2]], "integer")
    expect_type(imputed_clusters[[3]], "integer")

    expect_identical(sort(names(imputed_clusters[[1]])), sort(union(rownames(mat1), rownames(mat2))))

})





