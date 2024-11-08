

atac <- CreateMuscomicObject(
    type = "ATAC",
    mat_counts = mat_counts_atac_tumor,
    allele_counts = allele_counts_atac_tumor,
    features = peaks
)
rna <- CreateMuscomicObject(
    type = "RNA",
    mat_counts = mat_counts_rna_tumor,
    allele_counts = allele_counts_rna_tumor,
    features = genes
)

atac_ref <- CreateMuscomicObject(
    type = "ATAC",
    mat_counts = mat_counts_atac_ref,
    allele_counts = allele_counts_atac_ref,
    features = peaks
)
rna_ref <- CreateMuscomicObject(
    type = "RNA",
    mat_counts = mat_counts_rna_ref,
    allele_counts = allele_counts_rna_ref,
    features = genes
)

muscadet <- CreateMuscadetObject(
    omics = list(atac, rna),
    bulk_lrr = bulk_lrr,
    bulk.label = "WGS",
    genome = "hg38"
)
muscadet_ref <- CreateMuscadetObject(
    omics = list(atac_ref, rna_ref),
    genome = "hg38"
)



test_that("ComputeLogRatioATAC gives a correct object as output", {

    obj_atac <- computeLogRatioATAC(
        matTumor = matCounts(muscadet)$ATAC,
        matRef = matCounts(muscadet_ref)$ATAC,
        peaksCoord = coordFeatures(muscadet)$ATAC,
        genome = slot(muscadet, "genome"),
        minReads = 1,
        minPeaks = 10,
        quiet = TRUE
    )

    expect_identical(names(obj_atac), c("matTumor","matRef","params","coord"))
})

test_that("ComputeLogRatioRNA gives a correct object as output", {

    obj_rna <- computeLogRatioRNA(
        matTumor = matCounts(muscadet)$RNA,
        matRef = matCounts(muscadet_ref)$RNA,
        genesCoord = coordFeatures(muscadet)$RNA,
        genome = slot(muscadet, "genome"),
        refReads = 20,
        quiet = TRUE
    )

    expect_identical(names(obj_rna), c("matTumor","matRef","params","coord"))
})

test_that("ComputeLogRatioRNA gives a muscadet object as output", {

    obj_atac_all <- computeLogRatioATAC(
        matTumor = matCounts(muscadet)$ATAC,
        matRef = matCounts(muscadet_ref)$ATAC,
        peaksCoord = coordFeatures(muscadet)$ATAC,
        genome = slot(muscadet, "genome"),
        minReads = 1,
        minPeaks = 10,
        all_steps = TRUE,
        quiet = TRUE
    )
    expect_identical(names(obj_atac_all), c("step01","step02","step03","step04","step05","step07","step08","params","coord"))
    expect_identical(names(obj_atac_all$step08), c("matTumor","matRef","name"))
})

test_that("ComputeLogRatioRNA gives a muscadet object as output", {

    obj_rna_all <- computeLogRatioRNA(
        matTumor = matCounts(muscadet)$RNA,
        matRef = matCounts(muscadet_ref)$RNA,
        genesCoord = coordFeatures(muscadet)$RNA,
        genome = slot(muscadet, "genome"),
        refReads = 20,
        all_steps = TRUE,
        quiet = TRUE
    )
    expect_identical(names(obj_rna_all), c("step01","step02","step03","step04","step05","step06","step07","step08","params","coord"))
    expect_identical(names(obj_rna_all$step08), c("matTumor","matRef","name"))
})

test_that("ComputeLogRatio gives a correct muscadet object as output", {
    muscadet <- computeLogRatio(
        x = muscadet,
        reference = muscadet_ref,
        omic = "ATAC",
        method = "ATAC",
        minReads = 1,
        minPeaks = 10,
        quiet = TRUE
    )
    muscadet <- computeLogRatio(
        x = muscadet,
        reference = muscadet_ref,
        omic = "RNA",
        method = "RNA",
        refReads = 20,
        quiet = TRUE
    )
    expect_identical(as.character(class(muscadet)), "muscadet")
    expect_true(is.matrix(muscadet@omics[["ATAC"]]@coverage[["log.ratio"]]))
    expect_true(is.matrix(muscadet@omics[["RNA"]]@coverage[["log.ratio"]]))
})


test_that(
    "ComputeLogRatio run 2 times on the same omic with removal of raw counts after the 1st: error",
    {
        muscadet <- computeLogRatio(
            x = muscadet,
            reference = muscadet_ref,
            omic = "ATAC",
            method = "ATAC",
            minReads = 1,
            minPeaks = 10,
            quiet = TRUE
        )
        expect_error(
            muscadet <- computeLogRatio(
                x = muscadet,
                reference = muscadet_ref,
                omic = "ATAC",
                method = "ATAC",
                minReads = 1,
                minPeaks = 10,
                quiet = TRUE
            )
        )
    }
)




