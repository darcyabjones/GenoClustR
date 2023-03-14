test_that("Correlation matrix same", {
  load(testthat::test_path("data", "sclero_test.rda"))
  load(testthat::test_path("data", "sclero_test_corr.rda"))

  gene_names <- rownames(sclero_test)
  this_corr <- gene_corr(sclero_test, gene_names)
  expect_equal(sclero_test_corr, this_corr)
})


test_that("Ave correlation matrix same", {
  load(testthat::test_path("data", "sclero_test_corr.rda"))
  load(testthat::test_path("data", "sclero_test_ave_corr.rda"))

  this_ave_corr <- average_corr(sclero_test_corr)
  expect_equal(sclero_test_ave_corr, this_ave_corr)
})
