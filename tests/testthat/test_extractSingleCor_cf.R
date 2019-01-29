context('extractSingleCor.cf')

test_that('extract_symmetrized', {
  ## We generate some different correlators from the data that we already have.
  samplecf_boot <- bootstrap.cf(samplecf)
  cf1 <- mul.cf(samplecf_boot, 5.39)
  cf2 <- mul.cf(samplecf_boot, 1.48)
  cf3 <- mul.cf(samplecf_boot, 3.09)

  ## We glue them together and extract the second one again.
  corr <- c(cf1, cf2, cf3)
  ex2 <- extractSingleCor.cf(corr, 2)

  expect_equal(ex2$cf, cf2$cf)
  expect_equal(ex2$cf0, cf2$cf0)
  expect_equal(ex2$cf.tsboot$t0, cf2$cf.tsboot$t0)
  expect_equal(ex2$cf.tsboot$t, cf2$cf.tsboot$t)
})

test_that('extract_symmetrized', {
  ## We generate some different correlators from the data that we already have.
  ## `samplecf` is already symmetrized, so we need to let it forget about that.
  unsym <- samplecf
  unsym$Time <- ncol(unsym$cf)
  unsym$symmetrised <- FALSE
  samplecf_boot <- bootstrap.cf(unsym)
  cf1 <- mul.cf(samplecf_boot, 5.39)
  cf2 <- mul.cf(samplecf_boot, 1.48)
  cf3 <- mul.cf(samplecf_boot, 3.09)

  ## We glue them together and extract the second one again.
  corr <- c(cf1, cf2, cf3)
  ex2 <- extractSingleCor.cf(corr, 2)

  expect_equal(ex2$cf, cf2$cf)
  expect_equal(ex2$cf0, cf2$cf0)
  expect_equal(ex2$cf.tsboot$t0, cf2$cf.tsboot$t0)
  expect_equal(ex2$cf.tsboot$t, cf2$cf.tsboot$t)
})
