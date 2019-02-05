context('Remove Temporal')

test_that('equality', {
  ## Perform some fits to the sample cf such that we have energies to subtract.
  corr_boot <- bootstrap.cf(samplecf)

  fit1 <- matrixfit(corr_boot, 10, 20)
  fit2 <- matrixfit(corr_boot, 4, 9)

  corr_rt <- old_removeTemporal.cf(corr_boot, fit1, fit2)
  new_corr_rt <- removeTemporal.cf(corr_boot, fit1, fit2)

  expect_equal(corr_rt, new_corr_rt)
})
