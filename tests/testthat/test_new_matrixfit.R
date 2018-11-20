context('new_matrixfit')

samplecf_boot <- bootstrap.cf(samplecf, 100)

test_that('SimpleModel', {
  ## We perform the fit using the old `matrixfit` function.
  fit_old <- matrixfit(samplecf_boot, 10, 15, fit.method = 'lm')
  
  fit_new <- new_matrixfit(samplecf_boot, 10, 15, fit.method = 'lm')
  
  expect_equal(fit_old$t0, fit_new$t0, tolerance = 1e-4)
  expect_equal(fit_old$t[, 1:2], fit_new$t, tolerance = 1e-4)
  expect_equal(fit_old$dof, fit_new$dof)
  expect_equal(fit_old$chisqr, fit_new$chisqr)
  expect_equal(fit_old$Qval, fit_new$Qval, tolerance = 1e-4)
  
})