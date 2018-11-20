context('new_matrixfit')

samplecf_boot <- bootstrap.cf(samplecf, 100)

test_that('SimpleModel', {
  ## We perform the fit using the old `matrixfit` function.
  fit_old <- matrixfit(samplecf_boot, 10, 15)
  
  fit_new <- new_matrixfit(samplecf_boot, 10, 15)
})