context('new_matrixfit')

test_that('MatrixModel', {
  samplecf_boot <- bootstrap.cf(samplecf, 100)
  
  args <- list(samplecf_boot, 5, 10, fit.method = 'lm')
  fit_old <- do.call(matrixfit, args)
  fit_new <- do.call(new_matrixfit, args)
  
  expect_equal(fit_old$t0, fit_new$t0, tolerance = 1e-4)
  
  ## The `matrixfit` function gives the chi² as part of the bootstrap samples,
  ## we do not want this here. Therefore we need to cut.
  expect_equal(fit_old$t[, 1:2], fit_new$t, tolerance = 1e-4)
  expect_equal(fit_old$dof, fit_new$dof)
  expect_equal(fit_old$chisqr, fit_new$chisqr)
  expect_equal(fit_old$Qval, fit_new$Qval)
})

test_that('ShiftedModel', {
  samplecf_shited <- takeTimeDiff.cf(samplecf, 1)
  samplecf_boot <- bootstrap.cf(samplecf_shited, 100)
  
  args <- list(samplecf_boot, 5, 10, model = 'shifted', fit.method = 'lm', sym.vec = 'sinh')
  fit_old <- do.call(matrixfit, args)
  fit_new <- do.call(new_matrixfit, args)
  
  expect_equal(fit_old$t0, fit_new$t0, tolerance = 1e-4)
  
  ## The `matrixfit` function gives the chi² as part of the bootstrap samples,
  ## we do not want this here. Therefore we need to cut.
  expect_equal(fit_old$t[, 1:2], fit_new$t, tolerance = 1e-4)
  expect_equal(fit_old$dof, fit_new$dof)
  expect_equal(fit_old$chisqr, fit_new$chisqr)
  expect_equal(fit_old$Qval, fit_new$Qval)
})

test_that('TwoStateModel', {
  samplecf_boot <- bootstrap.cf(correlatormatrix, 100)
  gevp <- bootstrap.gevp(samplecf_boot, t0 = 4, element.order = 1:4)
  pc <- gevp2cf(gevp, 1)
  
  args <- list(pc, 2, 15, model = 'pc', fit.method = 'lm')
  fit_old <- do.call(matrixfit, args)
  fit_new <- do.call(new_matrixfit, args)
  
  expect_equal(fit_old$t0, fit_new$t0, tolerance = 1e-4)
  
  expect_gte(fit_new$t0[1], 0)
  expect_gte(fit_new$t0[2], 0)
  
  ## The `matrixfit` function gives the chi² as part of the bootstrap samples,
  ## we do not want this here. Therefore we need to cut.
  expect_equal(fit_old$t[, 1:3], fit_new$t)
  expect_equal(fit_old$dof, fit_new$dof)
  expect_equal(fit_old$chisqr, fit_new$chisqr)
  expect_equal(fit_old$Qval, fit_new$Qval)
})