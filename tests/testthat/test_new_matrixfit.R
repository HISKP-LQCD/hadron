context('new_matrixfit')

test_that('SingleModelJacobian', {
  time_extent <- 10
  m_size <- 1
  parlist <- make_parlist(m_size)
  par <- c(0.4, 2.3, 3.4)
  x <- 1:time_extent
  sym_vec <- 'cosh'
  neg_vec <- c(1)
  
  L <- diag(rep(1, time_extent))
  
  new_model <- SingleModel$new(time_extent, parlist, sym_vec, neg_vec, m_size)
  
  old_jac_vec <- dmatrixChi(
    par, x, 0, L, time_extent,
    parind = make_parind(parlist, time_extent),
    sign.vec = make_sign_vec(sym_vec, time_extent),
    ov.sign.vec = make_ov_sign_vec(neg_vec, time_extent))
  str(old_jac_vec)
  old_jac <- - matrix(old_jac_vec, ncol = length(par))
  new_jac <- new_model$prediction_jacobian(par, x)
  
  expect_equal(sum(is.na(old_jac)), 0)
  expect_equal(sum(is.na(new_jac)), 0)
  
  expect_equal(new_jac, old_jac)
})

test_that('SingleModel', {
  samplecf_boot <- bootstrap.cf(samplecf, 100)
  
  args <- list(samplecf_boot, 5, 10, fit.method = 'lm', model = 'single')
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

test_that('ShiftedModelPrediction', {
  time_extent <- 10
  delta_t <- 1
  m_size <- 1
  parlist <- make_parlist(m_size)
  par <- c(0.4, 2.3)
  x <- 1:time_extent
  sym_vec <- 'cosh'
  neg_vec <- c(1)
  
  L <- diag(rep(1, time_extent))
  
  new_model <- ShiftedModel$new(time_extent, parlist, sym_vec, neg_vec, m_size, delta_t)
  
  old_pred <- - c(matrixChi.shifted(par, x, 0, L, time_extent,
                                    parind = make_parind(parlist, time_extent),
                                    sign.vec = make_sign_vec(sym_vec, time_extent),
                                    ov.sign.vec = make_ov_sign_vec(neg_vec, time_extent),
                                    deltat = delta_t))
  new_pred <- new_model$prediction(par, x)
  expect_equal(new_pred, old_pred)
  
  old_jac <- - matrix(dmatrixChi.shifted(par, x, 0, L, time_extent, 
                                         parind = make_parind(parlist, time_extent),
                                         sign.vec = make_sign_vec(sym_vec, time_extent),
                                         ov.sign.vec = make_ov_sign_vec(neg_vec, time_extent),
                                         deltat = delta_t),
                      ncol = length(par))
  new_jac <- new_model$prediction_jacobian(par, x)
  expect_equal(old_jac, new_jac)
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
  
  expect_equal(fit_old$t0, fit_new$t0, tolerance = 1e-5)
  
  expect_gte(fit_new$t0[1], 0)
  expect_gte(fit_new$t0[2], 0)
  
  ## The `matrixfit` function gives the chi² as part of the bootstrap samples,
  ## we do not want this here. Therefore we need to cut.
  expect_equal(fit_old$t[, 1:3], fit_new$t)
  expect_equal(fit_old$dof, fit_new$dof)
  expect_equal(fit_old$chisqr, fit_new$chisqr)
  expect_equal(fit_old$Qval, fit_new$Qval)
})
