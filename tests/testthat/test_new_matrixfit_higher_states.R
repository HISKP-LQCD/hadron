context('new_matrixfit higher states')

test_that('two states', {
  corr_boot <- bootstrap.cf(samplecf, boot.R = 300)
  fit <- new_matrixfit(
    cf = corr_boot,
    t1 = 10,
    t2 = 20,
    useCov = FALSE,
    model = 'n_particles',
    higher_states = list(val = 0.1, boot = matrix(rnorm(300, 0.1, 0.01), ncol = 1), ampl = 1.0))
  
  expect_true(TRUE)
})

test_that('two states with cov', {
  corr_boot <- bootstrap.cf(samplecf, boot.R = 300)
  fit <- new_matrixfit(
    cf = corr_boot,
    t1 = 10,
    t2 = 20,
    useCov = TRUE,
    model = 'n_particles',
    higher_states = list(val = 0.1, boot = matrix(rnorm(300, 0.1, 0.01), ncol = 1), ampl = 1.0))
  
  expect_true(TRUE)
})