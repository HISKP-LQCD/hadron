context('bootstrap.nlsfit')

model <- function (par, x, ...) {
  par[1] + x * par[2]
}

bsamples_x <- parametric.bootstrap(300, 1:10, rep(0.1, 10))
bsamples_y <- parametric.bootstrap(300, 1:10, rep(0.1, 10))
bsamples_yx <- cbind(bsamples_y, bsamples_x)
cov_y <- cov(bsamples_y)
cov_yx <- cov(bsamples_yx)
x <- apply(bsamples_x, 2, mean)
y <- apply(bsamples_y, 2, mean)
mask <- 3:8

psamples <- rnorm(300, 0.1, 0.1)

test_that('y errors', {
  fit <- bootstrap.nlsfit(
    fn = model, 
    par.guess = c(0, 1),
    x = x,
    y = y,
    bsamples = bsamples_y,
    mask = mask)
  skip_on_os(os="windows")
  expect_true(TRUE)
})

test_that('y errors cov', {
  fit <- bootstrap.nlsfit(
    fn = model, 
    par.guess = c(0, 1),
    x = x,
    y = y,
    bsamples = bsamples_y,
    CovMatrix = cov_y,
    mask = mask)
  skip_on_os(os="windows")
  expect_true(TRUE)
})

test_that('xy errors', {
  fit <- bootstrap.nlsfit(
    fn = model, 
    par.guess = c(0, 1),
    x = x,
    y = y,
    bsamples = bsamples_yx,
    mask = mask)
  skip_on_os(os="windows")
  expect_true(TRUE)
})

test_that('xy errors cov', {
  fit <- bootstrap.nlsfit(
    fn = model, 
    par.guess = c(0, 1),
    x = x,
    y = y,
    bsamples = bsamples_yx,
    CovMatrix = cov_yx,
    mask = mask)
  skip_on_os(os="windows")
  expect_true(TRUE)
})


test_that('y errors with priors', {
  fit <- bootstrap.nlsfit(
    fn = model, 
    par.guess = c(0, 1),
    x = x,
    y = y,
    bsamples = bsamples_y,
    priors = list(param = 1, p = 0.1, psamples = psamples),
    mask = mask)
  
  expect_true(TRUE)
})

test_that('y errors cov with priors', {
  fit <- bootstrap.nlsfit(
    fn = model, 
    par.guess = c(0, 1),
    x = x,
    y = y,
    bsamples = bsamples_y,
    CovMatrix = cov_y,
    priors = list(param = 1, p = 0.1, psamples = psamples),
    mask = mask)
  
  expect_true(TRUE)
})

test_that('xy errors with priors', {
  fit <- bootstrap.nlsfit(
    fn = model, 
    par.guess = c(0, 1),
    x = x,
    y = y,
    bsamples = bsamples_yx,
    priors = list(param = 1, p = 0.1, psamples = psamples),
    mask = mask)
  
  expect_true(TRUE)
})

test_that('xy errors cov with priors', {
  fit <- bootstrap.nlsfit(
    fn = model, 
    par.guess = c(0, 1),
    x = x,
    y = y,
    bsamples = bsamples_yx,
    CovMatrix = cov_yx,
    priors = list(param = 1, p = 0.1, psamples = psamples),
    mask = mask)
  skip_on_os(os="windows")  
  expect_true(TRUE)
})
