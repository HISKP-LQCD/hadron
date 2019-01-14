context('parlist')

test_that('parlist_1', {
  corr_matrix_size <- 1
  actual <- make_parlist(corr_matrix_size)
  target <- array(c(1, 1), dim = c(2, 1))
  expect_equal(actual, target)
})

test_that('parlist_4', {
  corr_matrix_size <- 4
  actual <- make_parlist(corr_matrix_size)
  target <- array(c(1, 1, 1, 2, 2, 1, 2, 2), dim = c(2, 4))
  expect_equal(actual, target)
})

test_that('parind_1', {
  corr_matrix_size <- 1
  length_time <- 3
  parlist <- make_parlist(corr_matrix_size)
  actual <- make_parind(parlist, length_time)

  elements <- c(rep(c(2, 2), each = length_time))
  target <- array(elements, dim = c(length_time, 2))

  expect_equal(actual, target)
})

test_that('parind_4', {
  corr_matrix_size <- 4
  length_time <- 7
  parlist <- make_parlist(corr_matrix_size)
  actual <- make_parind(parlist, length_time, summands = 1)

  elements <- c(rep(c(2, 2, 3, 3, 2, 3, 2, 3), each = length_time))
  target <- array(elements, dim = c(corr_matrix_size * length_time, 2))
  expect_equal(actual, target)

  actual <- make_parind(parlist, length_time, summands = 2)
  target <- cbind(target, target + 2)
  expect_equal(actual, target)
})
