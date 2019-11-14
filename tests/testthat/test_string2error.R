context('string2error')

test_that('1', {
  expect_equal(string2error('3.10(2)'), c(3.10, 0.02))
  expect_equal(string2error('-3.10(2)'), c(-3.10, 0.02))
  
  expect_equal(string2error('310(2)'), c(310, 2))
  expect_equal(string2error('-310(2)'), c(-310, 2))
  
  expect_equal(string2error('0.001(2)'), c(0.001, 0.002))
  expect_equal(string2error('0.001(20)'), c(0.001, 0.020))
  expect_equal(string2error('0.001(200)'), c(0.001, 0.200))
  expect_equal(string2error('0.001(2000)'), c(0.001, 2.000))
  expect_equal(string2error('0.001(20000)'), c(0.001, 20.000))
})