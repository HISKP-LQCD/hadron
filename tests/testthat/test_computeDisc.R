context('computeDisc')

test_that('cross_vs_diagonal', {
  data(loopdata)
  X1 <- computeDisc(cf=loopdata, real=TRUE, subtract.vev=TRUE)
  X2 <- computeDisc(cf=loopdata, cf2=loopdata, real=TRUE, real2=TRUE, subtract.vev=TRUE, subtract.vev2=TRUE)
  expect_equal(sum(X1$cf-X2$cf), 0)

  X1 <- computeDisc(cf=loopdata, real=TRUE, subtract.vev=FALSE)
  X2 <- computeDisc(cf=loopdata, cf2=loopdata, real=TRUE, real2=TRUE, subtract.vev=FALSE, subtract.vev2=FALSE)
  expect_equal(sum(X1$cf-X2$cf), 0)
}
)
