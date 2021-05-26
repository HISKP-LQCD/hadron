context('bootstrap_rw.cf')

test_that('constant_rewfactor', {

  class(samplecf) <- append( class(samplecf), 'cf_indexed')
  samplecf$conf.index <- seq(0,1017)

  ## We multiply the reweighting factor with its inverse
  tmp <- multiply.rw(samplerw,samplerw_inverse)

  ## Perform reweighting with constant reweighting factors
  tmp1 <- bootstrap_rw.cf(cf=samplecf,rw=tmp,boot.R=1500)

  ## Perform bootstrapping withput any reweighting at all, with the same seed
  tmp2 <- bootstrap.cf(cf=samplecf,boot.R=1500)

  expect_equal( tmp1$cf0, tmp2$cf0 )
  expect_equal( tmp1$tsboot.se, tmp2$tsboot.se )
})
