context('tex_catwitherror')

test_that('small_error', {
    expect_equal(tex.catwitherror(1.23, 0.45, digits = 2, with.dollar = FALSE), '1.23(45)')
    expect_equal(tex.catwitherror(1.23, 0.450001, digits = 1, with.dollar = FALSE), '1.2(5)')
    expect_equal(tex.catwitherror(1.23, 0.001, digits = 1, with.dollar = FALSE), '1.230(1)')
    expect_equal(tex.catwitherror(123, 1, digits = 1, with.dollar = FALSE), '123(1)')
})

test_that('borderline_error', {
    expect_equal(tex.catwitherror(1.330563782105, 0.000966674080, digits = 1, with.dollar = FALSE), '1.331(1)')
})

test_that('another_borderline_error', {
    expect_equal(tex.catwitherror(0.031921636680, 0.000098, digits = 1, with.dollar=FALSE), '0.0319(1)')
})

test_that('even_nastier_borderline_error', {
    expect_equal(tex.catwitherror(1.330563782105, 0.996674080, digits = 2, with.dollar = FALSE), '1.3(10)')
})

test_that('scientific_notation', {
    expect_equal(tex.catwitherror(0.0008970, 0.0002106, with.dollar = FALSE, with.cdot = TRUE, digits = 2), '9.0(21)\\cdot 10^{-4}')
})

test_that('very_small_number', {
    expect_equal(tex.catwitherror(1.23e-20, 0.45e-20, digits = 2, with.dollar = FALSE, with.cdot = FALSE), '1.23(45)e-20')
    expect_equal(tex.catwitherror(1.23e-40, 0.45e-40, digits = 2), '$1.23(45)\\cdot 10^{-40}$')

    expect_equal(tex.catwitherror(1.23e-40, digits = 2, with.dollar = FALSE, with.cdot = FALSE), '1.2e-40')
})

test_that('same_error', {
    expect_equal(tex.catwitherror(0.00123, 0.00123, digits = 4, with.dollar = FALSE), '0.001230(1230)')
    expect_equal(tex.catwitherror(12.346, 12.346, digits = 4, with.dollar = FALSE), '12.35(1235)')
})

test_that('intermediate_error', {
    expect_equal(tex.catwitherror(175.2, 23.3, digits = 3, with.dollar = FALSE), '175.2(233)')
})

test_that('large_error', {
    expect_equal(tex.catwitherror(0.00123, 0.45, digits = 4, with.dollar = FALSE), '0.0012(4500)')
    expect_equal(tex.catwitherror(1.12345, 12.3, digits = 4, with.dollar = FALSE), '1.12(1230)')
})

test_that('similar_error', {
    expect_equal(tex.catwitherror(7492.8291130334482659, 1759.0859320695926726,
                                  digits = 3, with.dollar = FALSE), '7490(1760)')
})

test_that('no_error', {
    expect_equal(tex.catwitherror(0.00123, digits = 4, with.dollar = FALSE, flag = "#"), '0.001230')
})

test_that('zero_error', {
    expect_equal(tex.catwitherror(0.00000123, 0, digits = 4, with.dollar = FALSE, with.cdot = FALSE), '1.23(0)e-06')
    expect_equal(tex.catwitherror(0.00123, 0, digits = 4, with.dollar = FALSE), '0.00123(0)')
    expect_equal(tex.catwitherror(12.345, 0, digits = 3, with.dollar = FALSE), '12.3(0)')
})

test_that('zero_val_zero_err', {
    expect_equal(tex.catwitherror(0., 0., digits = 4, with.dollar = FALSE, flag="#"), '0.000(0)')
})

test_that('NA', {
    expect_equal(tex.catwitherror(NA, NA, digits = 4, with.dollar = FALSE), 'NA(  NA)')
})

test_that('vector', {
    for (i in 1:10) {
        x <- runif(1)
        dx <- runif(1)
        expect_equal(tex.catwitherror(x, dx, digits = 4),
                     tex.catwitherror(c(x, dx), digits = 4))
    }
})


test_that('vectorized_two_elements', {
    val <- runif(2)
    err <- runif(2)
    
    str_apply <- mapply(function (v, e) tex.catwitherror(v, e), val, err)
    str_vec <- tex.catwitherror(val, err)
    expect_equal(str_apply, str_vec)
})

test_that('vectorized_with_zero_error', {
    val <- runif(20)
    err <- c(runif(10), rep(0.0, 10))
    
    str_apply <- mapply(function (v, e) tex.catwitherror(v, e), val, err)
    str_vec <- tex.catwitherror(val, err)
    expect_equal(str_apply, str_vec)
})

test_that('vectorized_with_zero_error_at_front', {
    val <- runif(20)
    err <- c(rep(0.0, 10), runif(10))
    
    str_apply <- mapply(function (v, e) tex.catwitherror(v, e), val, err)
    str_vec <- tex.catwitherror(val, err)
    expect_equal(str_apply, str_vec)
})