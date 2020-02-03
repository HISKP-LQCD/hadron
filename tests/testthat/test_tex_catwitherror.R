context('tex_catwitherror')

test_that('small_error', {
    expect_equal(tex.catwitherror(1.23, 0.45, digits = 2, with.dollar = FALSE), '1.23(45)')
    expect_equal(tex.catwitherror(1.23, 0.450001, digits = 1, with.dollar = FALSE), '1.2(5)')
    expect_equal(tex.catwitherror(1.23, 0.001, digits = 1, with.dollar = FALSE), '1.230(1)')
    expect_equal(tex.catwitherror(123, 1, digits = 1, with.dollar = FALSE), '123(1)')
})

test_that('very_small_number', {
    expect_equal(tex.catwitherror(1.23e-20, 0.45e-20, digits = 2, with.dollar = FALSE), '1.23(45)e-20')
    expect_equal(tex.catwitherror(1.23e-40, 0.45e-40, digits = 2, with.dollar = TRUE), '$1.23(45)\\cdot 10^{-40}$')

    expect_equal(tex.catwitherror(1.23e-40, digits = 2, with.dollar = FALSE), '1.2e-40')
})

test_that('same_error', {
    expect_equal(tex.catwitherror(0.00123, 0.00123, digits = 4, with.dollar = FALSE), '0.001230(1230)')
    expect_equal(tex.catwitherror(12.345, 12.345, digits = 4, with.dollar = FALSE), '12.35(12.35)')
})

test_that('intermediate_error', {
    expect_equal(tex.catwitherror(175.2, 23.3, digits = 3, with.dollar = FALSE), '175.2(23.3)')
})

test_that('large_error', {
    expect_equal(tex.catwitherror(0.00123, 0.45, digits = 4, with.dollar = FALSE), '0.0012(0.45)')
    expect_equal(tex.catwitherror(1.12345, 12.3, digits = 4, with.dollar = FALSE), '1.12(12.3)')
})

test_that('similar_error', {
    expect_equal(tex.catwitherror(7492.8291130334482659, 1759.0859320695926726,
                                  digits = 3, with.dollar = FALSE), '7492.829(090)')
})

test_that('no_error', {
    expect_equal(tex.catwitherror(0.00123, digits = 4, with.dollar = FALSE), '0.001230')
})

test_that('zero_error', {
    expect_equal(tex.catwitherror(0.00123, 0, digits = 4, with.dollar = FALSE), '0.001230(0)')
    expect_equal(tex.catwitherror(12.345, 0, digits = 3, with.dollar = FALSE), '12.3(0)')
})

test_that('zero_val_zero_err', {
    expect_equal(tex.catwitherror(0., 0., digits = 4, with.dollar = FALSE), '0.000(0)')
})

test_that('NA', {
    expect_equal(tex.catwitherror(NA, NA, digits = 4, with.dollar = FALSE), 'NA(NA)')
})


test_that('vector', {
    for (i in 1:10) {
        x <- runif(1)
        dx <- runif(1)
        expect_equal(tex.catwitherror(x, dx, digits = 4),
                     tex.catwitherror(c(x, dx), digits = 4))
    }
})
