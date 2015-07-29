library(nrmisc)
context("Interleave function")

test_that("Output matches expected value", {
    expect_equal(interleave(rep(TRUE, 3), rep(FALSE, 3)), c(T, F, T, F, T, F))
    expect_equal(interleave(c(1,3), c(2,4)), 1:4)
    expect_equal(interleave(letters[1:2], letters[3:4]), c("a", "c", "b", "d"))
})

test_that("Erroneous input triggers error", {
    expect_error(interleave(data.frame(1:2), data.frame(3:4)), "must be a vector")
    expect_error(interleave(1:2, c("a", "b")), "must have the same mode")
    expect_error(interleave(1:2, 1), "must have the same length")
})


if (requireNamespace("lintr", quietly = TRUE)) {
    context("lints")
    test_that("Package style conforms to linters", {
        lintr::expect_lint_free()
    })
}
