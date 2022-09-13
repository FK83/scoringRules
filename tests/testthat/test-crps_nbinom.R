context("CRPS for negative binomial distribution")

if (requireNamespace("hypergeo", quietly = TRUE)) {
  test_that("computed values are correct", {
    expect_equal(crps_nbinom(         5,      100, prob = 0.95),   0.53215257681)
    expect_equal(crps_nbinom(1999797667, 200000.1, prob = 0.0001), 1045062.81311)
    expect_equal(crps_nbinom(        44,    200.1, prob = 0.82),   1.70859687251)
    expect_equal(crps_nbinom(        41,    200.1, prob = 0.83),   1.63833674721)
    expect_equal(crps_nbinom(       222,   2000.1, prob = 0.9),    3.66926727576)
    expect_equal(crps_nbinom(      2222,  20000.1, prob = 0.9),    11.6114812752)
  })
}