# test_that("identity and log1p give same MIC on toy data", {
#   df  <- data.frame(score  = ordered(rep(0:4, 8)),
#                     conc   = rep(0:4, 8),
#                     strain = factor(rep(c("A","B"), 20)))
#   fit1 <- ordinal::clm(score ~ strain * conc, data = df, method="optim")
#   fit2 <- ordinal::clm(score ~ strain * log1p(conc), data = df, method="optim")
#
#   res1 <- mic_solve(fit1,
#                     expand.grid(strain = levels(df$strain)),
#                     conc_name = "conc",
#                     transform_fun     = identity,
#                     inv_transform_fun = identity)
#
#   res2 <- mic_solve(fit2,
#                     expand.grid(strain = levels(df$strain)),
#                     conc_name = "conc")
#
#   expect_equal(res1$mic_estimates$MIC,
#                res2$mic_estimates$MIC,
#                tolerance = 1e-6)
# })



test_that("autoplot returns a ggplot", {
  skip_if_not_installed("ggplot2")
  df  <- data.frame(score = ordered(sample(0:4, 40, TRUE)),
                    conc  = runif(40, 0.1, 4),
                    strain = factor(sample(c("A","B"), 40, TRUE)))
  fit <- ordinal::clm(score ~ strain * conc, data = df)
  res <- mic_solve(fit, expand.grid(strain = levels(df$strain)),
                   conc_name = "conc", transform_fun = identity,
                   inv_transform_fun = identity)
  p   <- ggplot2::autoplot(res, type = "mic")
  expect_s3_class(p, "ggplot")
})
