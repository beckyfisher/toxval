# The hand-built curves below are chosen so each crossing is known exactly.
# `prob_vals = c(0.5, 0, 1)` reports the median and the extremes of two
# realisations, so every expected value can be read off the curves.
all_quantiles <- c(0.5, 0, 1)

# columns --------------------------------------------------------------------

test_that("summarise_metric returns a toxval tibble with the ecx columns", {
  x_vec <- 0:10
  pred <- compute_pred(list(rbind(10 - x_vec, 10 - 2 * x_vec)), x_vec)

  result <- summarise_metric(compute_ecx(pred), pred, "ecx")

  expect_s3_class(result, "toxval")
  expect_named(
    result,
    c(
      "metric",
      "direction",
      "estimate",
      "conf.low",
      "conf.high",
      "proportion",
      "type",
      "control",
      "reference",
      "asymptote"
    )
  )
  expect_identical(result$metric, c("ecx", "ecx"))
  expect_identical(result$direction, c("decreasing", "increasing"))
})

test_that("summarise_metric returns the nsec columns, with the ecnsec summary", {
  x_vec <- 0:10
  pred <- compute_pred(list(rbind(10 - x_vec, 10 - 2 * x_vec)), x_vec)

  result <- summarise_metric(compute_nsec(pred), pred, "nsec")

  expect_named(
    result,
    c(
      "metric",
      "direction",
      "estimate",
      "conf.low",
      "conf.high",
      "sig_val",
      "type",
      "anchor",
      "control",
      "reference",
      "asymptote",
      "ecnsec",
      "ecnsec.low",
      "ecnsec.high"
    )
  )
  expect_identical(result$metric, c("nsec", "nsec"))
})

test_that("summarise_metric puts the key columns after direction", {
  x_vec <- 0:10
  curve <- rbind(10 - x_vec, 10 - 2 * x_vec)

  pred <- compute_pred(
    list(curve, curve, curve, curve),
    x_vec,
    keys = data.frame(
      group = factor(c("A", "A", "B", "B")),
      response = factor(c("growth", "survival", "growth", "survival"))
    ),
    group_var = "site",
    multi_var = c("growth", "survival")
  )

  result <- summarise_metric(compute_ecx(pred), pred, "ecx")

  expect_identical(
    names(result)[1:4],
    c("metric", "direction", "group", "response")
  )
  expect_identical(result$group, factor(rep(c("A", "B"), each = 4)))
  expect_identical(
    result$response,
    factor(rep(rep(c("growth", "survival"), each = 2), times = 2))
  )
})

# the summary ----------------------------------------------------------------

test_that("summarise_metric takes quantiles of the crossings", {
  x_vec <- 0:10
  # EC10 references 9 for both; crossings at 1 and 0.5
  pred <- compute_pred(list(rbind(10 - x_vec, 10 - 2 * x_vec)), x_vec)

  result <- summarise_metric(
    compute_ecx(pred),
    pred,
    "ecx",
    prob_vals = all_quantiles
  )

  expect_equal(
    unlist(result[1, c("estimate", "conf.low", "conf.high")]),
    c(estimate = 0.75, conf.low = 0.5, conf.high = 1)
  )
})

test_that("summarise_metric reports the central control and reference", {
  x_vec <- 0:10
  pred <- compute_pred(list(rbind(10 - x_vec, 8 - x_vec)), x_vec)

  result <- summarise_metric(
    compute_ecx(pred),
    pred,
    "ecx",
    prob_vals = all_quantiles
  )

  expect_equal(result$control, c(9, 9))
  expect_equal(result$reference[1], median(c(9, 7.2)))
})

test_that("summarise_metric returns NA where no realisation crosses", {
  x_vec <- 0:10
  pred <- compute_pred(list(rbind(10 - x_vec, 10 - 2 * x_vec)), x_vec)

  expect_no_warning(
    result <- summarise_metric(compute_ecx(pred), pred, "ecx")
  )

  expect_identical(
    unlist(result[2, c("estimate", "conf.low", "conf.high")]),
    c(estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_)
  )
})

test_that("summarise_metric returns NA where a realisation past the reference shows no effect and the rest are censored", {
  x_vec <- 0:10
  # The increasing NSEC reference is the 99% quantile of the controls, 9.98.
  # The first control, 10, already lies above it, but that realisation falls,
  # away from the effect, and below the reference, so it is censored rather
  # than placed at x_control (T21); the other two never rise and are censored.
  # None crosses (T10).
  pred <- compute_pred(
    list(rbind(10 - x_vec, 9 - x_vec, 8 - x_vec)),
    x_vec,
    x_control = -1
  )

  expect_no_warning(
    result <- summarise_metric(compute_nsec(pred), pred, "nsec", draws = TRUE)
  )

  expect_identical(result$draws[[2]]$value.draw, c(Inf, Inf, Inf))
  expect_identical(
    unlist(result[2, c("estimate", "conf.low", "conf.high")]),
    c(estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_)
  )
})

test_that("summarise_metric places a realisation already past the reference at x_control (T21)", {
  x_vec <- 1:10
  # `direct` reference 5: the first crosses at 5, the second starts at 3.9
  pred <- compute_pred(
    list(rbind(10 - x_vec, 4 - x_vec / 10)),
    x_vec,
    x_control = 0
  )

  result <- summarise_metric(
    compute_ecx(pred, reference = 5, type = "direct"),
    pred,
    "ecx",
    prob_vals = all_quantiles
  )

  expect_equal(
    unlist(result[1, c("estimate", "conf.low", "conf.high")]),
    c(estimate = 2.5, conf.low = 0, conf.high = 5)
  )
})

# One realisation per case of T21, each past the decreasing NSEC reference at
# its control. The 26 `t21_background` realisations, `6 - x`, set that
# reference: with 9 controls at 4 and 26 at 6, the 0.25 quantile is 5.
# `10 - curves` mirrors all 35 onto the increasing side, where the 0.75
# quantile is 5 and the expected values are the same. `x_control` is -1, so a
# realisation placed there is distinct from any crossing.
t21_curves <- rbind(
  crosses = c(4, 6, 6, 6, 4, 4, 4, 4, 4, 4, 4),
  dips_then_crosses = c(4, 3, 6, 6, 6, 4, 4, 4, 4, 4, 4),
  towards = c(4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6),
  flat_then_towards = c(4, 4, 4, 3, 2, 1, 0, -1, -2, -3, -4),
  round_off_away_then_towards = c(4, 4 + 1e-12, 3, 2, 1, 0, -1, -2, -3, -4, -5),
  hormetic = c(4, 4.5, 4.8, 4.5, 3, 2, 1, 0, -1, -2, -3),
  away_never_leaves = c(4, 4.2, 4.4, 4.6, 4.8, 4.9, 4.9, 4.9, 4.9, 4.9, 4.9),
  towards_then_leaves = c(4, 3, 2, 2, 3, 4, 6, 7, 8, 9, 10),
  away_leaves = c(4, 5.5, 6, 7, 8, 9, 10, 10, 10, 10, 10)
)
t21_background <- matrix(6 - 0:10, nrow = 26, ncol = 11, byrow = TRUE)
t21_expected <- c(3.5, 4.5, -1, -1, -1, -1, -1, -1, Inf)

t21_pred <- function(curves) {
  n <- nrow(curves)
  compute_pred(
    list(curves),
    0:10,
    x_control = -1,
    asymptote = list(list(decreasing = rep(-10, n), increasing = rep(20, n)))
  )
}

t21_nsec <- function(curves) {
  pred <- t21_pred(curves)
  warnings <- character()
  result <- withCallingHandlers(
    summarise_metric(
      compute_nsec(pred, sig_val = 0.25),
      pred,
      "nsec",
      prob_vals = all_quantiles,
      draws = TRUE
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(result = result, warnings = warnings)
}

test_that("summarise_metric censors a decreasing NSEC realisation past the reference only where it changes away and leaves (T21)", {
  output <- t21_nsec(rbind(t21_curves, t21_background))

  expect_equal(output$result$reference[[1]], 5)
  expect_equal(output$result$draws[[1]]$value.draw[1:9], t21_expected)
  expect_match(output$warnings, "decreasing: 2.9% censored", all = FALSE)
})

test_that("summarise_metric censors an increasing NSEC realisation past the reference only where it changes away and leaves (T21)", {
  output <- t21_nsec(10 - rbind(t21_curves, t21_background))

  expect_equal(output$result$reference[[2]], 5)
  expect_equal(output$result$draws[[2]]$value.draw[1:9], t21_expected)
  expect_match(output$warnings, "increasing: 2.9% censored", all = FALSE)
})

test_that("summarise_metric leaves ecx under direct unchanged on the T21 curves (#74)", {
  # rule 2 of T21 applies to an NSEC only. A realisation past the supplied
  # level at the control is beyond it, so it stays at the control, the last
  # one included
  expected <- c(3.5, 4.5, -1, -1, -1, -1, -1, -1, -1)
  direct <- function(curves) {
    pred <- t21_pred(curves)
    suppressWarnings(summarise_metric(
      compute_ecx(pred, reference = 5, type = "direct"),
      pred,
      "ecx",
      draws = TRUE
    ))
  }

  expect_equal(direct(t21_curves)$draws[[1]]$value.draw, expected)
  expect_equal(direct(10 - t21_curves)$draws[[2]]$value.draw, expected)
})

test_that("summarise_metric leaves ecx under absolute, relative and range unchanged on the T21 curves (#74)", {
  # each reference comes from the realisation's own control, so no control is
  # past it, and every realisation takes its own crossing or is censored
  for (type in c("absolute", "relative", "range")) {
    for (curves in list(t21_curves, 10 - t21_curves)) {
      pred <- t21_pred(curves)
      rows <- compute_ecx(pred, type = type)

      result <- suppressWarnings(
        summarise_metric(rows, pred, "ecx", draws = TRUE)
      )

      expect_equal(
        lapply(result$draws, function(d) d$value.draw),
        lapply(rows$draws, function(x) replace(x, is.na(x), Inf))
      )
    }
  }
})

test_that("summarise_metric reports a censored quantile at the highest concentration (T21)", {
  x_vec <- 0:10
  # the second realisation reaches its EC10 reference of 9 only at x = 20
  pred <- compute_pred(list(rbind(10 - x_vec, 10 - x_vec / 20)), x_vec)

  expect_warning(
    result <- summarise_metric(
      compute_ecx(pred),
      pred,
      "ecx",
      prob_vals = all_quantiles
    ),
    regexp = "decreasing: 50.0% censored"
  )

  expect_equal(
    unlist(result[1, c("estimate", "conf.low", "conf.high")]),
    c(estimate = 10, conf.low = 1, conf.high = 10)
  )
})

test_that("summarise_metric does not warn when no reported quantile is censored", {
  x_vec <- 0:10
  # 1 of 80 realisations is censored, beyond the reported 2.5% to 97.5%
  slopes <- seq(0.2, 2, length.out = 79)
  curves <- rbind(
    t(vapply(slopes, function(a) 10 - a * x_vec, numeric(11))),
    10 - x_vec / 20
  )
  pred <- compute_pred(list(curves), x_vec)
  rows <- compute_ecx(pred)

  expect_no_warning(
    result <- summarise_metric(rows[1, ], pred, "ecx")
  )
  expect_lt(result$conf.high, 10)
})

test_that("summarise_metric excludes a realisation with no reference and warns (T21, T25)", {
  x_vec <- 0:10
  # the second control is negative, so its absolute reference is NA (T18)
  pred <- compute_pred(list(rbind(10 - x_vec, rep(-1, 11))), x_vec)

  expect_warning(
    result <- summarise_metric(
      compute_ecx(pred),
      pred,
      "ecx",
      prob_vals = all_quantiles
    ),
    regexp = "decreasing: 50.0% excluded"
  )

  expect_equal(
    unlist(result[1, c("estimate", "conf.low", "conf.high")]),
    c(estimate = 1, conf.low = 1, conf.high = 1)
  )
})

test_that("summarise_metric names the group in the warning", {
  x_vec <- 0:10
  pred <- compute_pred(
    list(rbind(10 - x_vec, 10 - x_vec), rbind(10 - x_vec, 10 - x_vec / 20)),
    x_vec,
    keys = data.frame(group = factor(c("A", "B"))),
    group_var = "site"
  )

  expect_warning(
    summarise_metric(
      compute_ecx(pred),
      pred,
      "ecx",
      prob_vals = all_quantiles
    ),
    regexp = "group B, decreasing: 50.0% censored"
  )
})

test_that("summarise_metric summarises control over the realisations not excluded", {
  x_vec <- 0:10
  # under `relative` the first realisation has no increasing asymptote and the
  # second no decreasing one, so each row excludes a different realisation
  pred <- compute_pred(
    list(rbind(10 - x_vec, 8 + x_vec)),
    x_vec,
    asymptote = list(list(decreasing = c(0, NA), increasing = c(NA, 20)))
  )

  result <- suppressWarnings(summarise_metric(
    compute_ecx(pred, type = "relative"),
    pred,
    "ecx"
  ))

  expect_identical(result$direction, c("decreasing", "increasing"))
  expect_equal(result$control, c(10, 8))
})

test_that("summarise_metric gives NA control where every realisation is excluded", {
  x_vec <- 0:10
  # under `relative` no realisation has an increasing asymptote
  pred <- compute_pred(
    list(rbind(10 - x_vec, 10 - 2 * x_vec)),
    x_vec,
    asymptote = list(list(decreasing = c(0, 0), increasing = c(NA_real_, NA_real_)))
  )

  result <- suppressWarnings(summarise_metric(
    compute_ecx(pred, type = "relative"),
    pred,
    "ecx"
  ))

  expect_identical(result$control[2], NA_real_)
  expect_identical(result$reference[2], NA_real_)
  expect_identical(result$estimate[2], NA_real_)
})

test_that("summarise_metric censors a realisation with no span under range", {
  x_vec <- 0:10
  # the second realisation is flat, so it never falls below its control (T38)
  pred <- compute_pred(list(rbind(10 - x_vec, rep(10, 11))), x_vec)

  warnings <- character()
  result <- withCallingHandlers(
    summarise_metric(
      compute_ecx(pred, type = "range"),
      pred,
      "ecx",
      draws = TRUE
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_identical(result$draws[[1]]$value.draw[[2]], Inf)
  expect_equal(result$reference[[1]], 9)
  expect_match(warnings, "decreasing: 50.0% censored", all = FALSE)
  expect_no_match(warnings, "% excluded")
})

test_that("summarise_metric keeps excluding a span of zero outside range", {
  x_vec <- 0:10
  # under `absolute` a control of 0 has no decreasing span
  pred <- compute_pred(list(rbind(10 - x_vec, 0 - x_vec)), x_vec)

  result <- suppressWarnings(summarise_metric(
    compute_ecx(pred),
    pred,
    "ecx",
    draws = TRUE
  ))

  expect_identical(result$draws[[1]]$value.draw[[2]], NA_real_)
})

test_that("summarise_metric sets control to NA under direct", {
  x_vec <- 0:10
  pred <- compute_pred(list(rbind(10 - x_vec, 10 - 2 * x_vec)), x_vec)

  result <- summarise_metric(
    compute_ecx(pred, reference = 5, type = "direct"),
    pred,
    "ecx"
  )

  expect_identical(result$control, c(NA_real_, NA_real_))
  expect_identical(result$reference, c(5, 5))
})

# xform ----------------------------------------------------------------------

test_that("summarise_metric applies xform before summarising", {
  x_vec <- 0:10
  pred <- compute_pred(list(rbind(10 - x_vec, 10 - 2 * x_vec)), x_vec)

  result <- summarise_metric(
    compute_ecx(pred),
    pred,
    "ecx",
    prob_vals = all_quantiles,
    xform = function(x) x^2
  )

  # the median of 1 and 0.25, not the square of the median of 1 and 0.5
  expect_equal(result$estimate[1], 0.625)
})

test_that("summarise_metric reports a censored quantile as xform of the highest concentration", {
  x_vec <- 0:10
  pred <- compute_pred(list(rbind(10 - x_vec, 10 - x_vec / 20)), x_vec)

  result <- suppressWarnings(summarise_metric(
    compute_ecx(pred),
    pred,
    "ecx",
    prob_vals = all_quantiles,
    xform = function(x) x * 10
  ))

  expect_equal(result$conf.high[1], 100)
})

# draws ----------------------------------------------------------------------

test_that("summarise_metric adds a draws list column when draws = TRUE", {
  x_vec <- 0:10
  pred <- compute_pred(list(rbind(10 - x_vec, 10 - x_vec / 20)), x_vec)

  result <- suppressWarnings(
    summarise_metric(compute_ecx(pred), pred, "ecx", draws = TRUE)
  )

  expect_identical(names(result)[ncol(result)], "draws")
  expect_equal(
    result$draws[[1]],
    dplyr::tibble(.draw = 1:2, value.draw = c(1, Inf))
  )
})

test_that("summarise_metric has no draws column by default", {
  x_vec <- 0:10
  pred <- compute_pred(list(rbind(10 - x_vec, 10 - 2 * x_vec)), x_vec)

  result <- summarise_metric(compute_ecx(pred), pred, "ecx")

  expect_false("draws" %in% names(result))
})

test_that("summarise_metric pairs ecnsec with the nsec realisations in draws", {
  x_vec <- 0:10
  pred <- compute_pred(list(rbind(10 - x_vec, 10 - 2 * x_vec)), x_vec)
  rows <- compute_nsec(pred)

  result <- summarise_metric(rows, pred, "nsec", draws = TRUE)

  expect_named(result$draws[[1]], c(".draw", "value.draw", "ecnsec.draw"))
  expect_identical(result$draws[[1]]$ecnsec.draw, rows$ecnsec[[1]])
})

test_that("summarise_metric draws unnest alongside the summary columns", {
  skip_if_not_installed("tidyr")
  x_vec <- 0:10
  pred <- compute_pred(list(rbind(10 - x_vec, 10 - 2 * x_vec)), x_vec)

  result <- summarise_metric(compute_nsec(pred), pred, "nsec", draws = TRUE)
  long <- tidyr::unnest(result, draws)

  expect_equal(nrow(long), 4L)
  expect_true(all(
    c("ecnsec", ".draw", "value.draw", "ecnsec.draw") %in% names(long)
  ))
})

# argument checks ------------------------------------------------------------

test_that("summarise_metric errors on prob_vals out of order", {
  x_vec <- 0:10
  pred <- compute_pred(list(rbind(10 - x_vec, 10 - 2 * x_vec)), x_vec)

  expect_error(
    summarise_metric(
      compute_ecx(pred),
      pred,
      "ecx",
      prob_vals = c(0.025, 0.5, 0.975)
    ),
    regexp = "`prob_vals` must be the central, lower and upper quantiles"
  )
})

test_that("summarise_metric errors when xform is not a function", {
  x_vec <- 0:10
  pred <- compute_pred(list(rbind(10 - x_vec, 10 - 2 * x_vec)), x_vec)

  expect_error(
    summarise_metric(compute_ecx(pred), pred, "ecx", xform = "log"),
    regexp = "`xform` must be a function"
  )
})

# fitted models --------------------------------------------------------------

test_that("summarise_metric matches today's decreasing EC10 on brms_model_1, less the excluded draws", {
  pred <- toxval_predict(brms_model_1, x_var = "x", resolution = 500)

  # the increasing direction is censored, and 8 of 1000 controls are negative
  result <- suppressWarnings(
    summarise_metric(compute_ecx(pred), pred, "ecx")
  )

  # `ecx(brms_model_1, x_var = "x", resolution = 500)` at a3f682f gives 0.8322880
  # (0.8195822 to 1.05). The per-draw values agree; the difference is the 8
  # draws with a negative control, whose reference is NA and which are excluded
  # (T18, T21) rather than given a crossing of a reference above the control.
  expect_equal(
    unlist(result[1, c("estimate", "conf.low", "conf.high")]),
    c(estimate = 0.832247, conf.low = 0.8197345, conf.high = 1.05),
    tolerance = 1e-6
  )
})

test_that("summarise_metric reports an absolute reference of 0.9 times the control on brms_model_1 (T34)", {
  pred <- toxval_predict(brms_model_1, x_var = "x", resolution = 500)

  # 8 of 1000 controls are negative and excluded from the row
  result <- suppressWarnings(
    summarise_metric(compute_ecx(pred), pred, "ecx")
  )

  expect_equal(result$reference[1], 0.9 * result$control[1])
})

test_that("summarise_metric warns for the censored increasing EC10 on brms_model_1", {
  pred <- toxval_predict(brms_model_1, x_var = "x", resolution = 500)

  expect_warning(
    summarise_metric(compute_ecx(pred), pred, "ecx"),
    regexp = "increasing: [0-9.]+% censored"
  )
})

test_that("summarise_metric summarises a drc NSEC from bootstrap realisations", {
  pred <- toxval_predict(
    nsec_drc_1,
    resolution = 200,
    n_boot = 200,
    seed = 42
  )

  result <- summarise_metric(compute_nsec(pred), pred, "nsec")

  expect_s3_class(result, "toxval")
  expect_true(is.finite(result$estimate[1]))
  expect_true(result$conf.low[1] <= result$estimate[1])
  expect_true(result$estimate[1] <= result$conf.high[1])
})
