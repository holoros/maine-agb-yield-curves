#!/usr/bin/env Rscript
# bakuzis_check.R  --  automated Bakuzis-matrix consistency assessment
# Usage: Rscript bakuzis_check.R trajectories.csv [form_factor=0.40] [qmd_rate=0.10]
# Input CSV (long, one row per group x age) needs columns:
#   site (productivity class, numeric or ordered), age, HT, TPH, BAPH, QMD
#   optional: origin (a second grouping), VOL (else VOL = BAPH*HT*form_factor)
# Output: a PASS/FLAG report and a recommendations block to stdout.
# Base R only. Falsifies structurally flawed projections; does not confirm accuracy.
#
# PATCHED September 7, 2026, against the defect found on the FVS-ACD linked survival
# projection (RT5-E). The shipped Reineke screen was
#     keep <- c(TRUE, diff(g$QMD) > 0.5)
# an absolute 0.5 cm threshold on the quadratic mean diameter increment between
# successive rows. That is a rate expressed as a difference, so it is a function of the
# projection step. It was written for five and ten year steps, where 0.5 cm is roughly
# 0.1 cm yr-1. Run on a one-year step it demands 0.5 cm yr-1, removes every point on any
# slowly growing trajectory, and the else branch then prints "not in a self-thinning
# regime" and adds no flag. The run returns FLAGS: 0, which reads as a clean pass on the
# single most diagnostic relation in the matrix, on an axis that was never evaluated.
# Two changes fix it and a third makes the failure unmissable.
#   1  The threshold is a RATE, qmd_rate cm yr-1 (default 0.10, the per-year equivalent of
#      the original 0.5 cm at a five year step), multiplied by the actual age step of each
#      trajectory. A step of any length now screens for the same biological condition.
#   2  A trajectory that loses stems but leaves too few points after the screen is FLAGGED
#      as not evaluated rather than passed over in silence. Only a trajectory that holds
#      its stem count is legitimately n/a, since Reineke has nothing to say about it.
#   3  The summary prints evaluated, not-evaluated and flagged counts beside FLAGS, so a
#      zero flag count cannot be read without the number of axes it rests on.

args <- commandArgs(trailingOnly = TRUE)
csv  <- if (length(args) >= 1) args[1] else stop("supply a trajectories CSV")
ff   <- if (length(args) >= 2) as.numeric(args[2]) else 0.40
qrate<- if (length(args) >= 3) as.numeric(args[3]) else 0.10   # cm yr-1, screen for self-thinning phase
d <- read.csv(csv, stringsAsFactors = FALSE)
if (!"origin" %in% names(d)) d$origin <- "all"
if (!"VOL"    %in% names(d)) d$VOL <- d$BAPH * d$HT * ff
d <- d[order(d$origin, d$site, d$age), ]
sites   <- sort(unique(d$site))
origins <- unique(d$origin)
flags <- character(0)
mono  <- function(x, inc = TRUE, tol = 1e-6) if (inc) all(diff(x) >= -tol) else all(diff(x) <= tol)
cat(strrep("=", 64), "\nBAKUZIS ASSESSMENT:", csv, "\n", strrep("=", 64), "\n", sep = "")

## 1 + 3 monotonicity (HT up, TPH down) -------------------------------------
cat("\n[1/3 Monotonicity] HT increasing, TPH decreasing with age:\n")
for (o in origins) for (s in sites) {
  g <- d[d$origin == o & d$site == s, ]
  okh <- mono(g$HT, TRUE); okt <- mono(g$TPH, FALSE)
  cat(sprintf("  %s site %s: HT_up=%s TPH_down=%s\n", o, s, okh, okt))
  if (!okh) flags <- c(flags, sprintf("HT not monotone up (%s site %s)", o, s))
  if (!okt) flags <- c(flags, sprintf("TPH not monotone down (%s site %s)", o, s))
}

## 1 + 7 site ordering (better site not overtaken in HT or VOL) -------------
cat("\n[1/7 Site ordering] HT and VOL ordered by site at each age:\n")
ageset <- sort(unique(d$age))
for (o in origins) {
  bh <- bv <- 0
  for (a in ageset) {
    h <- sapply(sites, function(s) { v <- d$HT [d$origin==o & d$site==s & d$age==a]; if(length(v)) v[1] else NA })
    v <- sapply(sites, function(s) { z <- d$VOL[d$origin==o & d$site==s & d$age==a]; if(length(z)) z[1] else NA })
    if (!any(is.na(h)) && is.unsorted(h)) bh <- bh + 1
    if (!any(is.na(v)) && is.unsorted(v)) bv <- bv + 1
  }
  cat(sprintf("  %s: HT inversions=%d/%d  VOL inversions=%d/%d\n", o, bh, length(ageset), bv, length(ageset)))
  if (bh) flags <- c(flags, sprintf("HT site-ordering inverted %dx (%s)", bh, o))
  if (bv) flags <- c(flags, sprintf("VOL site-ordering inverted %dx (%s) [check self-thinning/SDImax]", bv, o))
}

## 4 Reineke slope ----------------------------------------------------------
rein_eval <- rein_na <- rein_flag <- 0L
steps <- unlist(lapply(origins, function(o) unlist(lapply(sites, function(s) {
  a <- d$age[d$origin==o & d$site==s]; if (length(a) > 1) stats::median(diff(sort(a))) else NA_real_ }))))
step_med <- stats::median(steps, na.rm = TRUE)
cat(sprintf("\n[4 Reineke] log(TPH)~log(QMD) self-thinning slope (target ~ -1.605).\n  Median age step %.3g yr; self-thinning screen QMD advance > %.3g cm yr-1, i.e. > %.3g cm per step.\n",
            step_med, qrate, qrate*step_med))
for (o in origins) for (s in sites) {
  g <- d[d$origin==o & d$site==s, ]; g <- g[order(g$age), ]
  if (nrow(g) < 2) next
  stp <- stats::median(diff(g$age)); if (!is.finite(stp) || stp <= 0) stp <- step_med
  # self-thinning phase = points where QMD advances faster than qrate cm yr-1 (excludes capped QMD).
  # The threshold scales with the step, so the same biological condition is screened at any step length.
  keep <- c(TRUE, diff(g$QMD) > qrate*stp)
  gk <- g[keep, ]
  loses_stems <- (max(g$TPH) - min(g$TPH))/max(g$TPH) > 0.01   # 1% of stems is a real self-thinning phase
  if (nrow(gk) >= 3 && sd(log(gk$QMD)) > 1e-3) {
    sl <- coef(lm(log(TPH) ~ log(QMD), gk))[2]
    ok <- is.finite(sl) && sl > -2.2 && sl < -1.2
    rein_eval <- rein_eval + 1L
    cat(sprintf("  %s site %s: slope=%.2f on %d of %d points %s\n", o, s, sl, nrow(gk), nrow(g), ifelse(ok,"OK","FLAG")))
    if (!ok) { flags <- c(flags, sprintf("Reineke slope %.2f off (%s site %s)", sl, o, s)); rein_flag <- rein_flag + 1L }
  } else if (loses_stems) {
    # the trajectory does self-thin, so the axis is testable and the screen removed it. Fail loudly.
    rein_na <- rein_na + 1L
    cat(sprintf("  %s site %s: NOT EVALUATED, screen left %d of %d points while the stand loses %.1f%% of its stems - FLAG\n",
                o, s, nrow(gk), nrow(g), 100*(max(g$TPH)-min(g$TPH))/max(g$TPH)))
    flags <- c(flags, sprintf("Reineke NOT EVALUATED, screen left %d of %d points on a self-thinning trajectory (%s site %s); check the age step against qmd_rate",
                              nrow(gk), nrow(g), o, s))
  } else {
    rein_na <- rein_na + 1L
    cat(sprintf("  %s site %s: not in a self-thinning regime (stem count flat, QMD capped) - Reineke n/a\n", o, s))
  }
}

## 8 Eichhorn: volume-height site independence ------------------------------
cat("\n[8 Eichhorn] volume-height site independence (CV of VOL across sites at matched HT):\n")
for (o in origins) {
  hi <- max(sapply(sites, function(s) min(d$HT[d$origin==o & d$site==s])))
  hx <- min(sapply(sites, function(s) max(d$HT[d$origin==o & d$site==s])))
  if (hx > hi) {
    grid <- seq(hi, hx, length.out = 8)
    cv <- sapply(grid, function(hg) {
      vv <- sapply(sites, function(s) { g <- d[d$origin==o & d$site==s, ]; approx(g$HT, g$VOL, hg, rule=2)$y })
      sd(vv) / mean(vv) })
    m <- mean(cv); ok <- m < 0.15
    cat(sprintf("  %s: mean CV=%.1f%% %s\n", o, 100*m, ifelse(ok,"OK (site-independent)","FLAG (site-dependent V-H)")))
    if (!ok) flags <- c(flags, sprintf("Eichhorn: V-H site-dependent CV %.0f%% (%s)", 100*m, o))
  }
}

## summary ------------------------------------------------------------------
cat("\n", strrep("=", 64), "\nFLAGS: ", length(flags),
    sprintf("   [Reineke: %d trajectories evaluated, %d not evaluated, %d flagged; age step %.3g yr]\n",
            rein_eval, rein_na, rein_flag, step_med), sep = "")
if (length(flags)) for (f in flags) cat("  - ", f, "\n", sep = "") else cat("  none - projection passes all tested relations (falsification failed; not a confirmation of accuracy)\n")
if (rein_eval == 0) cat("  WARNING: the Reineke axis was never evaluated. A flag count of zero says nothing about self-thinning.\n")
cat(strrep("=", 64), "\nInterpret with references/bakuzis-method.md: map each flag to the responsible submodel and a refinement.\n", sep="")
