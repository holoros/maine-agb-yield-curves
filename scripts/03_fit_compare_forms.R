## Maine FIA biomass yield curves, step 3 (revision b): compare curve forms under two fitting strategies, select, and build
## generic AGB on stand age curves by forest type group and ecoregion section.
##
## Forms (multiplicative level a times shape f):
##   CR  Chapman-Richards, monotone          f = (1 - exp(-k t))^p
##   CRS CR with exponential senescence      f = (1 - exp(-k t))^p exp(-d t)               declines toward 0
##   CRD CR with bounded decline             f = (1 - exp(-k t))^p (1 - q (1 - exp(-d t)))  declines to a floor (1 - q)
##   HG  Hoerl gamma                         f = t^b exp(-c t)                               declines toward 0
## Strategies:
##   S1  remeasurement only: algebraic difference on every harvest free interval, arithmetic scale (red team must fix 1),
##       residual (B2 - B1 f(t2)/f(t1)) / (B1 sqrt(dt)); the level cancels, so shape comes only from how plots of known age change.
##   S2  joint: log B = log a + u_section + u_plot + log f(t) on harvest free plot visits (chronosequence plus repeated measures),
##       shape by profile maximum likelihood over lme4 fits, by type.
##   Level for both: log a = log B - log f(t), lmer with type fixed effects, type by section and plot random effects.
## Selection rule. Original rule (written before any fit): biological screen, then lowest CV RMSE of annual net AGB change on
##   intervals of 10 yr or more with midpoint age 80 yr or more, ties within 2% to fewer parameters. Revision b, made after a 4
##   resample test run showed every S1 shape unbounded (no peak, f(100)/f(300) of 0.1 to 0.3) and after the red team found the
##   magnitude screen was computed but not applied: (i) the magnitude screen is applied (statewide median curve must stay at or below
##   the type's observed 99.9th percentile AGB through 250 yr); (ii) a decline floor is added, B(250) at least 0.5 of the curve
##   maximum within 250 yr [ASSUMPTION: gross screen against collapse to zero; not a measured floor]; (iii) peak, if any, not
##   before 40 yr; (iv) strategy S2 is added as a candidate; ranking criterion unchanged, on the common subset of intervals every
##   candidate predicted. Equivalence: TOST on annual net AGB change, plot cluster bootstrap, intercept region +/-0.5 Mg ha-1 yr-1
##   absolute [ASSUMPTION: carbon accounting tolerance], slope region +/-0.25; smallest passing regions reported.
## Usage: Rscript 03_fit_compare_forms.R <out_dir> [nboot_curve=100] [nboot_eq=2000] [cores=24]
suppressMessages({ library(data.table); library(minpack.lm); library(lme4); library(ggplot2); library(patchwork); library(mgcv); library(parallel) })
a <- commandArgs(TRUE); OD <- a[1]; NBC <- if (length(a) > 1) as.integer(a[2]) else 100L
NBE <- if (length(a) > 2) as.integer(a[3]) else 2000L; NC <- if (length(a) > 3) as.integer(a[4]) else 24L
YD <- file.path(OD, "curves"); dir.create(YD, showWarnings = FALSE)
elog <- file.path(YD, "error_log.txt"); try_log <- function(expr, what) tryCatch(expr, error = function(e) {
  cat(sprintf("[%s] %s: %s\n", format(Sys.time()), what, conditionMessage(e)), file = elog, append = TRUE); NULL })
RNGkind("L'Ecuyer-CMRG"); set.seed(20261015)
FT <- c(SF = "Spruce-fir", MBB = "Maple-beech-birch", AB = "Aspen-birch", WRPH = "White-red pine-hemlock", OAK = "Oak and oak-pine",
        EAC = "Elm-ash-cottonwood")
I <- fread(file.path(OD, "me_intervals_DATA.csv")); V <- fread(file.path(OD, "me_plot_visits_DATA.csv"))
D <- I[hf == TRUE & !harvested1 & ftg %in% names(FT) & AGB1 >= 10 & AGB2 >= 0 & !is.na(STDAGE1) & dt >= 3]
D[, `:=`(t1 = pmax(STDAGE1, 1), t2 = pmax(STDAGE1, 1) + dt, dB = (AGB2 - AGB1) / dt)]
D[, `:=`(mid_age = t1 + dt / 2, horizon = fifelse(dt <= 7, "5", fifelse(dt <= 12, "10", "15-25")))]
cat(sprintf("intervals used: %d on %d plots; AGB2 == 0 kept: %d\n", nrow(D), uniqueN(D$plot_ix), D[AGB2 == 0, .N]))
VV <- V[ftg %in% names(FT) & !harvested & AGB > 0.5 & !is.na(STDAGE) & STDAGE >= 1]
VV[, cell := paste(ftg, section, sep = ":")]
cat(sprintf("visits for level and S2: %d on %d plots\n", nrow(VV), uniqueN(VV$plot_ix)))
mag <- V[ftg %in% names(FT) & !harvested, .(agb999 = quantile(AGB, .999), age99 = quantile(STDAGE, .99, na.rm = TRUE)), by = ftg]

F <- list(
  CR  = list(np = 2, st = c(lk = log(0.03), lp = log(2)),
             lf = function(t, th) exp(th[2]) * log1p(-exp(-exp(th[1]) * t))),
  CRS = list(np = 3, st = c(lk = log(0.04), lp = log(2), ld = log(0.002)),
             lf = function(t, th) exp(th[2]) * log1p(-exp(-exp(th[1]) * t)) - exp(th[3]) * t),
  CRD = list(np = 4, st = c(lk = log(0.04), lp = log(2), lq = qlogis(0.3), ld = log(0.02)),
             lf = function(t, th) exp(th[2]) * log1p(-exp(-exp(th[1]) * t)) + log1p(-plogis(th[3]) * (1 - exp(-exp(th[4]) * t)))),
  HG  = list(np = 2, st = c(lb = log(1.5), lc = log(0.01)),
             lf = function(t, th) exp(th[1]) * log(t) - exp(th[2]) * t))
nat <- function(m, th) switch(m, CR = c(k = exp(th[1]), p = exp(th[2])), CRS = c(k = exp(th[1]), p = exp(th[2]), d = exp(th[3])),
  CRD = c(k = exp(th[1]), p = exp(th[2]), q = plogis(th[3]), d = exp(th[4])), HG = c(b = exp(th[1]), c = exp(th[2])))

## S1: arithmetic ADA
fit_s1 <- function(d, fm, st = NULL) {
  obj <- function(th) { r <- (d$AGB2 - d$AGB1 * exp(fm$lf(d$t2, th) - fm$lf(d$t1, th))) / (d$AGB1 * sqrt(d$dt)); r[!is.finite(r)] <- 1e3; r }
  best <- NULL; starts <- c(list(fm$st), if (!is.null(st)) list(st), lapply(1:4, function(i) fm$st + rnorm(length(fm$st), 0, 0.5)))
  for (s0 in starts) { r <- tryCatch(nls.lm(s0, fn = obj, control = nls.lm.control(maxiter = 500)), error = function(e) NULL)
    if (!is.null(r) && all(is.finite(r$par)) && (is.null(best) || r$deviance < best$deviance)) best <- r }
  if (is.null(best)) return(NULL)
  n <- nrow(d); se <- tryCatch(sqrt(diag(chol2inv(chol(best$hessian))) * best$deviance / (n - length(best$par))), error = function(e) rep(NA_real_, length(best$par)))
  list(th = best$par, se = se, crit = best$deviance, n = n) }
## S2: joint mixed model on visits, profile ML over the shape parameters
ctrl <- lmerControl(calc.derivs = FALSE, check.conv.singular = "ignore")
dev_s2 <- function(th, v, fm) { la <- log(v$AGB) - fm$lf(v$STDAGE, th); if (any(!is.finite(la))) return(1e10)
  dd <- data.frame(la = la, section = v$section, plot_ix = v$plot_ix)
  m <- tryCatch(if (length(unique(dd$section)) > 1) lmer(la ~ 1 + (1 | section) + (1 | plot_ix), dd, REML = FALSE, control = ctrl) else
    lmer(la ~ 1 + (1 | plot_ix), dd, REML = FALSE, control = ctrl), error = function(e) NULL)
  if (is.null(m)) 1e10 else -2 * as.numeric(logLik(m)) }
fit_s2 <- function(v, fm, st = NULL, maxit = 400) {
  starts <- c(list(fm$st), if (!is.null(st)) list(st)); best <- NULL
  for (s0 in starts) { r <- tryCatch(optim(s0, dev_s2, v = v, fm = fm, method = "Nelder-Mead", control = list(maxit = maxit, reltol = 1e-9)), error = function(e) NULL)
    if (!is.null(r) && (is.null(best) || r$value < best$value)) best <- r }
  if (is.null(best)) return(NULL)
  list(th = best$par, se = rep(NA_real_, length(best$par)), crit = best$value, n = nrow(v), conv = best$convergence) }
fit_any <- function(strat, f, m, d, v, st = NULL, maxit = 400) if (strat == "S1") fit_s1(d[ftg == f], F[[m]], st) else fit_s2(v[ftg == f], F[[m]], st, maxit)

## level model, all types together
level_fit <- function(vv, thl) {
  vv <- copy(vv); vv[, la := log(AGB) - F[[thl[[.BY$ftg]]$form]]$lf(STDAGE, thl[[.BY$ftg]]$th), by = ftg]
  m <- lmer(la ~ 0 + ftg + (1 | cell) + (1 | plot_ix), data = vv, REML = TRUE, control = ctrl)
  fe <- fixef(m); re <- ranef(m)$cell; vc <- as.data.frame(VarCorr(m))
  list(fe = setNames(fe, sub("^ftg", "", names(fe))), re = setNames(re[, 1], rownames(re)), s_plot = vc$sdcor[vc$grp == "plot_ix"],
       s_cell = vc$sdcor[vc$grp == "cell"], s_res = vc$sdcor[vc$grp == "Residual"], resid = vv[, .(ftg, STDAGE, r = la - predict(m, re.form = ~ (1 | cell)))]) }
curve_of <- function(lv, thl, f, ce, tt) { u <- if (!is.null(ce) && ce %in% names(lv$re)) lv$re[[ce]] else 0
  exp(lv$fe[[f]] + u + F[[thl[[f]]$form]]$lf(tt, thl[[f]]$th)) }

## ---- 1. full data fits, all strategies x forms x types (parallel)
grid <- CJ(strat = c("S1", "S2"), ftg = names(FT), form = names(F))
full <- mclapply(seq_len(nrow(grid)), function(i) { g <- grid[i]; r <- try_log(fit_any(g$strat, g$ftg, g$form, D, VV), paste("full", g$strat, g$ftg, g$form))
  if (is.null(r)) NULL else c(as.list(g), list(fit = r)) }, mc.cores = NC)
full <- Filter(Negate(is.null), full)
P <- rbindlist(lapply(full, function(x) data.table(strat = x$strat, ftg = x$ftg, form = x$form, param = names(nat(x$form, x$fit$th)),
  est = nat(x$form, x$fit$th), est_unconstrained = x$fit$th, se_unconstrained = x$fit$se, crit = x$fit$crit, n = x$fit$n)))
thget <- function(s_, f_, m_) { x <- P[strat == s_ & ftg == f_ & form == m_]; if (!nrow(x)) NULL else x$est_unconstrained }
cat("\n== shape parameters\n"); print(P[, .(strat, ftg, form, param, est = signif(est, 4))])

## ---- 2. plot grouped 5 fold CV: interval projections (a) and generic curve predictions of held out visits (b)
plots <- unique(c(D$plot_ix, VV$plot_ix)); fold <- data.table(plot_ix = plots, fold = sample(rep(1:5, length.out = length(plots))))
D <- merge(D, fold, by = "plot_ix"); VV <- merge(VV, fold, by = "plot_ix")
cvg <- CJ(strat = c("S1", "S2"), ftg = names(FT), form = names(F), k = 1:5)
cvfits <- mclapply(seq_len(nrow(cvg)), function(i) { g <- cvg[i]
  r <- try_log(fit_any(g$strat, g$ftg, g$form, D[fold != g$k], VV[fold != g$k], st = thget(g$strat, g$ftg, g$form), maxit = 250), paste("cv", g$strat, g$ftg, g$form, g$k))
  if (is.null(r)) NULL else c(as.list(g), list(th = r$th)) }, mc.cores = NC)
cvfits <- Filter(Negate(is.null), cvfits)
O <- rbindlist(lapply(cvfits, function(x) { te <- D[fold == x$k & ftg == x$ftg]; fm <- F[[x$form]]
  te[, .(strat = x$strat, form = x$form, plot_ix, ftg, fold, dt, horizon, mid_age, t1, t2, AGB1, AGB2, dB,
         B2hat = AGB1 * exp(fm$lf(t2, x$th) - fm$lf(t1, x$th)))] }))
O[, dBhat := (B2hat - AGB1) / dt]
O[, key := paste(plot_ix, t1, t2, AGB1)]; ncand <- uniqueN(O[, .(strat, form)])
common <- O[, .N, by = key][N == ncand]$key; cat(sprintf("\ncommon subset: %d of %d intervals\n", length(common), uniqueN(O$key)))
O <- O[key %in% common]
O[, s := abs(AGB2 - B2hat) / (AGB1 * sqrt(dt))]
O[, q90 := { qq <- sapply(1:5, function(k) quantile(s[fold != k], .9)); qq[fold] }, by = .(strat, form, ftg)]
O[, covered := abs(AGB2 - B2hat) <= q90 * AGB1 * sqrt(dt)]
## (b) generic curve CV: level refit per fold and strategy x form (types share the level model)
Gcv <- rbindlist(mclapply(seq_len(nrow(CJ(strat = c("S1", "S2"), form = names(F), k = 1:5))), function(i) {
  g <- CJ(strat = c("S1", "S2"), form = names(F), k = 1:5)[i]
  thl <- setNames(lapply(names(FT), function(f) { x <- Filter(function(z) z$strat == g$strat && z$ftg == f && z$form == g$form && z$k == g$k, cvfits)
    if (!length(x)) NULL else list(form = g$form, th = x[[1]]$th) }), names(FT))
  if (any(sapply(thl, is.null))) return(NULL)
  lv <- try_log(level_fit(VV[fold != g$k], thl), paste("cvlevel", g$strat, g$form, g$k)); if (is.null(lv)) return(NULL)
  te <- VV[fold == g$k]; sdp <- sqrt(lv$s_plot^2 + lv$s_res^2)
  te[, med := curve_of(lv, thl, ftg[1], cell[1], STDAGE), by = .(ftg, cell)]
  te[, .(strat = g$strat, form = g$form, ftg, plot_ix, STDAGE, AGB, med, mean = med * exp(sdp^2 / 2),
         in90 = AGB >= med * exp(-1.645 * sdp) & AGB <= med * exp(1.645 * sdp))] }, mc.cores = NC))

## metrics and equivalence
met <- function(z) z[, .(n = .N, bias_dB = mean(dBhat - dB), rmse_dB = sqrt(mean((dBhat - dB)^2)), mae_dB = mean(abs(dBhat - dB)),
                       bias_B2 = mean(B2hat - AGB2), rmse_B2 = sqrt(mean((B2hat - AGB2)^2)), cover90 = mean(covered))]
M_all <- O[, met(.SD), by = .(ftg, strat, form)]; M_long <- O[dt >= 10 & mid_age >= 80, met(.SD), by = .(ftg, strat, form)]
M_h <- O[, met(.SD), by = .(ftg, strat, form, horizon)]
M_age <- O[, met(.SD), by = .(ftg, strat, form, ageclass = cut(mid_age, c(0, 40, 80, 120, 300), right = FALSE))]
Gm <- Gcv[, .(n = .N, bias_med = mean(med - AGB), bias_mean = mean(mean - AGB), rmse_mean = sqrt(mean((mean - AGB)^2)), cover90 = mean(in90)),
          by = .(ftg, strat, form, ageclass = cut(STDAGE, c(0, 40, 80, 120, 300), right = FALSE))]
Gall <- Gcv[, .(n = .N, bias_mean = mean(mean - AGB), rmse_mean = sqrt(mean((mean - AGB)^2)), cover90 = mean(in90)), by = .(ftg, strat, form)]
eq_boot <- function(d, nb) { ids <- unique(d$plot_ix); sp <- split(seq_len(nrow(d)), d$plot_ix)
  est <- function(z) { xb <- mean(z$dBhat); m <- lm(dB ~ I(dBhat - xb), z); c(coef(m)[1] - xb, coef(m)[2]) }
  e0 <- est(d); bs <- replicate(nb, tryCatch(est(d[unlist(sp[sample(as.character(ids), length(ids), TRUE)], use.names = FALSE)]), error = function(e) c(NA, NA)))
  ci0 <- quantile(bs[1, ], c(.05, .95), na.rm = TRUE); ci1 <- quantile(bs[2, ], c(.05, .95), na.rm = TRUE)
  data.table(int_dev = e0[1], int_lo = ci0[1], int_hi = ci0[2], slope = e0[2], slope_lo = ci1[1], slope_hi = ci1[2], nboot = nb,
             int_pass_0.5 = max(abs(ci0)) <= 0.5, slope_pass_0.25 = max(abs(ci1 - 1)) <= 0.25,
             int_min_region_abs = max(abs(ci0)), slope_min_region = max(abs(ci1 - 1))) }
eqg <- unique(O[, .(ftg, strat, form)])
EQ <- rbindlist(mclapply(seq_len(nrow(eqg)), function(i) { g <- eqg[i]; d <- O[ftg == g$ftg & strat == g$strat & form == g$form]
  cbind(g, eq_boot(d, NBE), eq_boot(d[dt >= 10], max(500, NBE %/% 4))[, setNames(.SD, paste0(names(.SD), "_long"))]) }, mc.cores = NC))

## ---- 3. level and curves for every candidate; screens; selection
TT <- 0:250
cand <- unique(P[, .(strat, form)])
lev <- list(); CT <- list(); resid_tab <- list()
for (i in seq_len(nrow(cand))) { g <- cand[i]
  thl <- setNames(lapply(names(FT), function(f) { th <- thget(g$strat, f, g$form); if (is.null(th)) NULL else list(form = g$form, th = th) }), names(FT))
  if (any(sapply(thl, is.null))) next
  lv <- try_log(level_fit(VV, thl), paste("level", g$strat, g$form)); if (is.null(lv)) next
  nm <- paste(g$strat, g$form); lev[[nm]] <- list(lv = lv, thl = thl)
  resid_tab[[nm]] <- lv$resid[, .(strat = g$strat, form = g$form, ftg, STDAGE, r)]
  for (f in names(FT)) CT[[length(CT) + 1]] <- data.table(strat = g$strat, form = g$form, ftg = f, age = TT[-1], AGB_med = curve_of(lv, thl, f, NULL, TT[-1]))
}
CT <- rbindlist(CT)
## early behavior screen (added after the first full run showed one candidate, AB S2 CRD, with a spike near age 1): AGB must rise
## monotonically from age 1 to its peak (tolerance 0.1% of the maximum) and start below 5% of the maximum at age 1
scr <- merge(CT[, .(max250 = max(AGB_med), at250 = AGB_med[age == 250], peak_age = age[which.max(AGB_med)],
                    early_ok = { pk <- which.max(AGB_med); all(diff(AGB_med[1:pk]) >= -0.001 * max(AGB_med)) && AGB_med[1] < 0.05 * max(AGB_med) }),
                 by = .(strat, form, ftg)], mag, by = "ftg")
scr[, `:=`(over999 = max250 > agb999, floor_ok = at250 >= 0.5 * max250, peak_ok = peak_age >= 40)][, bio_ok := !over999 & floor_ok & peak_ok & early_ok]
S <- merge(scr, M_long[, .(ftg, strat, form, n_long = n, rmse_long = rmse_dB, bias_long = bias_dB)], by = c("ftg", "strat", "form"), all.x = TRUE)
S <- merge(S, M_all[, .(ftg, strat, form, rmse_all = rmse_dB)], by = c("ftg", "strat", "form"), all.x = TRUE)
S <- merge(S, P[, .(np = .N), by = .(strat, ftg, form)], by = c("ftg", "strat", "form"))
S[, crit_basis := fifelse(!is.na(n_long) & n_long >= 30, "long_old", "all")][, crit := fifelse(crit_basis == "long_old", rmse_long, rmse_all)]
S[, best_crit := suppressWarnings(min(crit[bio_ok], na.rm = TRUE)), by = ftg][, within2 := bio_ok & crit <= best_crit * 1.02]
sel <- S[within2 == TRUE][order(ftg, np, crit)][, .SD[1], by = ftg][, .(ftg, strat, form, crit_basis)]
miss <- setdiff(names(FT), sel$ftg)
if (length(miss)) { cat("\nNO CANDIDATE PASSED THE SCREEN for:", miss, "; fallback = lowest criterion among curves under the 99.9th percentile, flagged\n")
  fb <- S[ftg %in% miss & !over999][order(ftg, crit)][, .SD[1], by = ftg][, .(ftg, strat, form, crit_basis = paste(crit_basis, "FALLBACK"))]
  sel <- rbind(sel, fb) }
cat("\n== screen\n"); print(S[order(ftg, strat, form), .(ftg, strat, form, max250 = round(max250), agb999 = round(agb999), peak_age, at250 = round(at250), early_ok, bio_ok, crit = round(crit, 3))])
cat("\n== selected\n"); print(sel)

## ---- 4. selected set: curves by section with joint plot bootstrap bands
thl_sel <- setNames(lapply(names(FT), function(f) { s <- sel[ftg == f]; if (!nrow(s)) NULL else list(form = s$form, strat = s$strat, th = thget(s$strat, f, s$form)) }), names(FT))
lv_sel <- level_fit(VV, thl_sel)
cells <- VV[, .(n_visits = .N, n_plots = uniqueN(plot_ix), age_q99 = quantile(STDAGE, .99), age_max = max(STDAGE)), by = .(ftg, section, cell)]
allcells <- rbind(cells[, .(ftg, section, cell, n_plots, age_q99)], mag[, .(ftg, section = "ALL", cell = paste0(ftg, ":ALL"), n_plots = NA_integer_, age_q99 = age99)])
sdp <- sqrt(lv_sel$s_plot^2 + lv_sel$s_res^2)
SC <- allcells[, .(age = TT[-1], AGB_med = curve_of(lv_sel, thl_sel, ftg, if (section == "ALL") NULL else cell, TT[-1])), by = .(ftg, section, cell, n_plots, age_q99)]
SC[, `:=`(AGB_mean = AGB_med * exp(sdp^2 / 2), AGB_p25 = AGB_med * exp(-0.674 * lv_sel$s_plot), AGB_p75 = AGB_med * exp(0.674 * lv_sel$s_plot),
          PI90_lo = AGB_med * exp(-1.645 * sdp), PI90_hi = AGB_med * exp(1.645 * sdp), extrapolated = age > age_q99)]
cl <- split(seq_len(nrow(D)), D$plot_ix); clv <- split(seq_len(nrow(VV)), VV$plot_ix); allp <- union(names(cl), names(clv))
bt <- mclapply(seq_len(NBC), function(b) { s <- sample(allp, length(allp), replace = TRUE)
  sd_ <- s[s %in% names(cl)]; Db <- D[unlist(cl[sd_], use.names = FALSE)]
  sv <- s[s %in% names(clv)]; Vb <- VV[unlist(clv[sv], use.names = FALSE)]
  Vb[, plot_ix := paste(plot_ix, rep(seq_along(sv), lengths(clv[sv])), sep = "_")]
  thl <- tryCatch(setNames(lapply(names(FT), function(f) { s0 <- thl_sel[[f]]
    r <- if (s0$strat == "S1") fit_s1(Db[ftg == f], F[[s0$form]], s0$th) else fit_s2(Vb[ftg == f], F[[s0$form]], s0$th, maxit = 200)
    list(form = s0$form, th = r$th) }), names(FT)), error = function(e) NULL); if (is.null(thl)) return(NULL)
  lv <- tryCatch(level_fit(Vb, thl), error = function(e) NULL); if (is.null(lv)) return(NULL)
  allcells[, .(b = b, age = TT[-1], AGB = curve_of(lv, thl, ftg, if (section == "ALL") NULL else cell, TT[-1])), by = cell] }, mc.cores = NC)
BT <- rbindlist(Filter(Negate(is.null), bt))
band <- BT[, .(lo95 = quantile(AGB, .025), hi95 = quantile(AGB, .975), nboot = .N), by = .(cell, age)]
SC <- merge(SC, band, by = c("cell", "age"), all.x = TRUE)
SC <- merge(SC, sel[, .(ftg, strat, form)], by = "ftg")
setorder(SC, ftg, section, age)

## outputs
fwrite(P, file.path(YD, "shape_params_DATA.csv")); fwrite(S, file.path(YD, "form_screen_selection_DATA.csv")); fwrite(sel, file.path(YD, "selected_DATA.csv"))
fwrite(M_all, file.path(YD, "cv_interval_all_DATA.csv")); fwrite(M_long, file.path(YD, "cv_interval_long_old_DATA.csv"))
fwrite(M_h, file.path(YD, "cv_interval_by_horizon_DATA.csv")); fwrite(M_age, file.path(YD, "cv_interval_by_age_DATA.csv"))
fwrite(Gm, file.path(YD, "cv_generic_by_age_DATA.csv")); fwrite(Gall, file.path(YD, "cv_generic_all_DATA.csv")); fwrite(EQ, file.path(YD, "equivalence_dB_DATA.csv"))
fwrite(CT, file.path(YD, "statewide_curves_allcandidates_DATA.csv")); fwrite(SC, file.path(YD, "generic_curves_selected_DATA.csv"))
fwrite(data.table(ftg = names(lv_sel$fe), mu = lv_sel$fe, s_plot = lv_sel$s_plot, s_cell = lv_sel$s_cell, s_res = lv_sel$s_res), file.path(YD, "level_selected_DATA.csv"))
fwrite(data.table(cell = names(lv_sel$re), u_cell = lv_sel$re), file.path(YD, "level_selected_cell_effects_DATA.csv")); fwrite(cells, file.path(YD, "cells_DATA.csv"))

## ---- 5. figures
theme_crsf <- function() theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2),
  strip.background = element_rect(fill = "grey95", colour = NA), legend.position = "bottom")
lab_agb <- expression("AGB (Mg ha"^-1*")"); fcol <- c(CR = "#1B9E77", CRS = "#D95F02", CRD = "#7570B3", HG = "#E7298A")
gref <- rbindlist(lapply(names(FT), function(f) { d <- VV[ftg == f]; m <- gam(log(AGB) ~ s(STDAGE, k = 6), data = d, method = "REML")
  g <- data.table(STDAGE = seq(2, quantile(d$STDAGE, .995), 2)); p <- predict(m, g, se.fit = TRUE)
  data.table(ftg = f, age = g$STDAGE, fit = exp(p$fit), lo = exp(p$fit - 1.96 * p$se.fit), hi = exp(p$fit + 1.96 * p$se.fit)) }))[, lab := FT[ftg]]
ct <- copy(CT)[, lab := FT[ftg]][, strategy := fifelse(strat == "S1", "S1 remeasurement only", "S2 joint")]
mg <- mag[, .(lab = FT[ftg], age99, agb999)]
pA <- ggplot() + geom_point(data = VV[, .(age = STDAGE, AGB, lab = FT[ftg])], aes(age, AGB), colour = "grey78", size = 0.1, alpha = 0.3) +
  geom_ribbon(data = gref, aes(age, ymin = lo, ymax = hi), fill = "grey50", alpha = 0.35) +
  geom_line(data = ct, aes(age, AGB_med, colour = form, linetype = strategy), linewidth = 0.5) + scale_colour_manual(values = fcol, name = "Form") +
  scale_linetype_manual(values = c("S1 remeasurement only" = 2, "S2 joint" = 1), name = NULL) +
  geom_vline(data = mg, aes(xintercept = age99), linetype = 3, linewidth = 0.3) + geom_hline(data = mg, aes(yintercept = agb999), linetype = 3, linewidth = 0.3, colour = "#B2182B") +
  facet_wrap(~lab, ncol = 3) + coord_cartesian(xlim = c(0, 250), ylim = c(0, max(mag$agb999) * 1.4)) + labs(x = "Stand age (yr)", y = lab_agb,
  caption = "Statewide median curves for every form and strategy. Grey band: chronosequence GAM, 95% confidence band. Dotted vertical: observed 99th percentile stand age;\ndotted red: observed 99.9th percentile AGB (magnitude screen).") + theme_crsf() + guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1))
try_log(ggsave(file.path(YD, "me-fia_agb-curve-forms-strategies.png"), pA, width = 8, height = 6, dpi = 300, device = ragg::agg_png), "fig forms")
sc <- SC[section != "ALL"][, lab := FT[ftg]]
pB <- ggplot(sc, aes(age, AGB_med, group = cell)) + geom_ribbon(aes(ymin = lo95, ymax = hi95, fill = section), alpha = 0.13) +
  geom_line(data = sc[extrapolated == FALSE], aes(colour = section), linewidth = 0.55) +
  geom_line(data = sc[extrapolated == TRUE], aes(colour = section), linewidth = 0.4, linetype = 2) +
  facet_wrap(~lab, ncol = 3) + labs(x = "Stand age (yr)", y = lab_agb, colour = "Ecoregion section", fill = "Ecoregion section",
  caption = sprintf("Selected form per type (see table). Band: 95%% plot cluster bootstrap interval on the median curve (%d resamples, shape and level refit on one plot draw).\nDashed: beyond the cell's 99th percentile observed stand age (extrapolation).", uniqueN(BT$b))) +
  theme_crsf() + guides(colour = guide_legend(nrow = 1))
try_log(ggsave(file.path(YD, "me-fia_agb-generic-curves-by-section.png"), pB, width = 8, height = 6, dpi = 300, device = ragg::agg_png), "fig generic")
gm <- copy(Gm)[, lab := FT[ftg]][, cand := paste(strat, form)]
pC <- ggplot(gm, aes(ageclass, bias_mean, colour = form, shape = strat, group = cand)) + geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_point(position = position_dodge(0.5)) + geom_line(position = position_dodge(0.5), linewidth = 0.25) + scale_colour_manual(values = fcol) +
  facet_wrap(~lab, ncol = 3, scales = "free_y") + labs(x = "Stand age class (yr)", y = expression("Generic curve CV bias, mean AGB (Mg ha"^-1*")"), shape = "Strategy",
  caption = "Plot grouped 5 fold cross validation of the generic curves (level and shape refit without the held out plots), predicting held out plot visits from age, type and section.") + theme_crsf()
try_log(ggsave(file.path(YD, "me-fia_cv-generic-bias-by-age.png"), pC, width = 8, height = 6, dpi = 300, device = ragg::agg_png), "fig cv")
rs <- rbindlist(resid_tab)[paste(strat, form) %in% paste(sel$strat, sel$form)]
rs <- merge(rs, sel[, .(ftg, strat, form)], by = c("ftg", "strat", "form"))[, lab := FT[ftg]]
pD <- ggplot(rs, aes(STDAGE, r)) + geom_point(size = 0.1, alpha = 0.2, colour = "grey50") + geom_smooth(method = "gam", formula = y ~ s(x, k = 6), colour = "#B2182B", linewidth = 0.5) +
  geom_hline(yintercept = 0, linewidth = 0.3) + facet_wrap(~lab, ncol = 3) + labs(x = "Stand age (yr)", y = "Level model residual (log scale)",
  caption = "Selected curves. A flat smooth means the fitted shape and the chronosequence agree; trend means they disagree (GAM with 95% confidence band).") + theme_crsf()
try_log(ggsave(file.path(YD, "me-fia_level-residuals-by-age.png"), pD, width = 8, height = 5.5, dpi = 300, device = ragg::agg_png), "fig resid")

cat("\n== CV interval, long horizon old stands\n"); print(M_long); cat("\n== CV interval, all\n"); print(M_all)
cat("\n== generic curve CV, all\n"); print(Gall); cat("\n== equivalence\n"); print(EQ)
cat("\n== coverage by horizon\n"); print(dcast(M_h, ftg + strat + form ~ horizon, value.var = "cover90"))
cat(sprintf("\nlevel selected: s_plot %.3f s_cell %.3f s_res %.3f; bootstrap resamples kept: %d\n", lv_sel$s_plot, lv_sel$s_cell, lv_sel$s_res, uniqueN(BT$b)))
