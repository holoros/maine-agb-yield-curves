## Maine FIA biomass yield curves, step 4: stress test of decline-capable forms and final capped curves.
## Aaron's direction (October 15, 2026): prefer a form that allows a decline; stress test and keep what works best; cap curves at the
## observed ceiling. Decline-capable forms: CRS, CRD, HG (CR kept as the monotone reference, not eligible as primary).
## Stress tests (each refit for both strategies S1 and S2, every type):
##   base     harvest free intervals and visits, recorded stand age (as step 3)
##   ud       strict undisturbed: no recorded disturbance of any kind on the interval or at the visit
##   agefix   plot consistent age: age0 = median over visits of (STDAGE - elapsed years since first visit), age = age0 + elapsed
##   noyoung  intervals and visits with stand age < 15 yr removed
##   salvage  bound on natural decline: intervals with a harvest that coincides with recorded natural disturbance are kept as losses
## Transfer: leave one ecoregion section out, generic curve predicts the held out section's visits from age.
## Uncertainty of the decline itself: plot cluster bootstrap (B = nboot) of the base fit, share of draws with a peak before 250 yr
##   and the distribution of B(200)/max and B(250)/max.
## Scoring rule (written before any step 4 output was seen), per type, among decline-capable candidates (form x strategy) that pass
##   the step 3 screens (early monotone rise, peak at or after 40 yr, B(250) at least 0.5 of the maximum, statewide median curve at
##   or below the observed 99.9th percentile through 250 yr) on the base fit:
##   score = rank(generic curve CV RMSE, step 3) + rank(LOSO RMSE) + rank(stress spread), where stress spread is the range of
##   B(200)/B(100) across the five data variants. Lowest score wins; ties go to fewer parameters, then to S2.
##   If no decline-capable candidate passes for a type, the best passing CR is used and the type is flagged.
## Cap (Aaron's choice): every cell curve is held flat once it reaches the type's observed 99th percentile AGB (pmin).
## Usage: Rscript 04_stress_decline.R <out_dir> [nboot=200] [cores=28]
suppressMessages({ library(data.table); library(minpack.lm); library(lme4); library(ggplot2); library(parallel) })
a <- commandArgs(TRUE); OD <- a[1]; NB <- if (length(a) > 1) as.integer(a[2]) else 200L; NC <- if (length(a) > 2) as.integer(a[3]) else 28L
SD <- file.path(OD, "stress"); dir.create(SD, showWarnings = FALSE)
elog <- file.path(SD, "error_log.txt"); try_log <- function(expr, what) tryCatch(expr, error = function(e) {
  cat(sprintf("[%s] %s: %s\n", format(Sys.time()), what, conditionMessage(e)), file = elog, append = TRUE); NULL })
RNGkind("L'Ecuyer-CMRG"); set.seed(4040)
FT <- c(SF = "Spruce-fir", MBB = "Maple-beech-birch", AB = "Aspen-birch", WRPH = "White-red pine-hemlock", OAK = "Oak and oak-pine",
        EAC = "Elm-ash-cottonwood")
I <- fread(file.path(OD, "me_intervals_DATA.csv")); V <- fread(file.path(OD, "me_plot_visits_DATA.csv"))
mag <- V[ftg %in% names(FT) & !harvested, .(agb99 = quantile(AGB, .99), agb999 = quantile(AGB, .999), age99 = quantile(STDAGE, .99, na.rm = TRUE)), by = ftg]

## ---- functions (identical to step 3)
F <- list(
  CR  = list(np = 2, st = c(lk = log(0.03), lp = log(2)),
             lf = function(t, th) exp(th[2]) * log1p(-exp(-exp(th[1]) * t))),
  CRS = list(np = 3, st = c(lk = log(0.04), lp = log(2), ld = log(0.002)),
             lf = function(t, th) exp(th[2]) * log1p(-exp(-exp(th[1]) * t)) - exp(th[3]) * t),
  CRD = list(np = 4, st = c(lk = log(0.04), lp = log(2), lq = qlogis(0.3), ld = log(0.02)),
             lf = function(t, th) exp(th[2]) * log1p(-exp(-exp(th[1]) * t)) + log1p(-plogis(th[3]) * (1 - exp(-exp(th[4]) * t)))),
  HG  = list(np = 2, st = c(lb = log(1.5), lc = log(0.01)),
             lf = function(t, th) exp(th[1]) * log(t) - exp(th[2]) * t))
fit_s1 <- function(d, fm, st = NULL) {
  obj <- function(th) { r <- (d$AGB2 - d$AGB1 * exp(fm$lf(d$t2, th) - fm$lf(d$t1, th))) / (d$AGB1 * sqrt(d$dt)); r[!is.finite(r)] <- 1e3; r }
  best <- NULL; starts <- c(list(fm$st), if (!is.null(st)) list(st), lapply(1:3, function(i) fm$st + rnorm(length(fm$st), 0, 0.5)))
  for (s0 in starts) { r <- tryCatch(nls.lm(s0, fn = obj, control = nls.lm.control(maxiter = 500)), error = function(e) NULL)
    if (!is.null(r) && all(is.finite(r$par)) && (is.null(best) || r$deviance < best$deviance)) best <- r }
  if (is.null(best)) NULL else list(th = best$par) }
ctrl <- lmerControl(calc.derivs = FALSE, check.conv.singular = "ignore")
dev_s2 <- function(th, v, fm) { la <- log(v$AGB) - fm$lf(v$STDAGE, th); if (any(!is.finite(la))) return(1e10)
  dd <- data.frame(la = la, section = v$section, plot_ix = v$plot_ix)
  m <- tryCatch(if (length(unique(dd$section)) > 1) lmer(la ~ 1 + (1 | section) + (1 | plot_ix), dd, REML = FALSE, control = ctrl) else
    lmer(la ~ 1 + (1 | plot_ix), dd, REML = FALSE, control = ctrl), error = function(e) NULL)
  if (is.null(m)) 1e10 else -2 * as.numeric(logLik(m)) }
fit_s2 <- function(v, fm, st = NULL, maxit = 300) {
  starts <- c(list(fm$st), if (!is.null(st)) list(st)); best <- NULL
  for (s0 in starts) { r <- tryCatch(optim(s0, dev_s2, v = v, fm = fm, method = "Nelder-Mead", control = list(maxit = maxit, reltol = 1e-9)), error = function(e) NULL)
    if (!is.null(r) && (is.null(best) || r$value < best$value)) best <- r }
  if (is.null(best)) NULL else list(th = best$par) }
fit_any <- function(strat, d, v, m, st = NULL, maxit = 300) if (strat == "S1") fit_s1(d, F[[m]], st) else fit_s2(v, F[[m]], st, maxit)
level_fit <- function(vv, thl) {
  vv <- copy(vv); vv[, la := log(AGB) - F[[thl[[.BY$ftg]]$form]]$lf(STDAGE, thl[[.BY$ftg]]$th), by = ftg]
  m <- lmer(la ~ 0 + ftg + (1 | cell) + (1 | plot_ix), data = vv, REML = TRUE, control = ctrl)
  fe <- fixef(m); re <- ranef(m)$cell; vc <- as.data.frame(VarCorr(m))
  list(fe = setNames(fe, sub("^ftg", "", names(fe))), re = setNames(re[, 1], rownames(re)), s_plot = vc$sdcor[vc$grp == "plot_ix"],
       s_cell = vc$sdcor[vc$grp == "cell"], s_res = vc$sdcor[vc$grp == "Residual"]) }
curve_of <- function(lv, thl, f, ce, tt) { u <- if (!is.null(ce) && ce %in% names(lv$re)) lv$re[[ce]] else 0
  exp(lv$fe[[f]] + u + F[[thl[[f]]$form]]$lf(tt, thl[[f]]$th)) }
shape_metrics <- function(m, th) { tt <- 1:250; y <- exp(F[[m]]$lf(tt, th)); mx <- max(y); pk <- which.max(y)
  data.table(peak_age = pk, r100 = y[100] / mx, r200 = y[200] / mx, r250 = y[250] / mx, g200_100 = y[200] / y[100],
             early_ok = all(diff(y[1:pk]) >= -0.001 * mx) && y[1] < 0.05 * mx) }

## ---- data variants
Vb <- V[ftg %in% names(FT) & !harvested & AGB > 0.5 & !is.na(STDAGE) & STDAGE >= 1]
Vb[, cell := paste(ftg, section, sep = ":")]
mkD <- function(Ii) { d <- Ii[ftg %in% names(FT) & AGB1 >= 10 & AGB2 >= 0 & !is.na(STDAGE1) & dt >= 3]
  d[, `:=`(t1 = pmax(STDAGE1, 1), t2 = pmax(STDAGE1, 1) + dt)]; d }
## plot consistent age
V[, first_myr := min(myr), by = plot_ix]; ag <- V[!is.na(STDAGE), .(age0 = median(STDAGE - (myr - first_myr))), by = plot_ix]
V <- merge(V, ag, by = "plot_ix", all.x = TRUE)
## salvage intervals: harvest after visit i that coincides with recorded natural disturbance (codes < 80) at the same visit
setorder(V, plot_ix, myr); V[, k := seq_len(.N), by = .(plot_ix, run)]
nat_h <- V[harvested & dist_any & !dist_human, .(plot_ix, run, k)]
Ix <- copy(I); Ix[, salv := FALSE]
if (nrow(nat_h)) { Ix[, id := .I]; jj <- Ix[nat_h, on = .(plot_ix, run), allow.cartesian = TRUE, nomatch = 0][k > i & k <= j, unique(id)]
  Ix[id %in% jj, salv := TRUE]; Ix[, id := NULL] }
cat(sprintf("salvage intervals added: %d\n", Ix[salv & !hf, .N]))
variants <- list(
  base    = list(D = mkD(I[hf == TRUE & !harvested1]), V = Vb),
  ud      = list(D = mkD(I[ud == TRUE & !harvested1]), V = Vb[dist_any == FALSE]),
  agefix  = list(D = { d <- merge(I[hf == TRUE & !harvested1], ag, by = "plot_ix")
                       d <- merge(d, unique(V[, .(plot_ix, first_myr)]), by = "plot_ix")
                       d[, STDAGE1 := age0 + (myr1 - first_myr)][STDAGE1 >= 1]; mkD(d) },
                 V = { v <- merge(Vb, unique(V[, .(plot_ix, first_myr, age0)]), by = "plot_ix"); v[, STDAGE := age0 + (myr - first_myr)][STDAGE >= 1] }),
  noyoung = list(D = mkD(I[hf == TRUE & !harvested1 & STDAGE1 >= 15]), V = Vb[STDAGE >= 15]),
  salvage = list(D = mkD(Ix[(hf == TRUE | salv == TRUE) & !harvested1]), V = Vb))
print(rbindlist(lapply(names(variants), function(n) data.table(variant = n, intervals = nrow(variants[[n]]$D), visits = nrow(variants[[n]]$V)))))

## ---- 1. stress fits: variant x strategy x form x type (shape only), parallel
g1 <- CJ(variant = names(variants), strat = c("S1", "S2"), form = names(F), ftg = names(FT))
st1 <- mclapply(seq_len(nrow(g1)), function(i) { g <- g1[i]; vv <- variants[[g$variant]]
  r <- try_log(fit_any(g$strat, vv$D[ftg == g$ftg], vv$V[ftg == g$ftg], g$form), paste("stress", g$variant, g$strat, g$form, g$ftg))
  if (is.null(r)) NULL else cbind(g, shape_metrics(g$form, r$th), data.table(th = list(r$th))) }, mc.cores = NC)
ST <- rbindlist(Filter(Negate(is.null), st1), fill = TRUE)
spread <- ST[, .(spread_g200_100 = max(g200_100) - min(g200_100), peak_min = min(peak_age), peak_max = max(peak_age),
                 r250_min = min(r250), r250_max = max(r250), n_var = .N), by = .(strat, form, ftg)]

## ---- 2. base level fits and screens for every candidate (statewide median curve with all-type level model)
base <- ST[variant == "base"]
TT <- 1:250; cand <- unique(base[, .(strat, form)]); scr <- list(); lvs <- list()
for (i in seq_len(nrow(cand))) { g <- cand[i]
  thl <- setNames(lapply(names(FT), function(f) { x <- base[strat == g$strat & form == g$form & ftg == f]; if (!nrow(x)) NULL else list(form = g$form, th = x$th[[1]]) }), names(FT))
  if (any(sapply(thl, is.null))) next
  lv <- try_log(level_fit(Vb, thl), paste("level", g$strat, g$form)); if (is.null(lv)) next
  lvs[[paste(g$strat, g$form)]] <- list(lv = lv, thl = thl)
  for (f in names(FT)) { y <- curve_of(lv, thl, f, NULL, TT)
    scr[[length(scr) + 1]] <- data.table(strat = g$strat, form = g$form, ftg = f, max250 = max(y), B100 = y[100], B200 = y[200], B250 = y[250]) } }
SCR <- merge(merge(rbindlist(scr), base[, .(strat, form, ftg, peak_age, r250, early_ok)], by = c("strat", "form", "ftg")), mag, by = "ftg")
SCR[, `:=`(over999 = max250 > agb999, floor_ok = r250 >= 0.5, peak_ok = peak_age >= 40)][, pass := early_ok & !over999 & floor_ok & peak_ok]
SCR[, decline_within_250 := peak_age < 250]

## ---- 3. leave one section out transfer (generic curve predicting held out section visits)
secs <- sort(unique(Vb$section)); g3 <- CJ(strat = c("S1", "S2"), form = names(F), hold = secs)
lo <- mclapply(seq_len(nrow(g3)), function(i) { g <- g3[i]; Dt <- variants$base$D[section != g$hold]; Vt <- Vb[section != g$hold]
  thl <- tryCatch(setNames(lapply(names(FT), function(f) { st <- base[strat == g$strat & form == g$form & ftg == f]$th
    r <- fit_any(g$strat, Dt[ftg == f], Vt[ftg == f], g$form, if (length(st)) st[[1]] else NULL, maxit = 200); list(form = g$form, th = r$th) }), names(FT)),
    error = function(e) NULL); if (is.null(thl)) return(NULL)
  lv <- tryCatch(level_fit(Vt, thl), error = function(e) NULL); if (is.null(lv)) return(NULL)
  te <- Vb[section == g$hold]; sdp <- sqrt(lv$s_plot^2 + lv$s_res^2 + lv$s_cell^2)
  te[, med := curve_of(lv, thl, .BY$ftg, NULL, STDAGE), by = ftg]
  te[, .(strat = g$strat, form = g$form, hold = g$hold, ftg, STDAGE, AGB, mean = med * exp(sdp^2 / 2),
         in90 = AGB >= med * exp(-1.645 * sdp) & AGB <= med * exp(1.645 * sdp))] }, mc.cores = NC)
LO <- rbindlist(Filter(Negate(is.null), lo))
LOm <- LO[, .(n = .N, loso_bias = mean(mean - AGB), loso_rmse = sqrt(mean((mean - AGB)^2)), loso_cover90 = mean(in90)), by = .(strat, form, ftg)]
LOsec <- LO[, .(n = .N, bias = mean(mean - AGB), rmse = sqrt(mean((mean - AGB)^2)), cover90 = mean(in90)), by = .(strat, form, ftg, hold)]

## ---- 4. scoring
G3 <- fread(file.path(OD, "curves", "cv_generic_all_DATA.csv"))[, .(strat, form, ftg, cv_rmse = rmse_mean, cv_bias = bias_mean, cv_cover90 = cover90)]
SC <- Reduce(function(x, y) merge(x, y, by = c("strat", "form", "ftg"), all.x = TRUE), list(SCR, spread, LOm, G3))
SC[, np := sapply(form, function(m) F[[m]]$np)]
SC[, eligible := pass & form != "CR"]
SC[eligible == TRUE, score := frank(cv_rmse) + frank(loso_rmse) + frank(spread_g200_100), by = ftg]
pick <- SC[eligible == TRUE][order(ftg, score, np, -as.integer(strat == "S2"))][, .SD[1], by = ftg][, .(ftg, strat, form, basis = "decline-capable")]
miss <- setdiff(names(FT), pick$ftg)
if (length(miss)) pick <- rbind(pick, SC[ftg %in% miss & pass == TRUE][order(ftg, cv_rmse)][, .SD[1], by = ftg][, .(ftg, strat, form, basis = "CR fallback (no decline-capable form passed)")])
cat("\n== pick\n"); print(pick)

## ---- 5. bootstrap the decline for the picks (base data), shapes and level refit on one plot draw
thl_pick <- setNames(lapply(names(FT), function(f) { p <- pick[ftg == f]; list(form = p$form, strat = p$strat, th = base[strat == p$strat & form == p$form & ftg == f]$th[[1]]) }), names(FT))
lv_pick <- level_fit(Vb, thl_pick)
D0 <- variants$base$D; cl <- split(seq_len(nrow(D0)), D0$plot_ix); clv <- split(seq_len(nrow(Vb)), Vb$plot_ix); allp <- union(names(cl), names(clv))
cells <- Vb[, .(n_plots = uniqueN(plot_ix), age_q99 = quantile(STDAGE, .99)), by = .(ftg, section, cell)]
allcells <- rbind(cells, mag[, .(ftg, section = "ALL", cell = paste0(ftg, ":ALL"), n_plots = NA_integer_, age_q99 = age99)])
bt <- mclapply(seq_len(NB), function(b) { s <- sample(allp, length(allp), replace = TRUE)
  sd_ <- s[s %in% names(cl)]; Db <- D0[unlist(cl[sd_], use.names = FALSE)]
  sv <- s[s %in% names(clv)]; Vv <- Vb[unlist(clv[sv], use.names = FALSE)]; Vv[, plot_ix := paste(plot_ix, rep(seq_along(sv), lengths(clv[sv])), sep = "_")]
  thl <- tryCatch(setNames(lapply(names(FT), function(f) { p <- thl_pick[[f]]
    r <- if (p$strat == "S1") fit_s1(Db[ftg == f], F[[p$form]], p$th) else fit_s2(Vv[ftg == f], F[[p$form]], p$th, maxit = 200)
    list(form = p$form, th = r$th) }), names(FT)), error = function(e) NULL); if (is.null(thl)) return(NULL)
  lv <- tryCatch(level_fit(Vv, thl), error = function(e) NULL); if (is.null(lv)) return(NULL)
  list(shape = rbindlist(lapply(names(FT), function(f) cbind(data.table(b = b, ftg = f), shape_metrics(thl[[f]]$form, thl[[f]]$th)))),
       curves = allcells[, .(b = b, age = TT, AGB = curve_of(lv, thl, ftg, if (section == "ALL") NULL else cell, TT)), by = cell]) }, mc.cores = NC)
bt <- Filter(Negate(is.null), bt)
BS <- rbindlist(lapply(bt, `[[`, "shape")); BC <- rbindlist(lapply(bt, `[[`, "curves"))
## draws whose refit shape fails the early monotone screen (degenerate spikes near age 1) are dropped per type, and counted
BC[, ftg := sub(":.*", "", cell)]; BC <- merge(BC, BS[, .(b, ftg, early_ok)], by = c("b", "ftg"))
drop_ct <- BS[, .(draws = .N, dropped_early = sum(!early_ok)), by = ftg]; cat("\n== bootstrap draws dropped by the early screen\n"); print(drop_ct)
BC <- BC[early_ok == TRUE]; BS <- BS[early_ok == TRUE]
BSs <- BS[, .(nboot = .N, p_decline_250 = mean(peak_age < 250), peak_q05 = quantile(peak_age, .05), peak_q50 = quantile(peak_age, .5),
              peak_q95 = quantile(peak_age, .95), r250_q05 = quantile(r250, .05), r250_q50 = quantile(r250, .5), r250_q95 = quantile(r250, .95)), by = ftg]
cat("\n== bootstrap decline\n"); print(BSs)

## ---- 6. final curves, capped at the type's observed 99th percentile AGB
sdp <- sqrt(lv_pick$s_plot^2 + lv_pick$s_res^2)
FC <- allcells[, .(age = TT, AGB_med_raw = curve_of(lv_pick, thl_pick, ftg, if (section == "ALL") NULL else cell, TT)), by = .(ftg, section, cell, n_plots, age_q99)]
band <- BC[, .(lo95_raw = quantile(AGB, .025), hi95_raw = quantile(AGB, .975), nboot = .N), by = .(cell, age)]
FC <- merge(FC, band, by = c("cell", "age"), all.x = TRUE); FC <- merge(FC, mag[, .(ftg, cap = agb99)], by = "ftg")
FC <- merge(FC, pick[, .(ftg, strat, form, basis)], by = "ftg")
FC[, `:=`(AGB_med = pmin(AGB_med_raw, cap), lo95 = pmin(lo95_raw, cap), hi95 = pmin(hi95_raw, cap), capped = AGB_med_raw > cap, extrapolated = age > age_q99)]
FC[, `:=`(AGB_mean = pmin(AGB_med_raw * exp(sdp^2 / 2), cap), PI90_lo = pmin(AGB_med_raw * exp(-1.645 * sdp), cap), PI90_hi = pmin(AGB_med_raw * exp(1.645 * sdp), cap))]
z <- unique(FC[, .(ftg, section, cell, n_plots, age_q99, cap, strat, form, basis)])[, `:=`(age = 0L, AGB_med_raw = 0, lo95_raw = 0, hi95_raw = 0, nboot = NA_integer_,
  AGB_med = 0, lo95 = 0, hi95 = 0, capped = FALSE, extrapolated = FALSE, AGB_mean = 0, PI90_lo = 0, PI90_hi = 0)]
FC <- rbind(z, FC, use.names = TRUE); setorder(FC, ftg, section, age)

## outputs
fwrite(ST[, !"th"], file.path(SD, "stress_shape_metrics_DATA.csv")); fwrite(spread, file.path(SD, "stress_spread_DATA.csv"))
fwrite(SCR, file.path(SD, "base_screen_DATA.csv")); fwrite(LOm, file.path(SD, "loso_DATA.csv")); fwrite(LOsec, file.path(SD, "loso_by_section_DATA.csv"))
fwrite(SC, file.path(SD, "scoring_DATA.csv")); fwrite(pick, file.path(SD, "pick_DATA.csv")); fwrite(merge(BSs, drop_ct, by = "ftg"), file.path(SD, "bootstrap_decline_DATA.csv"))
fwrite(FC, file.path(SD, "final_curves_capped_DATA.csv"))
fwrite(rbindlist(lapply(names(FT), function(f) { p <- thl_pick[[f]]; th <- p$th
  nm <- switch(p$form, CR = c("k", "p"), CRS = c("k", "p", "d"), CRD = c("k", "p", "q", "d"), HG = c("b", "c"))
  tr <- switch(p$form, CR = exp(th), CRS = exp(th), CRD = c(exp(th[1:2]), plogis(th[3]), exp(th[4])), HG = exp(th))
  data.table(ftg = f, strat = p$strat, form = p$form, param = nm, est = tr, mu = lv_pick$fe[[f]]) })), file.path(SD, "final_params_DATA.csv"))
fwrite(data.table(cell = names(lv_pick$re), u_cell = lv_pick$re, s_plot = lv_pick$s_plot, s_cell = lv_pick$s_cell, s_res = lv_pick$s_res), file.path(SD, "final_level_DATA.csv"))

## figures
theme_crsf <- function() theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2),
  strip.background = element_rect(fill = "grey95", colour = NA), legend.position = "bottom")
lab_agb <- expression("AGB (Mg ha"^-1*")")
fc <- FC[section != "ALL"][, lab := paste0(FT[ftg], " (", strat, " ", form, ")")]
p1 <- ggplot(fc, aes(age, AGB_med, group = cell)) + geom_ribbon(aes(ymin = lo95, ymax = hi95, fill = section), alpha = 0.13) +
  geom_line(data = fc[extrapolated == FALSE], aes(colour = section), linewidth = 0.55) + geom_line(data = fc[extrapolated == TRUE], aes(colour = section), linewidth = 0.4, linetype = 2) +
  geom_hline(data = unique(fc[, .(lab, cap)]), aes(yintercept = cap), linetype = 3, linewidth = 0.3) +
  facet_wrap(~lab, ncol = 3) + coord_cartesian(ylim = c(0, max(mag$agb99) * 1.05)) +
  labs(x = "Stand age (yr)", y = lab_agb, colour = "Ecoregion section", fill = "Ecoregion section",
       caption = sprintf("Final curves: decline-capable form per type chosen by the stress test score; capped at the type's observed 99th percentile AGB (dotted).\nBand: 95%% plot cluster bootstrap interval on the median curve (%d resamples). Dashed: beyond the cell's 99th percentile observed stand age.", uniqueN(BS$b))) +
  theme_crsf() + guides(colour = guide_legend(nrow = 1))
try_log(ggsave(file.path(SD, "me-fia_agb-final-curves-by-section.png"), p1, width = 8, height = 6, dpi = 300, device = ragg::agg_png), "fig final")
## stress figure is drawn by 05_stress_figure.R

cat("\n== base screen\n"); print(SCR[order(ftg, strat, form), .(ftg, strat, form, peak_age, r250 = round(r250, 2), max250 = round(max250), agb999 = round(agb999), pass)])
cat("\n== scoring (eligible)\n"); print(SC[eligible == TRUE][order(ftg, score), .(ftg, strat, form, cv_rmse = round(cv_rmse, 1), loso_rmse = round(loso_rmse, 1), spread = round(spread_g200_100, 3), score)])
cat("\n== stress spread\n"); print(spread[order(ftg, strat, form)])
cat(sprintf("\nfinal level: s_plot %.3f s_cell %.3f s_res %.3f; bootstrap kept %d\n", lv_pick$s_plot, lv_pick$s_cell, lv_pick$s_res, uniqueN(BS$b)))
