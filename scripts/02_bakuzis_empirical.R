## Maine FIA biomass yield curves, step 2: empirical Bakuzis assessment of observed long-term behavior.
## Inputs: out/me_plot_visits_DATA.csv and out/me_intervals_DATA.csv from step 1 (anonymised plot index only).
## Chronosequence (binned plot visits) plus every harvest free remeasurement interval. Plot cluster bootstrap for all intervals.
## Usage: Rscript 02_bakuzis_empirical.R <out_dir> <bakuzis_check_step-aware.R>
suppressMessages({ library(data.table); library(ggplot2); library(patchwork); library(quantreg); library(mgcv) })
a <- commandArgs(TRUE); OD <- a[1]; BK <- a[2]
FD <- file.path(OD, "bakuzis"); dir.create(FD, showWarnings = FALSE)
set.seed(4747); NB <- 1000
elog <- file.path(FD, "error_log.txt"); try_log <- function(expr, what) tryCatch(expr, error = function(e) {
  cat(sprintf("[%s] %s: %s\n", format(Sys.time()), what, conditionMessage(e)), file = elog, append = TRUE); NULL })
V <- fread(file.path(OD, "me_plot_visits_DATA.csv")); I <- fread(file.path(OD, "me_intervals_DATA.csv"))
FT <- c(SF = "Spruce-fir", MBB = "Maple-beech-birch", AB = "Aspen-birch", WRPH = "White-red pine-hemlock", OAK = "Oak and oak-pine",
        EAC = "Elm-ash-cottonwood")
V <- V[ftg %in% names(FT) & !harvested & !is.na(STDAGE) & !is.na(site) & TPH > 0 & is.finite(SDI)]
V[, site := factor(site, c("low", "mid", "high"))]
V[, agebin := pmin(floor(STDAGE / 10) * 10 + 5, 175)]
theme_crsf <- function() theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2),
  strip.background = element_rect(fill = "grey95", colour = NA), legend.position = "bottom")
sitecol <- c(low = "#B2182B", mid = "#888888", high = "#2166AC")
ftlab <- function(x) FT[x]

## magnitude reference (observed distribution)
mag <- V[, .(n = .N, AGB_q50 = median(AGB), AGB_q99 = quantile(AGB, .99), AGB_q999 = quantile(AGB, .999),
             BA_q99 = quantile(BAPH, .99), BA_q999 = quantile(BAPH, .999), SDI_q99 = quantile(SDI, .99), SDI_q999 = quantile(SDI, .999),
             age_q99 = quantile(STDAGE, .99), age_max = max(STDAGE)), by = ftg]
mag <- rbind(mag, V[, .(ftg = "ALL", n = .N, AGB_q50 = median(AGB), AGB_q99 = quantile(AGB, .99), AGB_q999 = quantile(AGB, .999),
  BA_q99 = quantile(BAPH, .99), BA_q999 = quantile(BAPH, .999), SDI_q99 = quantile(SDI, .99), SDI_q999 = quantile(SDI, .999),
  age_q99 = quantile(STDAGE, .99), age_max = max(STDAGE))])
fwrite(mag, file.path(FD, "magnitude_reference_DATA.csv"))

## binned trajectories (median per type x site x 10 yr age bin, n >= 10) and the step aware check
B <- V[, .(n = .N, HT = median(HT, na.rm = TRUE), TPH = median(TPH), BAPH = median(BAPH), QMD = median(QMD), VOL = median(VOL),
           AGB = median(AGB)), by = .(ftg, site, agebin)][n >= 10][order(ftg, site, agebin)]
fwrite(B, file.path(FD, "binned_trajectories_DATA.csv"))
bk <- B[, .(origin = ftg, site = as.integer(site), age = agebin, HT, TPH, BAPH, QMD, VOL)]
fwrite(bk, file.path(FD, "bakuzis_input_binned.csv"))
bkout <- try_log(system2("Rscript", c("--vanilla", BK, file.path(FD, "bakuzis_input_binned.csv"), "0.40", "0.10"), stdout = TRUE, stderr = TRUE), "bakuzis_check")
if (is.null(bkout)) bkout <- "bakuzis_check failed to run"
if (!is.null(attr(bkout, "status"))) bkout <- c(sprintf("bakuzis_check EXIT STATUS %s", attr(bkout, "status")), bkout)
writeLines(bkout, file.path(FD, "bakuzis_check_binned_REPORT.txt"))

## plot cluster bootstrap helper: resample whole plots
cl_boot <- function(dt, fun, nb = NB) { ids <- unique(dt$plot_ix); sp <- split(seq_len(nrow(dt)), dt$plot_ix)
  replicate(nb, { s <- sample(as.character(ids), length(ids), replace = TRUE)
    tryCatch(fun(dt[unlist(sp[s], use.names = FALSE)]), error = function(e) NA_real_) }) }

## [4] Reineke: upper boundary of log TPH on log QMD, tau 0.95 and 0.99, stems >= 12.7 cm (primary) and >= 2.54 cm
rein <- rbindlist(lapply(names(FT), function(f) rbindlist(lapply(c("ge12.7cm", "ge2.54cm"), function(basis) {
  d <- if (basis == "ge12.7cm") V[ftg == f & TPH13 > 0 & QMD13 > 12.7, .(plot_ix, lt = log(TPH13), lq = log(QMD13))] else
       V[ftg == f & TPH > 0, .(plot_ix, lt = log(TPH), lq = log(QMD))]
  if (nrow(d) < 100) return(NULL)
  rbindlist(lapply(c(0.95, 0.99), function(tau) {
    fit <- try_log(rq(lt ~ lq, tau = tau, data = d), "rq"); if (is.null(fit)) return(NULL)
    bs <- cl_boot(d, function(x) coef(rq(lt ~ lq, tau = tau, data = x))[2], nb = 400)
    data.table(ftg = f, basis = basis, tau = tau, n = nrow(d), slope = coef(fit)[2], lo = quantile(bs, .025), hi = quantile(bs, .975))
  })) }))))
## test against Reineke -1.605 with a stated tolerance of 0.30: FLAG when the 95% interval lies wholly outside [-1.905, -1.305]
rein[, verdict := fifelse(hi < -1.905 | lo > -1.305, "FLAG", "consistent within 0.30 of -1.605")]
fwrite(rein, file.path(FD, "reineke_boundary_DATA.csv"))

## [8] Eichhorn: log VOL on log HT with site tier; ratio high/low at median HT and CV across tiers
eich <- rbindlist(lapply(names(FT), function(f) try_log({
  d <- V[ftg == f & VOL > 1 & HT > 5][, site := droplevels(site)]; if (uniqueN(d$site) < 2 || nrow(d) < 60) return(NULL)
  m0 <- lm(log(VOL) ~ log(HT), d); m1 <- lm(log(VOL) ~ log(HT) * site, d)
  h <- median(d$HT); nd <- data.table(HT = h, site = factor(levels(d$site), levels(d$site)))
  pv <- exp(predict(m1, nd)); cv <- sd(pv) / mean(pv)
  bs <- cl_boot(d, function(x) { x <- copy(x)[, site := factor(as.character(site), levels(d$site))]; m <- lm(log(VOL) ~ log(HT) * site, x); p <- exp(predict(m, nd)); sd(p) / mean(p) }, nb = 400)
  data.table(ftg = f, n = nrow(d), r2_ht = summary(m0)$r.squared, r2_full = summary(m1)$r.squared, HT_med = h,
             cv_site = cv, cv_lo = quantile(bs, .025, na.rm = TRUE), cv_hi = quantile(bs, .975, na.rm = TRUE), verdict = fifelse(cv < 0.15, "PASS", "FLAG"))
}, paste("eichhorn", f))))
fwrite(eich, file.path(FD, "eichhorn_DATA.csv"))

## [1/7] site ordering of median HT and AGB by age bin, with plot bootstrap of the bin medians
so <- rbindlist(lapply(names(FT), function(f) rbindlist(lapply(c("HT", "AGB"), function(v) {
  d <- V[ftg == f, .(plot_ix, site, agebin, y = get(v))][!is.na(y)]
  bins <- d[, .N, by = .(agebin, site)][N >= 10][, .N, by = agebin][N == 3]$agebin
  rbindlist(lapply(bins, function(b) { x <- d[agebin == b]
    med <- x[, .(m = median(y)), by = site][order(site)]$m
    bs <- cl_boot(x, function(z) { r <- z[, .(m = median(y)), by = site]; setkey(r, site); r[c("low", "mid", "high")]$m }, nb = 300)
    lo <- apply(bs, 1, quantile, .05, na.rm = TRUE); hi <- apply(bs, 1, quantile, .95, na.rm = TRUE)
    inv <- is.unsorted(med); sig <- isTRUE(lo[1] > hi[2]) || isTRUE(lo[2] > hi[3]) || isTRUE(lo[1] > hi[3])
    data.table(ftg = f, var = v, agebin = b, low = med[1], mid = med[2], high = med[3], inverted = inv, inverted_beyond_90ci = sig)
  })) }))))
fwrite(so, file.path(FD, "site_ordering_DATA.csv"))

## [Sukachev] net stem change rate of trees >= 12.7 cm on harvest free consecutive intervals against relative density and site
H <- I[hf == TRUE & ftg %in% names(FT) & !is.na(site)]
## SDI and RD from stems >= 12.7 cm; SDImax as the 99th percentile of plot SDI13 by type (stated deviation from the brms plot level SDImax)
V[, SDI13 := fifelse(TPH13 > 0, TPH13 * (QMD13 / 25.4)^1.605, 0)]
sdimax <- V[, .(sdimax = quantile(SDI13, .99, na.rm = TRUE)), by = ftg]
H <- merge(H, sdimax, by = "ftg"); H[, QMD131 := sqrt(BA131 / pmax(TPH131, 1e-9) / 0.00007854)][, RD1 := fifelse(TPH131 > 0, TPH131 * (QMD131 / 25.4)^1.605, 0) / sdimax]
H[, site := factor(site, c("low", "mid", "high"))]
Hc <- H[consecutive == TRUE & TPH131 > 50]
Hc[, rstem := log(pmax(TPH132, 1) / TPH131) / dt]
suk <- rbindlist(lapply(names(FT), function(f) { d <- Hc[ftg == f]; if (nrow(d) < 80) return(NULL)
  m <- lm(rstem ~ RD1 + site, d); co <- summary(m)$coefficients
  bs <- cl_boot(d, function(x) coef(lm(rstem ~ RD1 + site, x))["RD1"], nb = 400)
  data.table(ftg = f, n = nrow(d), b_RD = co["RD1", 1], b_RD_lo = quantile(bs, .025, na.rm = TRUE), b_RD_hi = quantile(bs, .975, na.rm = TRUE),
             b_high_vs_low = if ("sitehigh" %in% rownames(co)) co["sitehigh", 1] else NA_real_) }))
suk[, verdict := fifelse(b_RD_hi < 0, "PASS (denser stands lose stems faster)", "FLAG or weak")]
fwrite(suk, file.path(FD, "sukachev_DATA.csv"))

## [yield density] AGB against TPH within age bins 35 to 85, Spearman
yd <- V[agebin %between% c(35, 85), .(n = .N, rho = cor(AGB, TPH, method = "spearman")), by = .(ftg, agebin)][n >= 30]
fwrite(yd, file.path(FD, "yield_density_DATA.csv"))

## LONG TERM BEHAVIOR: net AGB change on harvest free intervals (primary) and strict undisturbed (sensitivity) by age class
D <- I[ftg %in% names(FT) & hf == TRUE & !harvested1 & AGB1 > 5 & !is.na(STDAGE1)]
D[, `:=`(mid_age = STDAGE1 + dt / 2, dB = (AGB2 - AGB1) / dt, rgr = log(pmax(AGB2, 0.5) / AGB1) / dt)]
D[, ageclass := cut(mid_age, c(0, 20, 40, 60, 80, 100, 120, 250), right = FALSE)]
D[, horizon := fifelse(dt <= 7, "5 yr (one cycle)", fifelse(dt <= 12, "10 yr", "15 to 25 yr"))]
dec <- rbindlist(lapply(c("hf", "ud"), function(scr) { d <- if (scr == "hf") copy(D) else D[ud == TRUE]; d[, screen := scr]
  d[, { bs <- if (uniqueN(plot_ix) >= 15) cl_boot(.SD, function(z) c(mean(z$dB), mean(z$dB < 0)), nb = 500) else matrix(NA_real_, 2, 1)
        .(n = .N, plots = uniqueN(plot_ix), dB_mean = mean(dB), dB_lo = quantile(bs[1, ], .025, na.rm = TRUE), dB_hi = quantile(bs[1, ], .975, na.rm = TRUE),
          share_loss = mean(dB < 0), share_lo = quantile(bs[2, ], .025, na.rm = TRUE), share_hi = quantile(bs[2, ], .975, na.rm = TRUE),
          rgr_med = median(rgr)) }, by = .(screen, ftg, ageclass)] }))
setorder(dec, screen, ftg, ageclass); fwrite(dec, file.path(FD, "net_change_by_age_DATA.csv"))
dech <- D[, { bs <- if (uniqueN(plot_ix) >= 15) cl_boot(.SD, function(z) mean(z$dB), nb = 300) else NA_real_
  .(n = .N, dB_mean = mean(dB), lo = quantile(bs, .025, na.rm = TRUE), hi = quantile(bs, .975, na.rm = TRUE)) }, by = .(horizon, ageclass)]
setorder(dech, horizon, ageclass); fwrite(dech, file.path(FD, "net_change_by_age_horizon_DATA.csv"))

## GAM of annual net AGB change on mid age by type, plot random effect, consecutive and long intervals together, weights dt
gfits <- rbindlist(lapply(names(FT), function(f) { d <- D[ftg == f]; if (nrow(d) < 150) return(NULL)
  d[, pf := factor(plot_ix)]
  m <- try_log(gam(dB ~ s(mid_age, k = 6) + s(pf, bs = "re"), data = d, weights = dt / mean(dt), method = "REML"), paste("gam", f)); if (is.null(m)) return(NULL)
  g <- data.table(mid_age = seq(max(5, quantile(d$mid_age, .01)), quantile(d$mid_age, .995), length.out = 80), pf = d$pf[1])
  p <- predict(m, g, se.fit = TRUE, exclude = "s(pf)")
  data.table(ftg = f, mid_age = g$mid_age, fit = p$fit, lo = p$fit - 1.96 * p$se.fit, hi = p$fit + 1.96 * p$se.fit,
             age_zero_cross = { z <- g$mid_age[which(diff(sign(p$fit)) < 0)]; if (length(z)) z[1] else NA_real_ }) }))
fwrite(gfits, file.path(FD, "gam_net_change_DATA.csv"))

## figures -------------------------------------------------------------------------------------------------------------
lab_agb <- expression("AGB (Mg ha"^-1*")"); lab_dB <- expression("Net AGB change (Mg ha"^-1*" yr"^-1*")")
lab_ba <- expression("Basal area (m"^2*" ha"^-1*")"); lab_tph <- expression("Stems (ha"^-1*")"); lab_vol <- expression("Gross volume (m"^3*" ha"^-1*")")
labs_v <- list(age = "Stand age (yr)", HT = "Top height (m)", TPH = lab_tph, BAPH = lab_ba, QMD = "QMD (cm)", VOL = lab_vol, AGB = lab_agb)
mat_fig <- function(f) {
  vars <- c("age", "HT", "TPH", "BAPH", "QMD", "VOL", "AGB"); pl <- list()
  pts <- V[ftg == f][, age := STDAGE]; bb <- B[ftg == f][, age := agebin]
  for (r in 2:length(vars)) for (cc in 1:(length(vars) - 1)) {
    if (cc >= r) { pl[[length(pl) + 1]] <- plot_spacer(); next }
    x <- vars[cc]; y <- vars[r]
    g <- ggplot() + geom_point(data = pts, aes(.data[[x]], .data[[y]]), colour = "grey70", size = 0.15, alpha = 0.35) +
      geom_path(data = bb, aes(.data[[x]], .data[[y]], colour = site), linewidth = 0.6) + scale_colour_manual(values = sitecol, name = "Site tier") +
      labs(x = if (r == length(vars)) labs_v[[x]] else NULL, y = if (cc == 1) labs_v[[y]] else NULL) + theme_crsf() + theme(axis.text = element_text(size = 5))
    if (y == "TPH" || x == "TPH") g <- g + scale_y_continuous(trans = if (y == "TPH") "log10" else "identity")
    pl[[length(pl) + 1]] <- g }
  wrap_plots(pl, ncol = length(vars) - 1, guides = "collect") + plot_annotation(title = paste0(ftlab(f), ": observed Bakuzis matrix (Maine FIA, harvest free visits)"),
    caption = "Grey: plot visits. Lines: median by site tier and 10 yr stand age bin (n >= 10).") & theme(legend.position = "bottom")
}
for (f in names(FT)) try_log(ggsave(file.path(FD, sprintf("me-fia_bakuzis-matrix_%s.png", f)), mat_fig(f), width = 10, height = 9, dpi = 300, device = ragg::agg_png), paste("matrix", f))

gd <- copy(gfits)[, lab := ftlab(ftg)]; dd <- copy(dec[screen == "hf"])[, lab := ftlab(ftg)]
dd[, amid := c(10, 30, 50, 70, 90, 110, 150)[as.integer(ageclass)]]
p1 <- ggplot(gd, aes(mid_age, fit)) + geom_hline(yintercept = 0, linewidth = 0.3) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey80") +
  geom_line(linewidth = 0.6) + geom_pointrange(data = dd[n >= 15], aes(amid, dB_mean, ymin = dB_lo, ymax = dB_hi), colour = "#2166AC", size = 0.15) +
  facet_wrap(~lab, ncol = 3) + labs(x = "Stand age at interval midpoint (yr)", y = lab_dB,
  caption = "Line and band: GAM with plot random effect, 95% confidence band on the mean. Points: age class means with 95% plot cluster bootstrap intervals.\nAll harvest free remeasurement intervals (consecutive and longer).") + theme_crsf()
try_log(ggsave(file.path(FD, "me-fia_net-agb-change-by-age.png"), p1, width = 7.5, height = 5, dpi = 300, device = ragg::agg_png), "fig net change")
p2 <- ggplot(dd[n >= 15], aes(amid, share_loss)) + geom_pointrange(aes(ymin = share_lo, ymax = share_hi), size = 0.2) + geom_line(linewidth = 0.3) +
  facet_wrap(~lab, ncol = 3) + scale_y_continuous(labels = scales::percent) + labs(x = "Stand age at interval midpoint (yr)", y = "Intervals with net AGB loss",
  caption = "95% plot cluster bootstrap intervals; harvest free intervals.") + theme_crsf()
try_log(ggsave(file.path(FD, "me-fia_share-net-loss-by-age.png"), p2, width = 7.5, height = 5, dpi = 300, device = ragg::agg_png), "fig share")
rr <- V[TPH13 > 0 & QMD13 > 12.7][, lab := ftlab(ftg)]; rl <- rein[basis == "ge12.7cm" & tau == 0.99][, lab := ftlab(ftg)]
ic <- rbindlist(lapply(rl$ftg, function(f) { d <- V[ftg == f & TPH13 > 0 & QMD13 > 12.7]; co <- coef(rq(log(TPH13) ~ log(QMD13), tau = 0.99, data = d))
  data.table(ftg = f, lab = ftlab(f), q = exp(seq(log(13), log(max(d$QMD13)), length.out = 30)), b0 = co[1], b1 = co[2]) }))[, t := exp(b0 + b1 * log(q))]
p3 <- ggplot(rr, aes(QMD13, TPH13)) + geom_point(size = 0.2, alpha = 0.25, colour = "grey50") + geom_line(data = ic, aes(q, t), colour = "#B2182B") +
  scale_x_log10() + scale_y_log10() + facet_wrap(~lab, ncol = 3) + labs(x = "QMD of stems >= 12.7 cm (cm)", y = expression("Stems >= 12.7 cm (ha"^-1*")"),
  caption = "Red: 99th percentile quantile regression boundary.") + theme_crsf()
try_log(ggsave(file.path(FD, "me-fia_reineke-boundary.png"), p3, width = 7.5, height = 5, dpi = 300, device = ragg::agg_png), "fig reineke")

## summary to log
cat("\n== magnitude\n"); print(mag); cat("\n== reineke\n"); print(rein); cat("\n== eichhorn\n"); print(eich)
cat("\n== site ordering inversions\n"); print(so[, .(bins = .N, inverted = sum(inverted), beyond_ci = sum(inverted_beyond_90ci)), by = .(ftg, var)])
cat("\n== sukachev\n"); print(suk); cat("\n== yield density\n"); print(yd)
cat("\n== net change by age (hf)\n"); print(dec[screen == "hf"]); cat("\n== net change by age (ud)\n"); print(dec[screen == "ud"])
cat("\n== by horizon\n"); print(dech); cat("\n== GAM zero crossing age\n"); print(unique(gfits[, .(ftg, age_zero_cross)]))
cat("\n== bakuzis_check (binned)\n"); cat(bkout, sep = "\n")
