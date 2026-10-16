## Maine FIA biomass yield curves, step 6: benchmark against the prior Maine age based AGB curves (ME.AGB.R, supplied by Aaron on
## October 15, 2026). Prior form: AGB = b1 Age^b2 b3^Age (the Hoerl gamma, c = -ln b3) by EPA Level III ecoregion and FORTYPCD.
## EPA Level III is not in FIADB, so the plot's ECOSUBCD section is crosswalked [ASSUMPTION, approximate]: M211A -> 5.3.1;
## 211A, 211B, 211C, 211D -> 8.1.8; 221A -> 8.1.7. FORTYPCD values the prior code does not list use its default curve (flagged).
## Both curve sets are scored on the same harvest free plot visits (in sample for the new curves, which were fitted to these
## visits; the step 3 cross validated RMSE is reported beside them). Units: Mg ha-1 for both [ASSUMPTION for the prior code,
## stated "metric" by Aaron].
## Usage: Rscript 06_benchmark_prior.R <out_dir> <ME.AGB.R>
suppressMessages({ library(data.table); library(ggplot2) })
a <- commandArgs(TRUE); OD <- a[1]; source(a[2])
BD <- file.path(OD, "benchmark"); dir.create(BD, showWarnings = FALSE)
FT <- c(SF = "Spruce-fir", MBB = "Maple-beech-birch", AB = "Aspen-birch", WRPH = "White-red pine-hemlock", OAK = "Oak and oak-pine", EAC = "Elm-ash-cottonwood")
V <- fread(file.path(OD, "me_plot_visits_DATA.csv"))[ftg %in% names(FT) & !harvested & !is.na(STDAGE) & STDAGE >= 1 & AGB > 0]
V[, eco3 := fcase(section == "M211A", "5.3.1", section %in% c("211A", "211B", "211C", "211D"), "8.1.8", section == "221A", "8.1.7", default = NA_character_)]
## which ECO x FORTYPCD combinations the prior code lists (anything else falls to its default)
src <- readLines(a[2]); listed <- regmatches(src, regexpr("ECO=='[0-9.]+' & ForTyp==[0-9]+", src))
listed <- unique(sub("ECO=='([0-9.]+)' & ForTyp==([0-9]+)", "\\1_\\2", listed))
V[, prior_listed := paste(eco3, FORTYPCD, sep = "_") %in% listed]
V[, AGB_prior := mapply(function(e, f, t) ME.AGB.pred(e, f, t), eco3, FORTYPCD, STDAGE)]
FC <- fread(file.path(OD, "stress", "final_curves_capped_DATA.csv"))
V[, cell := paste(ftg, section, sep = ":")][, age_i := pmin(round(STDAGE), 250L)]
V <- merge(V, FC[, .(cell, age_i = age, AGB_new_med = AGB_med, AGB_new_mean = AGB_mean)], by = c("cell", "age_i"), all.x = TRUE)
V[, ageclass := cut(STDAGE, c(0, 40, 80, 120, 300), right = FALSE)]
met <- function(z) z[, .(n = .N, obs_mean = mean(AGB), prior_bias = mean(AGB_prior - AGB), prior_rmse = sqrt(mean((AGB_prior - AGB)^2)),
                       new_bias_mean = mean(AGB_new_mean - AGB), new_rmse_mean = sqrt(mean((AGB_new_mean - AGB)^2)),
                       new_bias_med = mean(AGB_new_med - AGB), share_prior_listed = mean(prior_listed))]
M1 <- V[, met(.SD), by = ftg]; M2 <- V[, met(.SD), by = .(ftg, ageclass)]; M0 <- met(V)[, ftg := "ALL"]
cv <- fread(file.path(OD, "curves", "cv_generic_all_DATA.csv")); pk <- fread(file.path(OD, "stress", "pick_DATA.csv"))
M1 <- merge(M1, merge(pk, cv, by = c("ftg", "strat", "form"))[, .(ftg, new_cv_rmse = rmse_mean)], by = "ftg", all.x = TRUE)
fwrite(rbind(M1, M0, fill = TRUE), file.path(BD, "benchmark_by_type_DATA.csv")); fwrite(M2, file.path(BD, "benchmark_by_type_age_DATA.csv"))
## prior curve shapes over 0 to 300 yr for the listed combinations in each of our type groups (biological screen)
pr <- rbindlist(lapply(listed, function(k) { e <- sub("_.*", "", k); f <- as.integer(sub(".*_", "", k)); t <- 1:300; y <- ME.AGB.pred(e, f, t)
  data.table(eco3 = e, FORTYPCD = f, peak_age = if (which.max(y) >= 300) NA_integer_ else which.max(y), B100 = y[100], B300 = y[300], r300 = y[300] / max(y)) }))
pr[, ftg := fcase(FORTYPCD %/% 10 * 10 == 120 | FORTYPCD == 167, "SF", FORTYPCD %/% 10 * 10 == 800, "MBB", FORTYPCD %/% 10 * 10 == 900, "AB",
                  FORTYPCD %/% 10 * 10 == 100, "WRPH", FORTYPCD %/% 100 %in% 4:5, "OAK", FORTYPCD %/% 100 == 7, "EAC", default = "OTH")]
fwrite(pr, file.path(BD, "prior_curve_shapes_DATA.csv"))
prs <- pr[, .(n_curves = .N, no_peak_300 = sum(is.na(peak_age)), collapse_r300_lt_0.25 = sum(r300 < 0.25, na.rm = TRUE),
              B100_min = min(B100), B100_max = max(B100)), by = ftg]
fwrite(prs, file.path(BD, "prior_curve_screen_DATA.csv"))
## figure: prior curves (thin, by FORTYPCD and EPA L3) against the new statewide type curve
theme_crsf <- function() theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2),
  strip.background = element_rect(fill = "grey95", colour = NA), legend.position = "bottom")
pc <- rbindlist(lapply(seq_len(nrow(pr)), function(i) data.table(ftg = pr$ftg[i], id = paste(pr$eco3[i], pr$FORTYPCD[i]), eco3 = pr$eco3[i], age = 1:250,
  AGB = ME.AGB.pred(pr$eco3[i], pr$FORTYPCD[i], 1:250))))[ftg %in% names(FT)][, lab := FT[ftg]]
nc <- FC[section == "ALL"][, lab := FT[ftg]]
bm <- V[, .(age = STDAGE, AGB, lab = FT[ftg])]
p <- ggplot() + geom_point(data = bm, aes(age, AGB), colour = "grey80", size = 0.1, alpha = 0.3) +
  geom_line(data = pc, aes(age, AGB, group = id, colour = eco3), linewidth = 0.3, alpha = 0.7) +
  geom_line(data = nc, aes(age, AGB_med), colour = "black", linewidth = 0.8) + geom_line(data = nc, aes(age, AGB_mean), colour = "black", linewidth = 0.5, linetype = 2) +
  facet_wrap(~lab, ncol = 3) + coord_cartesian(xlim = c(0, 250), ylim = c(0, 320)) + scale_colour_manual(values = c(`5.3.1` = "#1B9E77", `8.1.7` = "#D95F02", `8.1.8` = "#7570B3"), name = "Prior curve, EPA Level III") +
  labs(x = "Stand age (yr)", y = expression("AGB (Mg ha"^-1*")"),
       caption = "Colored: prior ME.AGB curves, one per FORTYPCD and EPA Level III ecoregion.\nBlack: new statewide curve (solid median, dashed mean). Grey: harvest free plot visits.") +
  theme_crsf() + theme(plot.caption = element_text(hjust = 0))
ggsave(file.path(BD, "me-fia_prior-vs-new-curves.png"), p, width = 8, height = 5.8, dpi = 300, device = ragg::agg_png)
cat("== by type\n"); print(rbind(M1, M0, fill = TRUE)); cat("\n== by type and age\n"); print(M2[order(ftg, ageclass)]); cat("\n== prior curve screen\n"); print(prs)
