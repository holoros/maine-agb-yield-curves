## Maine FIA biomass yield curves, step 7: graphical abstract (lead file of the Zenodo record).
## Three panels from deposited aggregate tables: (a) observed annual net AGB change by age with 95% plot cluster bootstrap intervals,
## (b) final capped statewide median curves with 95% bootstrap bands and extrapolation dashed, (c) bootstrap probability that each curve
## peaks before 250 yr. The probability panel is a bootstrap proportion over 200 resamples; no further interval is drawn on it.
## Usage: Rscript 07_graphical_abstract.R <out_dir> <png> [doi]
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
a <- commandArgs(TRUE); OD <- a[1]; PNG <- a[2]; DOI <- if (length(a) > 2) a[3] else ""
FT <- c(SF = "Spruce-fir", MBB = "Maple-beech-birch", AB = "Aspen-birch", WRPH = "White-red pine-hemlock", OAK = "Oak and oak-pine", EAC = "Elm-ash-cottonwood")
col <- c(SF = "#1a3d28", MBB = "#c5a55a", AB = "#7f9c6b", WRPH = "#3b6e8f", OAK = "#8c5a3c", EAC = "#9a9a9a")
th <- theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2),
  plot.title = element_text(face = "bold", size = 9), legend.position = "none")
dec <- fread(file.path(OD, "bakuzis", "net_change_by_age_DATA.csv"))[screen == "hf" & n >= 15 & !is.na(dB_lo)]
dec[, amid := c(10, 30, 50, 70, 90, 110, 150)[as.integer(factor(ageclass, levels = c("[0,20)", "[20,40)", "[40,60)", "[60,80)", "[80,100)", "[100,120)", "[120,250)")))]]
pA <- ggplot(dec, aes(amid, dB_mean, colour = ftg)) + geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_pointrange(aes(ymin = dB_lo, ymax = dB_hi), position = position_dodge(6), size = 0.25) + geom_line(position = position_dodge(6), linewidth = 0.4) +
  scale_colour_manual(values = col) + labs(title = "a. Observed net change", x = "Stand age (yr)", y = expression("Net AGB change (Mg ha"^-1*" yr"^-1*")")) + th
fc <- fread(file.path(OD, "stress", "final_curves_capped_DATA.csv"))[section == "ALL"]
pB <- ggplot(fc, aes(age, AGB_med, colour = ftg, fill = ftg)) + geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.12, colour = NA) +
  geom_line(data = fc[extrapolated == FALSE], linewidth = 0.8) + geom_line(data = fc[extrapolated == TRUE], linewidth = 0.6, linetype = 2) +
  scale_colour_manual(values = col, labels = FT) + scale_fill_manual(values = col, labels = FT) +
  labs(title = "b. Generic curves, 95% band", x = "Stand age (yr)", y = expression("AGB (Mg ha"^-1*")"), colour = NULL, fill = NULL) +
  th + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))
bs <- fread(file.path(OD, "stress", "bootstrap_decline_DATA.csv"))[, lab := FT[ftg]]
pC <- ggplot(bs, aes(reorder(lab, p_decline_250), p_decline_250, fill = ftg)) + geom_col(width = 0.65) + coord_flip(ylim = c(0, 1.15)) +
  geom_text(aes(label = sprintf("%.2f", p_decline_250)), hjust = -0.15, size = 2.6) + scale_fill_manual(values = col) +
  labs(title = "c. P(peak before 250 yr)", x = NULL, y = "Share of bootstrap refits") + scale_y_continuous(breaks = c(0, 0.5, 1)) + th
cap <- paste0("Maine FIA annual inventory 1999 to 2024, 15,500 harvest free remeasurement intervals on 2,383 plots.\nPanel a: 95% plot cluster bootstrap intervals. Panel b: dashed beyond the 99th percentile observed age; curves capped at the 99th percentile AGB.",
              if (nzchar(DOI)) paste0("\nData and code: https://doi.org/", DOI) else "")
p <- (pA | pB | pC) + plot_layout(widths = c(1, 1.15, 0.9)) +
  plot_annotation(title = "Maine forest types: long-term aboveground biomass behavior and generic yield curves",
                  caption = cap, theme = theme(plot.title = element_text(face = "bold", size = 12), plot.caption = element_text(hjust = 0, size = 7)))
ggsave(PNG, p, width = 8.4, height = 4.2, dpi = 300, device = ragg::agg_png, bg = "white")
cat("written", PNG, "\n")
