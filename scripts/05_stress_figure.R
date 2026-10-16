## Maine FIA biomass yield curves, step 5: stress test figure (redrawn; the step 4 version failed on a vector linetype in geom_hline).
## Usage: Rscript 05_stress_figure.R <out_dir>
suppressMessages({ library(data.table); library(ggplot2) })
OD <- commandArgs(TRUE)[1]; SD <- file.path(OD, "stress")
FT <- c(SF = "Spruce-fir", MBB = "Maple-beech-birch", AB = "Aspen-birch", WRPH = "White-red pine-hemlock", OAK = "Oak and oak-pine", EAC = "Elm-ash-cottonwood")
ST <- fread(file.path(SD, "stress_shape_metrics_DATA.csv"))[form != "CR"]
pk <- fread(file.path(SD, "pick_DATA.csv")); ST <- merge(ST, pk[, .(ftg, pick = paste(strat, form))], by = "ftg")
ST[, `:=`(lab = FT[ftg], cand = paste(strat, form), chosen = paste(strat, form) == pick)]
ST[, variant := factor(variant, c("base", "ud", "agefix", "noyoung", "salvage"))]
theme_crsf <- function() theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2),
  strip.background = element_rect(fill = "grey95", colour = NA), legend.position = "bottom")
p <- ggplot(ST, aes(variant, r250, colour = form, shape = strat)) +
  geom_hline(yintercept = 1, linewidth = 0.3) + geom_hline(yintercept = 0.5, linewidth = 0.3, linetype = 3) +
  geom_point(position = position_dodge(0.6), size = 1.5, alpha = 0.8) +
  geom_point(data = ST[chosen == TRUE], position = position_dodge(0.6), size = 3.2, shape = 1, colour = "black", inherit.aes = TRUE) +
  facet_wrap(~lab, ncol = 3) + scale_colour_manual(values = c(CRS = "#D95F02", CRD = "#7570B3", HG = "#E7298A")) +
  scale_y_continuous(limits = c(0, 1.02)) +
  labs(x = "Data variant", y = "B(250) / curve maximum", shape = "Strategy", colour = "Form",
       caption = paste0("Decline-capable forms only. 1 = no decline by 250 yr; dotted = 0.5 floor screen; circled = form chosen for the type.\n",
                        "Variants: base (harvest free), ud (strict undisturbed), agefix (plot consistent age),\nnoyoung (age 15 yr and older), salvage (natural disturbance salvage kept as loss).")) +
  theme_crsf() + theme(plot.caption = element_text(hjust = 0), axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(SD, "me-fia_stress-decline-r250.png"), p, width = 8, height = 5.8, dpi = 300, device = ragg::agg_png)
cat("figure written\n")
