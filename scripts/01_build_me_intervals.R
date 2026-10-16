## Maine FIA biomass yield curves, step 1: plot visit metrics and every valid remeasurement interval.
## Annual design only (INVYR >= 1999, DESIGNCD 1, 115, 116, 117). Single forested condition plots.
## Two interval screens: harvest free (primary; natural disturbance mortality kept as mortality) and strict
## undisturbed (sensitivity). LAT and LON are never read. Outputs carry an anonymised plot index only.
## Usage: Rscript 01_build_me_intervals.R <fia_dir> <out_dir>
suppressMessages(library(data.table))
a <- commandArgs(TRUE); FD <- a[1]; OD <- a[2]
dir.create(OD, showWarnings = FALSE, recursive = TRUE)
elog <- file.path(OD, "error_log.txt")
say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = file.path(OD, "build_log.txt"), append = TRUE) }
set.seed(20261015)

P <- fread(file.path(FD, "ME_PLOT.csv"), select = c("CN","PREV_PLT_CN","INVYR","MEASYEAR","MEASMON","DESIGNCD","KINDCD",
  "PLOT_STATUS_CD","REMPER","STATECD","UNITCD","COUNTYCD","PLOT"), colClasses = list(character = c("CN","PREV_PLT_CN")))
G <- fread(file.path(FD, "ME_PLOTGEOM.csv"), select = c("CN","ECOSUBCD"), colClasses = list(character = "CN"))
P <- merge(P, G, by = "CN", all.x = TRUE)
P <- P[INVYR >= 1999 & DESIGNCD %in% c(1, 115, 116, 117)]
P[, pid := paste(STATECD, UNITCD, COUNTYCD, PLOT, sep = "_")]
P[, eco := gsub(" ", "", ECOSUBCD)][, section := sub("[a-z]+$", "", eco)]
P[, myr := MEASYEAR + (fifelse(is.na(MEASMON), 7L, MEASMON) - 0.5) / 12]
say("annual plot visits: %d on %d plots", nrow(P), uniqueN(P$pid))

C <- fread(file.path(FD, "ME_COND.csv"), select = c("PLT_CN","CONDID","COND_STATUS_CD","CONDPROP_UNADJ","FORTYPCD","STDAGE",
  "SITECLCD","DSTRBCD1","DSTRBCD2","DSTRBCD3","TRTCD1","TRTCD2","TRTCD3","STDORGCD"), colClasses = list(character = "PLT_CN"))
C <- C[PLT_CN %in% P$CN]
## single condition forested visit: one condition, forested, proportion ~1
nc <- C[, .(ncond = .N, nfor = sum(COND_STATUS_CD == 1), propfor = sum(CONDPROP_UNADJ[COND_STATUS_CD == 1])), by = PLT_CN]
C1 <- C[PLT_CN %in% nc[ncond == 1 & nfor == 1 & propfor > 0.999]$PLT_CN]
dist_codes <- function(...) { x <- cbind(...); x[is.na(x)] <- 0; x }
C1[, trt := rowSums(dist_codes(TRTCD1, TRTCD2, TRTCD3) > 0) > 0]
## natural disturbance codes kept in the harvest free screen; 80 (human), 90+ (geologic, other) and 30s (fire) are recorded
d <- dist_codes(DSTRBCD1 = C1$DSTRBCD1, DSTRBCD2 = C1$DSTRBCD2, DSTRBCD3 = C1$DSTRBCD3)
C1[, dist_any := rowSums(d > 0) > 0][, dist_human := rowSums(d >= 80) > 0]
C1[, g10 := (FORTYPCD %/% 10) * 10]
C1[, ftg := fcase(g10 == 120, "SF", g10 == 800, "MBB", g10 == 900, "AB", g10 == 100, "WRPH",
  g10 %in% c(400, 500, 510, 520), "OAK", g10 == 700, "EAC", default = "OTH")]
C1[, site := fcase(SITECLCD %in% 1:4, "high", SITECLCD == 5, "mid", SITECLCD %in% 6:7, "low", default = NA_character_)]
say("single condition forested visits: %d", nrow(C1))

## trees: live >= 2.54 cm (1 in), cut trees flag removals
Tcols <- c("PLT_CN","STATUSCD","DIA","HT","TPA_UNADJ","DRYBIO_AG","VOLCFGRS","AGENTCD")
Tr <- fread(file.path(FD, "ME_TREE.csv"), select = Tcols, colClasses = list(character = "PLT_CN"))
Tr <- Tr[PLT_CN %in% C1$PLT_CN]
say("tree records on retained visits: %d", nrow(Tr))
## removals plus trees killed by silvicultural or land clearing activity and left on site (STATUSCD 2, AGENTCD 80)
cut_by <- Tr[STATUSCD == 3 | (STATUSCD == 2 & AGENTCD %in% 80:89), .(n_cut = .N), by = PLT_CN]
L <- Tr[STATUSCD == 1 & DIA >= 1 & !is.na(TPA_UNADJ)]
L[, `:=`(dcm = DIA * 2.54, tph = TPA_UNADJ * 2.47105, hm = HT * 0.3048)]
## per tree contributions; Mg ha-1 from lb ac-1: x 0.000453592 / 0.404686; m3 ha-1 from ft3 ac-1: x 0.0283168 x 2.47105
L[, `:=`(agb = DRYBIO_AG * TPA_UNADJ * 0.000453592 / 0.404686, ba = 0.00007854 * dcm^2 * tph,
         vol = fifelse(is.na(VOLCFGRS), 0, VOLCFGRS) * TPA_UNADJ * 0.0283168 * 2.47105)]
topht <- function(d, h, w) { o <- order(-d); d <- d[o]; h <- h[o]; w <- w[o]; cw <- cumsum(w); k <- cw <= 100 | seq_along(cw) == 1
  if (sum(w[k]) < 100 && any(!k)) { j <- which(!k)[1]; w2 <- c(w[k], 100 - sum(w[k])); h2 <- c(h[k], h[j]) } else { w2 <- w[k]; h2 <- h[k] }
  ok <- !is.na(h2); if (!any(ok)) return(NA_real_); sum(h2[ok] * w2[ok]) / sum(w2[ok]) }
S <- L[, .(AGB = sum(agb, na.rm = TRUE), BAPH = sum(ba), TPH = sum(tph), VOL = sum(vol),
           HT = topht(dcm, hm, tph), n_tree = .N,
           TPH13 = sum(tph[dcm >= 12.7]), BA13 = sum(ba[dcm >= 12.7]), agb_na = sum(is.na(agb))), by = PLT_CN]
S[, `:=`(QMD = sqrt(BAPH / TPH / 0.00007854), QMD13 = sqrt(BA13 / pmax(TPH13, 1e-9) / 0.00007854))]
S[, SDI := TPH * (QMD / 25.4)^1.605]
V <- merge(P[, .(CN, pid, INVYR, myr, section, eco, REMPER, PREV_PLT_CN)], C1[, .(PLT_CN, FORTYPCD, ftg, STDAGE, SITECLCD, site,
  STDORGCD, trt, dist_any, dist_human)], by.x = "CN", by.y = "PLT_CN")
V <- merge(V, S, by.x = "CN", by.y = "PLT_CN", all.x = TRUE)
V <- merge(V, cut_by, by.x = "CN", by.y = "PLT_CN", all.x = TRUE)
for (v in c("AGB","BAPH","TPH","VOL","TPH13","BA13","n_tree")) set(V, which(is.na(V[[v]])), v, 0)
V[is.na(n_cut), n_cut := 0L]
V[, `:=`(harvested = trt | n_cut > 0)]
V[STDAGE >= 999, STDAGE := NA]
setorder(V, pid, myr)
V[, visit := seq_len(.N), by = pid]
## anonymised plot index: plot identifiers never leave the server
key <- data.table(pid = unique(V$pid)); key[, plot_ix := sample(.N)]
V <- merge(V, key, by = "pid")
say("retained visits: %d on %d plots; AGB median %.1f, 99th %.1f, 99.9th %.1f Mg/ha", nrow(V), uniqueN(V$pid),
    median(V$AGB), quantile(V$AGB, .99), quantile(V$AGB, .999))

## every valid interval: all pairs i<j inside a run of consecutive visits (linked by PREV_PLT_CN) that are all single
## condition forested. Screen hf: no harvest at any visit after i up to j. Screen ud: also no disturbance of any kind.
V[, linked := shift(CN) == PREV_PLT_CN, by = pid]; V[is.na(linked), linked := FALSE]
V[, run := cumsum(!linked), by = pid]
pairs_one <- function(g) {
  n <- nrow(g); if (n < 2) return(NULL)
  ij <- CJ(i = seq_len(n), j = seq_len(n))[i < j]
  ij[, hf := mapply(function(i, j) !any(g$harvested[(i + 1):j]) && !any(g$dist_human[(i + 1):j]), i, j)]
  ij[, ud := mapply(function(i, j) !any(g$harvested[(i + 1):j]) && !any(g$dist_any[(i + 1):j]), i, j)]
  ij[, consecutive := j == i + 1]
  ij[, `:=`(CN1 = g$CN[i], CN2 = g$CN[j])]
  ij
}
I <- V[, pairs_one(.SD), by = .(pid, run)]
v1 <- V[, .(CN, plot_ix, section, ftg, site, SITECLCD, FORTYPCD, myr, STDAGE, AGB, BAPH, TPH, QMD, HT, VOL, SDI, TPH13, BA13, harvested1 = harvested)]
I <- merge(I, v1, by.x = "CN1", by.y = "CN")
vn <- c("myr","STDAGE","AGB","BAPH","TPH","QMD","HT","VOL","SDI","TPH13","BA13"); setnames(I, vn, paste0(vn, "1"))
v2 <- V[, .(CN, ftg2 = ftg, myr2 = myr, STDAGE2 = STDAGE, AGB2 = AGB, BAPH2 = BAPH, TPH2 = TPH, QMD2 = QMD, HT2 = HT, VOL2 = VOL, SDI2 = SDI, TPH132 = TPH13, BA132 = BA13)]
I <- merge(I, v2, by.x = "CN2", by.y = "CN")
I[, `:=`(dt = myr2 - myr1)][, age2_proj := STDAGE1 + dt]
say("all pairs: %d; harvest free: %d; strict undisturbed: %d; consecutive hf: %d", nrow(I), I[hf == TRUE, .N], I[ud == TRUE, .N],
    I[hf == TRUE & consecutive == TRUE, .N])
say("dt (hf): %s", paste(round(quantile(I[hf == TRUE]$dt, c(0, .25, .5, .75, .9, 1)), 1), collapse = " "))
say("recorded STDAGE2 minus projected age (hf, consecutive): median %.1f, IQR %.1f to %.1f, share |diff|>10: %.2f",
    I[hf & consecutive, median(STDAGE2 - age2_proj, na.rm = TRUE)], I[hf & consecutive, quantile(STDAGE2 - age2_proj, .25, na.rm = TRUE)],
    I[hf & consecutive, quantile(STDAGE2 - age2_proj, .75, na.rm = TRUE)], I[hf & consecutive, mean(abs(STDAGE2 - age2_proj) > 10, na.rm = TRUE)])

## server only key; outputs keep plot_ix only
fwrite(key, file.path(OD, "SERVERONLY_plot_key.csv"))
Vout <- V[, .(plot_ix, visit, run, INVYR, myr, section, ftg, FORTYPCD, site, SITECLCD, STDORGCD, STDAGE, AGB, BAPH, TPH, QMD, HT, VOL, SDI,
              TPH13, BA13, QMD13, n_tree, harvested, dist_any, dist_human)]
Iout <- I[, .(plot_ix, run, i, j, consecutive, hf, ud, harvested1, section, ftg, ftg2, site, SITECLCD, FORTYPCD, myr1, myr2, dt, STDAGE1, STDAGE2,
              age2_proj, AGB1, AGB2, BAPH1, BAPH2, TPH1, TPH2, QMD1, QMD2, HT1, HT2, VOL1, VOL2, SDI1, SDI2, TPH131, TPH132, BA131, BA132)]
fwrite(Vout, file.path(OD, "me_plot_visits_DATA.csv")); fwrite(Iout, file.path(OD, "me_intervals_DATA.csv"))
say("wrote %d visits, %d intervals", nrow(Vout), nrow(Iout))
print(Vout[, .N, by = .(ftg)][order(-N)]); print(Vout[, .N, by = .(section)][order(-N)]); print(Vout[, .N, by = .(site)])
print(Iout[hf == TRUE, .(n = .N, plots = uniqueN(plot_ix), dt_med = median(dt)), by = .(dtc = cut(dt, c(0, 6, 11, 16, 30)))][order(dtc)])
