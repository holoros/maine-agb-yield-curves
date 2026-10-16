# Maine FIA long-term AGB behavior and generic yield curves (October 15, 2026)

Compute: OSC Cardinal (Aaron named Cardinal because firebreather was down), `~/maine_agb_yc_20261015/`, R 4.5.2 (`module purge; module load gcc/12.3.0 R/4.5.2`).
Input: `~/fia_data/ME_{PLOT,PLOTGEOM,COND,TREE}.csv` (FIADB Maine, April 2026 pull). LAT and LON are never read.

| Step | Script | md5 | Job | Output |
|:--|:--|:--|:--|:--|
| 1 | 01_build_me_intervals.R | 17e462f5 | 15542964 | out_b/me_plot_visits_DATA.csv, me_intervals_DATA.csv (plot index anonymised; key SERVERONLY) |
| 2 | 02_bakuzis_empirical.R (+ bakuzis_check_step-aware.R, 297e2ef3) | 112acb17 | 15543196 | out_b/bakuzis/ |
| 3 | 03_fit_compare_forms.R | 6192e841 | 15543389 | out_b/curves/ |
| 4 | 04_stress_decline.R | 3a148208 | 15544016 | out_b/stress/ |
| 5 | 05_stress_figure.R | 42f98604 | 15544016 | out_b/stress/me-fia_stress-decline-r250.png |
| 6 | 06_benchmark_prior.R (+ ME.AGB.R, 802b0427) | dd582243 | 15544555 | out_b/benchmark/ |
| 7 | 07_graphical_abstract.R | 1e3928e9 | 15544848 | out_b/00_maine-agb-yield-curves_graphical_abstract.png |

Final aggregate bundle: `meyc_results_h.tgz` (md5 28d32a81), no plot identifiers (grep plot_ix = 0).

```bash
O=$PWD/out_b
Rscript --vanilla scripts/01_build_me_intervals.R ~/fia_data $O
Rscript --vanilla scripts/02_bakuzis_empirical.R $O scripts/bakuzis_check_step-aware.R
Rscript --vanilla scripts/03_fit_compare_forms.R $O 100 2000 28
Rscript --vanilla scripts/04_stress_decline.R $O 200 28
Rscript --vanilla scripts/05_stress_figure.R $O
Rscript --vanilla scripts/06_benchmark_prior.R $O scripts/ME.AGB.R
Rscript --vanilla scripts/07_graphical_abstract.R $O $O/00_maine-agb-yield-curves_graphical_abstract.png [doi]
```
Seeds: 20261015 (steps 1, 3), 4747 (step 2), 4040 (step 4); L'Ecuyer-CMRG streams under mclapply.

Sandbox builds (Python and Node, from the aggregate bundle): docs workbook by code/build_xlsx.py; memo by code/build_memo.js (rendered with LibreOffice, pages inspected).
ME.AGB.R (the earlier curves, supplied by Aaron) is not redistributed here; step 6 sources it on Cardinal.
