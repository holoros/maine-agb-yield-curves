# Generic aboveground biomass yield curves for Maine forest types and ecoregions from FIA remeasurements

Version 1.0.0 (October 15, 2026). Aaron R. Weiskittel, Center for Research on Sustainable Forests, University of Maine.

## What this is

Generic aboveground live tree biomass (AGB, Mg ha-1, FIADB DRYBIO_AG) yield curves for six Maine forest type groups and six USDA
Forest Service ecological sections, from 0 to 250 years. The record also holds an empirical Bakuzis assessment of observed long-term
stand behavior, a stress test of curve forms that allow a decline, and a benchmark against earlier Maine age based AGB curves.
Only aggregate tables, figures and code are included. No plot identifiers or coordinates appear in any file.

## Data

FIADB Maine annual inventory 1999 to 2024 (April 2026 extraction), single condition forested plots in the annual design. Every pair of
visits within a linked run of remeasurements was used, including first to last, giving 15,500 harvest free intervals on 2,383 of the 2,554 retained plots
(lengths up to 25.7 yr). Curve shapes used the 12,218 intervals on 2,194 plots that started at 10 Mg ha-1 or more and did not start
right after a harvest.

## Main results

- Observed net AGB change stays positive in the oldest well sampled age class of every major type: 1.41 Mg ha-1 yr-1 in spruce-fir
  (95% plot cluster bootstrap interval 0.98 to 1.92) and 1.48 in maple-beech-birch (0.87 to 2.05) beyond 120 years
  (data/net_change_by_age_DATA.csv).
- Decline-capable forms were fitted for every type. Spruce-fir is the only type that declines, peaking near 110 yr and falling to 82% of
  its age 100 value by 250 yr, with a bootstrap probability of a peak before 250 yr of 0.91 (data/bootstrap_decline_DATA.csv). The
  probability is 0.30 for elm-ash-cottonwood, 0.17 for maple-beech-birch and 0.04 or less for the others.
- All curves are capped at the type's observed 99th percentile AGB. Values beyond each cell's 99th percentile observed age are flagged
  as extrapolations.
- The earlier ME.AGB curves underpredict current FIADB AGB on the same harvest free visits by 29 Mg ha-1 on average (RMSE 52), against
  a mean bias of +8 and an RMSE of 46 Mg ha-1 for the new mean curves in sample (data/benchmark_by_type_DATA.csv).

## Files

Zenodo stores files without folders; the folder names below are the layout of the companion GitHub repository and the staging package. Every file name is unique.

| Path | Content |
|:--|:--|
| 00_maine-agb-yield-curves_graphical_abstract.png | Lead figure |
| data/final_curves_capped_DATA.csv | Final curves, annual 0 to 250 yr, by type and section and statewide (median, mean, 95% bootstrap interval, 90% prediction interval, cap and extrapolation flags) |
| data/final_params_DATA.csv, final_level_DATA.csv, pick_DATA.csv | Shape parameters, level effects, chosen form and strategy per type |
| data/scoring_DATA.csv, stress_shape_metrics_DATA.csv, bootstrap_decline_DATA.csv, loso_DATA.csv | Stress test, scoring, bootstrap decline, leave one section out transfer |
| data/*_DATA.csv (Bakuzis) | Magnitude reference, Reineke boundary, Eichhorn, Sukachev, site ordering, yield density, binned trajectories, net change by age and horizon |
| data/cv_*, equivalence_dB, form_screen_selection | Step 3 comparison of forms and strategies |
| data/benchmark_*, prior_curve_shapes | Benchmark against the earlier ME.AGB curves |
| data_dictionary.csv | Every column of every data file, with units |
| figures/ | Bakuzis matrices, net change, Reineke boundary, curve forms, final curves, stress test, benchmark |
| code/ | R scripts 01 to 07 (Cardinal), bakuzis_check_step-aware.R, build_xlsx.py, build_memo.js |
| docs/ | Workbook, technical memo (PDF), RUN.md (job ids and md5s), implementation notes |

## Curve equation

AGB(t) = min(cap, exp(mu_type + u_cell) x f(t)), where f is the chosen form (CR, CRS, CRD or HG; see data_dictionary.csv) with the
parameters in final_params_DATA.csv. The median curve is the typical plot; the mean curve multiplies by exp((s_plot^2 + s_res^2)/2)
and is the right input for population totals.

## Uncertainty

Curve bands are 95% plot cluster bootstrap confidence intervals on the median curve (200 resamples; shape and level refit on one plot
draw; draws failing the early rise screen dropped, 20 for spruce-fir and 4 for aspen-birch). PI90 is a 90% prediction interval for a
single plot. Generic curve cross validation coverage of PI90 is 0.86 to 0.98 by type. Equivalence tests on annual net change fail the
slope test for every candidate, so the curves are population curves and do not predict plot level change.

## Known limitations

- The oldest well sampled stands are 110 to 170 years old, so every plateau or decline beyond that is a model choice informed by shape.
- Remeasurement only shapes do not level off while chronosequence shapes do; the gap (about 1.1 to 1.4 Mg ha-1 yr-1 in net change for
  the four larger types) is not attributed.
- Aspen-birch uses the remeasurement only fit, a power function of age that rises to the cap; long projections should hand off to a
  successor type.
- The benchmark crosswalks ECOSUBCD sections to EPA Level III ecoregions approximately.
- Site ordering inverts for maple-beech-birch, which points to SITECLCD being a weak productivity tier for northern hardwoods.

## Restricted data audit

See RUNBOOK.md. LAT and LON were never read; plot identifiers stayed on the analysis server.

## Citation

Cite the concept DOI shown on the Zenodo record page. This version: https://doi.org/10.5281/zenodo.23406410 (see CITATION.cff).
