// CRSF technical memo: Maine FIA long-term behavior and generic AGB yield curves. Reads aggregate CSVs; writes docx.
// Usage: node build_memo.js <results_dir> <out.docx>
const fs = require("fs"), path = require("path");
const d = require("docx");
const { Document, Packer, Paragraph, TextRun, ImageRun, Table, TableRow, TableCell, WidthType, AlignmentType, BorderStyle,
        HeadingLevel, VerticalAlign, Footer, PageNumber, ShadingType } = d;
const R = process.argv[2], OUT = process.argv[3];
const splitCsv = l => { const o = []; let cur = "", q = false; for (const ch of l) { if (ch === '"') q = !q; else if (ch === "," && !q) { o.push(cur); cur = ""; } else cur += ch; } o.push(cur); return o; };
const csv = f => { const L = fs.readFileSync(path.join(R, f), "utf8").trim().split(/\r?\n/); const h = splitCsv(L[0]);
  return L.slice(1).map(l => { const v = splitCsv(l); const o = {}; h.forEach((k, i) => o[k] = v[i]); return o; }); };
const FT = { SF: "Spruce-fir", MBB: "Maple-beech-birch", AB: "Aspen-birch", WRPH: "White-red pine-hemlock", OAK: "Oak and oak-pine", EAC: "Elm-ash-cottonwood" };
const ORDER = ["SF", "MBB", "AB", "WRPH", "OAK", "EAC"];
const FONT = "Calibri", SZ = 22, TS = " ", MINUS = "−";

// rich text: ^{..} superscript, _{..} subscript, *..* italic, **..** bold
function runs(s, base = {}) {
  const out = []; const re = /(\*\*[^*]+\*\*|\*[^*]+\*|\^\{[^}]+\}|_\{[^}]+\})/g; let last = 0, m;
  while ((m = re.exec(s))) { if (m.index > last) out.push(new TextRun({ text: s.slice(last, m.index), font: FONT, size: SZ, ...base }));
    const t = m[0];
    if (t.startsWith("**")) out.push(new TextRun({ text: t.slice(2, -2), bold: true, font: FONT, size: SZ, ...base }));
    else if (t.startsWith("*")) out.push(new TextRun({ text: t.slice(1, -1), italics: true, font: FONT, size: SZ, ...base }));
    else if (t.startsWith("^")) out.push(new TextRun({ text: t.slice(2, -1).replace(/-/g, MINUS), superScript: true, font: FONT, size: SZ, ...base }));
    else out.push(new TextRun({ text: t.slice(2, -1), subScript: true, font: FONT, size: SZ, ...base }));
    last = m.index + t.length; }
  if (last < s.length) out.push(new TextRun({ text: s.slice(last), font: FONT, size: SZ, ...base }));
  return out; }
const P = (s, o = {}) => new Paragraph({ children: runs(s), spacing: { after: 140, line: 276 }, alignment: AlignmentType.JUSTIFIED, ...o });
const H = s => new Paragraph({ children: [new TextRun({ text: s, bold: true, font: FONT, size: 26, color: "1F3864" })], spacing: { before: 220, after: 100 }, heading: HeadingLevel.HEADING_1 });
const U = "Mg" + TS + "ha^{-1}", UY = "Mg" + TS + "ha^{-1}" + TS + "yr^{-1}";
const f0 = x => Math.round(+x).toString(), f1 = x => (+x).toFixed(1), f2 = x => (+x).toFixed(2);

// figures: caption below; tables: caption above
function fig(file, wIn, cap) { const buf = fs.readFileSync(path.join(R, file)); const dims = pngDims(buf); const w = wIn * 96, h = w * dims.h / dims.w;
  return [new Paragraph({ children: [new ImageRun({ data: buf, type: "png", transformation: { width: w, height: h } })], alignment: AlignmentType.CENTER, spacing: { before: 120, after: 60 } }),
          new Paragraph({ children: runs(cap, { size: 20 }), spacing: { after: 200 }, alignment: AlignmentType.JUSTIFIED })]; }
function pngDims(b) { return { w: b.readUInt32BE(16), h: b.readUInt32BE(20) }; }
const NONE = { style: BorderStyle.NONE, size: 0, color: "FFFFFF" }, RULE = { style: BorderStyle.SINGLE, size: 8, color: "000000" };
function table(cap, header, rows, widths) {
  const W = widths.reduce((a, b) => a + b, 0);
  const cell = (s, i, opts) => new TableCell({ width: { size: widths[i], type: WidthType.DXA }, verticalAlign: VerticalAlign.CENTER,
    margins: { top: 40, bottom: 40, left: 60, right: 60 }, borders: { left: NONE, right: NONE, top: opts.top ? RULE : NONE, bottom: opts.bottom ? RULE : NONE },
    children: [new Paragraph({ children: runs(String(s), { size: 18, bold: !!opts.bold }), alignment: AlignmentType.CENTER })] });
  const tr = [new TableRow({ tableHeader: true, cantSplit: true, children: header.map((h, i) => cell(h, i, { top: true, bottom: true, bold: true })) })];
  rows.forEach((r, k) => tr.push(new TableRow({ cantSplit: true, children: r.map((v, i) => cell(v, i, { bottom: k === rows.length - 1 })) })));
  return [new Paragraph({ children: runs(cap, { size: 20 }), spacing: { before: 160, after: 80 }, keepNext: true }),
          new Table({ width: { size: W, type: WidthType.DXA }, columnWidths: widths, borders: { top: NONE, bottom: NONE, left: NONE, right: NONE, insideHorizontal: NONE, insideVertical: NONE }, rows: tr }),
          new Paragraph({ children: [], spacing: { after: 160 } })]; }

// ---- data
const pick = csv("stress/pick_DATA.csv"), bs = csv("stress/bootstrap_decline_DATA.csv"), fc = csv("stress/final_curves_capped_DATA.csv");
const sc = csv("stress/scoring_DATA.csv"), mag = csv("bakuzis/magnitude_reference_DATA.csv"), dec = csv("bakuzis/net_change_by_age_DATA.csv");
const rein = csv("bakuzis/reineke_boundary_DATA.csv"), eich = csv("bakuzis/eichhorn_DATA.csv"), suk = csv("bakuzis/sukachev_DATA.csv");
const so = csv("bakuzis/site_ordering_DATA.csv"), yd = csv("bakuzis/yield_density_DATA.csv");
const by = (a, k) => Object.fromEntries(a.map(x => [x[k], x]));
const PK = by(pick, "ftg"), BS = by(bs, "ftg"), MAG = by(mag, "ftg");
const curve = (f, age) => fc.find(x => x.ftg === f && x.section === "ALL" && +x.age === age);
const SCs = f => sc.find(x => x.ftg === f && x.strat === PK[f].strat && x.form === PK[f].form);
const old = f => dec.find(x => x.screen === "hf" && x.ftg === f && x.ageclass.startsWith("[120"));
const old100 = f => dec.find(x => x.screen === "hf" && x.ftg === f && x.ageclass.startsWith("[100"));

// ---- Table 1: Bakuzis relations
const rr = f => rein.find(x => x.ftg === f && x.basis === "ge2.54cm" && x.tau === "0.99");
const r13 = f => rein.find(x => x.ftg === f && x.basis === "ge12.7cm" && x.tau === "0.99");
const soCount = (f, v) => { const a = so.filter(x => x.ftg === f && x.var === v); return a.length ? `${a.filter(x => x.inverted_beyond_90ci === "TRUE").length} of ${a.length}` : "n/a"; };
const t1rows = ORDER.map(f => { const e = eich.find(x => x.ftg === f), s = suk.find(x => x.ftg === f), r = rr(f), q = r13(f);
  return [FT[f], `${f2(r.slope).replace("-", MINUS)} (${f2(r.lo).replace("-", MINUS)} to ${f2(r.hi).replace("-", MINUS)})`,
          `${f2(q.slope).replace("-", MINUS)}`, e ? `${f0(100 * e.cv_site)}% ${e.verdict}` : "n/a",
          `${(+s.b_RD).toFixed(3).replace("-", MINUS)} PASS`, soCount(f, "AGB"), f0(MAG[f].AGB_q99)]; });

// ---- Table 2: final curves
const ci = (f, a) => { const c = curve(f, a); return `${f0(c.AGB_med)} (${f0(c.lo95)} to ${f0(c.hi95)})${c.extrapolated === "TRUE" ? "†" : ""}`; };
const t2rows = ORDER.map(f => { const b = BS[f], s = SCs(f);
  return [FT[f], `${PK[f].form} (${PK[f].strat})`, ci(f, 50), ci(f, 100), ci(f, 150), ci(f, 200), `${f2(b.p_decline_250)}`, f1(s.cv_rmse), `${f1(s.loso_rmse)} (${f2(s.loso_cover90)})`]; });

// ---- text
const T = {};
T.oldSF = old("SF"), T.oldMBB = old("MBB"), T.oldWRPH = old("WRPH");
const sfPk = +curve("SF", 100).AGB_med, sf250 = +curve("SF", 250).AGB_med;
const body = [];
body.push(new Paragraph({ children: [new TextRun({ text: "Long-term behavior and generic aboveground biomass yield curves for Maine forest types from FIA remeasurements", bold: true, font: FONT, size: 32, color: "1F3864" })], spacing: { after: 120 } }));
body.push(new Paragraph({ children: runs("CRSF technical memo  |  Aaron Weiskittel, Center for Research on Sustainable Forests, University of Maine  |  October 15, 2026  |  DRAFT", { size: 20, color: "595959" }), spacing: { after: 240 } }));

body.push(H("Bottom line"));
body.push(P(`Harvest free Maine FIA plots keep accumulating aboveground live tree biomass (AGB) well past 100 years. Net change slows with age but its mean stays positive through the oldest well sampled age class in every major forest type, at ${f1(T.oldSF.dB_mean)} ${UY} in spruce-fir (95% CI ${f1(T.oldSF.dB_lo)} to ${f1(T.oldSF.dB_hi)}) and ${f1(T.oldMBB.dB_mean)} in maple-beech-birch (${f1(T.oldMBB.dB_lo)} to ${f1(T.oldMBB.dB_hi)}) for stands older than 120 years, and ${f1(old100('WRPH').dB_mean)} in white-red pine-hemlock (${f1(old100('WRPH').dB_lo)} to ${f1(old100('WRPH').dB_hi)}) at 100 to 120 years, where its oldest class with enough plots for an interval ends. A model form that allows a decline was therefore fitted for every type and left to the data to decide whether a decline occurs. Spruce-fir is the only major type where it does, peaking near 110 years and falling to ${f0(100 * sf250 / sfPk)}% of its age 100 value by 250 years, with a bootstrap probability of a peak before 250 years of ${f2(BS.SF.p_decline_250)}. The other types level off or keep rising slowly within 250 years, and all curves are capped at the observed 99th percentile AGB of their type.`));
body.push(P(`The curves are generic, as they predict the typical AGB of a stand of a given type, ecoregion section and age, and they do not reproduce the trajectory of an individual plot. Held out plots are predicted from age with mean biases of ${MINUS}4 to +23 ${U} by type and 90% prediction intervals that cover 86 to 98% of observations, yet the same curves explain little of the plot to plot variation in net change, which is dominated by mortality events the curves cannot see. Everything beyond each cell's oldest observed stands, roughly 100 to 170 years, is an extrapolation and is marked as such in the tables.`));

body.push(H("Data and approach"));
body.push(P(`The analysis used the Maine annual inventory from 1999 to 2024 (FIADB, April 2026 extraction) on OSC Cardinal, restricted to single condition forested plots in the annual design. Plot AGB is the sum of tree DRYBIO_AG for live trees at least 2.54 cm dbh. Every pair of visits within a linked run of remeasurements was retained, not only consecutive pairs, giving 15,500 harvest free intervals on 2,383 of the 2,554 retained plots, with lengths up to 25.7 years (median 9.7). Curve shapes used the 12,218 intervals on 2,194 plots that started at 10 ${U} or more and did not start right after a harvest. An interval was dropped if any visit after its start recorded a treatment, a removed tree, or a tree killed by silviculture, and a second, strict screen also dropped any recorded natural disturbance. Forest type groups follow FORTYPCD and ecoregions follow the section of ECOSUBCD (211A, 211B, 211C, 211D, M211A and 221A). Plot identifiers stayed on the server and coordinates were never read.`));
body.push(P(`Four curve forms were compared, namely Chapman-Richards (CR), CR with exponential senescence (CRS), CR with a bounded decline to a floor (CRD), and the Hoerl gamma (HG). Each was fitted two ways. The remeasurement strategy (S1) estimates shape only from how plots of known age change over each interval, through an algebraic difference form in which the level cancels. The joint strategy (S2) estimates shape from the chronosequence and the repeated measures together in a mixed model with section and plot random effects. The level of each curve comes from a mixed model with a type fixed effect and type by section and plot random effects. Selection screens, the scoring rule and the stress tests were written into the scripts before their output was read.`));

body.push(H("Observed long-term behavior"));
body.push(P(`Table 1 summarizes the empirical Bakuzis assessment. The self-thinning boundary on all stems (99th percentile quantile regression) has slopes of ${MINUS}1.64 to ${MINUS}1.93 in the five larger types, consistent with Reineke's ${MINUS}1.605 within the stated tolerance of 0.30, while the elm-ash-cottonwood estimate is too imprecise to test. The boundary for stems of at least 12.7 cm is much shallower because ingrowth feeds that class, so it is not a closed self-thinning population. Eichhorn's rule holds in five of six types, the Sukachev effect holds in all six (denser stands lose stems faster), and yield-density behaves as expected, with denser stands carrying more AGB when young and no relation by 75 to 85 years. Site ordering of AGB inverts in maple-beech-birch in 3 of 12 age bins beyond the 90% interval, which points to SITECLCD being a weak productivity tier for northern hardwoods rather than to a dynamics failure.`));
body.push(...table("**Table 1.** Empirical Bakuzis relations from harvest free Maine FIA visits and intervals. Intervals are 95% plot cluster bootstrap intervals.",
  ["Forest type", "Reineke slope, all stems (95% CI)", "Slope, stems ≥12.7 cm", "Eichhorn site CV", "Sukachev slope on RD", "Site order inversions, AGB", "Observed 99th pct AGB (" + U + ")"],
  t1rows, [1500, 1500, 1050, 1150, 1150, 1100, 1110]));
body.push(P(`Net change by stand age is the direct evidence on decline (Figure 1). Across all harvest free intervals, mean net change falls from about 4 ${UY} at 20 to 40 years to 1.4 to 2.4 ${UY} beyond 100 years. The share of intervals with a net loss rises from under 5% in young stands to 20 to 25% beyond 100 years, so individual old stands do lose biomass, yet gains elsewhere outweigh those losses on average. The result holds on the strict undisturbed screen and across 5, 10 and 15 to 25 year horizons. Only elm-ash-cottonwood crosses zero, near 120 years, and that rests on fewer than 15 plots.`));
body.push(...fig("bakuzis/me-fia_net-agb-change-by-age.png", 6.3, `**Figure 1.** Annual net AGB change on harvest free remeasurement intervals by stand age at the interval midpoint. Lines and bands are GAMs with a plot random effect and 95% confidence bands; points are age class means with 95% plot cluster bootstrap intervals.`));

body.push(H("Curve forms, strategies and stress test"));
body.push(P(`The two strategies disagree, and the disagreement is the central technical finding (Figure 2). Remeasurement only shapes never level off, and in four of six types they exceed the observed 99.9th percentile AGB before 250 years. Joint shapes level off and predict held out plots from age with small bias, yet in the four larger types they underpredict observed net change by about 1.1 to 1.4 ${UY}. The chronosequence is flatter than present day growth implies. Partial harvest before 1999 in stands that now look old, the 1970s and 1980s spruce budworm outbreak, and an ambiguous stand age in multi-cohort stands can each produce that pattern, and this analysis cannot separate them.`));
body.push(...fig("curves/me-fia_agb-curve-forms-strategies_memo.png", 6.3, `**Figure 2.** Statewide median AGB curves for every form under the remeasurement only (S1, dashed) and joint (S2, solid) strategies. The grey band is a chronosequence GAM with 95% confidence band; dotted lines mark the observed 99th percentile stand age and the 99.9th percentile AGB. The aspen-birch CRD (S2) curve spikes near age 1 and was rejected by the early rise screen.`));
body.push(P(`Restricting the candidates to the three decline-capable forms, each was refitted under five data variants, namely the harvest free base, strict undisturbed, a plot consistent stand age, stands aged 15 years and older, and a salvage bound that keeps intervals where harvest followed natural disturbance. Candidates had to rise monotonically from age 1, peak no earlier than 40 years, keep at least half their maximum at 250 years and stay below the observed 99.9th percentile. Survivors were scored on generic curve cross validation RMSE, leave one section out RMSE and the spread of B(200)/B(100) across the five variants (Figure 4). In five of six types the chosen decline-capable form drove its decline parameter to near zero, as the data asked for no decline within 250 years. Spruce-fir kept a decline in four of five variants, with AGB at 250 years ranging from 0.43 to 0.91 of the peak, and the decline disappeared when stands younger than 15 years were removed. Bootstrap draws whose refit shape spiked near age 1 were dropped before the intervals were computed (20 of 200 for spruce-fir, 4 for aspen-birch).`));
body.push(...table("**Table 2.** Final generic AGB curves, statewide median (" + U + ") with 95% plot cluster bootstrap interval (200 resamples). † beyond the type's 99th percentile observed stand age. P(decline) is the bootstrap probability that the curve peaks before 250 years. CV and LOSO RMSE are in " + U + "; LOSO coverage of the 90% prediction interval is in parentheses.",
  ["Forest type", "Form (strategy)", "Age 50", "Age 100", "Age 150", "Age 200", "P(decline)", "CV RMSE", "LOSO RMSE (cover)"],
  t2rows, [1300, 950, 1000, 1000, 1000, 1000, 800, 800, 1210]));

body.push(H("Generic curves"));
body.push(P(`Figure 3 shows the final curves by section. Section effects are modest relative to plot to plot variation (log scale standard deviations of about 0.2 for type by section and 0.5 for plot), and sections with few plots are shrunk toward the type mean. The workbook gives annual values from 0 to 250 years for every type and section cell with median and mean curves, the bootstrap confidence interval, a 90% prediction interval for a single plot, and flags for extrapolation and capping. The median is the curve for a typical stand, while the mean is the right input when the curves are multiplied by area for a population total.`));
body.push(...fig("stress/me-fia_agb-final-curves-by-section.png", 6.3, `**Figure 3.** Final generic AGB curves by forest type and ecoregion section, with the selected form and strategy in each panel title. Bands are 95% plot cluster bootstrap intervals on the median curve; dashed segments lie beyond the cell's 99th percentile observed stand age; dotted lines are the caps at the type's observed 99th percentile AGB.`));
body.push(...fig("stress/me-fia_stress-decline-r250.png", 6.3, `**Figure 4.** Stress test of the decline. Each point is the ratio of AGB at 250 years to the curve maximum for one decline-capable form, strategy and data variant; circles mark the form chosen for the type.`));

body.push(H("Comparison with the earlier Maine curves"));
const BM = by(csv("benchmark/benchmark_by_type_DATA.csv"), "ftg");
body.push(P(`The earlier Maine age based curves (ME.AGB.R) use the same Hoerl gamma form, AGB = b_{1} Age^{b2} b_{3}^{Age}, fitted separately for 109 combinations of EPA Level III ecoregion and FORTYPCD. Scored on the same harvest free visits, with ECOSUBCD sections crosswalked approximately to the Level III ecoregions, they underpredict current FIADB AGB by ${f0(-BM.ALL.prior_bias)} ${U} on average (RMSE ${f0(BM.ALL.prior_rmse)}), and the shortfall grows with age to 49 to 75 ${U} beyond 120 years (Figure 5). The new mean curves have a bias of +${f0(BM.ALL.new_bias_mean)} and an RMSE of ${f0(BM.ALL.new_rmse_mean)} ${U} on these in sample visits and cross validated RMSE of 43 to 61. Of the 109 earlier curves, 30 do not peak within 300 years and 10 fall below a quarter of their maximum by 300 years, the type of cell to cell instability that pooling to forest type groups with shrinkage removes. The level offset may reflect a change in the biomass basis or in plot screening between the two efforts, which this comparison cannot resolve.`));
body.push(...fig("benchmark/me-fia_prior-vs-new-curves.png", 6.3, `**Figure 5.** Earlier ME.AGB curves (colored, one per FORTYPCD and EPA Level III ecoregion) against the new statewide curves (black, solid median and dashed mean) and harvest free plot visits (grey).`));
body.push(H("Limits and next steps"));
body.push(P(`Three limits matter most for 100 to 200 year projections. First, the oldest well sampled stands are 110 to 170 years old, so every decline or plateau beyond that is a model choice informed by shape and is not observed. Second, the aspen-birch curve came from the remeasurement only strategy because every joint candidate exceeded its magnitude screen, and it is effectively a power function of age (Hoerl gamma with its decline term at zero) that rises to the cap. A stand that stays aspen-birch for 200 years is not biologically realistic in Maine, so long projections of that type should hand off to a successor type once the stand ages past about 110 years. Third, the curves carry no climate, disturbance regime or management signal, so they describe the recent harvest free past and are not a forecast.`));
body.push(P(`The most useful next steps are an attribution of the chronosequence versus remeasurement gap, using stand origin and pre 1999 harvest history from the periodic inventories, a type transition matrix for aspen-birch and other short lived types, and a check of these curves against FVS-ACD projections of the same plots, which would test them against an independent growth model rather than against the data that built them.`));
body.push(P(`**Data and code.** Curve tables, diagnostics, figures and R scripts are archived on Zenodo (version DOI 10.5281/zenodo.23406410, reserved October 15, 2026, and in draft until published). The scripts will also be in the holoros/maine-agb-yield-curves repository.`));

const doc = new Document({ styles: { default: { document: { run: { font: FONT, size: SZ } } } },
  sections: [{ properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1300, bottom: 1300, left: 1300, right: 1300 } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [new TextRun({ children: [PageNumber.CURRENT], font: FONT, size: 18 })] })] }) },
    children: body }] });
Packer.toBuffer(doc).then(b => { fs.writeFileSync(OUT, b); console.log("wrote", OUT); });
