"""Rebuild main/supplementary tables and figures from the REPAIR outputs (FoCus Paper 1).
Inputs: staged copies of Repair/OutPW, Repair/OutExt, Repair/Outputs_FunctionalForm.
"""
import csv, math, re, os, statistics as st, json
from pathlib import Path

R = Path(os.environ.get("FOCUS_BASE_DIR", "path/to/FoCus")) / "Repair"
PW, EXT, FF = R / "OutPW", R / "OutExt", R / "Outputs_FunctionalForm"
OUT = Path(os.environ.get("TABLES_OUT_DIR", "tables_figures_output"))
OUT.mkdir(parents=True, exist_ok=True)

def rd(p, delim=";"):
    return list(csv.DictReader(open(p, encoding="utf-8-sig"), delimiter=delim))
def f(x):
    if x is None: return float("nan")
    x = str(x).strip().replace(",", ".")
    if x in ("", "NA", "NaN", "Inf", "-Inf"): return float("nan") if x in ("", "NA", "NaN") else float(x.replace("Inf", "inf"))
    return float(x)
def norm(s): return re.sub(r"[^a-z0-9]", "", s.lower())

# ---------- display names / N SNPs / DB / category from old S3 (unchanged per Abgleich) ----------
si = open(SCR / "si.md", encoding="utf-8").read()
s3 = si[si.index("## S3."):si.index("## S4.")]
old = {}
for line in s3.split("\n"):
    cells = re.split(r"\s{2,}", line.strip())
    if len(cells) >= 8 and re.match(r"\d[\d\s\xa0]* / \d+", cells[1]):
        old[norm(cells[0])] = dict(name=cells[0], db=cells[5], cat=cells[6], nsnp=int(re.sub(r"\D","",cells[7])))
NAME_FIX = {"Starch Sucrose Metabolism": "Starch/Sucrose Metabolism", "Fructose Mannose Metabolism": "Fructose/Mannose Metabolism",
            "Beta Oxidation": "Beta-Oxidation", "Glycine Serine Threonine Metabolism": "Glycine/Serine/Threonine Metabolism",
            "Valine Leucine Isoleucine Degradation": "Valine/Leucine/Isoleucine Degradation", "Cysteine Methionine Metabolism": "Cysteine/Methionine Metabolism",
            "Alanine Aspartate Glutamate Metabolism": "Alanine/Aspartate/Glutamate Metabolism", "Arginine Proline Metabolism": "Arginine/Proline Metabolism",
            "Glycolysis Gluconeogenesis": "Glycolysis/Gluconeogenesis", "Nicotinate Nicotinamide Metabolism": "Nicotinate/Nicotinamide Metabolism",
            "Mitochondrial Iron Sulfur Cluster Biogenesis": "Mitochondrial Iron-Sulfur Cluster Biogenesis"}
def meta(p):
    m = old.get(norm(p))
    if m is None: raise KeyError(p)
    m = dict(m); m["name"] = NAME_FIX.get(m["name"], m["name"])
    return m

track = rd(PW / "AllPathways_Results_Tracking_Table.csv")
T = {r["Pathway"]: r for r in track}
robust = [r["Pathway"] for r in track if r["Robust_Candidate"].upper() == "TRUE"]
dr2 = {r["Pathway"]: r for r in rd(PW / "Pathway_ALL_Multivariate_SNP_SumScore_TestSet_Regression_Results_DeltaR2_01_Main_LDpruned_All_Oriented.csv")}

def fmt(x, d=2):
    return f"{x:.{d}f}"
def ci(o, l, u, d=2): return f"{fmt(o,d)} ({fmt(l,d)}–{fmt(u,d)})"
def sci(x):
    if x == 0: return "0"
    if x < 1e-3:
        m, e = f"{x:.1e}".split("e")
        return f"{m} × 10^{int(e)}"
    return f"{x:.3g}" if x < 0.01 else f"{x:.3f}"

def write(name, header, rows):
    with open(OUT / name, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh, delimiter=";"); w.writerow(header); w.writerows(rows)
    return header, rows

tables = {}
# ---------- Table 2 ----------
rows = []
for p in sorted(robust, key=lambda p: -f(T[p]["Multivariate_OR"])):
    d = dr2[p]; m = meta(p)
    rows.append([m["name"], m["nsnp"], f'{d["N_train"]} / {d["N_test"]}',
                 fmt(f(d["OR_train"])), f'{fmt(f(d["CI_low_train"]))}–{fmt(f(d["CI_high_train"]))}',
                 fmt(f(d["OR_test"])), f'{fmt(f(d["CI_low_test"]))}–{fmt(f(d["CI_high_test"]))}',
                 sci(f(T[p]["Q_Value_BH"]))])
tables["Table 2"] = write("Table2_RobustCandidates.csv",
    ["Pathway", "N SNPs", "N (train / test)", "OR per SD (train)", "95% CI (train)", "OR per SD (test)", "95% CI (test)", "q-value (BH-FDR)"], rows)

# ---------- S2a flow ----------
flow = rd(PW / "AllPathways_PRISMA_Flow_Counts.csv")
lab = ["Total pathways analyzed", "With complete cached results (main model, sensitivity variants, test-set replication)",
       "FDR significance (BH-adjusted q < 0.05, training-set Sum-Score test)",
       "+ Sensitivity consistency (no direction reversal in LD-pruning or beta-/ridge-/elastic-net-weighted variants)",
       "+ Train/test direction consistency",
       "+ Effect-size threshold (training-set OR per SD ≥ 1.10 or ≤ 1/1.10) – final robust candidates"]
rows = [[lab[i], r["N_Remaining"], "–" if r["N_Excluded"] in ("NA", "") else r["N_Excluded"]] for i, r in enumerate(flow)]
tables["Table S2a"] = write("TableS2a_RobustCandidate_Flow.csv", ["Stage", "Pathways remaining", "Excluded at this stage"], rows)
flow_n = [int(r["N_Remaining"]) for r in flow]; flow_x = [r["N_Excluded"] for r in flow]

# ---------- S3 (all 54) ----------
def reason(r):
    if r["Robust_Candidate"].upper() == "TRUE": return "–"
    if r["Criterion_FDR"].upper() != "TRUE": return "not FDR-significant"
    if r["Criterion_Sensitivity"].upper() != "TRUE": return "direction reversal in sensitivity variant"
    if r["Criterion_TrainTest_Direction"].upper() != "TRUE": return "direction reversed in test set"
    if r["Criterion_EffectSize"].upper() != "TRUE": return "below effect-size threshold"
    return "?"
rows = []
for r in sorted(track, key=lambda r: f(r["Q_Value_BH"])):
    p = r["Pathway"]; m = meta(p); d = dr2.get(p)
    ntt = f'{d["N_train"]} / {d["N_test"]}' if d else "?"
    rows.append([m["name"], ntt, f'{f(r["Q_Value_BH"]):.3g}', "Yes" if r["Robust_Candidate"].upper() == "TRUE" else "No",
                 reason(r), m["db"], m["cat"], m["nsnp"]])
tables["Table S3"] = write("TableS3_All54_Pathways.csv",
    ["Pathway", "N (train / test)", "q-value (BH-FDR)", "Robust candidate", "Reason if not robust", "Database", "Category", "N SNPs"], rows)

# ---------- S4 functional form (single-covariate GAMs, LRT; everyday activity from all-covariate model as before) ----------
age = rd(FF / "Sensitivity_Age_Functional_Form_GAM_Summary.csv")[0]
allc = {r["Variable"]: r for r in rd(FF / "Sensitivity_AllCovariates_Functional_Form_GAM_Summary.csv")}
se = {r["Covariate"]: r for r in rd(FF / "Sensitivity_SportsEnergy_Batch_Summary.csv", ",")}
agebatch = rd(FF / "Age_OR_Before_After_Batch.csv", ",")
rows = [
 ["Age (years)", fmt(f(age["edf_smooth_age"])), sci(f(age["LRT_p_age"])),
  f'Non-linear; resolves under recruitment-arm adjustment (edf = 1.00; linear term p = {f(agebatch[1]["p"]):.3f})'],
 ["Sports activity (min/week)", fmt(f(se["sports_activity_min_week"]["edf_no_arm"])), sci(f(se["sports_activity_min_week"]["LRT_p_no_arm"])),
  f'Non-linear; persists under recruitment-arm adjustment (edf = {f(se["sports_activity_min_week"]["edf_with_arm"]):.2f}, LRT p = {f(se["sports_activity_min_week"]["LRT_p_with_arm"]):.3f})'],
 ["Energy intake (kcal/day)", fmt(f(se["energy_kcal"]["edf_no_arm"])), f'{f(se["energy_kcal"]["LRT_p_no_arm"]):.3f}',
  f'Modest non-linearity; linear under recruitment-arm adjustment (edf = {f(se["energy_kcal"]["edf_with_arm"]):.2f})'],
 ["Everyday activity (min/week)", fmt(f(allc["everyday_activity_min_week"]["edf"])), f'{f(allc["everyday_activity_min_week"]["p_value"]):.3f}',
  "No evidence of non-linearity"],
]
tables["Table S4"] = write("TableS4_FunctionalForm.csv", ["Covariate", "edf (GAM)", "LRT p-value", "Interpretation"], rows)

# ---------- S7 sensitivity variants (per-SD OR, training set) ----------
var = [("Main model", "01_Main_LDpruned_All_Oriented"), ("Unpruned", "02_Sensitivity_NoLDpruned_All_Oriented"),
       ("Beta-weighted", "05_Sensitivity_BetaWeighted_LDpruned_Beta_weighted"), ("Ridge-weighted", "05_Sensitivity_BetaWeighted_LDpruned_Ridge_weighted"),
       ("Elastic-net", "05_Sensitivity_BetaWeighted_LDpruned_Elasticnet_weighted"), ("Functional-VEP", "03_Sensitivity_Functional_VEP_LDpruned_Oriented"),
       ("LoF-Burden", "06_Sensitivity_LoF_VEP_LDpruned_Lof_burden")]
def score_row(p, suf):
    fn = PW / f"Pathway_{p}_Multivariate_SNP_SumScore_Regression_Results_{suf}.csv"
    if not fn.exists(): return None
    rr = rd(fn)
    x = [r for r in rr if r["term"] == "SNP_SumScore"]
    return x[0] if x else None
S7 = {}; rows = []
order = sorted(robust, key=lambda p: -f(T[p]["Multivariate_OR"]))
for p in order:
    row = [meta(p)["name"]]; S7[p] = {}
    for lab_, suf in var:
        x = score_row(p, suf)
        if x is None or not math.isfinite(f(x["OR_Score_Per_SD"])):
            row.append("n/a"); continue
        o, l, u, pv = f(x["OR_Score_Per_SD"]), f(x["OR_Score_Per_SD_CI_lower"]), f(x["OR_Score_Per_SD_CI_upper"]), f(x["p.value"])
        S7[p][lab_] = (o, l, u, pv)
        row.append(ci(o, l, u) + ("*" if pv < 0.05 else ""))
    rows.append(row)
tables["Table S7"] = write("TableS7_SensitivityVariants.csv", ["Pathway"] + [v[0] for v in var], rows)

# ---------- S8 covariates ----------
cov = [("SNP_SumScore", "SNP-Sum-Score (per allele)"), ("school_cert", "School-leaving certificate (low vs. reference)"),
       ("graduation", "Academic degree (none vs. reference)"), ("sample_age", "Age (per year)"),
       ("sports_activity", "Sports activity (per min/week)"), ("everyday_activity", "Everyday activity (per min/week)"),
       ("gender", "Sex (female vs. male)"), ("smoking_catformer", "Smoking: former vs. never"),
       ("smoking_catcurrent", "Smoking: current vs. never"), ("energy_kcal", "Energy intake (per kcal/day)")]
def stars(p): return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
def g3(x): return f"{x:.4f}" if abs(x - 1) < 0.01 else f"{x:.3g}"
rows = []; covstats = {k: [] for k, _ in cov}
for p in order:
    rr = rd(PW / f"Pathway_{p}_Multivariate_SNP_SumScore_Regression_Results_01_Main_LDpruned_All_Oriented.csv")
    row = [meta(p)["name"]]
    for key, _ in cov:
        x = [r for r in rr if r["term"].startswith(key)]
        if not x: row.append("–"); continue
        x = x[0]; o, l, u, pv = f(x["OR"]), f(x["OR_CI_lower"]), f(x["OR_CI_upper"]), f(x["p.value"])
        covstats[key].append((o, l, u, pv, x["term"]))
        row.append(f"{g3(o)} ({g3(l)}–{g3(u)}){stars(pv)}")
    rows.append(row)
tables["Table S8"] = write("TableS8_CovariateORs.csv", ["Pathway"] + [c[1] for c in cov], rows)

# ---------- S9 / S10 ----------
desc = rd(PW / "Q6_Descriptive_Tier_Summary_v2.csv")
tiername = {"Reference": "Reference (intronic/synonymous)", "Regulatory_UTR": "Regulatory/UTR", "Missense_Lower": "Missense, lower pathogenicity",
            "Missense_Higher": "Missense, higher pathogenicity", "LoF": "Loss-of-function²", "Missense_Pooled": "Missense, pooled³"}
rows = []
for r in desc:
    rows.append([tiername[r["Tier"]], r["N_total"], r["N_with_estimate"], f'{f(r["mean_abs_logOR"]):.3f} ({f(r["sd_abs_logOR"]):.3f})',
                 r["N_with_direction"], r["N_risk"], f'{100*f(r["risk_proportion"]):.1f}% ({100*f(r["risk_CI_lower"]):.1f}–{100*f(r["risk_CI_upper"]):.1f}%)',
                 f'{f(r["p_vs_50pct"]):.3f}'])
tables["Table S9"] = write("TableS9_SNPLevel_Descriptive.csv",
    ["Tier", "N SNPs", "N with estimate", "Mean |log OR| (SD)", "N with direction", "N risk-increasing", "Risk-increasing proportion (95% CI)", "p vs. 50%¹"], rows)

comb = rd(PW / "Q6_Results_Combined.csv"); maf = rd(PW / "Q6_MAF_Sensitivity_Results.csv")
tn = {"Missense_Higher": "Missense, higher pathogenicity", "Missense_Lower": "Missense, lower pathogenicity", "Regulatory_UTR": "Regulatory/UTR", "Missense_Pooled": "Missense, pooled", "Missense": "Missense, pooled"}
def sgn(x, d=3):
    if round(abs(x), d) == 0: return f"{0:.{d}f}"
    return ("+" if x >= 0 else "−") + f"{abs(x):.{d}f}"
def num(x, d=3):
    if round(abs(x), d) == 0: return f"{0:.{d}f}"
    return ("−" if x < 0 else "") + f"{abs(x):.{d}f}"
def pf(p): return f"{p:.4f}" if p < 0.01 else f"{p:.3f}" if p < 0.1 else f"{p:.2f}"
rows = []
for r in comb:
    if r["Hypothesis"] == "Effect-magnitude":
        rows.append(["Effect-magnitude (|log OR| scale)", tn[r["Tier"]], sgn(f(r["Estimate"])), f'[{num(f(r["CI_lower"]))}, {num(f(r["CI_upper"]))}]', pf(f(r["p_value"])), pf(f(r["q_value_BH"]))])
em = rows[:]
for r in maf:
    if r["Model"].startswith(("Split (b)", "Pooled (b)")) and not (r["Model"].startswith("Pooled") and r["Tier"] == "Regulatory_UTR"):
        rows.append(["Effect-magnitude, MAF-adjusted (|log OR| scale)", tn[r["Tier"]], sgn(f(r["Estimate"])), f'[{num(f(r["CI_lower"]))}, {num(f(r["CI_upper"]))}]', pf(f(r["p_value"])), "–"])
for r in comb:
    if r["Hypothesis"] != "Effect-magnitude":
        rows.append(["Risk-direction-concentration (OR scale)", tn[r["Tier"]], f'OR {f(r["OR"]):.2f}', f'[{f(r["OR_CI_lower"]):.2f}, {f(r["OR_CI_upper"]):.2f}]', pf(f(r["p_value"])), pf(f(r["q_value_BH"]))])
tables["Table S10"] = write("TableS10_SNPLevel_Models.csv", ["Hypothesis", "Tier", "Estimate", "95% CI", "p-value", "q-value (BH)"], rows)
zrows = [r for r in maf if r["Model"].startswith("Split (d)")]

# ---------- S11 permutation ----------
perm = rd(PW / "AllPathways_Permutation_Test_Results.csv")
rows = [[meta(r["pathway"])["name"], r["n_snps_pathway"], f'{f(r["observed_stat"]):.2f}', f'{f(r["null_mean"]):.2f} ± {f(r["null_sd"]):.2f}', f'{f(r["p_value"]):.3f}'] for r in perm]
tables["Table S11"] = write("TableS11_Permutation.csv", ["Pathway", "N SNPs", "Observed |log(OR per SD)|", "Null distribution, mean ± SD", "Competitive p-value"], rows)

# ---------- key numbers for text ----------
K = {}
K["n_fdr"] = flow_n[2]; K["n_robust"] = len(robust)
allOR = [f(r["Multivariate_OR"]) for r in track]
K["OR_all_range"] = (min(allOR), max(allOR)); K["OR_all_median"] = st.median(allOR)
rob_or = [(f(dr2[p]["OR_train"]), f(dr2[p]["CI_low_train"]), f(dr2[p]["CI_high_train"])) for p in robust]
lo = min(rob_or); hi = max(rob_or); K["OR_robust_min"] = lo; K["OR_robust_max"] = hi
te = [f(dr2[p]["OR_test"]) for p in robust]; K["OR_test_range"] = (min(te), max(te))
K["test_ci_all_include_1"] = all(f(dr2[p]["CI_low_test"]) < 1 < f(dr2[p]["CI_high_test"]) for p in robust)
# AUC from master tables
aucs = []
for p in T:
    mf = PW / f"Pathway_{p}_Multivariate_Master_Table_Results_All_Configurations.csv"
    if mf.exists():
        for r in rd(mf):
            if r["Analysis"] == "main_ld_oriented": aucs.append(f(r["AUC"]))
K["AUC_robust_range"] = (min(aucs), max(aucs), st.median(aucs)) if aucs else None
K["AUC_n"] = len(aucs)
base = [f(dr2[p]["R2_test_covonly"]) for p in robust]; dtr = [f(dr2[p]["Delta_R2_train"]) for p in robust]
dte = [f(dr2[p]["Delta_R2_test"]) for p in robust]
rel = [f(dr2[p]["Delta_R2_test"]) / f(dr2[p]["R2_test_covonly"]) for p in robust]
f2 = [f(dr2[p]["Delta_R2_test"]) / (1 - f(dr2[p]["R2_test_full"])) for p in robust]
K["R2_base_test"] = (st.mean(base), min(base), max(base)); K["dR2_train_mean"] = st.mean(dtr)
K["dR2_test"] = (st.mean(dte), st.median(dte), min(dte), max(dte)); K["rel"] = (st.mean(rel), st.median(rel), min(rel), max(rel)); K["f2_mean"] = st.mean(f2)
# shrinkage functional vs main
ret = [math.log(S7[p]["Functional-VEP"][0]) / math.log(S7[p]["Main model"][0]) for p in robust if "Functional-VEP" in S7[p]]
K["functional_retained"] = (st.mean(ret), st.stdev(ret)); K["functional_shrinkage"] = (1 - st.mean(ret), st.stdev(ret))
K["functional_sig"] = sum(1 for p in robust if "Functional-VEP" in S7[p] and S7[p]["Functional-VEP"][3] < 0.05)
K["functional_same_dir"] = sum(1 for p in robust if "Functional-VEP" in S7[p] and S7[p]["Functional-VEP"][0] > 1)
lm = [math.log(S7[p]["Main model"][0]) for p in robust]; lf = [math.log(S7[p]["Functional-VEP"][0]) for p in robust]
ns = [meta(p)["nsnp"] for p in robust]
def corr(a, b):
    ma, mb = st.mean(a), st.mean(b); return sum((x-ma)*(y-mb) for x, y in zip(a, b)) / math.sqrt(sum((x-ma)**2 for x in a) * sum((y-mb)**2 for y in b))
K["r_func_main"] = corr(lm, lf); K["r_main_nsnp"] = corr(lm, ns); K["r_func_nsnp"] = corr(lf, ns)
lof = [p for p in robust if "LoF-Burden" in S7[p]]
K["lof_n"] = len(lof); K["lof_sig"] = sum(S7[p]["LoF-Burden"][3] < 0.05 for p in lof); K["lof_concord"] = sum(S7[p]["LoF-Burden"][0] > 1 for p in lof)
K["unpruned_identical"] = all(abs(S7[p]["Unpruned"][0] - S7[p]["Main model"][0]) < 1e-9 for p in robust)
K["en_n"] = sum("Elastic-net" in S7[p] for p in robust)
K["beta_ridge_all_sig_risk"] = all(S7[p][v][3] < 0.05 and S7[p][v][0] > 1 for p in robust for v in ("Unpruned", "Beta-weighted", "Ridge-weighted") if v in S7[p])
# covariates
def rng(key): xs = [c[0] for c in covstats[key]]; return (min(xs), max(xs))
def nsig(key): return sum(c[3] < 0.05 for c in covstats[key])
for k in ("SNP_SumScore", "school_cert", "sample_age", "sports_activity", "everyday_activity", "energy_kcal", "graduation", "gender"):
    K["cov_" + k] = (rng(k), nsig(k), len(covstats[k]))
K["energy_sig_pathways"] = [meta(p)["name"] for p, c in zip(order, covstats["energy_kcal"]) if c[3] < 0.05]
# permutation
pv = [f(r["p_value"]) for r in perm]; K["perm_p"] = (min(pv), max(pv))
# external
ext = rd(EXT / "AllPathways_ExternalOrientation_Score_Results_Stage1.csv")
exr = [r for r in ext if r["Pathway"] in robust and r["Status"] == "OK"]
K["ext_testable"] = [meta(r["Pathway"])["name"] for r in exr]
K["ext_or"] = [(meta(r["Pathway"])["name"], f(r["Oriented_OR"]), f(r["Oriented_CI_lower"]), f(r["Oriented_CI_upper"]), f(r["Oriented_Pvalue"]),
               f(r["BetaWeighted_OR"]), f(r["BetaWeighted_CI_lower"]), f(r["BetaWeighted_CI_upper"]), f(r["BetaWeighted_Pvalue"])) for r in exr]
ep = rd(EXT / "AllPathways_External_Permutation_Test_Results.csv")
K["ext_perm"] = [(meta(r["pathway"])["name"], f(r["external_or_ref"]), f(r["p_value"])) for r in ep]
# jaccard pairs
cl = rd(PW / "AllPathways_Gene_Overlap_Jaccard_RobustCandidates_Clusters.csv")
from collections import defaultdict
cc = defaultdict(list)
for r in cl: cc[r["Cluster_ID"]].append(meta(r["Pathway"])["name"])
K["jaccard_clusters"] = len(cc); K["jaccard_pairs"] = [v for v in cc.values() if len(v) > 1]
# batch sensitivity
bd = []
for p in robust:
    for r in rd(PW / f"Pathway_{p}_Multivariate_Master_Table_Results_All_Configurations.csv"):
        if r["Analysis"] == "main_ld_oriented" and r["Delta_percent_Batch"] not in ("NA", ""):
            bd.append(abs(f(r["Delta_percent_Batch"])))
K["batch_absdelta_pct"] = (st.mean(bd), max(bd)) if bd else None
K["Q6"] = {"higher": [r for r in comb if r["Tier"] == "Missense_Higher" and r["Hypothesis"] == "Effect-magnitude"][0],
           "zrows": zrows}
json.dump(K, open(OUT.parent / "key_numbers.json", "w"), indent=1, default=str)

# ---------- Excel workbook ----------
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment
wb = Workbook(); wb.remove(wb.active)
for name, (hdr, rws) in tables.items():
    ws = wb.create_sheet(name.replace("Table ", "T"))
    ws.append([name + " (REPAIR data, generated 2026-11-03)"]); ws["A1"].font = Font(bold=True)
    ws.append(hdr)
    for c in ws[2]: c.font = Font(bold=True)
    for r in rws: ws.append(r)
    for col in ws.columns:
        ws.column_dimensions[col[1].column_letter].width = min(48, max(10, max(len(str(c.value or "")) for c in col[1:]) + 2))
wb.save(OUT / "Repair_Tables_All.xlsx")
print("tables:", list(tables))

# ---------- figures ----------
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
plt.rcParams["font.family"] = "DejaVu Sans"
# Figure 3 forest
names = [meta(p)["name"] for p in order]; o = [f(dr2[p]["OR_train"]) for p in order]
l = [f(dr2[p]["CI_low_train"]) for p in order]; u = [f(dr2[p]["CI_high_train"]) for p in order]
fig, ax = plt.subplots(figsize=(8.5, 0.42 * len(order) + 1.6))
y = list(range(len(order)))[::-1]
ax.errorbar(o, y, xerr=[[a - b for a, b in zip(o, l)], [c - a for a, c in zip(o, u)]], fmt="o", color="#B22222", ecolor="#B22222", capsize=3, ms=6, lw=1.6)
ax.axvline(1, ls="--", color="grey", lw=1)
ax.set_xscale("log"); ax.set_xticks([1, 2, 3, 5, 7, 10]); ax.set_xticklabels(["1", "2", "3", "5", "7", "10"])
ax.set_yticks(y); ax.set_yticklabels(names); ax.grid(axis="x", color="#e5e5e5"); ax.set_axisbelow(True)
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.set_xlabel("Odds ratio per SD of pathway-level SNP sum score, 95% CI (log scale)")
ax.set_title("Robust candidate pathways — multivariate association with obesity")
plt.tight_layout()
for ext_ in ("png", "pdf", "tiff"):
    plt.savefig(OUT / f"Figure3_ForestPlot_RobustCandidates_REPAIR.{ext_}", dpi=300, facecolor="white")
plt.close()
# Figure 2 flow
stages = [f"{flow_n[0]} pathways selected\n(KEGG-first, two-zone model)", f"{flow_n[1]} pathways with complete\ncached association results",
          f"{flow_n[2]} pathways\nq < 0.05 (BH-FDR)", f"{flow_n[3]} pathways\nno significant effect reversal in\ncore sensitivity analyses",
          f"{flow_n[4]} pathways\ntest-set direction = training direction", f"{flow_n[5]} robust candidate pathways\nOR ≥ 1.10 per SD (risk direction)"]
fig, ax = plt.subplots(figsize=(6.8, 11)); ax.set_xlim(0, 100); ax.set_ylim(0, 112); ax.axis("off")
ax.text(50, 109, "Robust-candidate selection flow", ha="center", fontsize=14)
H = 12; gap = 5.5; top = 102
for i, s in enumerate(stages):
    yb = top - (i + 1) * H - i * gap
    ax.add_patch(FancyBboxPatch((14, yb), 64, H, boxstyle="round,pad=0.25,rounding_size=0.8", fc="#EEF3F7" if i < 2 else "#F5F5F5", ec="#222222", lw=1.4))
    ax.text(46, yb + H / 2, s, ha="center", va="center", fontsize=10.5, linespacing=1.35)
    if i < len(stages) - 1:
        ax.add_patch(FancyArrowPatch((46, yb - 0.3), (46, yb - gap + 0.3), arrowstyle="-|>", mutation_scale=15, lw=1.6, color="#222222"))
        if i >= 1:
            ym = yb - gap / 2
            ax.plot([46, 80], [ym, ym], color="#B22222", lw=1.2)
            ax.text(81, ym, f"{flow_x[i+1]} excluded", va="center", fontsize=9.5, color="#B22222")
for ext_ in ("png", "pdf", "tiff"):
    plt.savefig(OUT / f"Figure2_RobustCandidate_Flow_REPAIR.{ext_}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print("done")
