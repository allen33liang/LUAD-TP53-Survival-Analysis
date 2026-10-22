"""
03_survival_analysis.py  -  ANALYSES 2-4: TP53 mutation type and overall survival

  Analysis 2  Overall survival in each fine-grained TP53 group (are groups big enough?)
  Analysis 3  Kaplan-Meier curves + log-rank tests
                3A wildtype vs mutant
                3B wildtype vs missense vs truncating/splice
  Analysis 4  Cox proportional-hazards models:  Survival ~ TP53 group + Stage
  Checks      proportional-hazards test, sensitivity analyses

Inputs:   data/processed/phase1_final_table.csv, data/processed/tp53_mutation_records.csv
Outputs:  results/figures/*.png, results/tables/*.csv, results/survival_results.json

Run from the project folder:   python src/03_survival_analysis.py
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import logrank_test, multivariate_logrank_test, proportional_hazard_test
from lifelines.utils import median_survival_times

PROJECT_DIR = Path(__file__).resolve().parents[1]
PROCESSED_DIR = PROJECT_DIR / "data" / "processed"
FIG_DIR = PROJECT_DIR / "results" / "figures"
TAB_DIR = PROJECT_DIR / "results" / "tables"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TAB_DIR.mkdir(parents=True, exist_ok=True)

COLORS = {"Wildtype": "#2c6fbb", "Mutant": "#c62828", "Missense": "#e07020",
          "Truncating/Splice": "#c62828", "Nonsense": "#c62828", "Frameshift": "#1b9e77",
          "Splice": "#7b3fa0", "Other/mixed": "#8a8a8a",
          "I": "#1b7837", "II": "#5aae61", "III": "#e08214", "IV": "#b2182b"}
THREE_GROUPS = ["Wildtype", "Missense", "Truncating/Splice"]
SHORT = {"Truncating/Splice": "Trunc/Splice"}
TRUNCATING = {"Nonsense", "Frameshift", "Splice"}

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False})
RESULTS = {}   # everything the report quotes is collected here and saved as JSON


# ============================================================ data preparation
def three_group(types: str) -> str:
    """Pre-specified grouping (decided after looking at group SIZES, not at survival):
    Missense-only | Truncating/Splice-only (nonsense, frameshift, splice site) | Other/mixed."""
    if types == "Wildtype":
        return "Wildtype"
    parts = set(types.split(";"))
    if parts == {"Missense"}:
        return "Missense"
    if parts <= TRUNCATING:
        return "Truncating/Splice"
    return "Other/mixed"


def fine_group(types: str) -> str:
    parts = set(types.split(";"))
    if types == "Wildtype":
        return "Wildtype"
    if len(parts) == 1 and parts <= {"Missense"} | TRUNCATING:
        return parts.pop()
    return "Other/mixed"


def stage_group(raw) -> str:
    """'Stage IIIA' / 'STAGE IIIA' -> 'III'. Missing, [Discrepancy] etc. -> NaN."""
    if pd.isna(raw):
        return np.nan
    found = pd.Series([str(raw).upper()]).str.extract(r"STAGE\s+(IV|III|II|I)")[0].iloc[0]
    return found if isinstance(found, str) else np.nan


def load_data() -> pd.DataFrame:
    df = pd.read_csv(PROCESSED_DIR / "phase1_final_table.csv")
    n_all = len(df)
    df = df.dropna(subset=["OS_MONTHS", "OS_STATUS"]).copy()
    df["time"] = df["OS_MONTHS"]
    df["event"] = df["OS_STATUS"].str.startswith("1").astype(int)   # 1:DECEASED -> 1
    df["group3"] = df["TP53_mutation_types_found"].map(three_group)
    df["fine_group"] = df["TP53_mutation_types_found"].map(fine_group)
    df["stage"] = df["AJCC_PATHOLOGIC_TUMOR_STAGE"].map(stage_group)
    RESULTS["n_all"] = n_all
    RESULTS["n_analysed"] = len(df)
    RESULTS["n_deaths"] = int(df["event"].sum())
    print(f"Patients with survival data: {len(df)} of {n_all} ({df['event'].sum()} deaths)")
    return df


# ============================================================ helpers
def km_fit(df, label):
    return KaplanMeierFitter().fit(df["time"], df["event"], label=label)


def km_summary(df, label):
    """n, deaths, median OS (95% CI) and 3-year OS (95% CI) for one group."""
    kmf = km_fit(df, label)
    med_ci = median_survival_times(kmf.confidence_interval_)
    ci_3y = kmf.confidence_interval_survival_function_.loc[:36].iloc[-1]
    return {"group": label, "n": len(df), "deaths": int(df["event"].sum()),
            "median_OS_months": float(kmf.median_survival_time_),
            "median_CI_low": float(med_ci.iloc[0, 0]), "median_CI_high": float(med_ci.iloc[0, 1]),
            "OS_3yr_pct": round(100 * float(kmf.predict(36)), 1),
            "OS_3yr_CI_low": round(100 * float(ci_3y.iloc[0]), 1),
            "OS_3yr_CI_high": round(100 * float(ci_3y.iloc[1]), 1)}


def fit_cox(df, covariates, label):
    """Cox model with explicit 0/1 dummy columns (reference groups are left out)."""
    data = df[["time", "event"] + covariates]
    cph = CoxPHFitter().fit(data, duration_col="time", event_col="event")
    s = cph.summary
    rows = [{"model": label, "n": len(data), "events": int(data["event"].sum()), "term": term,
             "HR": s.loc[term, "exp(coef)"], "CI_low": s.loc[term, "exp(coef) lower 95%"],
             "CI_high": s.loc[term, "exp(coef) upper 95%"], "p": s.loc[term, "p"]} for term in covariates]
    return cph, data, pd.DataFrame(rows)


def add_dummies(df):
    df = df.copy()
    df["Missense_vs_WT"] = (df["group3"] == "Missense").astype(int)
    df["TruncSplice_vs_WT"] = (df["group3"] == "Truncating/Splice").astype(int)
    df["Mutant_vs_WT"] = (df["TP53_status"] == "Mutant").astype(int)
    for st in ["II", "III", "IV"]:                       # Stage I is the reference
        df[f"Stage_{st}_vs_I"] = (df["stage"] == st).astype(int)
    df["Age_per_10yr"] = (df["AGE"] - df["AGE"].mean()) / 10
    df["Male_vs_Female"] = (df["SEX"] == "Male").astype(int)
    return df


def fmt_p(p):
    return "<0.001" if p < 0.001 else f"{p:.3f}"


def p_text(p):
    """'p < 0.001' or 'p = 0.043' for figure labels."""
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}"


# ============================================================ descriptive table
def patient_characteristics(df):
    """Table 1: who is in each group? (checks whether stage/age differ between groups)"""
    rows = []
    groups = [("All patients", df)] + [(g, df[df["group3"] == g]) for g in THREE_GROUPS + ["Other/mixed"]]
    for name, d in groups:
        st = d["stage"].value_counts()
        rows.append({
            "group": name, "n": len(d),
            "age_median": round(d["AGE"].median(), 0),
            "age_IQR": f"{d['AGE'].quantile(.25):.0f}-{d['AGE'].quantile(.75):.0f}",
            "male_pct": round(100 * (d["SEX"] == "Male").mean(), 1),
            **{f"stage_{s}": int(st.get(s, 0)) for s in ["I", "II", "III", "IV"]},
            "stage_unknown": int(d["stage"].isna().sum()),
            "stage_III_IV_pct": round(100 * d["stage"].isin(["III", "IV"]).sum() / d["stage"].notna().sum(), 1),
            "deaths": int(d["event"].sum()),
        })
    table = pd.DataFrame(rows)
    table.to_csv(TAB_DIR / "patient_characteristics.csv", index=False)

    # chi-square test: is stage distribution different between the 3 groups?
    from scipy.stats import chi2_contingency, kruskal
    sub = df[df["group3"].isin(THREE_GROUPS)]
    ct = pd.crosstab(sub["group3"], sub["stage"])
    chi2, p_stage, _, _ = chi2_contingency(ct)
    p_age = kruskal(*[sub.loc[sub.group3 == g, "AGE"].dropna() for g in THREE_GROUPS]).pvalue
    RESULTS["stage_by_group_chi2_p"] = p_stage
    RESULTS["age_by_group_kruskal_p"] = p_age
    rf = KaplanMeierFitter().fit(df["time"], 1 - df["event"])           # reverse KM
    RESULTS["median_followup_months"] = float(rf.median_survival_time_)
    print("\nTable 1 - patient characteristics\n", table.to_string(index=False))
    print(f"Stage distribution across 3 groups: chi-square p = {p_stage:.3f}; age: Kruskal-Wallis p = {p_age:.3f}")
    print(f"Median follow-up (reverse Kaplan-Meier): {rf.median_survival_time_:.1f} months")
    return table


# ============================================================ Analysis 2
def analysis2_fine_groups(df):
    order = ["Wildtype", "Missense", "Nonsense", "Frameshift", "Splice", "Other/mixed"]
    table = pd.DataFrame([km_summary(df[df["fine_group"] == g], g) for g in order])
    table.to_csv(TAB_DIR / "os_by_fine_group.csv", index=False)
    RESULTS["fine_groups"] = table.to_dict(orient="records")
    print("\nAnalysis 2 - OS by fine-grained TP53 group\n", table.round(1).to_string(index=False))

    # Supplementary figure: why small groups are unstable
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    for g in order:
        d = df[df["fine_group"] == g]
        km_fit(d, f"{g} (n={len(d)})").plot_survival_function(
            ax=ax, ci_show=False, color=COLORS[g], lw=1.8, show_censors=True,
            censor_styles={"ms": 4, "marker": "|"})
    ax.set_xlim(0, 150); ax.set_ylim(0, 1.02)
    ax.set_xlabel("Overall survival (months)"); ax.set_ylabel("Survival probability")
    ax.set_title("Overall survival by fine-grained TP53 mutation type", loc="left")
    ax.legend(frameon=False, fontsize=8.5, loc="lower left")
    ax.grid(alpha=0.25)
    fig.savefig(FIG_DIR / "figS1_km_fine_groups.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return table


# ============================================================ Analysis 3
def analysis3_km(df):
    # --- 3A: wildtype vs mutant
    wt, mut = df[df.TP53_status == "Wildtype"], df[df.TP53_status == "Mutant"]
    lr = logrank_test(wt.time, mut.time, wt.event, mut.event)
    RESULTS["3A"] = {"Wildtype": km_summary(wt, "Wildtype"), "Mutant": km_summary(mut, "Mutant"),
                     "logrank_chi2": lr.test_statistic, "logrank_p": lr.p_value}
    print(f"\nAnalysis 3A - log-rank wildtype vs mutant: chi2 = {lr.test_statistic:.2f}, p = {lr.p_value:.3f}")

    # --- 3B: three groups
    d3 = df[df.group3.isin(THREE_GROUPS)]
    mlr = multivariate_logrank_test(d3.time, d3.group3, d3.event)
    pairs = {}
    for a, b in [("Wildtype", "Missense"), ("Wildtype", "Truncating/Splice"), ("Missense", "Truncating/Splice")]:
        x, y = d3[d3.group3 == a], d3[d3.group3 == b]
        pairs[f"{a} vs {b}"] = logrank_test(x.time, y.time, x.event, y.event).p_value
    RESULTS["3B"] = {"groups": {g: km_summary(d3[d3.group3 == g], g) for g in THREE_GROUPS},
                     "logrank_chi2": mlr.test_statistic, "logrank_p": mlr.p_value, "pairwise_p": pairs}
    print(f"Analysis 3B - log-rank 3 groups: chi2 = {mlr.test_statistic:.2f}, p = {mlr.p_value:.3f}")
    for k, v in pairs.items():
        print(f"   pairwise {k}: p = {v:.3f} (unadjusted for multiple testing)")

    km_table = pd.DataFrame([RESULTS["3A"]["Wildtype"], RESULTS["3A"]["Mutant"]] +
                            [RESULTS["3B"]["groups"][g] for g in ["Missense", "Truncating/Splice"]])
    km_table.to_csv(TAB_DIR / "km_summary.csv", index=False)

    # --- Figure 3
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.4))
    panels = [(axes[0], [("Wildtype", wt), ("Mutant", mut)], "A. TP53 wildtype vs mutant", lr.p_value,
               {"Wildtype": "TP53 wildtype", "Mutant": "TP53 mutant"}),
              (axes[1], [(g, d3[d3.group3 == g]) for g in THREE_GROUPS], "B. By TP53 mutation type",
               mlr.p_value, {g: g for g in THREE_GROUPS})]
    for ax, groups, title, p, names in panels:
        fitters, legend_names = [], []
        for g, d in groups:
            k = km_fit(d, SHORT.get(g, g))           # short label = row name in the at-risk table
            k.plot_survival_function(ax=ax, ci_show=False, color=COLORS[g], lw=2, show_censors=True,
                                     censor_styles={"ms": 5, "marker": "|"})
            fitters.append(k)
            legend_names.append(f"{names[g]} (n={len(d)})")
        ax.set_xlim(0, 150); ax.set_ylim(0, 1.02)
        ax.set_xlabel("Overall survival (months)"); ax.set_ylabel("Survival probability")
        ax.set_title(title, loc="left", fontsize=12)
        ax.text(0.97, 0.95, f"log-rank p = {p:.3f}", transform=ax.transAxes, ha="right", fontsize=10.5)
        ax.legend(ax.get_legend_handles_labels()[0], legend_names, frameon=False, loc="lower left", fontsize=9.5)
        ax.grid(alpha=0.25)
        ax.set_xticks(range(0, 151, 30))
        add_at_risk_counts(*fitters, ax=ax, rows_to_show=["At risk"], fontsize=8.5)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig3_km_tp53.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return d3


# ============================================================ Analysis 4
def analysis4_cox(df):
    df = add_dummies(df)
    stage_cols = ["Stage_II_vs_I", "Stage_III_vs_I", "Stage_IV_vs_I"]

    # Stage alone: is stage a strong prognostic factor in this cohort?
    ds = df.dropna(subset=["stage"])
    stage_lr = multivariate_logrank_test(ds.time, ds.stage, ds.event)
    RESULTS["stage_logrank_p"] = stage_lr.p_value
    RESULTS["stage_km"] = {s: km_summary(ds[ds.stage == s], s) for s in ["I", "II", "III", "IV"]}

    tables, models = [], {}
    # ---- binary comparison: all models on the SAME patients (stage known) so they are comparable
    b = ds
    _, _, t = fit_cox(b, ["Mutant_vs_WT"], "Binary: unadjusted"); tables.append(t)
    cph_b, data_b, t = fit_cox(b, ["Mutant_vs_WT"] + stage_cols, "Binary: + stage"); tables.append(t)
    b2 = b.dropna(subset=["AGE"])
    _, _, t = fit_cox(b2, ["Mutant_vs_WT"] + stage_cols + ["Age_per_10yr", "Male_vs_Female"],
                      "Binary: + stage + age + sex"); tables.append(t)

    # ---- three groups (other/mixed set aside)
    g = ds[ds.group3.isin(THREE_GROUPS)]
    tp53_cols = ["Missense_vs_WT", "TruncSplice_vs_WT"]
    _, _, t = fit_cox(g, tp53_cols, "3-group: unadjusted"); tables.append(t)
    cph_g, data_g, t = fit_cox(g, tp53_cols + stage_cols, "3-group: + stage"); tables.append(t)
    models["main"] = (cph_g, data_g)
    g2 = g.dropna(subset=["AGE"])
    _, _, t = fit_cox(g2, tp53_cols + stage_cols + ["Age_per_10yr", "Male_vs_Female"],
                      "3-group: + stage + age + sex"); tables.append(t)

    cox = pd.concat(tables, ignore_index=True)
    cox.to_csv(TAB_DIR / "cox_models.csv", index=False)
    RESULTS["cox"] = cox.to_dict(orient="records")
    print("\nAnalysis 4 - Cox models")
    for _, r in cox.iterrows():
        print(f"  {r.model:32s} {r.term:20s} HR={r.HR:.2f} ({r.CI_low:.2f}-{r.CI_high:.2f}) p={fmt_p(r.p)} n={r.n}")

    # ---- proportional-hazards assumption (Schoenfeld residual test)
    ph = {}
    for name, (cph, data) in {"binary_stage": (cph_b, data_b), "3group_stage": (cph_g, data_g)}.items():
        res = proportional_hazard_test(cph, data, time_transform="rank")
        ph[name] = {term: float(p) for term, p in zip(res.summary.index, res.summary["p"])}
    RESULTS["ph_test"] = ph
    print("\nProportional-hazards test p-values (p < 0.05 = possible violation):")
    for k, v in ph.items():
        print("  ", k, {t: round(p, 3) for t, p in v.items()})

    # ---- Figure 4: stage KM + forest plot
    fig = plt.figure(figsize=(14, 5.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 0.8], wspace=0.62)
    ax = fig.add_subplot(gs[0])
    fitters, legend_names = [], []
    for s in ["I", "II", "III", "IV"]:
        d = ds[ds.stage == s]
        legend_names.append(f"Stage {s} (n={len(d)})")
        k = km_fit(d, f"Stage {s}")
        k.plot_survival_function(ax=ax, ci_show=False, color=COLORS[s], lw=2, show_censors=True,
                                 censor_styles={"ms": 5, "marker": "|"})
        fitters.append(k)
    ax.set_xlim(0, 150); ax.set_ylim(0, 1.02); ax.set_xticks(range(0, 151, 30))
    ax.set_xlabel("Overall survival (months)"); ax.set_ylabel("Survival probability")
    ax.set_title("A. Overall survival by AJCC stage", loc="left", fontsize=12)
    ax.text(0.97, 0.95, f"log-rank {p_text(stage_lr.p_value)}", transform=ax.transAxes, ha="right")
    ax.legend(ax.get_legend_handles_labels()[0], legend_names, frameon=False, loc="lower left", fontsize=9.5)
    ax.grid(alpha=0.25)
    add_at_risk_counts(*fitters, ax=ax, rows_to_show=["At risk"], fontsize=8.5)

    ax = fig.add_subplot(gs[1])
    main = cox[cox.model == "3-group: + stage"].reset_index(drop=True)
    unadj = cox[cox.model == "3-group: unadjusted"].set_index("term")
    labels = {"Missense_vs_WT": "Missense vs wildtype", "TruncSplice_vs_WT": "Truncating/Splice vs wildtype",
              "Stage_II_vs_I": "Stage II vs I", "Stage_III_vs_I": "Stage III vs I", "Stage_IV_vs_I": "Stage IV vs I"}
    ys = np.arange(len(main))[::-1]
    for y, r in zip(ys, main.itertuples()):
        col = COLORS["Missense"] if r.term == "Missense_vs_WT" else (
            COLORS["Truncating/Splice"] if r.term == "TruncSplice_vs_WT" else "#444444")
        ax.plot([r.CI_low, r.CI_high], [y, y], color=col, lw=2)
        ax.plot(r.HR, y, "s", color=col, ms=8)
        if r.term in unadj.index:   # show the unadjusted estimate for comparison
            u = unadj.loc[r.term]
            ax.plot(u.HR, y + 0.25, "o", mfc="white", mec=col, ms=6)
            ax.plot([u.CI_low, u.CI_high], [y + 0.25, y + 0.25], color=col, lw=1, ls=":")
        ax.text(1.03, y, f"{r.HR:.2f} ({r.CI_low:.2f}-{r.CI_high:.2f})   {p_text(r.p)}",
                va="center", fontsize=9, transform=ax.get_yaxis_transform())
    ax.axvline(1, color="#888888", ls="--", lw=1)
    ax.set_xscale("log"); ax.set_xlim(0.5, 8); ax.set_xticks([0.5, 1, 2, 4, 8]); ax.set_xticklabels(["0.5", "1", "2", "4", "8"])
    ax.set_yticks(ys); ax.set_yticklabels([labels[t] for t in main.term])
    ax.set_xlabel("Hazard ratio (log scale, 95% CI)")
    ax.set_title(f"B. Cox model: TP53 type + stage\n(n={int(main.n[0])}, {int(main.events[0])} deaths)", loc="left", fontsize=12)
    ax.text(0.0, -0.17, "■ adjusted for stage      ○ unadjusted (same patients)", transform=ax.transAxes,
            ha="left", fontsize=8.5, color="#444444")
    ax.text(1.03, len(main) - 0.4, "HR (95% CI)", fontsize=9, weight="bold", transform=ax.get_yaxis_transform())
    fig.savefig(FIG_DIR / "fig4_stage_and_cox.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return cox


# ============================================================ sensitivity analyses
def sensitivity(df):
    """Does the answer change if the mutation classification rules change?"""
    rec = pd.read_csv(PROCESSED_DIR / "tp53_mutation_records.csv")

    def regroup(records):
        types = records.groupby("PATIENT_ID")["TP53_mutation_type"].agg(lambda t: ";".join(sorted(set(t))))
        out = df.drop(columns=["group3"]).copy()
        out["types"] = out["PATIENT_ID"].map(types).fillna("Wildtype")
        out["group3"] = out["types"].map(three_group)
        out["TP53_status"] = np.where(out["types"] == "Wildtype", "Wildtype", "Mutant")
        return out

    variants = {
        "Main analysis (Silent excluded; splice-region kept separate)": df,
        "S1: original Phase 1 rules (Silent counted as 'Other')":
            regroup(rec.assign(TP53_mutation_type=rec.TP53_mutation_type.where(rec.protein_altering, "Other"))),
        "S2: exon-boundary splice-region variants counted as Splice":
            regroup(rec[rec.protein_altering].assign(
                TP53_mutation_type=rec.TP53_mutation_type.replace({"Splice_Region": "Splice"}))),
    }
    rows = []
    for name, d in variants.items():
        d = add_dummies(d)
        g = d[d.group3.isin(THREE_GROUPS)]
        mlr = multivariate_logrank_test(g.time, g.group3, g.event)
        wt, mut = d[d.TP53_status == "Wildtype"], d[d.TP53_status == "Mutant"]
        p2 = logrank_test(wt.time, mut.time, wt.event, mut.event).p_value
        gs = g.dropna(subset=["stage"])
        _, _, t = fit_cox(gs, ["Missense_vs_WT", "TruncSplice_vs_WT", "Stage_II_vs_I", "Stage_III_vs_I",
                               "Stage_IV_vs_I"], name)
        t = t.set_index("term")
        rows.append({"analysis": name,
                     "n_wt/mis/trunc": f"{(g.group3 == 'Wildtype').sum()}/{(g.group3 == 'Missense').sum()}/"
                                       f"{(g.group3 == 'Truncating/Splice').sum()}",
                     "logrank_p_WT_vs_mut": p2, "logrank_p_3group": mlr.p_value,
                     "HR_missense_adj": t.loc["Missense_vs_WT", "HR"], "p_missense_adj": t.loc["Missense_vs_WT", "p"],
                     "HR_trunc_adj": t.loc["TruncSplice_vs_WT", "HR"], "p_trunc_adj": t.loc["TruncSplice_vs_WT", "p"]})
    # S3: stratified Cox - each stage gets its own baseline hazard, so the model no longer
    # assumes that the stage effect is constant over time (stage III failed the PH test)
    d = add_dummies(df)
    g = d[d.group3.isin(THREE_GROUPS)].dropna(subset=["stage"])
    cph = CoxPHFitter().fit(g[["time", "event", "Missense_vs_WT", "TruncSplice_vs_WT", "stage"]],
                            duration_col="time", event_col="event", strata=["stage"])
    sm = cph.summary
    rows.append({"analysis": "S3: main classification, Cox stratified by stage",
                 "n_wt/mis/trunc": rows[0]["n_wt/mis/trunc"],
                 "logrank_p_WT_vs_mut": np.nan, "logrank_p_3group": np.nan,
                 "HR_missense_adj": sm.loc["Missense_vs_WT", "exp(coef)"], "p_missense_adj": sm.loc["Missense_vs_WT", "p"],
                 "HR_trunc_adj": sm.loc["TruncSplice_vs_WT", "exp(coef)"], "p_trunc_adj": sm.loc["TruncSplice_vs_WT", "p"]})
    table = pd.DataFrame(rows)
    table.to_csv(TAB_DIR / "sensitivity_analyses.csv", index=False)
    RESULTS["sensitivity"] = table.to_dict(orient="records")
    print("\nSensitivity analyses\n", table.round(3).to_string(index=False))


def main():
    df = load_data()
    patient_characteristics(df)
    analysis2_fine_groups(df)
    analysis3_km(df)
    analysis4_cox(df)
    sensitivity(df)
    with open(PROJECT_DIR / "results" / "survival_results.json", "w") as fh:
        json.dump(RESULTS, fh, indent=2, default=float)
    print(f"\nAll results saved to {PROJECT_DIR / 'results'}")


if __name__ == "__main__":
    main()
