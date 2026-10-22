"""
02_mutation_landscape.py  -  ANALYSIS 1: how common are TP53 mutations, and what kinds?

Inputs:   data/processed/phase1_final_table.csv
          data/processed/tp53_mutation_records.csv
Outputs:  results/figures/fig1_study_flow.png
          results/figures/fig2_tp53_landscape.png
          results/tables/mutation_record_counts.csv
          results/tables/patient_category_counts.csv

Run from the project folder:   python src/02_mutation_landscape.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle

PROJECT_DIR = Path(__file__).resolve().parents[1]
PROCESSED_DIR = PROJECT_DIR / "data" / "processed"
FIG_DIR = PROJECT_DIR / "results" / "figures"
TAB_DIR = PROJECT_DIR / "results" / "tables"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TAB_DIR.mkdir(parents=True, exist_ok=True)

COLORS = {"Wildtype": "#2c6fbb", "Missense": "#e07020", "Nonsense": "#c62828",
          "Frameshift": "#1b9e77", "Splice": "#7b3fa0", "Other/mixed": "#8a8a8a",
          "Splice_Region": "#b07cc6", "In_frame_Indel": "#4a9a6a"}
TRUNCATING = {"Nonsense", "Frameshift", "Splice"}

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False})


def patient_category(types: str) -> str:
    """Fine-grained patient category used for Analysis 1 and 2."""
    parts = set(types.split(";"))
    if types == "Wildtype":
        return "Wildtype"
    if len(parts) == 1 and parts <= {"Missense", "Nonsense", "Frameshift", "Splice"}:
        return parts.pop()
    return "Other/mixed"   # in-frame indels, splice-region only, or >1 mutation class


def mutation_tables(final: pd.DataFrame, records: pd.DataFrame):
    """Record-level and patient-level counts (records != patients!)."""
    rec = (records["Variant_Classification"].value_counts()
           .rename_axis("Variant_Classification").reset_index(name="n_records"))
    rec["percent_of_records"] = (100 * rec["n_records"] / rec["n_records"].sum()).round(1)
    rec.to_csv(TAB_DIR / "mutation_record_counts.csv", index=False)

    final = final.copy()
    final["category"] = final["TP53_mutation_types_found"].map(patient_category)
    order = ["Wildtype", "Missense", "Nonsense", "Frameshift", "Splice", "Other/mixed"]
    pat = final["category"].value_counts().reindex(order).rename_axis("category").reset_index(name="n_patients")
    pat["percent_of_patients"] = (100 * pat["n_patients"] / len(final)).round(1)
    pat.to_csv(TAB_DIR / "patient_category_counts.csv", index=False)

    n_mut = (final["TP53_status"] == "Mutant").sum()
    print(f"TP53-mutant patients: {n_mut}/{len(final)} ({100 * n_mut / len(final):.1f}%)")
    print("\nRecord level:\n", rec.to_string(index=False))
    print("\nPatient level:\n", pat.to_string(index=False))
    return rec, pat


def plot_landscape(pat: pd.DataFrame, records: pd.DataFrame):
    fig = plt.figure(figsize=(11, 4.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.9], wspace=0.25)

    # --- A: patients per category
    ax = fig.add_subplot(gs[0])
    bars = ax.bar(pat["category"], pat["n_patients"],
                  color=[COLORS[c] for c in pat["category"]], edgecolor="white")
    for b, n, p in zip(bars, pat["n_patients"], pat["percent_of_patients"]):
        ax.text(b.get_x() + b.get_width() / 2, n + 4, f"{n}\n({p:.0f}%)", ha="center", va="bottom", fontsize=8.5)
    ax.set_ylabel("Number of patients")
    ax.set_ylim(0, pat["n_patients"].max() * 1.22)
    ax.set_title("A. TP53 status of 566 LUAD patients", loc="left", fontsize=11)
    ax.tick_params(axis="x", rotation=35)
    for lbl in ax.get_xticklabels():
        lbl.set_ha("right")

    # --- B: lollipop plot along the p53 protein
    ax = fig.add_subplot(gs[1])
    rec = records[records["protein_altering"]].copy()
    rec["pos"] = pd.to_numeric(rec["Protein_position"].astype(str).str.extract(r"(\d+)")[0], errors="coerce")
    rec = rec.dropna(subset=["pos"])
    rec["cls"] = rec["TP53_mutation_type"].where(rec["TP53_mutation_type"].isin(
        ["Missense"] + sorted(TRUNCATING)), "Other")
    rec["plot_cls"] = rec["cls"].replace({"Nonsense": "Truncating/Splice", "Frameshift": "Truncating/Splice",
                                          "Splice": "Truncating/Splice"})
    counts = rec.groupby(["pos", "plot_cls"]).size().reset_index(name="n")
    pcol = {"Missense": COLORS["Missense"], "Truncating/Splice": COLORS["Nonsense"], "Other": "#8a8a8a"}
    for cls, sub in counts.groupby("plot_cls"):
        ax.vlines(sub["pos"], 0, sub["n"], color="#bbbbbb", lw=0.8, zorder=1)
        ax.scatter(sub["pos"], sub["n"], s=22 + 6 * sub["n"], color=pcol[cls], edgecolor="white",
                   lw=0.5, zorder=2, label=f"{cls} ({sub['n'].sum()})")
    domains = [(1, 61, "TAD", "#d9e6f5"), (64, 92, "PRD", "#e8e8e8"),
               (94, 292, "DNA-binding domain", "#fde3cf"), (323, 356, "TET", "#d8efe0")]
    y0 = -1.6
    ax.add_patch(Rectangle((1, y0 - 0.3), 392, 0.6, fc="#eeeeee", ec="none", zorder=0))
    for s, e, name, col in domains:
        ax.add_patch(Rectangle((s, y0 - 0.6), e - s, 1.2, fc=col, ec="#666666", lw=0.6, zorder=1))
        ax.text((s + e) / 2, y0, name, ha="center", va="center", fontsize=8)
    # label the most frequently hit residues (residue letter taken from HGVSp_Short, e.g. p.R273L -> R273)
    top = counts.sort_values("n", ascending=False).head(5)
    for pos, cls, n in top[["pos", "plot_cls", "n"]].itertuples(index=False):
        here = rec[(rec["pos"] == pos) & (rec["plot_cls"] == cls)]
        if here["TP53_mutation_type"].eq("Splice_Region").any():
            label = f"c.375 splice\n(exon 4 end)"
        else:
            change = here["HGVSp_Short"].astype(str).str.extract(r"p\.([A-Z]\d+)")[0].mode()
            label = change.iloc[0] if len(change) else f"{int(pos)}"
        dx = {249: -10, 273: 10, 125: -24, 158: 8}.get(int(pos), 0)
        ax.annotate(label, (pos, n), xytext=(dx, 7), textcoords="offset points",
                    ha="center", va="bottom", fontsize=7.5)
    ax.set_xlim(0, 395)
    ax.set_ylim(y0 - 1, counts["n"].max() + 5.5)
    ax.set_xlabel("p53 amino-acid position (1-393)")
    ax.set_ylabel("Mutation records")
    ax.set_yticks(range(0, int(counts["n"].max()) + 2, 2))
    ax.spines["left"].set_bounds(0, counts["n"].max() + 1)
    ax.legend(frameon=False, fontsize=8.5, loc="upper left", title="Mutation class (records)", title_fontsize=8.5)
    ax.set_title("B. Where TP53 mutations fall on the p53 protein", loc="left", fontsize=11)

    fig.savefig(FIG_DIR / "fig2_tp53_landscape.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_study_flow(final: pd.DataFrame):
    """CONSORT-style diagram of how many patients enter each analysis."""
    final = final.copy()
    final["category"] = final["TP53_mutation_types_found"].map(patient_category)
    has_os = final[["OS_MONTHS", "OS_STATUS"]].notna().all(axis=1)
    a = final[has_os]
    stage_ok = a["AJCC_PATHOLOGIC_TUMOR_STAGE"].astype(str).str.upper().str.contains(r"STAGE I")
    n_wt = (a["TP53_status"] == "Wildtype").sum()
    n_mut = (a["TP53_status"] == "Mutant").sum()
    n_mis = (a["category"] == "Missense").sum()
    n_tr = a["category"].isin(["Nonsense", "Frameshift", "Splice"]).sum()
    n_oth = (a["category"] == "Other/mixed").sum()

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

    def box(x, y, w, h, text, fc="#f4f7fb", ec="#2c6fbb", bold=False):
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.4,rounding_size=1.5",
                                    fc=fc, ec=ec, lw=1.2))
        ax.text(x, y, text, ha="center", va="center", fontsize=9.5, weight="bold" if bold else "normal")

    def arrow(x1, y1, x2, y2):
        ax.annotate("", (x2, y2), (x1, y1), arrowprops=dict(arrowstyle="-|>", color="#555555", lw=1))

    def note(y, text, x=50):
        ax.text(x, y, text, ha="center", va="center", fontsize=8.8, style="italic", color="#444444")

    n_deaths = int(a["OS_STATUS"].str.startswith("1").sum())
    box(50, 94, 60, 8, f"TCGA LUAD PanCancer Atlas (cBioPortal)\n{len(final)} patients, all with whole-exome mutation data", bold=True)
    arrow(50, 89.6, 50, 79.2)
    ax.plot([50, 67], [84.5, 84.5], color="#555555", lw=1)
    box(82, 84.5, 29, 6.5, f"Excluded: {(~has_os).sum()} patients\nwithout overall survival data", fc="#fbf1f1", ec="#c62828")
    box(50, 74.5, 40, 7.5, f"{len(a)} patients analysed\n({n_deaths} deaths)", bold=True)
    arrow(50, 70.4, 27, 63.4); arrow(50, 70.4, 73, 63.4)
    box(27, 58.5, 32, 8, f"TP53 wildtype\nn = {n_wt}", fc="#e8f0fa")
    box(73, 58.5, 32, 8, f"TP53 mutant\nn = {n_mut}", fc="#fdeee4", ec="#e07020")
    note(48.5, f"Analysis 3A:\nwildtype ({n_wt}) vs mutant ({n_mut})", x=22)
    arrow(73, 54.2, 49, 42.6); arrow(73, 54.2, 72, 42.6); arrow(73, 54.2, 92, 42.6)
    box(49, 38, 19, 8, f"Missense only\nn = {n_mis}", fc="#fdeee4", ec="#e07020")
    box(71.5, 38, 21, 8, f"Truncating/Splice\nonly, n = {n_tr}", fc="#fbe9e9", ec="#c62828")
    box(92, 38, 13, 8, f"Other/mixed\nn = {n_oth}", fc="#f1f1f1", ec="#8a8a8a")
    note(27, f"Analysis 3B: wildtype ({n_wt}) vs missense ({n_mis}) vs truncating/splice ({n_tr})\n"
               f"Other/mixed ({n_oth}) set aside: too small and biologically mixed")
    box(50, 13, 78, 8, f"Analysis 4 (stage-adjusted Cox model): {int(stage_ok.sum())} patients with known AJCC stage\n"
                      f"({len(a) - int(stage_ok.sum())} with missing or discrepant stage excluded)", fc="#f6f6f6", ec="#555555")
    fig.savefig(FIG_DIR / "fig1_study_flow.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    final = pd.read_csv(PROCESSED_DIR / "phase1_final_table.csv")
    records = pd.read_csv(PROCESSED_DIR / "tp53_mutation_records.csv")
    _, pat = mutation_tables(final, records)
    plot_landscape(pat, records)
    plot_study_flow(final)
    print(f"\nFigures saved to {FIG_DIR}")


if __name__ == "__main__":
    main()
