"""
01_build_phase1_table.py  -  PHASE 1: build the patient-level analysis table

Goal: one row per patient with
    PATIENT_ID | survival (OS_MONTHS, OS_STATUS) | AGE | SEX | tumour stage
    | TP53_status | TP53 mutation type(s)

Inputs (download from cBioPortal - see data/raw/README.md):
    data/raw/luad_tcga_pan_can_atlas_2018/data_mutations.txt
    data/raw/luad_tcga_pan_can_atlas_2018/data_clinical_patient.txt

Outputs:
    data/processed/phase1_final_table.csv     (one row per patient)
    data/processed/tp53_mutation_records.csv  (one row per TP53 mutation record)

Run from the project folder:   python src/01_build_phase1_table.py
"""
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------- paths
PROJECT_DIR = Path(__file__).resolve().parents[1]          # the repo folder
RAW_DIR = PROJECT_DIR / "data" / "raw" / "luad_tcga_pan_can_atlas_2018"
MUTATIONS_FILE = RAW_DIR / "data_mutations.txt"
CLINICAL_FILE = RAW_DIR / "data_clinical_patient.txt"

PROCESSED_DIR = PROJECT_DIR / "data" / "processed"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- settings
# Variant classes that do NOT change the p53 protein sequence.
# These are dropped before deciding whether a patient is TP53-mutant.
NON_PROTEIN_ALTERING = {"Silent"}

# Raw MAF Variant_Classification  ->  simplified category
CLASSIFICATION_MAPPING = {
    "Missense_Mutation": "Missense",
    "Nonsense_Mutation": "Nonsense",
    "Frame_Shift_Del": "Frameshift",
    "Frame_Shift_Ins": "Frameshift",
    "Splice_Site": "Splice",
    "In_Frame_Del": "In_frame_Indel",
    "In_Frame_Ins": "In_frame_Indel",
    "Splice_Region": "Splice_Region",   # kept separate so it is visible (see report 3.6)
}

# Clinical columns to keep (stage is needed for the stage-adjusted Cox model)
REQUIRED_CLINICAL = ["PATIENT_ID", "OS_MONTHS", "OS_STATUS"]
OPTIONAL_CLINICAL = ["AGE", "SEX", "AJCC_PATHOLOGIC_TUMOR_STAGE"]


# ---------------------------------------------------------------- loading
def load_table(path: Path, label: str) -> pd.DataFrame:
    """Read a cBioPortal tab-separated file (lines starting with # are metadata)."""
    print(f"Loading {label} from: {path}")
    table = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    print(f"  -> {table.shape[0]} rows, {table.shape[1]} columns")
    return table


# ---------------------------------------------------------------- TP53 steps
def filter_tp53(mutations: pd.DataFrame) -> pd.DataFrame:
    """Keep only TP53 mutation records."""
    required = ["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification", "Variant_Type"]
    missing = [c for c in required if c not in mutations.columns]
    if missing:
        raise ValueError(f"Missing required columns in mutation data: {missing}")

    tp53 = mutations[mutations["Hugo_Symbol"] == "TP53"].copy()
    print(f"\nTP53 mutation records: {tp53.shape[0]}")
    print(f"TP53 unique sample barcodes: {tp53['Tumor_Sample_Barcode'].nunique()}")
    print("\nTP53 Variant_Classification breakdown (raw records):")
    print(tp53["Variant_Classification"].value_counts().to_string())
    return tp53


def add_patient_id(table: pd.DataFrame) -> pd.DataFrame:
    """TCGA sample barcode -> patient ID: keep the first three parts.
    e.g. TCGA-05-4382-01  ->  TCGA-05-4382"""
    table = table.copy()
    table["PATIENT_ID"] = table["Tumor_Sample_Barcode"].str.split("-").str[:3].str.join("-")
    return table


def classify_records(tp53: pd.DataFrame) -> pd.DataFrame:
    """Simplify mutation types and flag records that do not change the protein."""
    tp53 = tp53.copy()
    tp53["TP53_mutation_type"] = (
        tp53["Variant_Classification"].map(CLASSIFICATION_MAPPING).fillna("Other")
    )
    tp53["protein_altering"] = ~tp53["Variant_Classification"].isin(NON_PROTEIN_ALTERING)

    n_dropped = (~tp53["protein_altering"]).sum()
    print(f"\nNon-protein-altering records excluded from patient status: {n_dropped} "
          f"({', '.join(sorted(NON_PROTEIN_ALTERING))})")
    print("Simplified TP53 mutation types (protein-altering records):")
    print(tp53.loc[tp53["protein_altering"], "TP53_mutation_type"].value_counts().to_string())
    return tp53


def build_patient_level_tp53(tp53: pd.DataFrame) -> pd.DataFrame:
    """Collapse mutation RECORDS into one row per PATIENT."""
    altering = tp53[tp53["protein_altering"]]
    n_samples = altering["Tumor_Sample_Barcode"].nunique()
    n_patients = altering["PATIENT_ID"].nunique()
    print(f"\n{n_samples} TP53-mutant samples map to {n_patients} unique patients")
    if n_samples != n_patients:
        print("  NOTE: some patients have >1 tumour sample - check before analysis.")

    patient_level = (
        altering.groupby("PATIENT_ID")
        .agg(
            TP53_mutation_types_found=("TP53_mutation_type", lambda t: ";".join(sorted(set(t)))),
            TP53_n_mutations=("TP53_mutation_type", "size"),
            TP53_protein_changes=("HGVSp_Short", lambda p: ";".join(p.astype(str))),
        )
        .reset_index()
    )
    patient_level["TP53_status"] = "Mutant"
    return patient_level


# ---------------------------------------------------------------- merge
def build_final_table(patient_level: pd.DataFrame, clinical: pd.DataFrame) -> pd.DataFrame:
    """Merge TP53 information onto the clinical table (one row per patient).
    Every patient in this study was sequenced (case list: 566/566 samples with
    mutation data), so a patient with no TP53 record is TP53 wildtype."""
    missing = [c for c in REQUIRED_CLINICAL if c not in clinical.columns]
    if missing:
        raise ValueError(f"Missing required columns in clinical data: {missing}")
    keep = REQUIRED_CLINICAL + [c for c in OPTIONAL_CLINICAL if c in clinical.columns]
    if "AJCC_PATHOLOGIC_TUMOR_STAGE" not in keep:
        print("WARNING: no AJCC_PATHOLOGIC_TUMOR_STAGE column - stage-adjusted model will not run.")

    final = clinical[keep].merge(patient_level, on="PATIENT_ID", how="left")
    final["TP53_status"] = final["TP53_status"].fillna("Wildtype")
    final["TP53_mutation_types_found"] = final["TP53_mutation_types_found"].fillna("Wildtype")
    final["TP53_n_mutations"] = final["TP53_n_mutations"].fillna(0).astype(int)
    final["TP53_protein_changes"] = final["TP53_protein_changes"].fillna("")

    if final["PATIENT_ID"].duplicated().any():
        raise ValueError("Duplicate PATIENT_IDs after merge - check the input files.")

    print(f"\nFinal Phase-1 table: {final.shape[0]} patients")
    print(final["TP53_status"].value_counts().to_string())
    return final


# ---------------------------------------------------------------- QC
def run_qc(mutations, tp53, clinical, final):
    """Phase 1 completion checklist (from the project guide, section 19)."""
    def tick(ok):
        return "YES" if ok else "   "

    has_os = final[["OS_MONTHS", "OS_STATUS"]].notna().all(axis=1)
    print("\n" + "=" * 60)
    print("PHASE 1 QC CHECKLIST")
    print("=" * 60)
    print(f"[{tick(len(mutations) > 0)}] data_mutations.txt loaded ({len(mutations)} rows)")
    print(f"[{tick(len(clinical) > 0)}] data_clinical_patient.txt loaded ({len(clinical)} rows)")
    print(f"[{tick(len(tp53) > 0)}] TP53 mutation records extracted ({len(tp53)} rows)")
    print(f"[{tick(tp53['Tumor_Sample_Barcode'].nunique() > 0)}] TP53 sample IDs identified "
          f"({tp53['Tumor_Sample_Barcode'].nunique()} samples)")
    print(f"[{tick(final['TP53_mutation_types_found'].nunique() > 1)}] TP53 mutation types classified")
    print(f"[{tick(final['PATIENT_ID'].is_unique)}] One row per patient (no duplicates)")
    print(f"[{tick('OS_MONTHS' in final.columns)}] OS_MONTHS present")
    print(f"[{tick('OS_STATUS' in final.columns)}] OS_STATUS present")
    print(f"[{tick('AJCC_PATHOLOGIC_TUMOR_STAGE' in final.columns)}] Tumour stage present")
    print(f"    Patients with complete survival data: {has_os.sum()} of {len(final)}")

    print("\nPatients per mutation-type combination:")
    counts = final["TP53_mutation_types_found"].value_counts()
    print(counts.to_string())
    small = counts[counts < 10]
    if not small.empty:
        print("\nWARNING: these groups have < 10 patients - merge categories before "
              "running survival statistics:")
        print(small.to_string())


# ---------------------------------------------------------------- main
def main():
    for f in (MUTATIONS_FILE, CLINICAL_FILE):
        if not f.exists():
            raise SystemExit(f"Missing raw data file: {f}\n"
                             "Download the study from cBioPortal first - see data/raw/README.md")
    mutations = load_table(MUTATIONS_FILE, "mutation data")
    clinical = load_table(CLINICAL_FILE, "clinical data")

    tp53 = filter_tp53(mutations)
    tp53 = add_patient_id(tp53)
    tp53 = classify_records(tp53)

    patient_level = build_patient_level_tp53(tp53)
    final = build_final_table(patient_level, clinical)
    run_qc(mutations, tp53, clinical, final)

    record_cols = ["PATIENT_ID", "Tumor_Sample_Barcode", "Variant_Classification", "Variant_Type",
                   "TP53_mutation_type", "protein_altering", "HGVSc", "HGVSp_Short", "Protein_position"]
    tp53[[c for c in record_cols if c in tp53.columns]].to_csv(
        PROCESSED_DIR / "tp53_mutation_records.csv", index=False)
    final.to_csv(PROCESSED_DIR / "phase1_final_table.csv", index=False)
    print(f"\nSaved: {PROCESSED_DIR / 'phase1_final_table.csv'}")
    print(f"Saved: {PROCESSED_DIR / 'tp53_mutation_records.csv'}")


if __name__ == "__main__":
    main()
