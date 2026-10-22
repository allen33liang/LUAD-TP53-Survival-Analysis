# Raw data (not tracked in git)

The scripts expect the cBioPortal study folder here:

```
data/raw/luad_tcga_pan_can_atlas_2018/
├── data_mutations.txt
└── data_clinical_patient.txt
```

## How to download

1. Go to https://www.cbioportal.org/datasets
2. Search for **Lung Adenocarcinoma (TCGA, PanCancer Atlas)**. The study ID is `luad_tcga_pan_can_atlas_2018`.
3. Click the download icon to get the study archive.
4. Extract it, then move the `luad_tcga_pan_can_atlas_2018` folder into `data/raw/`.

Only `data_mutations.txt` and `data_clinical_patient.txt` are used. Don't open `data_mutations.txt` in Excel, because it can corrupt large tab-separated files. Use a text editor or pandas.

The study's `case_lists/cases_sequenced.txt` lists 566 samples with mutation data, which equals all 566 patients. This is why patients without a TP53 record can be classed as wildtype.
