# VCF QC Pipeline

A minimal, **teaching-friendly** command-line tool to inspect, summarise, and troubleshoot VCF files - now with **one-click MultiQC HTML reports**.

**Key features**
- Counts SNPs, INDELs, transitions, transversions; computes **Ti/Tv**.
- Checks header integrity (contigs, INFO, FORMAT fields).
- Flags malformed rows & multiallelic sites.
- Outputs a clean **CSV summary** for regression checks.
- **NEW:** Standalone **MultiQC** builder:
  - Per-sample **table section** (custom content).
  - Optional **General Statistics** panel (headline KPIs).
- Ships with a **privacy-safe synthetic VCF** for demo use.

---

## Why this matters (for labs and teaching)

- **Pipeline gatekeeper** - catch header/reference mismatches, malformed rows, or odd Ti/Tv *before* heavy annotation/aggregation.
- **Regression checks** - compare CSVs across releases; flag shifts in SNP/INDEL counts or Ti/Tv.
- **Teaching tool** - learn VCF anatomy, variant composition, and QC signals.
- **Privacy safe** - bundled with `examples/demo_sample.vcf.gz`.

---

## Repository layout

```
.
├── README.md
├── docs/
├── examples/
│   └── demo_sample.vcf.gz      # demo VCF (bgzipped)
├── inputs/
├── results/
│   ├── qc_logs/                # run logs (from the QC runner)
│   ├── qc_summary.csv          # CSV summary (from the QC runner)
│   ├── multiqc_cc/             # MultiQC custom-content tables (generated)
│   └── multiqc/                # MultiQC HTML + assets (generated)
├── scripts/
│   ├── run_vcf_qc.py           # QC runner → writes results/qc_summary.csv
│   ├── utils.py
│   ├── plot_af_dp.py           # optional plotting helper
│   ├── make_vcf_qc_multiqc.sh  # NEW: VCF → MultiQC TSV (+ HTML)
│   └── report_vcf_qc.sh        # NEW: convenience wrapper for the above
└── step3_seed.sh
```

---

## Quick Start

## Setup

### Clone the repository

```bash
git clone https://github.com/samson-olofinsae/vcf_qc_pipeline.git
cd vcf_qc_pipeline
```


### 1) Run QC on a VCF (CSV summary)
```bash
# Your own data
python scripts/run_vcf_qc.py --vcf <your_sample>.vcf.gz --out results/qc_summary.csv

# Try the bundled demo
python scripts/run_vcf_qc.py --vcf examples/demo_sample.vcf.gz --out results/qc_summary.csv
```

### 2) Build the MultiQC report (table + General Statistics)

**Option A — wrapper (auto-detects `examples/` or `results/vcf/`):**
```bash
bash scripts/report_vcf_qc.sh
# or specify the folder that contains your .vcf.gz files:
bash scripts/report_vcf_qc.sh examples
# e.g., for real outputs placed under results/vcf/:
bash scripts/report_vcf_qc.sh results/vcf
```

**Option B — direct builder (explicit paths; writes both TSV + HTML):**
```bash
# MAKE_GS=1 also emits a numeric-only “General Statistics” TSV
MAKE_GS=1 scripts/make_vcf_qc_multiqc.sh <VCF_DIR> results/multiqc_cc/vcf_qc_summary_mqc.tsv results/multiqc/vcf_qc.html
```

**Open the HTML**
```bash
# Linux:
xdg-open results/multiqc/vcf_qc.html

# WSL:
wslview results/multiqc/vcf_qc.html
```

Generated files:
- `results/multiqc_cc/vcf_qc_summary_mqc.tsv` - MultiQC **table** (custom content).
- `results/multiqc_cc/vcf_qc_summary_gs_mqc.tsv` - **General Statistics** KPIs (numeric-only; written when `MAKE_GS=1`).
- `results/multiqc/vcf_qc.html` — HTML report (requires `multiqc` installed).

To overwrite an existing HTML:
```bash
multiqc results -o results/multiqc -n vcf_qc.html -f
```

---

## CSV column definitions (from `run_vcf_qc.py`)

- **samples** - number of sample columns (0 for site-only VCFs).
- **variants** - total variant rows parsed (malformed rows are skipped and counted in `warnings`).
- **snps** - REF length = 1 and all ALT alleles length = 1.
- **indels** - at least one ALT length differs from REF.
- **ti** - transitions (A↔G, C↔T).
- **tv** - transversions (other single-base substitutions).
- **titv** - `ti / tv` (this runner reports `0.0` if `tv=0`).
- **multiallelic_sites** - rows with multiple ALT alleles (e.g., `A,T`).
- **contigs / info_fields / format_fields** - header counts.
- **warnings** - malformed rows skipped.

---

## MultiQC table (what it shows, per sample)

**Section:** *VCF QC - variant summary*  
Columns produced by `make_vcf_qc_multiqc.sh` (via `bcftools`):

- **total_vars** - total variant records.
- **snvs / indels** - counts by class.
- **pass_rate** - PASS / total (`bcftools view -f PASS`).
- **titv** - transitions/transversions (`bcftools stats`).
- **het_hom** - heterozygous : homozygous genotype ratio (`bcftools stats`); blank if hom=0 or no GTs.
- **med_DP**, **med_QUAL** - medians from `bcftools query` (may be blank if tags absent).
- **multi_allelic_pct** - % of sites that are **not** biallelic (computed with `-m2 -M2`).
- **filtered_pct** - % non-PASS = `100 × (1 − pass_rate)`.

> Sample names are inferred from filenames; the builder strips common suffixes like `_final_variants`, `_variants`, `.filtered`, `.hardfiltered` for readability.

---

## Quick interpretation tips

- **Ti/Tv ratio**:  
  Germline WGS ≈ 2.0–2.2 • Exome/targeted ≈ 2.6–3.3 • Somatic often lower/more variable.  
  Much lower than expected can indicate technical artefacts or permissive filtering.
- **warnings > 0** (CSV) - malformed rows → validate with `bcftools`.
- **multiallelic_sites** - normal; some tools require splitting beforehand.

---

## Example output (demo file)

Running `examples/demo_sample.vcf.gz` yields (CSV excerpt):
```
samples=1
variants=3
snps=2
indels=1
ti=2
tv=0
titv=0.0
multiallelic_sites=1
contigs=2
info_fields=1
format_fields=1
warnings=0
```
The MultiQC report summarises the same VCF with:
- `total_vars=3`, `snvs=2`, `indels=1`
- `pass_rate=1.0`, `filtered_pct=0.0`
- `titv=2.0`
- `het_hom` blank (no hom calls / no GTs)
- `med_DP≈20`, `med_QUAL='.'` if QUAL absent
- `multi_allelic_pct≈33.3%` (1/3 records multiallelic)

---

## Requirements

- **bcftools** ≥ 1.10 (`view`, `query`, `stats`)
- **Python** 3.8+ (for the runner)
- *(Optional)* **MultiQC** for HTML: `pip install multiqc` (or `mamba/conda install -c bioconda multiqc`)

---

## Troubleshooting

- **“bcftools not found”** → install and ensure it’s on `PATH`.
- **“No .vcf.gz files found”** → pass the correct folder (`examples/` for demo, or `results/vcf/` for real outputs).
- **Blank `het_hom`** → no homozygous genotypes or no GT fields.
- **`.` in `med_QUAL`** → QUAL absent in VCF (not an error).
- **HTML won’t open in image viewers** → use a browser (`xdg-open`, `wslview`).

---

## License
MIT - © 2025 Samson Olofinsae.

## Contact
- Issues: [Report a bug or request a feature](https://github.com/samson-olofinsae/vcf_qc_pipeline/issues)
- Maintainer: https://github.com/samson-olofinsae
- Email: olofinsae.samson@gmail.com

