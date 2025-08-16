# VCF QC Pipeline

An educational command-line tool to inspect, summarise, and troubleshoot VCF files.
Computes basic variant statistics (SNP/INDEL counts, Ti/Tv), inspects headers, and generates a CSV summary.


## Why this matters (for labs & teaching)

- **Pipeline gate:** catch header/reference mismatches, malformed lines, or Ti/Tv anomalies *before* spending compute on annotation/aggregation.
- **Regression checks:** diff the CSV across pipeline releases; alert on drift in variant counts or Ti/Tv.
- **Teaching value:** a minimal, dependency-light script that helps new analysts understand VCF anatomy and QC signals (SNP/INDEL balance, Ti/Tv, multiallelic handling).
- **Privacy-safe demo:** ships with a tiny **synthetic VCF** (`examples/demo_sample.vcf.gz`); no PHI or real patient data.

---

## Quick Start

```bash
# run on your VCF (.vcf or .vcf.gz)
python scripts/run_vcf_qc.py --vcf <your_sample>.vcf.gz --out results/qc_summary.csv



### Try it now (with the bundled example)

```bash
python scripts/run_vcf_qc.py --vcf examples/demo_sample.vcf.gz --out results/qc_summary.csv



## CSV column definitions (what each field means)

- **samples** — number of *sample* columns present in the VCF (columns after `FORMAT` in the `#CHROM` header line). If your VCF has only site-level data and no samples, this will be 0.
- **variants** — count of variant records parsed (non-header lines). Malformed rows are skipped and counted under `warnings`.
- **snps** — single-nucleotide polymorphisms where `REF` length = 1 and **all** `ALT` alleles are length = 1.
- **indels** — insertions/deletions where any `ALT` allele length ≠ `REF` length.
- **ti** — number of **transitions** among SNPs (A↔G or C↔T changes).
- **tv** — number of **transversions** among SNPs (all other single-base changes: A↔C, A↔T, C↔G, G↔T).
- **titv** — transition/transversion ratio = `ti / tv`.  
  - In this training tool, if `tv = 0`, the value is `0.0` (you can change it to `Inf` if you prefer; see comment in `run_vcf_qc.py`).
- **multiallelic_sites** — sites where the `ALT` field contains multiple alleles separated by commas (e.g., `A,T`).
- **contigs** — number of `##contig=...` definitions seen in the header (declared reference contigs).
- **info_fields** — number of `##INFO=...` definitions (site annotation fields defined in the header).
- **format_fields** — number of `##FORMAT=...` definitions (per-sample genotype/format fields defined in the header).
- **warnings** — lines skipped because they were malformed (e.g., fewer than 8 tab-separated VCF columns).

### Quick interpretation tips
- A **higher Ti/Tv** than 1 is typical for real SNP data; very low Ti/Tv may indicate technical artefacts or parsing issues.
- A **non-zero `warnings`** count usually means there are rows with missing/short columns; consider validating the VCF.
- Many **multiallelic sites** can be normal in real datasets; just note they’re counted once per line.


### Example (from the included demo)
`demo_sample.vcf.gz` → `samples=1`, `variants=3`, `snps=2`, `indels=1`, `ti=2`, `tv=0`, `titv=0.0`, `multiallelic_sites=1`, `contigs=2`, `info_fields=1`, `format_fields=1`, `warnings=0`.


## Real-world application: how to use these QC checks

This tool gives you fast, lightweight signals about a VCF’s health. Here’s how each metric in the CSV helps you make decisions in practice:

- **samples**  
  Confirms how many sample columns exist.  
  *Use it to*: catch missing sample columns (e.g., site-only VCFs) or unexpected multi-sample inputs.

- **variants / snps / indels**  
  Basic site counts and composition. Most germline callsets are SNP-heavy; a sudden drop in total variants or a big shift in SNP:INDEL balance often means a pipeline/config change.  
  *Use it to*: compare runs/releases; spot regressions (e.g., new filters removing too many INDELs).

- **ti / tv / titv**  
  Transition (A↔G, C↔T) vs transversion counts and their ratio.  
  *Typical expectations (human, germline)*:  
  - Whole-genome: Ti/Tv ~ **2.0–2.2**  
  - Exome/targeted: Ti/Tv ~ **2.6–3.3** (capture & filters matter)  
  *Somatic/tumour callsets* often show lower or more variable Ti/Tv.  
  *Use it to*: flag potential artefacts or overly permissive/strict filters.  
  *Note*: With very few SNPs (or when `tv=0`), Ti/Tv is not very informative.

- **multiallelic_sites**  
  Lines with multiple ALT alleles (e.g., `A,T`). Normal in many callsets, but some tools prefer split biallelics.  
  *Use it to*: decide whether to normalise/split before downstream steps.  
  *Tip*: `bcftools norm -m -any in.vcf.gz -Oz -o out.vcf.gz`

- **contigs / info_fields / format_fields**  
  Header integrity checks.  
  - **contigs** should match your reference build; a mismatch (e.g., missing `chr` prefixes) can break tools.  
  - **info_fields** / **format_fields** should define fields you actually use in the body.  
  *Use it to*: detect reference build mismatches and incomplete headers early.  
  *Tips*:  
  - Rename chromosomes: `bcftools annotate --rename-chrs mapping.txt in.vcf.gz -Oz -o out.vcf.gz`  
  - Normalise REF/ALT against the FASTA: `bcftools norm -f ref.fa -Oz -o out.vcf.gz`

- **warnings**  
  Count of malformed lines skipped (e.g., fewer than 8 tab-separated columns).  
  *Use it to*: decide if you need to validate/repair the file before analysis.  
  *Tips*:  
  - Quick validation: `bcftools view -H in.vcf.gz >/dev/null` (non-zero exit on serious parse errors)  
  - Normalise/sanitise: `bcftools norm` (see above)

### Practical workflows this enables

- **Gatekeeping before expensive steps**  
  Fail fast if headers don’t match the reference, Ti/Tv is wildly off, or warnings > 0. Avoid wasting time on annotation/aggregation that will fail later.

- **Regression checks in pipelines**  
  Keep a small “golden” sample set. After code/config changes, compare CSV metrics to previous runs; alert on large deviations (e.g., ±10–20% in variant counts or Ti/Tv shifts).

- **Outlier detection across cohorts**  
  Aggregate per-sample CSVs and plot distributions. Outliers in `variants`, `snps/indels`, or `titv` often identify sample swaps, contamination, or filter drift.

- **Reference build hygiene**  
  If `contigs=0` or clearly mismatched, fix headers or re-lift before downstream tools (liftover/rename as appropriate).

- **Pre-normalisation decisions**  
  Many `multiallelic_sites`? Split/normalise now if your next tool expects biallelic inputs.

### Limitations (by design)

This QC is intentionally lightweight: it doesn’t compute depth/allele-balance histograms, HWE, or sample relatedness. For deeper diagnostics, pair it with tools like:
- `bcftools stats` (+ `plot-vcfstats`) for rich callset summaries
- `peddy` or cohort QA tools for ancestry/sex/relatedness (germline projects)
- caller-specific metrics (e.g., GATK/Picard `CollectVariantCallingMetrics`)

> Rule of thumb: **consistency beats absolutes**. Baseline these metrics on a known-good run; flag anything that drifts unexpectedly given your organism, assay (WGS/WES/panel), and filters.
