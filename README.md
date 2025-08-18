# VCF QC Pipeline

A minimal **command-line tool** to inspect, summarise, and troubleshoot VCF files.  
Designed for both **teaching** and **real-world pipeline checks**.

**Key features**
- Counts SNPs, INDELs, transitions, and transversions  
- Computes transition/transversion (Ti/Tv) ratio  
- Checks header integrity (contigs, INFO, FORMAT fields)  
- Flags malformed rows and multiallelic sites  
- Outputs results as a clean **CSV summary**  
- Ships with a **privacy-safe synthetic VCF** for demo use  

---

## Why this matters (for labs and teaching)

- **Pipeline gatekeeper** - catch header/reference mismatches, malformed rows, or odd Ti/Tv ratios *before* running heavy annotation and aggregation.  
- **Regression checks** - compare CSV outputs across pipeline releases; flag shifts in SNP/INDEL counts or Ti/Tv.  
- **Teaching tool** - helps students learn VCF anatomy, variant composition, and QC signals.  
- **Privacy safe** - bundled with `examples/demo_sample.vcf.gz`, a synthetic file with no patient data.  

---

## Quick Start

```bash
# Run QC on your VCF (.vcf or .vcf.gz)
python scripts/run_vcf_qc.py --vcf <your_sample>.vcf.gz --out results/qc_summary.csv
```

### Try it now (with the bundled example)

```bash
python scripts/run_vcf_qc.py --vcf examples/demo_sample.vcf.gz --out results/qc_summary.csv
```

---

## CSV column definitions

- **samples** - number of sample columns in the VCF (after `FORMAT` in the `#CHROM` header). Site-only VCFs will report 0.  
- **variants** - total number of variant rows parsed. Malformed rows are skipped and counted in `warnings`.  
- **snps** - rows where REF length = 1 and all ALT alleles are length = 1.  
- **indels** - rows with at least one ALT allele of different length than REF.  
- **ti** - number of transition SNPs (A<->G or C<->T).  
- **tv** - number of transversion SNPs (all other single-base substitutions).  
- **titv** - ratio = `ti / tv`. If `tv=0`, this pipeline reports `0.0` (configurable).  
- **multiallelic_sites** - rows with multiple ALT alleles (comma-separated, e.g., `A,T`).  
- **contigs** - number of `##contig=...` lines declared in the header.  
- **info_fields** - number of `##INFO=...` field definitions.  
- **format_fields** - number of `##FORMAT=...` field definitions.  
- **warnings** - number of malformed rows skipped (e.g., fewer than 8 tab-separated columns).  

---

## Quick interpretation tips

- **Ti/Tv ratio**:  
  - Germline WGS - ~2.0 to 2.2  
  - Exome/targeted - ~2.6 to 3.3  
  - Somatic/tumour callsets - often lower or more variable  
  A Ti/Tv much lower than expected suggests technical artefacts or overly permissive filtering.  

- **warnings > 0** - indicates malformed rows. Use bcftools to validate.  
- **multiallelic_sites** - normal in many datasets but some tools require splitting.  

---

## Example output (demo file)

Running `demo_sample.vcf.gz` produces:

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

---

## Real-world applications

This tool provides **fast health signals** for VCF files. Use cases include:

- **Gatekeeping before annotation**  
  Stop early if headers do not match reference, Ti/Tv is unrealistic, or warnings > 0.  

- **Regression testing**  
  Keep a small golden dataset. After pipeline updates, compare new QC metrics to expected values. Flag large deviations (for example, ±20 percent in variant counts).  

- **Outlier detection across cohorts**  
  Aggregate CSVs across samples and plot distributions. Outliers in `variants`, `titv`, or SNP:INDEL ratios may indicate contamination or sample swaps.  

- **Reference build hygiene**  
  If contigs=0 or contigs do not match, fix headers or rename chromosomes before downstream tools.  
  - Rename chromosomes:  
    ```bash
    bcftools annotate --rename-chrs mapping.txt in.vcf.gz -Oz -o out.vcf.gz
    ```  
  - Normalise REF/ALT against FASTA:  
    ```bash
    bcftools norm -f ref.fa -Oz -o out.vcf.gz
    ```  

- **Pre-normalisation decisions**  
  Many multiallelic sites? Split to biallelic before annotation:  
  ```bash
  bcftools norm -m -any in.vcf.gz -Oz -o out.vcf.gz
  ```  

---

## Limitations (by design)

This QC tool is intentionally lightweight. It does **not** compute:  
- Depth or allele-balance histograms  
- Hardy-Weinberg equilibrium  
- Relatedness or ancestry metrics  

For deeper QC use alongside:  
- `bcftools stats` (+ `plot-vcfstats`)  
- `peddy` for ancestry, sex, relatedness  
- Caller-specific metrics (for example, GATK/Picard `CollectVariantCallingMetrics`)  

> Rule of thumb: **Consistency beats absolutes**. Baseline your QC metrics on a known-good run, then flag unexpected drift.  

---

## Citation

If you use or adapt this tool in teaching or pipelines, please cite:  
**Samson Olofinsae. VCF QC Pipeline. GitHub 2025.**  
<https://github.com/samson-olofinsae/vcf_qc_pipeline>

---

## Contact

Questions or suggestions?  
- Open an issue on GitHub  
- Connect via [GitHub profile](https://github.com/samson-olofinsae)  
