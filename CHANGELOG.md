# Changelog

## v0.1.0: Initial public release
- Lightweight VCF QC runner → **CSV** metrics (Ti/Tv, SNP/INDEL counts, header checks, warnings).
- **NEW: MultiQC reporting** via `scripts/make_vcf_qc_multiqc.sh` + `scripts/report_vcf_qc.sh`:
  - Custom table section (PASS rate, Ti/Tv, het:hom, med DP/QUAL, multi-allelic %, etc.)
  - Optional **General Statistics** TSV (set `MAKE_GS=1`).
- Demo VCF in `examples/` and updated README with run instructions.
- `.gitignore` refined to keep real **inputs/** and **outputs/** local; demo artifacts under **results/**.
- **CITATION.cff** added (for DOI on release).
