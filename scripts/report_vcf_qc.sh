#!/usr/bin/env bash
# Usage: scripts/report_vcf_qc.sh [VCF_DIR] [OUT_TSV] [OUT_HTML]
set -euo pipefail
VCF_DIR="${1:-examples/}"
OUT_TSV="${2:-results/multiqc_cc/vcf_qc_summary_mqc.tsv}"
OUT_HTML="${3:-results/multiqc/vcf_qc.html}"
MAKE_GS=1 scripts/make_vcf_qc_multiqc.sh "$VCF_DIR" "$OUT_TSV" "$OUT_HTML"


