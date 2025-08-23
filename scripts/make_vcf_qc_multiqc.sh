#!/usr/bin/env bash
# Build MultiQC custom-content for VCF QC across a folder of .vcf.gz files (bcftools required).
# Usage:
#   scripts/make_vcf_qc_multiqc.sh [VCF_DIR] [OUT_TSV] [OUT_HTML]
# Defaults:
#   VCF_DIR=results/vcf
#   OUT_TSV=results/multiqc_cc/vcf_qc_summary_mqc.tsv
#   OUT_HTML=results/multiqc/vcf_qc.html

set -euo pipefail

VCF_DIR="${1:-results/vcf}"
OUT_TSV="${2:-results/multiqc_cc/vcf_qc_summary_mqc.tsv}"
OUT_HTML="${3:-results/multiqc/vcf_qc.html}"
OUT_GS="${OUT_TSV/_mqc.tsv/_gs_mqc.tsv}"  # numeric-only generalstats (optional)

command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not found in PATH." >&2; exit 127; }

mkdir -p "$(dirname "$OUT_TSV")" "$(dirname "$OUT_HTML")"

# Table section
cat > "$OUT_TSV" <<'EOF'
# id: vcf_qc_summary
# section_name: VCF QC — variant summary
# description: Key VCF metrics per sample (bcftools-derived; PASS only for pass_rate; Ti/Tv and Het/Hom from bcftools stats)
# plot_type: table
# file_format: tsv
Sample	total_vars	snvs	indels	pass_rate	titv	het_hom	med_DP	med_QUAL	multi_allelic_pct	filtered_pct
EOF

# Optional: numeric “General Statistics” panel (enable by setting MAKE_GS=1 when calling)
if [[ "${MAKE_GS:-0}" == "1" ]]; then
  cat > "$OUT_GS" <<'EOF'
# id: vcf_qc_summary_gs
# section_name: VCF QC — summary KPIs
# plot_type: generalstats
# file_format: tsv
Sample	total_vars	snvs	indels	pass_rate	titv	med_DP	med_QUAL	multi_allelic_pct	filtered_pct
EOF
fi

shopt -s nullglob
found_any=0

# helper: median from stdin
med() { awk '{a[NR]=$1} END{if(NR){if(NR%2)print a[int((NR+1)/2)]; else printf "%.2f", (a[NR/2]+a[NR/2+1])/2}}'; }

for VCF in "$VCF_DIR"/*.vcf.gz; do
  found_any=1
  base="$(basename "$VCF")"
  name="${base%.vcf.gz}"
  # strip common suffixes to get a clean sample id; tweak if your filenames differ
  sample="${name%_final_variants}"
  sample="${sample%_variants}"
  sample="${sample%.filtered}"
  sample="${sample%.hardfiltered}"

  total=$(bcftools view -H "$VCF" | wc -l || echo 0)
  snvs=$(bcftools view -v snps   -H "$VCF" | wc -l || echo 0)
  indels=$(bcftools view -v indels -H "$VCF" | wc -l || echo 0)

  pass=$(bcftools view -f PASS -H "$VCF" | wc -l || echo 0)
  pass_rate="0.0000"; [[ "$total" -gt 0 ]] && pass_rate=$(awk -v p="$pass" -v t="$total" 'BEGIN{printf "%.4f", p/t}')

  stats="$(bcftools stats "$VCF")"
  titv=$(awk -F'\t' '/^TSTV\t/ {r=$5} END{if(r=="")r=0; print r+0}' <<<"$stats")
  het=$(awk -F'\t' '/^PSC\t/ {het+=$8} END{print het+0}' <<<"$stats")
  hom=$(awk -F'\t' '/^PSC\t/ {hom+=$7} END{print hom+0}' <<<"$stats")
  het_hom=""; [[ "${hom:-0}" -gt 0 ]] && het_hom=$(awk -v h="$het" -v o="$hom" 'BEGIN{printf "%.3f", (o>0? h/o : 0)}')

  med_DP=$(bcftools query -f '%DP\n'   "$VCF" 2>/dev/null | sort -n | med || true)
  med_QUAL=$(bcftools query -f '%QUAL\n' "$VCF" 2>/dev/null | sort -n | med || true)

  total_bi=$(bcftools view -m2 -M2 -H "$VCF" | wc -l || echo 0)
  multi_pct="0.00"; [[ "$total" -gt 0 ]] && multi_pct=$(awk -v t="$total" -v b="$total_bi" 'BEGIN{printf "%.2f", ((t-b)/t)*100}')

  filtered_pct=$(awk -v p="$pass_rate" 'BEGIN{printf "%.2f", (1-p)*100}')

  printf "%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$sample" "$total" "$snvs" "$indels" "$pass_rate" "$titv" "${het_hom:-}" \
    "${med_DP:-}" "${med_QUAL:-}" "$multi_pct" "$filtered_pct" >> "$OUT_TSV"

  if [[ "${MAKE_GS:-0}" == "1" ]]; then
    printf "%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$sample" "$total" "$snvs" "$indels" "$pass_rate" "$titv" "${med_DP:-}" "${med_QUAL:-}" "$multi_pct" "$filtered_pct" >> "$OUT_GS"
  fi
done

[[ "$found_any" -eq 0 ]] && echo "WARNING: No .vcf.gz files found in $VCF_DIR" >&2

# Optional: build HTML
if command -v multiqc >/dev/null 2>&1; then
  multiqc results -o "$(dirname "$OUT_HTML")" -n "$(basename "$OUT_HTML")" >/dev/null 2>&1 || true
  echo "MultiQC report: $OUT_HTML"
else
  echo "Tip: pip install multiqc  # to render HTML"
fi
