#!/usr/bin/env python3

# USAGE:
# run on your VCF (.vcf or .vcf.gz)
#python scripts/run_vcf_qc.py --vcf inputs/<your_sample>.vcf.gz --out results/qc_summary.csv



import argparse, csv, os, sys
from datetime import datetime
from utils import open_maybe_gzip, variant_iter, is_snp, is_indel, transitions_and_transversions

def log(msg, log_path):
    stamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    line = f'[{stamp}] {msg}'
    print(line)
    with open(log_path, 'a') as fh:
        fh.write(line + '\n')

def summarise_vcf(vcf_path, log_path):
    stats = {'file': os.path.basename(vcf_path),'samples':0,'variants':0,'snps':0,'indels':0,
             'ti':0,'tv':0,'titv':0.0,'multiallelic_sites':0,'contigs':0,'info_fields':0,'format_fields':0,'warnings':0}

    # parse headers
    with open_maybe_gzip(vcf_path) as fh:
        for raw in fh:
            line = raw.strip()
            if line.startswith('##contig'):
                stats['contigs'] += 1
            if line.startswith('##INFO='):
                stats['info_fields'] += 1
            if line.startswith('##FORMAT='):
                stats['format_fields'] += 1
            if line.startswith('#CHROM'):
                fields = line.split('\t')
                # samples start at column 10 (index 9); if none present, set 0
                stats['samples'] = max(0, len(fields) - 9) if len(fields) > 9 else 0
                break

    # variant lines
    for parts in variant_iter(vcf_path):
        if parts is None:
            stats['warnings'] += 1
            log('Warning: malformed variant line skipped', log_path)
            continue
        ref, alt = parts[3], parts[4]
        stats['variants'] += 1
        if ',' in alt:
            stats['multiallelic_sites'] += 1
        if is_snp(ref, alt):
            stats['snps'] += 1
            ti, tv = transitions_and_transversions(ref, alt)
            stats['ti'] += ti; stats['tv'] += tv
        elif is_indel(ref, alt):
            stats['indels'] += 1

    stats['titv'] = round(stats['ti']/stats['tv'], 3) if stats['tv'] > 0 else 0.0
    return stats

def main():
    ap = argparse.ArgumentParser(description='Educational VCF QC summariser (pure-Python, gzip-friendly).')
    ap.add_argument('--vcf', required=True, help='Path to VCF (.vcf or .vcf.gz)')
    ap.add_argument('--out', required=True, help='Output CSV path (will append if exists)')
    ap.add_argument('--logdir', default='results/qc_logs', help='Directory for log files')
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out) or '.', exist_ok=True)
    os.makedirs(args.logdir, exist_ok=True)

    log_path = os.path.join(args.logdir, f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    log(f"Starting VCF QC on {args.vcf}", log_path)

    if not os.path.exists(args.vcf):
        log(f"ERROR: VCF not found: {args.vcf}", log_path)
        sys.exit(1)

    try:
        stats = summarise_vcf(args.vcf, log_path)
    except Exception as e:
        log(f"ERROR during parsing: {e}", log_path)
        sys.exit(2)

    header = ['file','samples','variants','snps','indels','ti','tv','titv','multiallelic_sites','contigs','info_fields','format_fields','warnings']
    write_header = not os.path.exists(args.out) or os.path.getsize(args.out) == 0
    with open(args.out, 'a', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=header)
        if write_header:
            w.writeheader()
        w.writerow(stats)

    log(f"Completed. Variants={stats['variants']} SNPs={stats['snps']} INDELs={stats['indels']} Ti/Tv={stats['titv']}", log_path)
    log(f"Summary appended to {args.out}", log_path)

if __name__ == '__main__':
    main()
