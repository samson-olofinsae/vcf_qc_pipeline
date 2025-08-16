#!/usr/bin/env python3
import argparse

def main():
    ap = argparse.ArgumentParser(description='Placeholder for AF/DP plotting.')
    ap.add_argument('--csv', required=False, help='QC summary CSV (future use)')
    args = ap.parse_args()
    print('Plotting module is a placeholder. Extend with cyvcf2/pandas/matplotlib as needed.')

if __name__ == '__main__':
    main()
