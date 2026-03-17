#!/usr/bin/env python3

import argparse
import pandas as pd
import os

def read_tsv(file_path, column_name):
    df = pd.read_csv(
        file_path,
        sep='\t',
        usecols=['clade_name', column_name],
        dtype={column_name: float},
        comment='#'
    )
    df.set_index('clade_name', inplace=True)
    return df

def merge_data(files, column_name):
    merged_df = None

    if column_name == 'ANI':
        column_name = 'ANI (if strain-level)'
    elif column_name == 'Coverage':
        column_name = 'Coverage (if strain-level)'

    for file in files:
        # Use filename as sample name to avoid duplicates
        sample_name = os.path.basename(file).replace('.tsv', '')
        df = read_tsv(file, column_name)
        df.rename(columns={column_name: sample_name}, inplace=True)
        if merged_df is None:
            merged_df = df
        else:
            merged_df = merged_df.join(df, how='outer')

    merged_df.dropna(how="all", inplace=True)
    merged_df.fillna(0, inplace=True)
    return merged_df

def main():
    parser = argparse.ArgumentParser(description="Merge TSV files and filter by sequence_abundance")
    parser.add_argument("files", nargs='+', help="Input TSV files")
    parser.add_argument("-c", "--column", default="sequence_abundance",
                        help="Column to merge (default: sequence_abundance)")
    parser.add_argument("-o", "--output", required=True, help="Output file path")
    parser.add_argument("--min_abundance", type=float, default=0,
                        help="Minimum sequence_abundance to keep (default: 0)")

    args = parser.parse_args()

    merged_df = merge_data(args.files, args.column)

    # Filter by minimum sequence_abundance if requested
    if args.min_abundance > 0:
        merged_df = merged_df[merged_df >= args.min_abundance].dropna(how='all')

    merged_df.to_csv(args.output, sep='\t')
    print(f"Merged and filtered data written to {args.output}")

if __name__ == "__main__":
    main()
