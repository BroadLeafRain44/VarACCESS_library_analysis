#!/usr/bin/env python3
"""
SAM stream helpers for filtering BAM reads by sequence length.

Subcommands:
  count   — read-length histogram from samtools view (no header required)
  filter  — keep only reads whose sequence length is in --lengths (header kept)
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter


def count_lengths(infile):
    length_counts = Counter()
    for line in infile:
        fields = line.strip().split("\t")
        if len(fields) >= 10:
            length_counts[len(fields[9])] += 1

    print("Read length distribution:")
    print("Length (bp)\tCount")
    for length in sorted(length_counts.keys()):
        print(f"{length}\t{length_counts[length]}")

    total_reads = sum(length_counts.values())
    print(f"\nTotal reads: {total_reads}")


def filter_lengths(infile, outfile, lengths):
    target = set(lengths)
    reads_kept = 0
    reads_removed = 0

    for line in infile:
        if line.startswith("@"):
            outfile.write(line)
            continue

        fields = line.strip().split("\t")
        if len(fields) >= 10:
            read_length = len(fields[9])
            if read_length in target:
                outfile.write(line)
                reads_kept += 1
            else:
                reads_removed += 1

    lengths_label = " and ".join(f"{L} bp" for L in sorted(target))
    sys.stderr.write(f"Reads kept ({lengths_label}): {reads_kept}\n")
    sys.stderr.write(f"Reads removed: {reads_removed}\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("count", help="Print read-length counts from SAM alignments on stdin")

    p_filter = sub.add_parser("filter", help="Keep only reads with given sequence lengths")
    p_filter.add_argument(
        "--lengths",
        nargs="+",
        type=int,
        required=True,
        help="Allowed read lengths in bp (e.g. 189 190)",
    )

    args = parser.parse_args(argv)

    if args.command == "count":
        count_lengths(sys.stdin)
    elif args.command == "filter":
        filter_lengths(sys.stdin, sys.stdout, args.lengths)


if __name__ == "__main__":
    main()
