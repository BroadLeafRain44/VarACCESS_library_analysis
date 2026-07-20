#!/usr/bin/env python3
"""Write the base complement of sequences in a FASTQ file (headers/qualities unchanged)."""

from __future__ import annotations

import argparse
import os
import sys

COMPLEMENT_MAP = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "a": "t",
    "t": "a",
    "c": "g",
    "g": "c",
    "N": "N",
    "n": "n",
}


def fastq_to_complement(infile, outfile, max_reads=None):
    infile = os.path.expanduser(infile)
    outfile = os.path.expanduser(outfile)

    out_dir = os.path.dirname(outfile)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(infile, "r", buffering=1024) as f_in, open(outfile, "w", buffering=1024) as f_out:
        line_count = 0
        read_count = 0

        for line in f_in:
            line_count += 1

            # Line 2 of every 4-line block is the DNA sequence
            if line_count % 4 == 2:
                read_count += 1
                seq_line = line.rstrip("\n\r")
                for base in seq_line:
                    f_out.write(COMPLEMENT_MAP.get(base, base))
                f_out.write("\n")

                if read_count % 100000 == 0:
                    print(f"Processed {read_count:,} reads...")
            else:
                f_out.write(line)

            if max_reads is not None and read_count >= max_reads and line_count % 4 == 0:
                break

    print(f"Done! Converted {read_count:,} reads to complement. Saved to {outfile}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_fastq", help="Input FASTQ path")
    parser.add_argument("output_fastq", help="Output complement FASTQ path")
    parser.add_argument(
        "--max-reads",
        type=int,
        default=None,
        help="Optional limit on number of reads to convert",
    )
    args = parser.parse_args(argv)
    fastq_to_complement(args.input_fastq, args.output_fastq, max_reads=args.max_reads)


if __name__ == "__main__":
    main()
