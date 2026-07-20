#!/usr/bin/env python3
"""
Compare reads to reference sequences and calculate mismatch percentage,
including edit rates at TC and GA motifs (C preceded by T, G followed by A).

Outputs a TSV file with Reference_ID, Mismatch_Percentage, TC/GA motif edit rates,
and other C->T/G->A statistics.

Usage in Python:
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "motifs", "08_compare_reads_to_reference_with_motifs.py"
    )
    motifs = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(motifs)
    results, motif_info = motifs.compare_reads_with_motifs('sample.sorted.bam', 'reference.fasta')

    # Or save directly to TSV
    motifs.compare_reads_with_motifs('sample.sorted.bam', 'reference.fasta', output='results.tsv')

To run the sample batch from a Python console (avoids copy-paste IndentationError):
    exec(open('08_compare_reads_to_reference_with_motifs.py').read())
    # Run this from the directory containing the script, or use the full path
"""

import pysam
import argparse
import sys


def compare_reads_with_motifs(
    bam_file,
    reference_fasta,
    min_quality=20,
    output=None,
    verbose=True,
    max_references=None,
    ignore_first_n_bases=0,
):
    """
    Compare reads to reference and calculate mismatch % plus TC/GA motif edit rates.

    Args:
        bam_file: Path to sorted BAM file (must be indexed)
        reference_fasta: Path to reference FASTA file
        min_quality: Minimum base quality score (default: 20)
        output: Optional output TSV path (default: None)
        verbose: Print progress (default: True)
        max_references: Max references to process (default: None)
        ignore_first_n_bases: Ignore first N bases of each read (default: 0)

    Returns:
        (results_dict, motif_info_dict) where results has mismatch_pct per ref,
        motif_info has TC/GA stats per ref
    """
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
        fasta = pysam.FastaFile(reference_fasta)
    except Exception as e:
        error_msg = f"Error opening files: {e}"
        if verbose:
            print(error_msg, file=sys.stderr)
        raise IOError(error_msg) from e

    if verbose:
        print(f"Processing {bam_file}...", file=sys.stderr)

    results = {}
    motif_info = {}

    references = list(fasta.references)
    if max_references is not None:
        references = references[:max_references]
        if verbose:
            print(f"Limiting to first {max_references} reference sequences", file=sys.stderr)

    for ref_name in references:
        ref_seq = fasta.fetch(ref_name).upper()
        ref_length = len(ref_seq)

        # All C and G positions
        c_positions = set(i for i, base in enumerate(ref_seq) if base == "C")
        g_positions = set(i for i, base in enumerate(ref_seq) if base == "G")

        # TC motif: C preceded by T (pos-1 == T)
        tc_positions = set(pos for pos in c_positions if pos > 0 and ref_seq[pos - 1] == "T")
        # GA motif: G followed by A (pos+1 == A)
        ga_positions = set(pos for pos in g_positions if pos + 1 < ref_length and ref_seq[pos + 1] == "A")

        c_position_edits = {pos: {"C": 0, "T": 0} for pos in c_positions}
        g_position_edits = {pos: {"G": 0, "A": 0} for pos in g_positions}

        total_bases = 0
        mismatches = 0
        read_count = 0

        try:
            for read in bam.fetch(ref_name):
                if read.is_unmapped:
                    continue
                read_count += 1
                aligned_pairs = read.get_aligned_pairs(matches_only=False, with_seq=False)
                if aligned_pairs is None:
                    continue

                for query_pos, ref_pos in aligned_pairs:
                    if query_pos is None or ref_pos is None:
                        continue
                    if ref_pos >= ref_length:
                        continue
                    if ignore_first_n_bases and query_pos < ignore_first_n_bases:
                        continue
                    if read.query_qualities and read.query_qualities[query_pos] < min_quality:
                        continue

                    query_base = read.query_sequence[query_pos].upper()
                    ref_base_char = ref_seq[ref_pos]
                    total_bases += 1
                    if query_base != ref_base_char:
                        mismatches += 1

                    if ref_pos in c_positions:
                        if query_base == "T":
                            c_position_edits[ref_pos]["T"] += 1
                        elif query_base == "C":
                            c_position_edits[ref_pos]["C"] += 1
                    if ref_pos in g_positions:
                        if query_base == "A":
                            g_position_edits[ref_pos]["A"] += 1
                        elif query_base == "G":
                            g_position_edits[ref_pos]["G"] += 1
        except Exception as e:
            if verbose:
                print(f"Warning: Error processing reads for {ref_name}: {e}", file=sys.stderr)
            continue

        if total_bases == 0:
            continue

        mismatch_pct = (mismatches / total_bases) * 100
        results[ref_name] = mismatch_pct

        # --- All C / G stats ---
        total_c_obs = total_t_obs = 0
        total_g_obs = total_a_obs = 0
        edited_c_positions = edited_g_positions = 0
        for pos, counts in c_position_edits.items():
            t = counts["C"] + counts["T"]
            if t > 0:
                total_c_obs += counts["C"]
                total_t_obs += counts["T"]
                if (counts["T"] / t) >= 0.20:
                    edited_c_positions += 1
        for pos, counts in g_position_edits.items():
            t = counts["G"] + counts["A"]
            if t > 0:
                total_g_obs += counts["G"]
                total_a_obs += counts["A"]
                if (counts["A"] / t) >= 0.20:
                    edited_g_positions += 1

        total_c_t = total_c_obs + total_t_obs
        total_g_a = total_g_obs + total_a_obs
        c_to_t_edit_rate = (total_t_obs / total_c_t * 100) if total_c_t > 0 else 0.0
        g_to_a_edit_rate = (total_a_obs / total_g_a * 100) if total_g_a > 0 else 0.0
        c_positions_edited_pct = (edited_c_positions / len(c_positions) * 100) if c_positions else 0.0
        g_positions_edited_pct = (edited_g_positions / len(g_positions) * 100) if g_positions else 0.0

        # --- TC motif: C preceded by T ---
        tc_c_obs = tc_t_obs = 0
        tc_edited_positions = 0
        for pos in tc_positions:
            counts = c_position_edits[pos]
            t = counts["C"] + counts["T"]
            if t > 0:
                tc_c_obs += counts["C"]
                tc_t_obs += counts["T"]
                if (counts["T"] / t) >= 0.20:
                    tc_edited_positions += 1
        tc_total = tc_c_obs + tc_t_obs
        tc_edit_rate = (tc_t_obs / tc_total * 100) if tc_total > 0 else 0.0
        tc_positions_edited_pct = (tc_edited_positions / len(tc_positions) * 100) if tc_positions else 0.0

        # --- GA motif: G followed by A ---
        ga_g_obs = ga_a_obs = 0
        ga_edited_positions = 0
        for pos in ga_positions:
            counts = g_position_edits[pos]
            t = counts["G"] + counts["A"]
            if t > 0:
                ga_g_obs += counts["G"]
                ga_a_obs += counts["A"]
                if (counts["A"] / t) >= 0.20:
                    ga_edited_positions += 1
        ga_total = ga_g_obs + ga_a_obs
        ga_edit_rate = (ga_a_obs / ga_total * 100) if ga_total > 0 else 0.0
        ga_positions_edited_pct = (ga_edited_positions / len(ga_positions) * 100) if ga_positions else 0.0

        # Combined TC/GA edit rate: (GA_A_Obs + TC_T_Obs) / (TC_Total + GA_Total) * 100
        combined_tc_ga_total = tc_total + ga_total
        tc_ga_edit_rate = (
            (ga_a_obs + tc_t_obs) / combined_tc_ga_total * 100
            if combined_tc_ga_total > 0
            else 0.0
        )

        motif_info[ref_name] = {
            "mismatch_percentage": mismatch_pct,
            "c_to_t_edit_rate": c_to_t_edit_rate,
            "c_positions_edited_percentage": c_positions_edited_pct,
            "edited_c_positions": edited_c_positions,
            "total_c_positions": len(c_positions),
            "total_c_observations": total_c_obs,
            "total_t_observations": total_t_obs,
            "total_c_t_observations": total_c_t,
            "g_to_a_edit_rate": g_to_a_edit_rate,
            "g_positions_edited_percentage": g_positions_edited_pct,
            "edited_g_positions": edited_g_positions,
            "total_g_positions": len(g_positions),
            "total_g_observations": total_g_obs,
            "total_a_observations": total_a_obs,
            "total_g_a_observations": total_g_a,
            "tc_edit_rate": tc_edit_rate,
            "tc_positions_edited_percentage": tc_positions_edited_pct,
            "tc_edited_positions": tc_edited_positions,
            "total_tc_positions": len(tc_positions),
            "tc_c_observations": tc_c_obs,
            "tc_t_observations": tc_t_obs,
            "tc_total_observations": tc_total,
            "ga_edit_rate": ga_edit_rate,
            "ga_positions_edited_percentage": ga_positions_edited_pct,
            "ga_edited_positions": ga_edited_positions,
            "total_ga_positions": len(ga_positions),
            "ga_g_observations": ga_g_obs,
            "ga_a_observations": ga_a_obs,
            "ga_total_observations": ga_total,
            "tc_ga_edit_rate": tc_ga_edit_rate,
            "read_count": read_count,
        }

        if verbose:
            print(
                f"  {ref_name}: {mismatch_pct:.2f}% mismatch | "
                f"C->T: {c_to_t_edit_rate:.2f}% (TC motif: {tc_edit_rate:.2f}%) | "
                f"G->A: {g_to_a_edit_rate:.2f}% (GA motif: {ga_edit_rate:.2f}%)",
                file=sys.stderr,
            )

    bam.close()
    fasta.close()

    if output:
        with open(output, "w") as f:
            f.write(
                "Reference_ID\tMismatch_Percentage\t"
                "C_to_T_Edit_Rate\tC_Positions_Edited_Pct\tEdited_C_Positions\tTotal_C_Positions\t"
                "Total_C_Observations\tTotal_T_Observations\tTotal_C+T_Observations\t"
                "TC_Edit_Rate\tTC_Positions_Edited_Pct\tTC_Edited_Positions\tTotal_TC_Positions\t"
                "TC_C_Observations\tTC_T_Observations\tTC_Total_Observations\t"
                "G_to_A_Edit_Rate\tG_Positions_Edited_Pct\tEdited_G_Positions\tTotal_G_Positions\t"
                "Total_G_Observations\tTotal_A_Observations\tTotal_G+A_Observations\t"
                "GA_Edit_Rate\tGA_Positions_Edited_Pct\tGA_Edited_Positions\tTotal_GA_Positions\t"
                "GA_G_Observations\tGA_A_Observations\tGA_Total_Observations\t"
                "TC_GA_Edit_Rate\t"
                "Read_Count\n"
            )
            for ref_name in sorted(results.keys()):
                info = motif_info[ref_name]
                f.write(
                    f"{ref_name}\t{results[ref_name]:.4f}\t"
                    f"{info['c_to_t_edit_rate']:.4f}\t{info['c_positions_edited_percentage']:.4f}\t"
                    f"{info['edited_c_positions']}\t{info['total_c_positions']}\t"
                    f"{info['total_c_observations']}\t{info['total_t_observations']}\t{info['total_c_t_observations']}\t"
                    f"{info['tc_edit_rate']:.4f}\t{info['tc_positions_edited_percentage']:.4f}\t"
                    f"{info['tc_edited_positions']}\t{info['total_tc_positions']}\t"
                    f"{info['tc_c_observations']}\t{info['tc_t_observations']}\t{info['tc_total_observations']}\t"
                    f"{info['g_to_a_edit_rate']:.4f}\t{info['g_positions_edited_percentage']:.4f}\t"
                    f"{info['edited_g_positions']}\t{info['total_g_positions']}\t"
                    f"{info['total_g_observations']}\t{info['total_a_observations']}\t{info['total_g_a_observations']}\t"
                    f"{info['ga_edit_rate']:.4f}\t{info['ga_positions_edited_percentage']:.4f}\t"
                    f"{info['ga_edited_positions']}\t{info['total_ga_positions']}\t"
                    f"{info['ga_g_observations']}\t{info['ga_a_observations']}\t{info['ga_total_observations']}\t"
                    f"{info['tc_ga_edit_rate']:.4f}\t"
                    f"{info['read_count']}\n"
                )
        if verbose:
            print(f"Results written to {output}", file=sys.stderr)

    return results, motif_info


def main():
    parser = argparse.ArgumentParser(
        description="Compare reads to reference, report mismatch % and TC/AG motif edit rates",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("bam_file", help="Path to sorted BAM file (must be indexed)")
    parser.add_argument("reference_fasta", help="Path to reference FASTA file")
    parser.add_argument("-q", "--min-quality", type=int, default=20, help="Min base quality (default: 20)")
    parser.add_argument("-o", "--output", type=str, default=None, help="Output TSV path")
    parser.add_argument("-n", "--max-references", type=int, default=None, help="Max references to process")
    parser.add_argument(
        "--ignore-first-n-bases",
        type=int,
        default=0,
        help="Ignore first N bases per read (default: 0)",
    )
    args = parser.parse_args()

    results, motif_info = compare_reads_with_motifs(
        args.bam_file,
        args.reference_fasta,
        min_quality=args.min_quality,
        output=args.output,
        verbose=True,
        max_references=args.max_references,
        ignore_first_n_bases=args.ignore_first_n_bases,
    )

    if args.output is None and results:
        print(
            "Reference_ID\tMismatch_Pct\tC_to_T_Rate\tTC_Edit_Rate\tG_to_A_Rate\tGA_Edit_Rate\t"
            "TC_GA_Edit_Rate\t"
            "Edited_TC\tTotal_TC\tEdited_GA\tTotal_GA\tRead_Count"
        )
        for ref_name in sorted(results.keys()):
            info = motif_info[ref_name]
            print(
                f"{ref_name}\t{results[ref_name]:.4f}\t"
                f"{info['c_to_t_edit_rate']:.4f}\t{info['tc_edit_rate']:.4f}\t"
                f"{info['g_to_a_edit_rate']:.4f}\t{info['ga_edit_rate']:.4f}\t"
                f"{info['tc_ga_edit_rate']:.4f}\t"
                f"{info['tc_edited_positions']}\t{info['total_tc_positions']}\t"
                f"{info['ga_edited_positions']}\t{info['total_ga_positions']}\t"
                f"{info['read_count']}"
            )


def run_all_samples():
    """Run compare_reads_with_motifs on all sample and ctrl BAMs. Ignores first 9 bases (barcode)."""
    # EDIT ME: set reference FASTAs and BAM/output paths for your project
    main_reference_fasta = 'reference/YOUR_REFERENCE.fasta'
    complement_reference_fasta = 'reference/YOUR_REFERENCE_complement.fasta'

    # sample1 main and complement
    results_sample1_main, motif_info_sample1_main = compare_reads_with_motifs(
        'output/sample1.filtered_189_190bp.sorted.bam',
        main_reference_fasta,
        output='results_sample1_main_motifs.tsv',
        max_references=None,
        ignore_first_n_bases=9,
    )
    results_sample1_complement, motif_info_sample1_complement = compare_reads_with_motifs(
        'output/sample1.filtered_189_190bp.complement.sorted.bam',
        complement_reference_fasta,
        output='results_sample1_complement_motifs.tsv',
        max_references=None,
        ignore_first_n_bases=9,
    )

    # sample2 main and complement
    results_sample2_main, motif_info_sample2_main = compare_reads_with_motifs(
        'output/sample2.filtered_189_190bp.sorted.bam',
        main_reference_fasta,
        output='results_sample2_main_motifs.tsv',
        max_references=None,
        ignore_first_n_bases=9,
    )
    results_sample2_complement, motif_info_sample2_complement = compare_reads_with_motifs(
        'output/sample2.filtered_189_190bp.complement.sorted.bam',
        complement_reference_fasta,
        output='results_sample2_complement_motifs.tsv',
        max_references=None,
        ignore_first_n_bases=9,
    )

    # sample3 main and complement
    results_sample3_main, motif_info_sample3_main = compare_reads_with_motifs(
        'output/sample3.filtered_189_190bp.sorted.bam',
        main_reference_fasta,
        output='results_sample3_main_motifs.tsv',
        max_references=None,
        ignore_first_n_bases=9,
    )
    results_sample3_complement, motif_info_sample3_complement = compare_reads_with_motifs(
        'output/sample3.filtered_189_190bp.complement.sorted.bam',
        complement_reference_fasta,
        output='results_sample3_complement_motifs.tsv',
        max_references=None,
        ignore_first_n_bases=9,
    )

    # ctrl1 main and complement
    results_ctrl1_main, motif_info_ctrl1_main = compare_reads_with_motifs(
        'output/ctrl1.filtered_189_190bp.sorted.bam',
        main_reference_fasta,
        output='results_ctrl1_main_motifs.tsv',
        max_references=None,
        ignore_first_n_bases=9,
    )
    results_ctrl1_complement, motif_info_ctrl1_complement = compare_reads_with_motifs(
        'output/ctrl1.filtered_189_190bp.complement.sorted.bam',
        complement_reference_fasta,
        output='results_ctrl1_complement_motifs.tsv',
        max_references=None,
        ignore_first_n_bases=9,
    )


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main()
    else:
        run_all_samples()
