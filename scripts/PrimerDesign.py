#!/usr/bin/env python3
import argparse
import datetime
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import primer3

class PrimerDesignWorkflow:

    def __init__(self, representative_fasta, output_dir="results", amplicon_size=400, n_pools=2, min_base_freq=0.2):
        self.representative_fasta = representative_fasta
        self.output_dir = output_dir
        self.amplicon_size = amplicon_size
        self.n_pools = n_pools
        self.min_base_freq = min_base_freq

    def check_dependencies(self):
        result = subprocess.run(["which", "primalscheme3"], capture_output=True)
        if result.returncode != 0:
            print("\nMissing dependency: primalscheme3")
            print("Install with: pip install primalscheme3")
            return False
        return True

    def count_sequences(self):
        with open(self.representative_fasta, encoding="utf-8") as f:
            return sum(1 for line in f if line.startswith(">"))

    def run_primalscheme(self, msa_file, output_dir):
        print(f"\nDesigning primers using PrimalScheme3…")
        print(f"Amplicon size: {self.amplicon_size} bp")
        print(f"Pools: {self.n_pools}")
        print(f"Minimum base frequency: {self.min_base_freq}")
        cmd = [
            "primalscheme3", "scheme-create",
            "--msa", msa_file,
            "--output", str(output_dir),
            "--amplicon-size", str(self.amplicon_size),
            "--n-pools", str(self.n_pools),
            "--min-base-freq", str(self.min_base_freq),
            "--force",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"PrimalScheme3 failed: {result.stderr}")
            if "No valid primerpairs found" in result.stderr:
                print(
                    "\nHint: this usually means the input sequences are too diverse for any single "
                    "primer site to work across all of them (e.g. mixed subclades, or real biological "
                    "length/indel variation). Try raising --min-base-freq further (e.g. 0.5) to let PrimalScheme3 "
                    "ignore rare variant bases when designing primers, rather than requiring 100% coverage."
                )
            return False
        return True

    # We no longer compare primers against PrimalScheme3's single chosen reference (PrimerSequenceMismatch.py
    # checks against every sequence instead), so its reference-specific outputs (reference.fasta,
    # plot.html, primer.html, config.json, work/*) aren't useful to keep — just the bed files
    # (still the real primer coordinates/sequences) and the representative FASTA that was the input.
    def cleanup_scheme_output(self, output_dir, file_prefix):
        representative_name = Path(self.representative_fasta).name
        reference_name = f"{file_prefix}_reference.fasta"

        for path in output_dir.rglob("*"):
            if not path.is_file():
                continue
            if path.suffix == ".bed" or path.name == representative_name:
                continue
            path.unlink()

        # flatten the representative FASTA up to output_dir (PrimalScheme3 may nest it, e.g. in
        # work/), renaming it to "<name>_reference.fasta" — inside the primer design output it now
        # serves as the reference set primers were designed against (and PrimerSequenceMismatch.py checks against)
        for path in output_dir.rglob(representative_name):
            path.rename(output_dir / reference_name)

        # remove now-empty subdirectories
        for path in sorted(output_dir.rglob("*"), reverse=True):
            if path.is_dir() and not any(path.iterdir()):
                path.rmdir()

        # add a column-name header to each bed file (as a "#" comment, so it stays valid
        # BED format and anything that already skips "#" lines — e.g. PrimerSequenceMismatch.py — is unaffected)
        self.add_bed_header(output_dir / "primer.bed",
                             ["chrom", "start", "end", "name", "pool", "strand", "sequence", "pc"])
        self.add_bed_header(output_dir / "amplicon.bed",
                             ["chrom", "start", "end", "name", "pool"])
        self.add_bed_header(output_dir / "primertrim.amplicon.bed",
                             ["chrom", "start", "end", "name", "pool"])

        # prefix bed filenames with the scheme name so files from different schemes are
        # identifiable at a glance once collected outside their per-scheme output folder
        for name in ("primer.bed", "amplicon.bed", "primertrim.amplicon.bed"):
            path = output_dir / name
            if path.exists():
                path.rename(output_dir / f"{file_prefix}_{name}")

        return output_dir / reference_name

    @staticmethod
    def add_bed_header(bed_path, columns):
        if not bed_path.exists():
            return
        with open(bed_path, encoding="utf-8") as f:
            lines = f.readlines()
        comment_lines = [l for l in lines if l.startswith("#")]
        data_lines = [l for l in lines if not l.startswith("#")]
        with open(bed_path, "w", encoding="utf-8") as f:
            f.writelines(comment_lines)
            f.write("# " + "\t".join(columns) + "\n")
            f.writelines(data_lines)

    @staticmethod
    def count_bed_rows(bed_path):
        if not bed_path.exists():
            return 0
        with open(bed_path, encoding="utf-8") as f:
            return sum(1 for line in f if line.strip() and not line.startswith("#"))

    # a lone folder-name token stripped from a primer name (e.g. "<hash>_15_LEFT") groups every
    # variant PrimalScheme3 generated at that amplicon+strand position — the trailing "_<n>" is
    # the variant index, dropped here so variants of the same position count as one group
    @staticmethod
    def strip_variant(name):
        prefix, _, suffix = name.rpartition("_")
        return prefix if suffix.isdigit() else name

    @staticmethod
    def write_primer_info(primer_bed_path, info_tsv_path, start_time, command, total_pairs):
        rows = []
        with open(primer_bed_path, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                name, pool, seq = parts[3], parts[4], parts[6]
                length = len(seq)
                gc = sum(1 for base in seq.upper() if base in "GC") / length * 100
                tm = primer3.calc_tm(seq)
                rows.append((name, pool, seq, length, gc, tm))

        variant_groups = defaultdict(int)
        for name, *_ in rows:
            variant_groups[PrimerDesignWorkflow.strip_variant(name)] += 1
        variant_counts = list(variant_groups.values())

        lengths = [r[3] for r in rows]
        gcs = [r[4] for r in rows]
        tms = [r[5] for r in rows]

        with open(info_tsv_path, "w", encoding="utf-8") as f:
            f.write(f"Time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Command: {command}\n\n")

            f.write("name\tpool\tseq\tlength\t%gc\ttm\n")
            for name, pool, seq, length, gc, tm in rows:
                f.write(f"{name}\t{pool}\t{seq}\t{length}\t{gc:.2f}\t{tm:.2f}\n")

            f.write(f"\nTotal primers: {len(rows)}\n")
            f.write(f"Total primer pairs: {total_pairs}\n")

            f.write(f"\nMinimum alternative primers: {min(variant_counts) if variant_counts else 0}\n")
            f.write(f"Maximum alternative primers: {max(variant_counts) if variant_counts else 0}\n")

            f.write(f"\nSize range: {min(lengths)}-{max(lengths)} bp\n")
            f.write(f"Average size: {sum(lengths) / len(lengths):.0f} bp\n")

            f.write(f"\n%gc range: {min(gcs):.2f}-{max(gcs):.2f}\n")
            f.write(f"Average %gc: {sum(gcs) / len(gcs):.2f}\n")

            f.write(f"\nTm range: {min(tms):.2f}-{max(tms):.2f}\n")
            f.write(f"Average Tm: {sum(tms) / len(tms):.2f}\n")

    def run(self):
        start_time = datetime.datetime.now()
        command = " ".join([os.path.basename(sys.argv[0])] + sys.argv[1:])

        stem = Path(self.representative_fasta).stem
        if stem.endswith("_representative_aligned"):
            file_prefix = stem[:-len("_representative_aligned")]
        elif stem.endswith("_representative"):
            file_prefix = stem[:-len("_representative")]
        else:
            file_prefix = stem
        folder_name = f"{file_prefix}_primers"
        output_dir = Path(self.output_dir) / folder_name

        if not os.path.exists(self.representative_fasta):
            print(f"\n✗ Input file not found: {self.representative_fasta}")
            sys.exit(1)

        seq_count = self.count_sequences()
        if seq_count < 2:
            print(f"\n✗ Not enough sequences ({seq_count}) for primer design. Skipping.")
            sys.exit(0)

        if not self.check_dependencies():
            sys.exit(1)

        output_dir.mkdir(parents=True, exist_ok=True)

        if not self.run_primalscheme(self.representative_fasta, output_dir):
            sys.exit(1)

        reference_fasta = self.cleanup_scheme_output(output_dir, file_prefix)

        primer_bed = output_dir / f"{file_prefix}_primer.bed"
        amplicon_bed = output_dir / f"{file_prefix}_amplicon.bed"
        trimmed_bed = output_dir / f"{file_prefix}_primertrim.amplicon.bed"

        primer_count = self.count_bed_rows(primer_bed)
        amplicon_count = self.count_bed_rows(amplicon_bed)
        trimmed_count = self.count_bed_rows(trimmed_bed)

        primer_info_tsv = output_dir / f"{file_prefix}_summary.tsv"
        self.write_primer_info(primer_bed, primer_info_tsv, start_time, command, amplicon_count)

        print(f"✓ Saved {seq_count} reference sequences to: {reference_fasta}")
        print(f"✓ Saved {primer_count} primers to: {primer_bed}")
        print(f"✓ Saved {amplicon_count} amplicons to: {amplicon_bed}")
        print(f"✓ Saved {trimmed_count} trimmed amplicons to: {trimmed_bed}")
        print(f"✓ Saved summary (length/%GC/Tm) to: {primer_info_tsv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Design primers from representative sequences using PrimalScheme3")
    parser.add_argument(
        "-i", "--input", required=True,
        help="Representative FASTA file produced by Treemmer.py",
    )
    parser.add_argument(
        "-o", "--output-dir", default="results",
        help="Directory where the output folder will be saved (default: results/)",
    )
    parser.add_argument(
        "--amplicon-size", type=int, default=400,
        help="Amplicon size in bp (default: 400)",
    )
    parser.add_argument(
        "--n-pools", type=int, default=2,
        help="Number of sequencing pools (default: 2)",
    )
    parser.add_argument(
        "--min-base-freq", type=float, default=0.2,
        help="Minimum frequency (0-1) a variant base must reach across the input sequences to be considered "
             "when designing a degenerate primer (default: 0.2)",
    )
    args = parser.parse_args()

    workflow = PrimerDesignWorkflow(args.input, args.output_dir, args.amplicon_size, args.n_pools, args.min_base_freq)
    workflow.run()