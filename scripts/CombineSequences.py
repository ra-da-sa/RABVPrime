#!/usr/bin/env python3
import argparse
import datetime
import os
import sys
from pathlib import Path

def read_fasta_records(path):
    records = []
    header, parts = None, []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(parts)))
                header, parts = line[1:].strip().split()[0], []
            elif line:
                parts.append(line)
        if header is not None:
            records.append((header, "".join(parts)))
    return records

# header-driven, case-insensitive — returns {accession: (country, subclade, year, host)}
def read_metadata(path):
    lookup = {}
    if not os.path.exists(path):
        return lookup
    with open(path, encoding="utf-8") as f:
        header_lower = [h.lower() for h in f.readline().rstrip("\n").split("\t")]
        try:
            id_idx       = header_lower.index("accession")
            country_idx  = header_lower.index("country")
            subclade_idx = header_lower.index("subclade")
            year_idx     = header_lower.index("year")
            host_idx     = header_lower.index("host")
        except ValueError:
            return lookup
        for line in f:
            row = line.rstrip("\n").split("\t")
            if len(row) > max(id_idx, country_idx, subclade_idx, year_idx, host_idx):
                lookup[row[id_idx]] = (row[country_idx], row[subclade_idx], row[year_idx], row[host_idx])
    return lookup

def folder_stem(dir_path):
    name = Path(dir_path).name
    return name[:-len("_sequences")] if name.endswith("_sequences") else name

def run(retrieved_dir, submitted_dir, output_dir):
    retrieved_dir = Path(retrieved_dir)
    submitted_dir = Path(submitted_dir)

    retrieved_fasta = retrieved_dir / f"{retrieved_dir.name}_sequences.fasta"
    submitted_fasta = submitted_dir / f"{submitted_dir.name}_sequences.fasta"
    retrieved_meta  = retrieved_dir / f"{retrieved_dir.name}_metadata.tsv"
    submitted_meta  = submitted_dir / f"{submitted_dir.name}_metadata.tsv"

    for path in (retrieved_fasta, submitted_fasta):
        if not path.exists():
            print(f"✗ Input file not found: {path}")
            sys.exit(1)

    retrieved_records = read_fasta_records(retrieved_fasta)
    submitted_records = read_fasta_records(submitted_fasta)

    # both inputs must already share one alignment coordinate system (RetrieveSequences.py and
    # SubmitSequences.py both align against the same reference panel) — if they don't, combining
    # would silently misalign every sequence rather than fail loudly, so check first
    lengths = {len(seq) for _, seq in retrieved_records + submitted_records}
    if len(lengths) > 1:
        print(f"✗ Alignment length mismatch between inputs: {sorted(lengths)}")
        print("  Both inputs must share one coordinate system — re-run RetrieveSequences.py/"
              "SubmitSequences.py against the same database so they align against the same "
              "reference panel.")
        sys.exit(1)

    retrieved_meta_lookup = read_metadata(retrieved_meta)
    submitted_meta_lookup = read_metadata(submitted_meta)

    # combine, de-duplicating by accession — the same accession can legitimately appear in both
    # (e.g. a submitted sequence that SubmitSequences.py found already exists in the DB, and that
    # DB accession also happens to match the retrieved filter); keep the retrieved copy since its
    # metadata is always fully populated from the database
    combined = {seq_id: seq for seq_id, seq in submitted_records}
    overlap = [seq_id for seq_id, _ in retrieved_records if seq_id in combined]
    for seq_id, seq in retrieved_records:
        combined[seq_id] = seq

    # No <4-sequences guard here — QualityControl.py (the next pipeline step in every use case)
    # already checks the post-QC count against Treemmer's minimum, so checking again here would
    # just be a redundant, earlier-but-less-accurate copy of the same check.

    retrieved_stem = folder_stem(retrieved_dir)
    submitted_stem = folder_stem(submitted_dir)
    combined_stem  = f"{retrieved_stem}_{submitted_stem}"

    folder_name  = f"{combined_stem}_sequences"
    run_dir      = Path(output_dir) / folder_name
    run_dir.mkdir(parents=True, exist_ok=True)
    fasta_path   = run_dir / f"{folder_name}_sequences.fasta"
    tsv_path     = run_dir / f"{folder_name}_metadata.tsv"
    summary_path = run_dir / f"{folder_name}_summary.txt"

    with open(fasta_path, "w", encoding="utf-8") as f:
        for seq_id, seq in combined.items():
            f.write(f">{seq_id}\n{seq}\n")

    with open(tsv_path, "w", encoding="utf-8") as f:
        f.write("Accession\tCountry\tSubclade\tYear\tHost\n")
        for seq_id in combined:
            country, subclade, year, host = (
                retrieved_meta_lookup.get(seq_id)
                or submitted_meta_lookup.get(seq_id)
                or ("", "", "", "")
            )
            f.write(f"{seq_id}\t{country}\t{subclade}\t{year}\t{host}\n")

    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(f"Time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("\nCombined from:\n")
        f.write(f"  retrieved: {retrieved_dir} ({len(retrieved_records)} sequences)\n")
        f.write(f"  submitted: {submitted_dir} ({len(submitted_records)} sequences)\n")
        if overlap:
            f.write(f"\n{len(overlap)} accession(s) present in both inputs (kept the retrieved "
                    f"copy): {', '.join(overlap)}\n")
        f.write(f"\nTotal combined sequences: {len(combined)}\n")
        f.write(f"Alignment length: {next(iter(lengths))} bp\n")

    print(f"✓ Combined {len(retrieved_records)} retrieved + {len(submitted_records)} submitted "
          f"= {len(combined)} unique sequence(s)")
    if overlap:
        print(f"  ({len(overlap)} accession(s) present in both — kept the retrieved copy)")
    print(f"✓ Saved combined sequences to: {fasta_path}")
    print(f"✓ Saved metadata to: {tsv_path}")
    print(f"✓ Saved summary to: {summary_path}")

    return run_dir

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Combine a RetrieveSequences.py output folder and a SubmitSequences.py "
                     "output folder into one unified <name>_sequences/ folder for QualityControl.py "
                     "— both inputs must share the same alignment coordinate system."
    )
    parser.add_argument(
        "--retrieved", required=True,
        help="RetrieveSequences.py output folder",
    )
    parser.add_argument(
        "--submitted", required=True,
        help="SubmitSequences.py output folder",
    )
    parser.add_argument(
        "-o", "--output_dir", default="results",
        help="Directory where the combined output folder will be saved (default: results/)",
    )
    args = parser.parse_args()

    for label, path in (("--retrieved", args.retrieved), ("--submitted", args.submitted)):
        if not os.path.isdir(path):
            print(f"✗ Not a directory: {path} ({label})")
            sys.exit(1)

    run(args.retrieved, args.submitted, args.output_dir)
