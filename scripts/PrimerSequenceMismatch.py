import argparse
import csv
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from datetime import datetime
from pathlib import Path

# -----------------------------
# 1. NUCLEOTIDE DEFINITIONS
# -----------------------------
ntCharToSubChars = {
    "A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"],
    "R": ["A", "G"], "Y": ["C", "T"], "K": ["G", "T"],
    "M": ["A", "C"], "S": ["C", "G"], "W": ["A", "T"],
    "B": ["C", "G", "T"], "D": ["A", "G", "T"],
    "H": ["A", "C", "T"], "V": ["A", "C", "G"],
    "N": ["A", "C", "G", "T"]
}

knownMismatchIssues = [
    {"ppId": "RdRP_SARSr-P1", "mismatch": "R15480C"},
    {"ppId": "RdRP_SARSr-P1", "mismatch": "T15489A"},
]

def build_nt_char_to_allowed(nt_char_to_subchars):
    result = {}
    for nt_char, sub_chars in nt_char_to_subchars.items():
        if len(sub_chars) == 1:
            result[nt_char] = [nt_char]
        else:
            full = list(sub_chars)
            for other_char, other_subs in nt_char_to_subchars.items():
                if set(other_subs).issubset(set(sub_chars)) and other_char not in full:
                    full.append(other_char)
            result[nt_char] = full
    return result

ntCharToAllowed = build_nt_char_to_allowed(ntCharToSubChars)

# -----------------------------
# 2. SEQUENCE MODIFICATIONS
# -----------------------------
def apply_mod(seq, mod):
    if mod["type"] == "replaceSingle":
        i = mod["loc"] - 1
        return seq[:i] + mod["replacement"] + seq[i + 1:]

    elif mod["type"] == "replaceRegion":
        start = mod["locStart"] - 1
        end = mod["locEnd"]
        return seq[:start] + mod["replacement"] + seq[end:]

    elif mod["type"] == "insertRegion":
        loc = mod["afterLoc"]
        return seq[:loc] + mod["insertion"] + seq[loc:]

    return seq


# -----------------------------
# 3. MISMATCH ANALYSIS
# -----------------------------
def is_known_mismatch(pp_id, mismatch_string):
    return any(
        k["ppId"] == pp_id and k["mismatch"] == mismatch_string
        for k in knownMismatchIssues
    )

def analyze_primer(primer_id, primer_seq, reference_seq, start_pos):
    """
    Compare primer against reference sequence
    """
    issues = []

    for i, base in enumerate(primer_seq):
        ref_pos = start_pos + i
        ref_base = reference_seq[ref_pos - 1]

        allowed = ntCharToAllowed[base]

        if ref_base not in allowed:
            mismatch_string = f"{base}{ref_pos}{ref_base}"
            known = is_known_mismatch(primer_id, mismatch_string)

            issues.append({
                "type": "mismatch",
                "primerBase": base,
                "refBase": ref_base,
                "position": ref_pos,
                "known": known
            })

    return issues

# -----------------------------
# 4. REPORT GENERATION
# -----------------------------
def todays_date():
    return datetime.now().strftime("%d/%m/%Y")

def generate_report(primer_id, primer_seq, reference_seq, start_pos):
    issues = analyze_primer(primer_id, primer_seq, reference_seq, start_pos)

    report = {
        "primer": primer_id,
        "sequence": primer_seq,
        "issues": issues,
        "numIssues": len(issues),
        "date": todays_date()
    }

    return report

# -----------------------------
# 5. CLI: check primers (from PrimerDesign.py's primer.bed) against a FASTA of sequences
# -----------------------------
def read_fasta(path):
    seqs = {}
    header, parts = None, []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(parts)
                header, parts = line[1:].strip().split()[0], []
            elif line:
                parts.append(line)
        if header is not None:
            seqs[header] = "".join(parts)
    return seqs

def read_primer_bed(path):
    primers = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            primers.append({"name": parts[3], "pool": parts[4], "strand": parts[5], "seq": parts[6]})
    return primers

# locates every primer's best-matching site in target_seq with a SINGLE BLAST call (all primers
# as one multi-fasta query — far fewer subprocess spawns than one call per primer), handling
# strand automatically. Returns {primer_name: {"qseq":..., "sseq":..., "pident":...}}
def blast_locate_all(primers, target_header, target_seq):
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as qf:
        for p in primers:
            qf.write(f">{p['name']}\n{p['seq']}\n")
        query_path = qf.name
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as sf:
        sf.write(f">{target_header}\n{target_seq.replace('-', '')}\n")
        subject_path = sf.name
    try:
        result = subprocess.run(
            ["blastn", "-task", "blastn-short", "-query", query_path, "-subject", subject_path,
             "-outfmt", "6 qseqid qseq sseq pident", "-dust", "no",
             "-max_target_seqs", "1", "-max_hsps", "1"],
            capture_output=True, text=True,
        )
        hits = {}
        for line in result.stdout.strip().splitlines():
            if not line.strip():
                continue
            qseqid, qseq, sseq, pident = line.split("\t")
            if qseqid not in hits:  # first row per primer = best HSP (BLAST sorts by score)
                hits[qseqid] = {"qseq": qseq, "sseq": sseq, "pident": float(pident)}
        return hits
    finally:
        os.unlink(query_path)
        os.unlink(subject_path)

# compares the BLAST-aligned primer/target substrings base-by-base, IUPAC-ambiguity aware.
# position is 1-based and counted along the primer (query) itself, not the reference genome.
def compare_aligned(qseq, sseq):
    issues = []
    for i, (q, s) in enumerate(zip(qseq, sseq), start=1):
        if q == "-" or s == "-":
            issues.append({"type": "indel", "position": i, "primerBase": q, "refBase": s})
            continue
        allowed = ntCharToAllowed.get(q.upper(), [q.upper()])
        if s.upper() not in allowed:
            issues.append({"type": "mismatch", "position": i, "primerBase": q, "refBase": s})
    return issues

def safe_filename(header):
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in header)

# PrimalScheme3 names primers "<scheme-hash>_<n>_LEFT_<variant>" / "<scheme-hash>_<n>_RIGHT_<variant>",
# but an externally-supplied primer.bed (e.g. checking someone else's scheme) may omit the
# trailing "_<variant>" entirely (e.g. "rabv_ea_1_LEFT") — strip whichever form is present to get
# the amplicon a primer belongs to
def amplicon_name(primer_name):
    for tag in ("_LEFT_", "_RIGHT_"):
        idx = primer_name.find(tag)
        if idx != -1:
            return primer_name[:idx]
    for tag in ("_LEFT", "_RIGHT"):
        if primer_name.endswith(tag):
            return primer_name[:-len(tag)]
    return primer_name

# same "with or without a trailing _<variant>" flexibility as amplicon_name() above
def is_left_primer(primer_name):
    return "_LEFT_" in primer_name or primer_name.endswith("_LEFT")

def is_right_primer(primer_name):
    return "_RIGHT_" in primer_name or primer_name.endswith("_RIGHT")

# the scheme-hash prefix is constant across one primer.bed, so the trailing number alone
# already uniquely identifies the amplicon within this report
def amplicon_number(primer_name):
    return amplicon_name(primer_name).rsplit("_", 1)[-1]

# mismatches within this many bases of the primer's 3' end get flagged by end in addition to their
# variant type, since they're the ones that most affect extension (matches the lab's R script:
# position >= primer_length - 5; there is no equivalent 5' check, by design, matching that script)
TERMINAL_WINDOW = 5

# every applicable label is reported together (e.g. an insertion that also falls in the 3' window
# is "insertion,3' mismatch", not just one or the other) — mirrors the lab's R script exactly,
# including labeling a plain, non-terminal base mismatch as "snp"
def classify_issue(issue, primer_length):
    position = issue["position"]
    labels = []
    if issue["type"] == "indel":
        # primerBase == "-": the primer is missing a base the target has (primer deletion).
        # refBase == "-": the primer has a base the target doesn't (primer insertion).
        labels.append("deletion" if issue["primerBase"] == "-" else "insertion")
    if position >= primer_length - TERMINAL_WINDOW:
        labels.append("3' mismatch")
    if not labels:
        return "snp"
    return ",".join(labels)

# dataset-wide coverage table — one row per amplicon, identical across every accession's file.
# a sequence counts as "covered" by an amplicon only if at least one forward AND at least one
# reverse primer variant for it is clean against that sequence (PrimalScheme3 generates several
# variants per side specifically so any one of them working is enough).
def calculate_coverage(primers, accessions, primer_issue_accessions):
    amplicons = defaultdict(list)
    for primer in primers:
        amplicons[amplicon_name(primer["name"])].append(primer)

    total = len(accessions)
    coverage_rows = []
    for amp_id, amp_primers in amplicons.items():
        forward = [p for p in amp_primers if is_left_primer(p["name"])]
        reverse = [p for p in amp_primers if is_right_primer(p["name"])]

        covered = 0
        for acc in accessions:
            fwd_clean = any(acc not in primer_issue_accessions.get(p["name"], ()) for p in forward)
            rev_clean = any(acc not in primer_issue_accessions.get(p["name"], ()) for p in reverse)
            if fwd_clean and rev_clean:
                covered += 1

        coverage_rows.append({
            "amplicon": amplicon_number(amp_primers[0]["name"]),
            "forward": "/".join(p["seq"] for p in forward) if forward else "-",
            "reverse": "/".join(p["seq"] for p in reverse) if reverse else "-",
            "coverage": f"{covered / total * 100:.0f}%",
        })
    coverage_rows.sort(key=lambda r: int(r["amplicon"]) if r["amplicon"].isdigit() else r["amplicon"])
    return coverage_rows

def run(bed_path, fasta_path, output_dir):
    primers = read_primer_bed(bed_path)
    sequences = read_fasta(fasta_path)

    if not primers:
        print(f"✗ No primers found in: {bed_path}")
        sys.exit(1)
    if not sequences:
        print(f"✗ No sequences found in: {fasta_path}")
        sys.exit(1)

    if subprocess.run(["which", "blastn"], capture_output=True).returncode != 0:
        print("✗ blastn not found — install BLAST+ (brew install blast).")
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Comparing {len(primers)} primers against {len(sequences)} reference sequences using BLAST...\n")
    print("Assessing primer mismatches...\n")

    detail_fieldnames = ["amplicon", "primer", "position_in_primer", "ref_base", "primer_base", "types"]
    coverage_fieldnames = ["amplicon", "forward", "reverse", "coverage"]

    # per primer, the set of accessions it had *any* issue against — used to build the coverage
    # table below (a primer can have several mismatch rows in one accession; that's still only
    # one accession the primer doesn't cleanly work on)
    primer_issue_accessions = defaultdict(set)
    per_accession_rows = {}

    for target_header, target_seq in sequences.items():
        hits = blast_locate_all(primers, target_header, target_seq)

        rows = []
        for primer in primers:
            hit = hits.get(primer["name"])
            if not hit:
                rows.append({
                    "amplicon": amplicon_number(primer["name"]), "primer": primer["name"],
                    "position_in_primer": "-", "ref_base": "no alignment found", "primer_base": "-", "types": "-",
                })
                primer_issue_accessions[primer["name"]].add(target_header)
                continue

            if hit["pident"] >= 100.0:
                continue  # perfect match — not flagged, skip from the report

            issues = compare_aligned(hit["qseq"], hit["sseq"])
            if not issues:
                continue

            primer_issue_accessions[primer["name"]].add(target_header)
            for issue in issues:
                rows.append({
                    "amplicon": amplicon_number(primer["name"]), "primer": primer["name"],
                    "position_in_primer": issue["position"],
                    "ref_base": issue["refBase"], "primer_base": issue["primerBase"],
                    "types": classify_issue(issue, len(primer["seq"])),
                })

        per_accession_rows[target_header] = rows

    print("Calculating amplicon coverage...")
    coverage_rows = calculate_coverage(primers, list(sequences.keys()), primer_issue_accessions)

    for target_header, rows in per_accession_rows.items():
        accession_path = output_dir / f"{safe_filename(target_header)}_mismatches.tsv"
        with open(accession_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=coverage_fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(coverage_rows)
            f.write("\n")
            writer = csv.DictWriter(f, fieldnames=detail_fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    print(f"\n✓ Saved {len(sequences)} mismatch reports (one per reference sequence) to: {output_dir.resolve()}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Evaluate primers for mismatches against each sequence provided in a FASTA file. "
                     "An individual TSV file is generated for each accession, containing primer coverage "
                     "and detailed information on primer-sequence mismatches."
    )
    parser.add_argument("-b", "--bed", dest="primer_bed", required=True,
                         help="Primer BED file containing the designed primers, either generated by "
                              "PrimerDesign.py or provided by the user.")
    parser.add_argument("-f", "--fasta", default=None,
                         help="FASTA file of sequences against which the primers are evaluated. "
                              "(default: <name>_qc/<name>_qc_sequences.fasta)")
    parser.add_argument("-o", "--output-dir", default=None,
                         help="Directory where the per-accession mismatch TSV files will be saved. "
                              "(default: <name>_mismatches/)")
    args = parser.parse_args()

    bed_dir = Path(args.primer_bed).parent
    name = bed_dir.name[:-len("_primers")] if bed_dir.name.endswith("_primers") else bed_dir.name

    if args.fasta:
        fasta_path = args.fasta
    else:
        fasta_path = bed_dir.parent / f"{name}_qc" / f"{name}_qc_sequences.fasta"
        if not fasta_path.exists():
            print(f"✗ No post-QC FASTA found at: {fasta_path}")
            print("  Run QualityControl.py first, or pass one explicitly with -f/--fasta.")
            sys.exit(1)
        print(f"Using post-QC FASTA: {fasta_path}")

    output_dir = Path(args.output_dir) if args.output_dir else bed_dir.parent / f"{name}_mismatches"

    run(args.primer_bed, fasta_path, output_dir)
