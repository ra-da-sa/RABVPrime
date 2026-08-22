#!/usr/bin/env python3
"""Single entry point for all 6 use cases — infers which pipeline to run from which inputs you
provide (--country and/or --fasta) and whether you're designing new primers and then evaluating
them (--design-and-evaluate) or only checking an existing/external set (--evaluate), then
dispatches to the matching run_pipeline*.sh script.
"""
import argparse
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent / "scripts"

# matches sanitize() in RetrieveSequences.py / SubmitSequences.py / etc.: lowercase, collapse
# whitespace to "_", "/" to "_"
def sanitize(text):
    return "_".join(part for part in text.strip().lower().replace("/", "_").split() if part)

# matches the PRIMER_BED_BASENAME stripping in run_pipeline_verify.sh / run_pipeline_verify_combine.sh
def primer_bed_stem(primer_bed_path):
    name = Path(primer_bed_path).name
    if name.endswith(".bed"):
        name = name[: -len(".bed")]
    if name.endswith(".primer"):
        name = name[: -len(".primer")]
    return sanitize(name)

# (mode, has_country, has_fasta) -> (use case number, script name)
USE_CASES = {
    ("design",   True,  False): (1, "run_pipeline.sh"),
    ("design",   False, True):  (2, "run_pipeline_submit.sh"),
    ("design",   True,  True):  (3, "run_pipeline_combine.sh"),
    ("evaluate", True,  False): (4, "run_pipeline_verify.sh"),
    ("evaluate", False, True):  (5, "run_pipeline_verify_submit.sh"),
    ("evaluate", True,  True):  (6, "run_pipeline_verify_combine.sh"),
}

# builds the positional args (everything except output_dir, which the caller appends) for the
# underlying run_pipeline*.sh script, in the exact order each one expects
def build_positional_args(use_case, args):
    if use_case == 1:
        return [args.country, args.database]
    if use_case == 2:
        return [args.fasta, args.database]
    if use_case == 3:
        return [args.country, args.fasta, args.database]
    if use_case == 4:
        return [args.country, args.primers, args.database]
    if use_case == 5:
        # run_pipeline_verify_submit.sh takes an <output_name> instead of a database (it never
        # touches RABV-DB) — derive one the same way UC4/UC6 derive their own output folder names
        output_name = f"{sanitize(Path(args.fasta).stem)}_{primer_bed_stem(args.primers)}"
        return [args.fasta, args.primers, output_name]
    if use_case == 6:
        return [args.country, args.fasta, args.primers, args.database]
    raise AssertionError(f"unhandled use case: {use_case}")

def main():
    parser = argparse.ArgumentParser(
        description="Run any of the 6 pipeline use cases through one entry point.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  Use case 1 (design + evaluate, DB only):
    RABVPrime.py --design-and-evaluate --country Morocco --database rabv_V-gDB_10042026.db

  Use case 2 (design + evaluate, own FASTA only):
    RABVPrime.py --design-and-evaluate --fasta EA_2024.fasta --database rabv_V-gDB_10042026.db

  Use case 3 (design + evaluate, DB + own FASTA):
    RABVPrime.py --design-and-evaluate --country Morocco --fasta EA_2024.fasta --database rabv_V-gDB_10042026.db

  Use case 4 (evaluate only, DB only):
    RABVPrime.py --evaluate --country Morocco --primers rabv_ea.primer.bed --database rabv_V-gDB_10042026.db

  Use case 5 (evaluate only, own FASTA only):
    RABVPrime.py --evaluate --fasta rabv_ea.reference.fasta --primers rabv_ea.primer.bed

  Use case 6 (evaluate only, DB + own FASTA):
    RABVPrime.py --evaluate --country Morocco --fasta rabv_ea.reference.fasta --primers rabv_ea.primer.bed --database rabv_V-gDB_10042026.db
""",
    )

    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument("--design-and-evaluate", action="store_true",
                             help="Design new primers, then evaluate them against the same sequences")
    mode_group.add_argument("--evaluate", action="store_true",
                             help="Evaluate an existing/external primer set only — no design step")

    parser.add_argument("-c", "--country", default=None, help="Country name to retrieve from RABV-DB (e.g. 'Morocco')")
    parser.add_argument("--fasta", default=None, help="Your own FASTA file of sequences")
    parser.add_argument("--primers", default=None, help="Existing primer.bed to check (required with --evaluate)")
    parser.add_argument("-db", "--database", default=None, help="Path to the RABV-DB SQLite database file")
    parser.add_argument("-o", "--output_dir", default="results", help="Output directory (default: results/)")
    parser.add_argument("--dry-run", action="store_true", help="Print which use case and command would run, without running it")

    args = parser.parse_args()
    mode = "design" if args.design_and_evaluate else "evaluate"

    if not args.country and not args.fasta:
        parser.error("at least one of --country or --fasta is required")
    if mode == "design" and args.primers:
        parser.error("--primers is only used with --evaluate — --design-and-evaluate generates its own primers")
    if mode == "evaluate" and not args.primers:
        parser.error("--evaluate requires --primers <primer.bed>")

    use_case, script_name = USE_CASES[(mode, bool(args.country), bool(args.fasta))]

    # use case 5 (evaluate, own FASTA only) is the one pipeline that never touches RABV-DB
    if use_case != 5 and not args.database:
        parser.error(f"--database is required for use case {use_case} ({script_name})")

    positional_args = build_positional_args(use_case, args) + [args.output_dir]
    script_path = SCRIPT_DIR / script_name

    print(f"→ Use case {use_case}: {script_name} {' '.join(positional_args)}", flush=True)
    if args.dry_run:
        return

    result = subprocess.run(["bash", str(script_path), *positional_args])
    sys.exit(result.returncode)

if __name__ == "__main__":
    main()
