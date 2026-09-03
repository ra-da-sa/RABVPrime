#!/usr/bin/env python3
"""Single entry point for all 6 use cases — infers which pipeline to run from which inputs you
provide (--country and/or --subclade, and/or --fasta) and whether you're designing new primers
and then evaluating them (--design-and-evaluate) or only checking an existing/external set
(--evaluate), then runs that pipeline's steps directly.

--country and --subclade are both optional database-side filters and can be used alone or
together (e.g. only sequences that are both from Morocco AND subclade AF1a) — this mirrors what
RetrieveSequences.py itself already supports.
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

# matches RetrieveSequences.py's own stem formula: sanitized countries, then sanitized
# subclades, joined with "_"
def db_stem(countries, subclades):
    parts = [sanitize(c) for c in (countries or [])] + [sanitize(s) for s in (subclades or [])]
    return "_".join(parts)

# matches the PRIMER_BED_BASENAME stripping previously done in run_pipeline_verify.sh /
# run_pipeline_verify_combine.sh
def primer_bed_stem(primer_bed_path):
    name = Path(primer_bed_path).name
    if name.endswith(".bed"):
        name = name[: -len(".bed")]
    if name.endswith(".primer"):
        name = name[: -len(".primer")]
    return sanitize(name)

# (mode, has_db_filter, has_fasta) -> use case number
USE_CASES = {
    ("design",   True,  False): 1,
    ("design",   False, True):  2,
    ("design",   True,  True):  3,
    ("evaluate", True,  False): 4,
    ("evaluate", False, True):  5,
    ("evaluate", True,  True):  6,
}

class Pipeline:
    def __init__(self, use_case, args):
        self.use_case = use_case
        self.args = args
        self.step_num = 0
        self.total_steps = {1: 6, 2: 6, 3: 8, 4: 4, 5: 2, 6: 5}[use_case]

    def announce(self, label):
        self.step_num += 1
        print(f"\n=== {self.step_num}/{self.total_steps} {label} ===")

    def run_script(self, script_name, script_args):
        command = [sys.executable, str(SCRIPT_DIR / script_name), *[str(a) for a in script_args]]
        print(" ".join(command))
        result = subprocess.run(command)
        if result.returncode != 0:
            sys.exit(result.returncode)

    def require_file(self, path, step_label):
        if not Path(path).is_file():
            print(f"✗ Expected output not found: {path}")
            print(f"  Pipeline stopped after: {step_label}")
            sys.exit(1)

    def require_dir(self, path, step_label):
        if not Path(path).is_dir():
            print(f"✗ Expected output not found: {path}")
            print(f"  Pipeline stopped after: {step_label}")
            sys.exit(1)

    # shared tail of every use case: sequence-level then population-level mismatch analysis
    # against an already-designed/-provided primer.bed
    def run_mismatch_steps(self, primer_bed, qc_fasta, mismatches_dir,
                            countries=None, subclades=None, fasta_name=None):
        self.announce("Sequence-level mismatches")
        script_args = ["-b", primer_bed]
        if qc_fasta is not None:
            script_args += ["-f", qc_fasta, "-o", mismatches_dir]
        self.run_script("PrimerSequenceMismatch.py", script_args)
        self.require_dir(mismatches_dir, "PrimerSequenceMismatch.py")

        self.announce("Population-level mismatches")
        script_args = ["-i", mismatches_dir, "-b", primer_bed]
        if countries:
            script_args += ["--country", *countries]
        if subclades:
            script_args += ["--subclade", *subclades]
        if fasta_name:
            script_args += ["--fasta-name", fasta_name]
        self.run_script("PrimerPopulationMismatch.py", script_args)
        print(f"\n✓ Pipeline complete — see {mismatches_dir}/ for final outputs.")

    # shared retrieve -> QC -> Treemmer -> PrimerDesign chain used by UC1/2/3 (design mode)
    def run_design_chain(self, seq_fasta, stem):
        self.announce("Quality control")
        self.run_script("QualityControl.py", ["-i", seq_fasta, "-o", self.args.output_dir])
        qc_fasta = f"{self.args.output_dir}/{stem}_qc/{stem}_qc_sequences.fasta"
        self.require_file(qc_fasta, "QualityControl.py")

        self.announce("Treemmer")
        self.run_script("Treemmer.py", ["-i", qc_fasta, "-o", self.args.output_dir])
        representative_fasta = f"{self.args.output_dir}/{stem}_treemmer/{stem}_representative.fasta"
        self.require_file(representative_fasta, "Treemmer.py")

        self.announce("Primer design")
        self.run_script("PrimerDesign.py", ["-i", representative_fasta, "-o", self.args.output_dir])
        primer_bed = f"{self.args.output_dir}/{stem}_primers/{stem}_primer.bed"
        self.require_file(primer_bed, "PrimerDesign.py")
        return primer_bed

    def retrieve(self, stem):
        self.announce("Retrieve sequences")
        script_args = ["-db", self.args.database, "-o", self.args.output_dir]
        if self.args.country:
            script_args += ["-c", *self.args.country]
        if self.args.subclade:
            script_args += ["-s", *self.args.subclade]
        self.run_script("RetrieveSequences.py", script_args)
        seq_fasta = f"{self.args.output_dir}/{stem}_sequences/{stem}_sequences_sequences.fasta"
        self.require_file(seq_fasta, "RetrieveSequences.py")
        return seq_fasta

    # use case 1: design + evaluate, database filter (country and/or subclade) only
    def uc1(self):
        stem = db_stem(self.args.country, self.args.subclade)
        seq_fasta = self.retrieve(stem)
        primer_bed = self.run_design_chain(seq_fasta, stem)
        mismatches_dir = f"{self.args.output_dir}/{stem}_mismatches"
        self.run_mismatch_steps(primer_bed, None, mismatches_dir,
                                 countries=self.args.country, subclades=self.args.subclade)

    # use case 2: design + evaluate, own FASTA only
    def uc2(self):
        stem = sanitize(Path(self.args.fasta).stem)
        self.announce("Submit sequences")
        self.run_script("SubmitSequences.py",
                         ["-i", self.args.fasta, "-db", self.args.database, "-o", self.args.output_dir])
        seq_fasta = f"{self.args.output_dir}/{stem}_sequences/{stem}_sequences_sequences.fasta"
        self.require_file(seq_fasta, "SubmitSequences.py")

        primer_bed = self.run_design_chain(seq_fasta, stem)
        mismatches_dir = f"{self.args.output_dir}/{stem}_mismatches"
        self.run_mismatch_steps(primer_bed, None, mismatches_dir,
                                 fasta_name=Path(self.args.fasta).name)

    # use case 3: design + evaluate, database filter (country and/or subclade) + own FASTA
    def uc3(self):
        country_stem = db_stem(self.args.country, self.args.subclade)
        self.retrieve(country_stem)
        retrieved_dir = f"{self.args.output_dir}/{country_stem}_sequences"

        submit_stem = sanitize(Path(self.args.fasta).stem)
        self.announce("Submit sequences")
        self.run_script("SubmitSequences.py",
                         ["-i", self.args.fasta, "-db", self.args.database, "-o", self.args.output_dir])
        submitted_dir = f"{self.args.output_dir}/{submit_stem}_sequences"
        self.require_file(f"{submitted_dir}/{submit_stem}_sequences_sequences.fasta", "SubmitSequences.py")

        combined_stem = f"{country_stem}_{submit_stem}"
        self.announce("Combine sequences")
        self.run_script("CombineSequences.py",
                         ["--retrieved", retrieved_dir, "--submitted", submitted_dir, "-o", self.args.output_dir])
        seq_fasta = f"{self.args.output_dir}/{combined_stem}_sequences/{combined_stem}_sequences_sequences.fasta"
        self.require_file(seq_fasta, "CombineSequences.py")

        primer_bed = self.run_design_chain(seq_fasta, combined_stem)
        mismatches_dir = f"{self.args.output_dir}/{combined_stem}_mismatches"
        self.run_mismatch_steps(primer_bed, None, mismatches_dir,
                                 countries=self.args.country, subclades=self.args.subclade,
                                 fasta_name=Path(self.args.fasta).name)

    # use case 4: evaluate only, database filter (country and/or subclade) only
    def uc4(self):
        stem = db_stem(self.args.country, self.args.subclade)
        self.retrieve(stem)
        self.announce("Quality control")
        self.run_script("QualityControl.py",
                         ["-i", f"{self.args.output_dir}/{stem}_sequences/{stem}_sequences_sequences.fasta",
                          "-o", self.args.output_dir])
        qc_fasta = f"{self.args.output_dir}/{stem}_qc/{stem}_qc_sequences.fasta"
        self.require_file(qc_fasta, "QualityControl.py")

        output_name = f"{stem}_{primer_bed_stem(self.args.primers)}"
        mismatches_dir = f"{self.args.output_dir}/{output_name}_mismatches"
        self.run_mismatch_steps(self.args.primers, qc_fasta, mismatches_dir,
                                 countries=self.args.country, subclades=self.args.subclade)

    # use case 5: evaluate only, own FASTA only — never touches the database
    def uc5(self):
        output_name = f"{sanitize(Path(self.args.fasta).stem)}_{primer_bed_stem(self.args.primers)}"
        mismatches_dir = f"{self.args.output_dir}/{output_name}_mismatches"
        self.run_mismatch_steps(self.args.primers, self.args.fasta, mismatches_dir,
                                 fasta_name=Path(self.args.fasta).name)

    # use case 6: evaluate only, database filter (country and/or subclade) + own FASTA
    def uc6(self):
        country_stem = db_stem(self.args.country, self.args.subclade)
        self.retrieve(country_stem)
        self.announce("Quality control")
        self.run_script("QualityControl.py",
                         ["-i", f"{self.args.output_dir}/{country_stem}_sequences/{country_stem}_sequences_sequences.fasta",
                          "-o", self.args.output_dir])
        qc_fasta = f"{self.args.output_dir}/{country_stem}_qc/{country_stem}_qc_sequences.fasta"
        self.require_file(qc_fasta, "QualityControl.py")

        own_stem = sanitize(Path(self.args.fasta).stem)
        combined_stem = f"{country_stem}_{own_stem}"
        self.announce("Combine QC-passed sequences with user-provided sequences")
        merged_dir = Path(self.args.output_dir) / f"{combined_stem}_verify_sequences"
        merged_fasta = merged_dir / f"{combined_stem}_verify_sequences.fasta"
        merged_dir.mkdir(parents=True, exist_ok=True)
        with open(merged_fasta, "w", encoding="utf-8") as out:
            out.write(Path(qc_fasta).read_text(encoding="utf-8"))
            out.write(Path(self.args.fasta).read_text(encoding="utf-8"))
        print(f"✓ Saved combined sequences to: {merged_fasta}")

        output_name = f"{combined_stem}_{primer_bed_stem(self.args.primers)}"
        mismatches_dir = f"{self.args.output_dir}/{output_name}_mismatches"
        self.run_mismatch_steps(self.args.primers, str(merged_fasta), mismatches_dir,
                                 countries=self.args.country, subclades=self.args.subclade,
                                 fasta_name=Path(self.args.fasta).name)

    def run(self):
        {1: self.uc1, 2: self.uc2, 3: self.uc3, 4: self.uc4, 5: self.uc5, 6: self.uc6}[self.use_case]()

def main():
    parser = argparse.ArgumentParser(
        description="Run any of the 6 pipeline use cases through one entry point.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  Use case 1 (design + evaluate, DB only, by country):
    RABVPrime.py --design-and-evaluate --country Morocco --database rabv_V-gDB_10042026.db

  Use case 1 (design + evaluate, DB only, by subclade):
    RABVPrime.py --design-and-evaluate --subclade AF1a --database rabv_V-gDB_10042026.db

  Use case 2 (design + evaluate, own FASTA only):
    RABVPrime.py --design-and-evaluate --fasta EA_2024.fasta --database rabv_V-gDB_10042026.db

  Use case 3 (design + evaluate, DB + own FASTA):
    RABVPrime.py --design-and-evaluate --country Morocco --fasta EA_2024.fasta --database rabv_V-gDB_10042026.db

  Use case 4 (evaluate only, DB only, by subclade):
    RABVPrime.py --evaluate --subclade AF1a --primers rabv_ea.primer.bed --database rabv_V-gDB_10042026.db

  Use case 5 (evaluate only, own FASTA only):
    RABVPrime.py --evaluate --fasta rabv_ea.reference.fasta --primers rabv_ea.primer.bed

  Use case 6 (evaluate only, DB + own FASTA, by country and subclade together):
    RABVPrime.py --evaluate --country Morocco --subclade AF1a --fasta rabv_ea.reference.fasta --primers rabv_ea.primer.bed --database rabv_V-gDB_10042026.db
""",
    )

    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument("--design-and-evaluate", action="store_true",
                             help="Design new primers, then evaluate them against the same sequences")
    mode_group.add_argument("--evaluate", action="store_true",
                             help="Evaluate an existing/external primer set only — no design step")

    parser.add_argument("-c", "--country", nargs="+", default=None,
                         help="One or more country names to retrieve from RABV-DB (e.g. 'Morocco'). "
                              "Can be combined with --subclade.")
    parser.add_argument("-s", "--subclade", nargs="+", default=None,
                         help="One or more subclade names to retrieve from RABV-DB (e.g. 'AF1a'). "
                              "Can be combined with --country.")
    parser.add_argument("--fasta", default=None, help="Your own FASTA file of sequences")
    parser.add_argument("--primers", default=None, help="Existing primer.bed to check (required with --evaluate)")
    parser.add_argument("-db", "--database", default=None, help="Path to the RABV-DB SQLite database file")
    parser.add_argument("-o", "--output_dir", default="results", help="Output directory (default: results/)")

    args = parser.parse_args()
    mode = "design" if args.design_and_evaluate else "evaluate"

    if not args.country and not args.subclade and not args.fasta:
        parser.error("at least one of --country, --subclade, or --fasta is required")
    if mode == "design" and args.primers:
        parser.error("--primers is only used with --evaluate — --design-and-evaluate generates its own primers")
    if mode == "evaluate" and not args.primers:
        parser.error("--evaluate requires --primers <primer.bed>")

    has_db_filter = bool(args.country) or bool(args.subclade)
    use_case = USE_CASES[(mode, has_db_filter, bool(args.fasta))]

    # use case 5 (evaluate, own FASTA only) is the one pipeline that never touches RABV-DB
    if use_case != 5 and not args.database:
        parser.error(f"--database is required for use case {use_case}")

    print(f"→ Use case {use_case}", flush=True)
    Pipeline(use_case, args).run()

if __name__ == "__main__":
    main()
