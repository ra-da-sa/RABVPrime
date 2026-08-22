#!/usr/bin/env python3
import argparse
import datetime
import os
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path

# Same convention as RetrieveSequences.py's MIN_GENOME_LENGTH: 90% of the full genome (11932 bp).
# Submissions shorter than this are partial sequences and are excluded before alignment.
MIN_SEQUENCE_LENGTH = 10739

class SubmitSequences:

    def __init__(self, input_fasta, db_path, output_dir="results", min_length=MIN_SEQUENCE_LENGTH,
                 rebuild_references=False, country=None, subclade=None):
        self.input_fasta             = input_fasta
        self.db_path                 = db_path
        self.output_dir              = output_dir
        self.min_length              = min_length
        self.rebuild_references = rebuild_references
        # user-asserted labels for this submission (not validated/filtered against the database —
        # SubmitSequences.py never filters by country/subclade, this is purely what the user says
        # this batch represents), shown in the console output and pipeline log as-is
        self.country                 = country
        self.subclade                = subclade
        self.conn                    = None

    def connect_db(self):
        if not os.path.exists(self.db_path):
            raise FileNotFoundError(f"Database file not found: {self.db_path}")
        self.conn = sqlite3.connect(self.db_path)

    @staticmethod
    def _read_fasta_records(lines):
        records = []
        header, seq_parts = None, []
        for line in lines:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header, seq_parts = line[1:].strip(), []
            elif line:
                seq_parts.append(line)
        if header is not None:
            records.append((header, "".join(seq_parts)))
        return records

    @staticmethod
    def sanitize_id(header):
        accession = header.split()[0]
        return "".join(c if c.isalnum() or c in "._-" else "_" for c in accession)

    @staticmethod
    def sanitize(text):
        return "_".join(part for part in text.strip().lower().replace("/", "_").split() if part)

    # returns (usable records, sequences skipped for being too short, sequences skipped for being empty);
    # also de-duplicates sanitized IDs so two differently-formatted headers never collide on disk.
    # Duplicate sequence content within the submission is intentionally NOT removed here — QC.py's
    # dedup_sequences() (seqkit rmdup -s) already does that, against the combined/QC'd set, so
    # doing it here too would just be redundant.
    def load_submitted_sequences(self):
        with open(self.input_fasta, encoding="utf-8") as f:
            raw_records = self._read_fasta_records(f)
        if not raw_records:
            print(f"✗ No sequences found in: {self.input_fasta}")
            sys.exit(1)

        records, seen, skipped_short, skipped_empty = [], set(), [], []
        for header, seq in raw_records:
            seq_id   = self.sanitize_id(header)
            ungapped = seq.replace("-", "").replace(".", "").upper()
            if not ungapped:
                skipped_empty.append(seq_id)
                continue
            if len(ungapped) < self.min_length:
                skipped_short.append((seq_id, len(ungapped)))
                continue

            unique_id, n = seq_id, 2
            while unique_id in seen:
                unique_id = f"{seq_id}_{n}"
                n += 1
            seen.add(unique_id)
            records.append((unique_id, ungapped))
        return records, skipped_short, skipped_empty

    # the reference panel is every accession meta_data marks as 'reference' or 'master' — together
    # they span the full clade diversity (Cosmopolitan, Africa, Asian, Arctic, Bats, RAC-SK, ...), so
    # aligning against all of them (rather than one single reference) gives a much better anchor for
    # a submitted sequence of unknown/divergent lineage
    def fetch_reference_records(self):
        cursor = self.conn.cursor()
        cursor.execute(
            "SELECT s.header, s.sequence"
            " FROM meta_data m JOIN sequences s ON s.header = m.locus"
            " WHERE m.accession_type IN ('reference', 'master')"
            " ORDER BY m.accession_type, m.locus"
        )
        records = cursor.fetchall()
        if not records:
            print("✗ No accessions marked accession_type='reference'/'master' were found in the database.")
            sys.exit(1)
        return records

    def check_dependencies(self):
        return [tool for tool in ("mafft",) if subprocess.run(["which", tool], capture_output=True).returncode != 0]

    # builds the multi-sequence alignment of the reference panel (expensive — tens of minutes for
    # hundreds of near-full-genome sequences) and caches it next to the database so later runs, or
    # runs against a different input FASTA, reuse it instead of rebuilding from scratch
    def reference_panel_cache_path(self):
        db_stem = Path(self.db_path).stem
        return Path(self.db_path).parent / f"{db_stem}_reference_panel.fasta"

    def build_reference_panel(self, cache_path):
        reference_records = self.fetch_reference_records()
        print(f"Building reference panel alignment from {len(reference_records)} reference/master accessions "
              "(this can take several minutes, but only needs to happen once)...")

        with tempfile.TemporaryDirectory() as tmp:
            panel_input = os.path.join(tmp, "panel.fasta")
            with open(panel_input, "w", encoding="utf-8") as f:
                for header, seq in reference_records:
                    f.write(f">{header}\n{seq}\n")

            result = subprocess.run(
                ["mafft", "--retree", "1", "--thread", "-1", panel_input],
                capture_output=True, text=True,
            )
            if result.returncode != 0:
                print(f"✗ MAFFT failed while building the reference panel alignment: {result.stderr.strip()}")
                sys.exit(1)

            with open(cache_path, "w", encoding="utf-8") as f:
                f.write(result.stdout)

    def resolve_reference_panel(self):
        cache_path = self.reference_panel_cache_path()
        if self.rebuild_references or not cache_path.exists():
            self.build_reference_panel(cache_path)
        else:
            print(f"Using cached reference panel alignment: {cache_path}")
        panel_records = self._read_fasta_records(cache_path.read_text(encoding="utf-8").splitlines())
        return str(cache_path), len(panel_records), len(panel_records[0][1]) if panel_records else 0

    # aligns submitted sequences against the reference panel with `mafft --add --keeplength`, which
    # keeps the output length pinned to the panel's existing coordinate system (any insertions relative
    # to it are dropped, not shifted in), so results line up column-for-column with the reference panel
    def run_mafft(self, panel_path, records, aligned_path):
        with tempfile.TemporaryDirectory() as tmp:
            new_path = os.path.join(tmp, "submitted.fasta")
            with open(new_path, "w", encoding="utf-8") as f:
                for seq_id, seq in records:
                    f.write(f">{seq_id}\n{seq}\n")

            result = subprocess.run(
                ["mafft", "--add", new_path, "--keeplength", "--thread", "-1", panel_path],
                capture_output=True, text=True,
            )
            if result.returncode != 0:
                print(f"✗ MAFFT failed: {result.stderr.strip()}")
                sys.exit(1)

            aligned = self._read_fasta_records(result.stdout.splitlines())

        submitted_ids = {seq_id for seq_id, _ in records}
        aligned_lookup = {header: seq for header, seq in aligned if header in submitted_ids}
        with open(aligned_path, "w", encoding="utf-8") as f:
            for seq_id, _ in records:
                f.write(f">{seq_id}\n{aligned_lookup[seq_id]}\n")

    # pulls every existing sequence from the database once, so the exact-duplicate check
    # (raw string equality) doesn't have to query the database per-submission
    def fetch_all_db_sequences(self):
        print(f"Pulling existing sequences from the database ({os.path.abspath(self.db_path)}) "
              "to check whether your submissions already exist...")
        cursor = self.conn.cursor()
        cursor.execute("SELECT header, sequence FROM sequences")
        return cursor.fetchall()

    # exact check: does a submitted sequence's content (not just its accession/header) already
    # exist verbatim under some accession in the database? Comparison ignores case and gap characters.
    def find_exact_duplicates(self, records, db_records):
        by_sequence = {}
        for header, seq in db_records:
            normalized = seq.replace("-", "").replace(".", "").upper()
            by_sequence.setdefault(normalized, []).append(header)

        duplicates = {}
        for seq_id, seq in records:
            hits = by_sequence.get(seq)
            if hits:
                duplicates[seq_id] = hits
        return duplicates

    # returns the Bats/RAC-SK clade (or "bat" host) label if any of the given accessions is
    # associated with one, else None — used to decide whether an exact-duplicate submission
    # must be excluded from further analysis
    def clade_label_if_bat_or_rac(self, accessions):
        cursor = self.conn.cursor()
        for accession in accessions:
            cursor.execute("SELECT EPA_major_clade, host FROM meta_data WHERE locus = ?", (accession,))
            row = cursor.fetchone()
            if not row:
                continue
            major_clade, host = row
            if major_clade in ("Bats", "RAC-SK") or (host and "bat" in host.lower()):
                return major_clade or host
        return None

    # A submission that's an exact duplicate of an existing database sequence already has real
    # Country/Subclade/Year/Host on file under that accession — pull it in rather than leaving it
    # blank. Genuinely new submissions still get blank/None metadata (the DB has nothing to offer,
    # and the user supplied no metadata of their own). Row shape matches QualityControl.py's
    # build_metadata (seq_id, country, subclade, year, host) so the same summary-table logic applies.
    def build_metadata_rows(self, records, exact_duplicates):
        cursor = self.conn.cursor()
        rows = []
        for seq_id, _ in records:
            country = subclade = year = host = None
            hits = exact_duplicates.get(seq_id)
            if hits:
                cursor.execute(
                    "SELECT country, EPA_minor_clade, collection_year, host FROM meta_data WHERE locus = ?",
                    (hits[0],),
                )
                row = cursor.fetchone()
                if row:
                    country, subclade, year, host = row
            rows.append((seq_id, country, subclade, year, host))
        return rows

    # Columns match QualityControl.py's expected schema so its output can be QC'd unchanged.
    def write_metadata(self, rows, tsv_path):
        with open(tsv_path, "w", encoding="utf-8") as f:
            f.write("Accession\tCountry\tSubclade\tYear\tHost\n")
            for seq_id, country, subclade, year, host in rows:
                f.write(f"{seq_id}\t{country or ''}\t{subclade or ''}\t{year or ''}\t{host or ''}\n")
        print(f"✓ Saved metadata to: {tsv_path}")

    def write_summary(self, start_time, command, initial_count, skipped_short, skipped_empty,
                       exact_duplicates, excluded_bat_rac, records, rows,
                       aligned_path, tsv_path, summary_path):
        total = len(records)

        def col_w(d, header):
            return max(max((len(k) for k in d), default=0), len(header)) + 2

        def write_table(f, title, counts, col_label, sort_key):
            f.write(f"\n{title}\n")
            w = col_w(counts, col_label)
            f.write(f"{col_label:<{w}}{'Count':>8}{'%':>8}\n")
            f.write("-" * (w + 16) + "\n")
            for k, v in sorted(counts.items(), key=sort_key):
                f.write(f"{k:<{w}}{v:>8}{v/total*100:>7.1f}%\n")

        with open(summary_path, "w", encoding="utf-8") as f:
            f.write(f"Time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Command: {command}\n")

            f.write(f"\nDetected {initial_count} submitted sequences.\n")

            f.write(f"\nFiltering sequences with length below the threshold ({self.min_length:,} bp)...\n")
            if skipped_short:
                acc_w = max((len(seq_id) for seq_id, _ in skipped_short), default=9)
                f.write(f"{'':<{acc_w}}{'Length (bp)':>14}\n")
                for seq_id, length in skipped_short:
                    f.write(f"{seq_id:<{acc_w}}{length:>14}\n")
            f.write(f"Removed {len(skipped_short)} sequences with length below the threshold.\n")

            f.write("\nComparing submitted sequences against the database...\n")
            if exact_duplicates:
                for seq_id, hits in exact_duplicates.items():
                    note = f" — Bats/RAC-SK ({excluded_bat_rac[seq_id]}), excluded from further analysis" \
                        if seq_id in excluded_bat_rac else ""
                    f.write(f"  {seq_id} is identical to {', '.join(hits)}{note}\n")
                remaining_duplicate_count = len(exact_duplicates) - len(excluded_bat_rac)
                f.write(f"INFO: {remaining_duplicate_count} submitted sequence(s) are already present in the "
                        f"database. Existing metadata will be retrieved and written to the output metadata file.\n")
                if excluded_bat_rac:
                    f.write(f"INFO: {len(excluded_bat_rac)} of these are identical to a Bats/RAC-SK database "
                            f"sequence and are excluded from further analysis.\n")
            else:
                f.write("INFO: No submitted sequence is an exact duplicate of an existing database sequence.\n")

            f.write("\nAligning retained sequences with reference sequences using MAFFT...\n")

            f.write(f"\nTotal sequences removed: {initial_count - total}\n")

            f.write(f"\nTotal sequences : {total}\n")
            if records:
                lengths = [len(seq) for _, seq in records]
                f.write(f"Sequence length range: {min(lengths)}-{max(lengths)} bp\n")
                f.write(f"Average sequence length: {sum(lengths) / len(lengths):.0f} bp\n")

                country_counts, subclade_counts, host_counts, year_counts = {}, {}, {}, {}
                for _, country, subclade, year, host in rows:
                    country_counts[country or "None"]  = country_counts.get(country or "None", 0) + 1
                    subclade_counts[subclade or "None"] = subclade_counts.get(subclade or "None", 0) + 1
                    host_counts[host or "None"]         = host_counts.get(host or "None", 0) + 1
                    year_counts[year or "None"]         = year_counts.get(year or "None", 0) + 1

                write_table(f, "Country Summary", country_counts, "Country", lambda x: -x[1])
                write_table(f, "Subclade Summary", subclade_counts, "Subclade", lambda x: -x[1])

                years = [int(y) for y in year_counts if y.isdigit()]
                f.write("\nTemporal Summary\n")
                f.write(f"Year range: {min(years)}-{max(years)}\n" if years else "Year range: -\n")
                w = col_w(year_counts, "Year")
                f.write(f"{'Year':<{w}}{'Count':>8}{'%':>8}\n")
                f.write("-" * (w + 16) + "\n")
                for k, v in sorted(year_counts.items(), key=lambda x: (x[0] == "None", x[0])):
                    f.write(f"{k:<{w}}{v:>8}{v/total*100:>7.1f}%\n")

                write_table(f, "Host Summary", host_counts, "Host", lambda x: -x[1])

            f.write("\nResults written to:\n")
            f.write(f"\t✓ Saved {total} sequences to: {aligned_path}\n")
            f.write(f"\t✓ Saved metadata to: {tsv_path}\n")
            f.write(f"\t✓ Saved summary to: {summary_path}\n")

        print(f"✓ Saved summary to: {summary_path}")

    def run(self):
        start_time = datetime.datetime.now()
        command    = " ".join([os.path.basename(sys.argv[0])] + sys.argv[1:])

        if not os.path.exists(self.input_fasta):
            print(f"✗ Input file not found: {self.input_fasta}")
            sys.exit(1)

        missing = self.check_dependencies()
        if missing:
            print(f"✗ Missing dependency: {', '.join(missing)} (brew install mafft)")
            sys.exit(1)

        try:
            self.connect_db()
        except FileNotFoundError as error:
            print(error)
            sys.exit(1)

        print(f"Using database: {os.path.abspath(self.db_path)}\n")

        records, skipped_short, skipped_empty = self.load_submitted_sequences()
        initial_count = len(records) + len(skipped_short) + len(skipped_empty)
        if not records:
            print("✗ No sequences remained after length filtering.")
            sys.exit(0)

        panel_path, panel_size, panel_length = self.resolve_reference_panel()
        print(f"Reference panel: {panel_size} accessions, {panel_length} aligned columns")

        print()
        db_records = self.fetch_all_db_sequences()

        exact_duplicates = self.find_exact_duplicates(records, db_records)
        excluded_bat_rac = {}
        if exact_duplicates:
            print(f"Note: {len(exact_duplicates)} submitted sequence(s) already exist in the database "
                  "(identical sequence content, not just a similar accession):")
            for seq_id, hits in exact_duplicates.items():
                label = self.clade_label_if_bat_or_rac(hits)
                if label:
                    excluded_bat_rac[seq_id] = label
                    print(f"  ✗ {seq_id} is identical to existing accession(s) {', '.join(hits)} "
                          f"({label}) — contains Bat/RAC-SK sequence content and will be excluded "
                          "from further analysis.")
                else:
                    print(f"  ✓ {seq_id} is identical to existing accession(s): {', '.join(hits)}")
        else:
            print("No submitted sequence is an exact duplicate of an existing database sequence.")

        if excluded_bat_rac:
            records = [(seq_id, seq) for seq_id, seq in records if seq_id not in excluded_bat_rac]
            if not records:
                print("\n✗ All submitted sequences were excluded (identical to a Bat/RAC-SK database sequence).")
                sys.exit(0)

        # No <4-sequences guard here — QualityControl.py (the next pipeline step in every use
        # case) already checks the post-QC count against Treemmer's minimum, so checking again
        # here would just be a redundant, earlier-but-less-accurate copy of the same check.

        print(f"\nAligning {len(records)} submitted sequence(s) to the reference panel using MAFFT...")

        stem        = self.sanitize(Path(self.input_fasta).stem)
        folder_name = f"{stem}_sequences"
        run_dir     = os.path.join(self.output_dir, folder_name)
        os.makedirs(run_dir, exist_ok=True)

        aligned_path = os.path.join(run_dir, f"{folder_name}_sequences.fasta")
        tsv_path     = os.path.join(run_dir, f"{folder_name}_metadata.tsv")
        summary_path = os.path.join(run_dir, f"{folder_name}_summary.txt")

        self.run_mafft(panel_path, records, aligned_path)
        print(f"✓ Saved {len(records)} aligned sequences to: {aligned_path}")

        rows = self.build_metadata_rows(records, exact_duplicates)
        self.write_metadata(rows, tsv_path)
        self.write_summary(start_time, command, initial_count, skipped_short, skipped_empty,
                            exact_duplicates, excluded_bat_rac, records, rows,
                            aligned_path, tsv_path, summary_path)

        print(f"\nNext step: python scripts/QualityControl.py -i {aligned_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Align sequences submitted by the user in a FASTA file to the reference "
                     "sequences in the database (accession_type IN ('reference', 'master')) using MAFFT."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="FASTA file of sequences submitted by the user",
    )
    parser.add_argument(
        "-db", "--database",
        required=True,
        help="Path to the SQLite database file",
    )
    parser.add_argument(
        "-o", "--output_dir",
        default="results",
        help="Directory where output files will be saved (default: results/)",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=MIN_SEQUENCE_LENGTH,
        help=f"Minimum ungapped sequence length (bp) used to determine whether a sequence is "
             f"included (default: {MIN_SEQUENCE_LENGTH} bp)",
    )
    parser.add_argument(
        "--rebuild-references",
        action="store_true",
        help="Rebuild the cached reference panel alignment from the database instead of reusing "
             "the cached copy alongside the database file (needed after the database is updated)",
    )
    parser.add_argument(
        "--country", default=None,
        help="Country label for the submitted sequences (e.g. Kenya)",
    )
    parser.add_argument(
        "--subclade", default=None,
        help="Subclade label for the submitted sequences (e.g. AF1b)",
    )
    args = parser.parse_args()

    SubmitSequences(args.input, args.database, args.output_dir, args.min_length, args.rebuild_references,
                     args.country, args.subclade).run()
