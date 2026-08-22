#!/usr/bin/env python3
import argparse
import datetime
import os
import sqlite3
import sys

from SubmitSequences import SubmitSequences

class SequenceExtractor:

    # subclade_input/country_input: None, or a non-empty list of user-typed names (argparse's
    # -s/--subclade and -c/--country both use nargs='+', so even a single value arrives as a
    # one-item list)
    def __init__(self, subclade_input, country_input, db_path, output_dir, rebuild_references=False):
        self.subclade_input     = subclade_input
        self.country_input      = country_input
        self.db_path            = db_path
        self.output_dir         = output_dir
        self.rebuild_references = rebuild_references
        self.conn = None

    def connect_db(self):
        if not os.path.exists(self.db_path):
            raise FileNotFoundError(f"Database file not found: {self.db_path}")
        self.conn = sqlite3.connect(self.db_path)

    def find_subclades(self, user_input):
        cursor = self.conn.cursor()
        cursor.execute(
            "SELECT DISTINCT EPA_minor_clade FROM meta_data WHERE EPA_minor_clade = ?",
            (user_input,),
        )
        rows = cursor.fetchall()
        if rows:
            return [r[0] for r in rows]
        cursor.execute(
            "SELECT DISTINCT EPA_minor_clade FROM meta_data WHERE lower(EPA_minor_clade) = ?",
            (user_input.lower(),),
        )
        rows = cursor.fetchall()
        if rows:
            return [r[0] for r in rows]
        cursor.execute(
            "SELECT DISTINCT EPA_minor_clade FROM meta_data WHERE lower(EPA_minor_clade) LIKE ?",
            (f"%{user_input.lower()}%",),
        )
        return [r[0] for r in cursor.fetchall()]

    def find_countries(self, user_input):
        cursor = self.conn.cursor()
        cursor.execute(
            "SELECT DISTINCT country FROM meta_data WHERE country = ?",
            (user_input,),
        )
        rows = cursor.fetchall()
        if rows:
            return [r[0] for r in rows]
        cursor.execute(
            "SELECT DISTINCT country FROM meta_data WHERE lower(country) = ?",
            (user_input.lower(),),
        )
        rows = cursor.fetchall()
        if rows:
            return [r[0] for r in rows]
        cursor.execute(
            "SELECT DISTINCT country FROM meta_data WHERE lower(country) LIKE ?",
            (f"%{user_input.lower()}%",),
        )
        return [r[0] for r in cursor.fetchall()]

    def choose_from_list(self, candidates, label):
        if not candidates:
            return None
        if len(candidates) == 1:
            return candidates[0]
        print(f"Multiple matching {label}s found:")
        for i, name in enumerate(candidates, 1):
            print(f"  {i}. {name}")
        while True:
            sel = input(f"Select a {label} by number: ").strip()
            if not sel.isdigit():
                print("Please enter a number.")
                continue
            sel = int(sel)
            if 1 <= sel <= len(candidates):
                return candidates[sel - 1]
            print(f"Please select a number between 1 and {len(candidates)}.")

    # resolves each requested subclade name to its canonical EPA_minor_clade value (prompting to
    # disambiguate any that match more than one), in the order given, de-duplicated in case two
    # inputs (or a repeated one) resolve to the same canonical name
    def resolve_subclades(self):
        resolved = []
        for user_input in self.subclade_input:
            candidates = self.find_subclades(user_input.strip())
            if not candidates:
                print(f"No subclade found matching '{user_input}'.")
                sys.exit(1)
            name = self.choose_from_list(candidates, "subclade")
            if name not in resolved:
                resolved.append(name)
        return resolved

    # resolves each requested country name to its canonical value (prompting to disambiguate any
    # that match more than one), in the order given, de-duplicated in case two inputs (or a
    # repeated one) resolve to the same canonical name
    def resolve_countries(self):
        resolved = []
        for user_input in self.country_input:
            candidates = self.find_countries(user_input.strip())
            if not candidates:
                print(f"No country found matching '{user_input}'.")
                sys.exit(1)
            name = self.choose_from_list(candidates, "country")
            if name not in resolved:
                resolved.append(name)
        return resolved

    # Minimum length = 90% of full genome (11932 bp)
    MIN_GENOME_LENGTH = 10739

    # record tuple: (locus, major_clade, minor_clade, country, host, collection_year, sequence).
    # Sourced from the raw "sequences" table (not the DB's precomputed "sequence_alignment"),
    # because that precomputed alignment was built externally with Nextalign and lands in a
    # different coordinate system (12005 cols) than the MAFFT-against-reference-panel alignment
    # SubmitSequences.py computes (12285 cols) — pulling raw sequence and aligning it ourselves
    # via the same MAFFT step (see run()) keeps retrieved and submitted sequences compatible.
    # subclades/countries: list of one or more canonical values, or falsy for "any"
    def query_sequences(self, subclades, countries):
        cursor = self.conn.cursor()
        if subclades and countries:
            subclade_placeholders = ",".join("?" for _ in subclades)
            country_placeholders = ",".join("?" for _ in countries)
            cursor.execute(
                "SELECT m.locus, m.EPA_major_clade, m.EPA_minor_clade, m.country, m.host, m.collection_year, s.sequence"
                " FROM meta_data m"
                " JOIN sequences s ON s.header = m.locus"
                f" WHERE m.EPA_minor_clade IN ({subclade_placeholders}) AND m.country IN ({country_placeholders})"
                " AND m.EPA_major_clade NOT IN ('Bats', 'RAC-SK')"
                " AND m.EPA_major_all NOT LIKE '%Bats%'"
                " AND m.EPA_major_all NOT LIKE '%RAC-SK%'"
                " AND CAST(m.length AS REAL) >= ?"
                " ORDER BY m.locus",
                (*subclades, *countries, self.MIN_GENOME_LENGTH),
            )
        elif subclades:
            placeholders = ",".join("?" for _ in subclades)
            cursor.execute(
                "SELECT m.locus, m.EPA_major_clade, m.EPA_minor_clade, m.country, m.host, m.collection_year, s.sequence"
                " FROM meta_data m"
                " JOIN sequences s ON s.header = m.locus"
                f" WHERE m.EPA_minor_clade IN ({placeholders})"
                " AND m.EPA_major_clade NOT IN ('Bats', 'RAC-SK')"
                " AND m.EPA_major_all NOT LIKE '%Bats%'"
                " AND m.EPA_major_all NOT LIKE '%RAC-SK%'"
                " AND CAST(m.length AS REAL) >= ?"
                " ORDER BY m.locus",
                (*subclades, self.MIN_GENOME_LENGTH),
            )
        else:
            placeholders = ",".join("?" for _ in countries)
            cursor.execute(
                "SELECT m.locus, m.EPA_major_clade, m.EPA_minor_clade, m.country, m.host, m.collection_year, s.sequence"
                " FROM meta_data m"
                " JOIN sequences s ON s.header = m.locus"
                f" WHERE m.country IN ({placeholders})"
                " AND m.EPA_major_clade NOT IN ('Bats', 'RAC-SK')"
                " AND m.EPA_major_all NOT LIKE '%Bats%'"
                " AND m.EPA_major_all NOT LIKE '%RAC-SK%'"
                " AND m.EPA_minor_clade IS NOT NULL"
                " AND CAST(m.length AS REAL) >= ?"
                " ORDER BY m.locus",
                (*countries, self.MIN_GENOME_LENGTH),
            )
        return cursor.fetchall()

    # 1. metadata file
    def write_tsv(self, records, tsv_path):
        with open(tsv_path, "w", encoding="utf-8") as f:
            f.write("Accession\tCountry\tSubclade\tYear\tHost\n")
            for r in records:
                accession, minor_clade, country, host, year = r[0], r[2], r[3], r[4], r[5]
                f.write(f"{accession}\t{country or ''}\t{minor_clade or ''}\t{year or ''}\t{host or ''}\n")
        print(f"✓ Saved metadata to: {tsv_path}")

    # 2. summary (run transcript + summary statistics + output file list)
    def write_summary(self, start_time, command, records, remark, fasta_path, tsv_path, summary_path):
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

            f.write(f"\nTotal sequences retrieved: {total}\n")

            if records:
                lengths = [len(r[-1].replace('-', '').replace('.', '')) for r in records]
                f.write(f"Sequence length range: {min(lengths)}-{max(lengths)} bp\n")
                f.write(f"Average sequence length: {sum(lengths)/len(lengths):.0f} bp\n")

            if remark:
                f.write(f"\n{remark}\n")

            if records:
                country_counts, subclade_counts, host_counts, year_counts = {}, {}, {}, {}
                for r in records:
                    country_counts[r[3] or "None"]  = country_counts.get(r[3] or "None", 0) + 1
                    subclade_counts[r[2] or "None"] = subclade_counts.get(r[2] or "None", 0) + 1
                    host_counts[r[4] or "None"]     = host_counts.get(r[4] or "None", 0) + 1
                    year_counts[r[5] or "None"]     = year_counts.get(r[5] or "None", 0) + 1

                write_table(f, "Country Summary", country_counts, "Country", lambda x: -x[1])
                write_table(f, "Subclade Summary", subclade_counts, "Subclade", lambda x: -x[1])

                years = [int(y) for y in year_counts if y.isdigit()]
                f.write("\nTemporal Summary\n")
                f.write(f"Year range: {min(years)}-{max(years)}\n" if years else "Year range: -\n")
                w = col_w(year_counts, "Year")
                f.write(f"{'Year':<{w}}{'Count':>8}\n")
                f.write("-" * (w + 8) + "\n")
                for k, v in sorted(year_counts.items(), key=lambda x: (x[0] == "None", x[0])):
                    f.write(f"{k:<{w}}{v:>8}\n")

                write_table(f, "Host Summary", host_counts, "Host", lambda x: -x[1])

            f.write("\nResults written to:\n")
            if fasta_path:
                f.write(f"\t✓ Saved {total} sequences to: {fasta_path}\n")
            if tsv_path:
                f.write(f"\t✓ Saved metadata to: {tsv_path}\n")
            f.write(f"\t✓ Saved summary to: {summary_path}\n")

        print(f"✓ Saved summary to: {summary_path}")

    @staticmethod
    def sanitize(text):
        return "_".join(
            part for part in text.strip().lower().replace("/", "_").split() if part
        )

    def run(self):
        start_time = datetime.datetime.now()
        command = " ".join([os.path.basename(sys.argv[0])] + sys.argv[1:])

        print(f"Using database: {os.path.abspath(self.db_path)}\n")

        # connect to the database
        try:
            self.connect_db()
        except FileNotFoundError as error:
            print(error)
            sys.exit(1)

        # find the country/countries and/or subclade(s) specificed by the user
        subclades = self.resolve_subclades() if self.subclade_input else []
        countries = self.resolve_countries() if self.country_input  else []
        # joined form used everywhere a single string is needed (folder stem, messages, the
        # pipeline log) — "+" rather than ", " so it stays filesystem-safe once sanitized
        subclade_label = "+".join(subclades) if subclades else None
        country_label = "+".join(countries) if countries else None

        print("Retrieving sequences from the database...")

        # retrieve the sequences from the country/countries and/or subclade(s) specified by the user
        records = self.query_sequences(subclades, countries)

        # output folder stem — the folder itself is only created on a successful run
        parts = []
        parts.extend(self.sanitize(country) for country in countries)
        parts.extend(self.sanitize(subclade) for subclade in subclades)
        stem = "_".join(parts)

        if not records:
            if countries and subclades:
                filter_desc = f"{country_label} (country) and {subclade_label} (subclade)"
            else:
                filter_label = "country" if countries else "subclade"
                filter_desc = f"{country_label or subclade_label} ({filter_label})"
            print(f"\nNo sequences were found for {filter_desc}.")
            print("Please select a different country and/or subclade.")
            sys.exit(0)

        # No <4-sequences guard here — QualityControl.py (the next pipeline step in every use
        # case) already checks the post-QC count against Treemmer's minimum, so checking again
        # here would just be a redundant, earlier-but-less-accurate copy of the same check (QC's
        # count reflects dedup/N-filtering too, this one doesn't).

        # align the raw retrieved sequences against the same reference panel SubmitSequences.py
        # uses, via the same MAFFT step, so retrieved and submitted sequences always share one
        # coordinate system (the database's own precomputed alignment does not — see query_sequences)
        submitter = SubmitSequences(None, self.db_path, self.output_dir,
                                     rebuild_references=self.rebuild_references)
        missing = submitter.check_dependencies()
        if missing:
            print(f"✗ Missing dependency: {', '.join(missing)} (brew install mafft)")
            sys.exit(1)
        submitter.connect_db()
        panel_path, panel_size, panel_length = submitter.resolve_reference_panel()
        print(f"Reference panel: {panel_size} accessions, {panel_length} aligned columns")

        # output — folder is only created once we know the run succeeded
        # files are named after their containing folder (e.g. morocco_sequences/morocco_sequences_metadata.tsv)
        # so they don't collide in basename with QualityControl.py's/Treemmer.py's same-named outputs
        folder_name = f"{stem}_sequences"
        run_dir = os.path.join(self.output_dir, folder_name)
        os.makedirs(run_dir, exist_ok=True)
        summary_path = os.path.join(run_dir, f"{folder_name}_summary.txt")
        fasta_path   = os.path.join(run_dir, f"{folder_name}_sequences.fasta")
        tsv_path     = os.path.join(run_dir, f"{folder_name}_metadata.tsv")

        print(f"\nAligning {len(records)} sequence(s) to the reference panel using MAFFT...")
        submitter.run_mafft(panel_path, [(r[0], r[-1]) for r in records], fasta_path)

        # re-read the aligned sequences so metadata/summary stats reflect what was actually
        # written (mafft --keeplength can drop insertions relative to the panel)
        with open(fasta_path, encoding="utf-8") as f:
            aligned_lookup = dict(SubmitSequences._read_fasta_records(f.read().splitlines()))
        records = [r[:-1] + (aligned_lookup[r[0]],) for r in records]

        # write metadata and summary files
        # (QualityControl.py writes its own copies into a separate <stem>_QC/ folder)
        print(f"✓ Saved {len(records)} sequences to: {fasta_path}")
        self.write_tsv(records, tsv_path)
        self.write_summary(start_time, command, records, None, fasta_path, tsv_path, summary_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Retrieve sequences from a selected country and/or subclade in the database. "
                     "Align the retrieved sequences to the reference panel (accession_type IN "
                     "('reference', 'master')) using MAFFT."
    )
    parser.add_argument(
        "-c", "--country",
        nargs="+",
        default=None,
        help="One or more country names, space-separated (e.g. 'China' or 'Laos Vietnam') — "
             "sequences matching any of them are retrieved and combined into one output",
    )
    parser.add_argument(
        "-s", "--subclade",
        nargs="+",
        default=None,
        help="One or more subclade names, space-separated (e.g. 'CA1' or 'SEA1a SEA1b') — "
             "sequences matching any of them are retrieved and combined into one output",
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
        "--rebuild-references",
        action="store_true",
        help="Rebuild the cached reference panel alignment from the database instead of reusing "
             "the cached copy alongside the database file (needed after the database is updated)",
    )
    args = parser.parse_args()

    if not args.subclade and not args.country:
        parser.error("At least one of -s/--subclade or -c/--country must be specified.")

    SequenceExtractor(args.subclade, args.country, args.database, args.output_dir,
                       args.rebuild_references).run()
