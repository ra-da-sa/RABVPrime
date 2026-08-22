#!/usr/bin/env python3
import argparse
import datetime
import os
import shutil
import subprocess
import sys
from pathlib import Path

class QualityControl:

    def __init__(self, fasta_path, output_dir, max_n_percent=20.0, max_gap_percent=20.0):
        self.fasta_path           = fasta_path
        self.output_dir           = output_dir
        self.max_n_percent        = max_n_percent
        self.max_gap_percent      = max_gap_percent

    def read_accessions(self):
        accessions = []
        with open(self.fasta_path, encoding="utf-8") as f:
            for line in f:
                if line.startswith(">"):
                    accessions.append(line[1:].strip().split()[0])
        return accessions

    def count_sequences(self):
        with open(self.fasta_path, encoding="utf-8") as f:
            return sum(1 for line in f if line.startswith(">"))

    def trim_alignment(self):
        if subprocess.run(['which', 'trimal'], capture_output=True).returncode != 0:
            print("Warning: trimal not found — skipping gap trimming (brew install trimal).")
            return "trimal not found — gap trimming skipped."

        tmp_path = self.fasta_path + ".tmp"
        result = subprocess.run(
            ['trimal', '-in', self.fasta_path, '-out', tmp_path, '-noallgaps'],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            print(f"Warning: trimal failed — keeping untrimmed alignment ({result.stderr.strip()}).")
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
            return f"trimal failed — keeping untrimmed alignment ({result.stderr.strip()})."

        os.replace(tmp_path, self.fasta_path)
        return None

    # returns (per-removal transcript lines, warning-or-None); rewrites self.fasta_path in place
    def dedup_sequences(self):
        if subprocess.run(['which', 'seqkit'], capture_output=True).returncode != 0:
            print("Warning: seqkit not found — skipping duplicate sequence removal (brew install seqkit).")
            return [], "seqkit not found — duplicate sequence removal skipped."

        tmp_path     = self.fasta_path + ".tmp"
        dupinfo_path = self.fasta_path + ".dupinfo.tmp"
        result = subprocess.run(
            ['seqkit', 'rmdup', '-s', self.fasta_path, '-o', tmp_path, '-D', dupinfo_path],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            print(f"Warning: seqkit rmdup failed — keeping all sequences ({result.stderr.strip()}).")
            for path in (tmp_path, dupinfo_path):
                if os.path.exists(path):
                    os.remove(path)
            return [], f"seqkit rmdup failed — keeping all sequences ({result.stderr.strip()})."

        lines = []
        if os.path.exists(dupinfo_path):
            with open(dupinfo_path, encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    _, id_list = line.split("\t", 1)
                    ids = [i.strip() for i in id_list.split(",")]
                    retained, removed = ids[0], ids[1:]
                    for accession in removed:
                        lines.append(
                            f"{accession} is identical to {retained}; retaining {retained} for subsequent analysis."
                        )
            os.remove(dupinfo_path)

        os.replace(tmp_path, self.fasta_path)
        return lines, None

    def _read_fasta_records(self):
        records = []
        header, seq_parts = None, []
        with open(self.fasta_path, encoding="utf-8") as f:
            for line in f:
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

    # calculation only, same shape as calculate_gap_content() below — does not remove anything.
    # Removal happens separately, in remove_high_n_content_sequences(), using this result.
    # returns a dict: table (accession, n_pct, status) sorted worst-first, plus min/max/avg
    def calculate_n_content(self):
        records = self._read_fasta_records()

        table = []
        for header, seq in records:
            accession = header.split()[0]
            ungapped  = seq.replace("-", "").replace(".", "")
            n_count   = ungapped.upper().count("N")
            n_pct     = (n_count / len(ungapped) * 100) if ungapped else 100.0
            status    = "passed" if n_pct <= self.max_n_percent else "failed"
            table.append((accession, n_pct, status))

        pcts = [pct for _, pct, _ in table]
        return {
            "table": sorted(table, key=lambda x: -x[1]),
            "min":   min(pcts) if pcts else 0.0,
            "max":   max(pcts) if pcts else 0.0,
            "avg":   sum(pcts) / len(pcts) if pcts else 0.0,
        }

    # removes sequences calculate_n_content() marked "failed", rewriting fasta_path in place —
    # the one QC check that actually deletes sequences rather than just flagging them.
    # returns (excluded_count, retained_count)
    def remove_high_n_content_sequences(self, n_result):
        failed = {accession for accession, _, status in n_result["table"] if status == "failed"}

        records = self._read_fasta_records()
        kept = [(header, seq) for header, seq in records if header.split()[0] not in failed]

        with open(self.fasta_path, "w", encoding="utf-8") as f:
            for header, seq in kept:
                f.write(f">{header}\n{seq}\n")

        return len(records) - len(kept), len(kept)

    # informational only — flags high-gap sequences but does not remove them, unlike
    # remove_high_n_content_sequences
    def calculate_gap_content(self):
        records = self._read_fasta_records()

        table = []
        for header, seq in records:
            accession = header.split()[0]
            length    = len(seq)
            gap_count = seq.count("-") + seq.count(".")
            gap_pct   = (gap_count / length * 100) if length else 0.0
            status    = "passed" if gap_pct <= self.max_gap_percent else "flagged"
            table.append((accession, gap_pct, status))

        pcts    = [pct for _, pct, _ in table]
        flagged = sum(1 for _, _, status in table if status == "flagged")
        return {
            "table":   sorted(table, key=lambda x: -x[1]),
            "flagged": flagged,
            "min":     min(pcts) if pcts else 0.0,
            "max":     max(pcts) if pcts else 0.0,
            "avg":     sum(pcts) / len(pcts) if pcts else 0.0,
        }

    def compute_length_stats(self):
        lengths = [len(seq.replace("-", "").replace(".", "")) for _, seq in self._read_fasta_records()]
        if not lengths:
            return {"min": 0, "max": 0, "avg": 0.0}
        return {"min": min(lengths), "max": max(lengths), "avg": sum(lengths) / len(lengths)}

    # reads RetrieveSequences.py's metadata.tsv, keyed by accession, header-driven (case-insensitive)
    def read_source_metadata(self, tsv_path):
        if not os.path.exists(tsv_path):
            raise FileNotFoundError(
                f"Source metadata file not found: {tsv_path}\n"
                "Expected it alongside the input FASTA, written by RetrieveSequences.py."
            )

        lookup = {}
        with open(tsv_path, encoding="utf-8") as f:
            header = f.readline().rstrip("\n").split("\t")
            header_lower = [h.lower() for h in header]
            try:
                id_idx       = header_lower.index("accession")
                country_idx  = header_lower.index("country")
                subclade_idx = header_lower.index("subclade")
            except ValueError:
                raise ValueError(
                    f"Source metadata file is missing required columns "
                    f"(accession, country, subclade): {tsv_path}"
                )
            year_idx = header_lower.index("year") if "year" in header_lower else None
            host_idx = header_lower.index("host") if "host" in header_lower else None

            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(id_idx, country_idx, subclade_idx):
                    continue
                accession = parts[id_idx]
                country   = parts[country_idx] or None
                subclade  = parts[subclade_idx] or None
                host      = parts[host_idx] if host_idx is not None and len(parts) > host_idx else None
                year      = parts[year_idx] if year_idx is not None and len(parts) > year_idx else None
                lookup[accession] = (country, subclade, host or None, year or None)
        return lookup

    # row tuple: (accession, country, subclade, host, collection_year)
    def build_metadata(self, accessions, source_lookup):
        rows = []
        for accession in accessions:
            country, subclade, host, year = source_lookup.get(accession, (None, None, None, None))
            rows.append((accession, country, subclade, host, year))
        return rows

    def write_metadata(self, rows, tsv_path):
        with open(tsv_path, "w", encoding="utf-8") as f:
            f.write("Accession\tCountry\tSubclade\tYear\tHost\n")
            for accession, country, subclade, host, year in rows:
                f.write(f"{accession}\t{country or ''}\t{subclade or ''}\t{year or ''}\t{host or ''}\n")
        print(f"✓ Saved metadata to: {tsv_path}")

    def summarize_label(self, values):
        distinct = set(values)
        return next(iter(distinct)) if len(distinct) == 1 else "Multiple"

    @staticmethod
    def sanitize(text):
        return "_".join(
            part for part in text.strip().lower().replace("/", "_").split() if part
        )

    # Data content can be homogeneous (or diverse) by coincidence (e.g. every Moroccan record
    # happens to be subclade AF1a, or a user-submitted FASTA happens to span many countries)
    # regardless of whether that dimension was ever actually filtered on. The folder/file name
    # (built from sanitize(country)[_sanitize(subclade)]) is the reliable signal for which
    # dimension(s) were actually filtered on — if it doesn't match either one (e.g. a
    # SubmitSequences.py run, which never filters by country/subclade at all), report "None"
    # for both rather than surfacing what's incidentally true of the data.
    def resolve_filter_labels(self, file_prefix, country_label, subclade_label):
        country_token  = self.sanitize(country_label) if country_label not in (None, "None", "Multiple") else None
        subclade_token = self.sanitize(subclade_label) if subclade_label not in (None, "None", "Multiple") else None

        if country_token and subclade_token and file_prefix == f"{country_token}_{subclade_token}":
            return country_label, subclade_label
        if country_token and file_prefix == country_token:
            return country_label, "None"
        if subclade_token and file_prefix == subclade_token:
            return "None", subclade_label
        return "None", "None"

    def write_summary(self, start_time, command, trim_warning, dedup_lines, dedup_warning, dedup_retained,
                       gap_result, n_result, length_stats, records, initial_count, fasta_path, tsv_path, summary_path):
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

            f.write("\nRemoving duplicate sequences...\n")
            if dedup_warning:
                f.write(f"{dedup_warning}\n")
            else:
                for line in dedup_lines:
                    f.write(f"{line}\n")
                f.write(f"Removed {len(dedup_lines)} sequences; retaining {dedup_retained} unique sequences.\n")

            f.write("\nTrimming gaps from the alignment...\n")
            if trim_warning:
                f.write(f"{trim_warning}\n")

            f.write("\nCalculating gap content...\n")
            gap_acc_w = max((len(acc) for acc, _, _ in gap_result["table"]), default=9)
            f.write(f"{'':<{gap_acc_w}}{'Gap content':>14}  {'Status':<8}\n")
            for accession, pct, status in gap_result["table"]:
                f.write(f"{accession:<{gap_acc_w}}{pct:>13.2f}%  {status}\n")
            f.write(f"Minimum gap content: {gap_result['min']:.2f}%\n")
            f.write(f"Maximum gap content: {gap_result['max']:.2f}%\n")
            f.write(f"Average gap content: {gap_result['avg']:.2f}%\n")
            if gap_result["flagged"]:
                f.write(f"WARNING: {gap_result['flagged']} sequences have high gap content (>{self.max_gap_percent:.0f}%).\n")
            else:
                f.write(f"INFO: No sequences have high gap content (>{self.max_gap_percent:.0f}%).\n")

            f.write(f"\nRemoving sequences with N-content above the threshold ({self.max_n_percent:.0f}%)...\n")
            acc_w = max((len(acc) for acc, _, _ in n_result["table"]), default=9)
            f.write(f"{'':<{acc_w}}{'N-content':>12}  {'Status':<6}\n")
            for accession, pct, status in n_result["table"]:
                f.write(f"{accession:<{acc_w}}{pct:>11.2f}%  {status}\n")
            f.write(f"Minimum N-content: {n_result['min']:.2f}%\n")
            f.write(f"Maximum N-content: {n_result['max']:.2f}%\n")
            f.write(f"Average N-content: {n_result['avg']:.2f}%\n")
            f.write(f"Removed {n_result['excluded']} sequences with N-content above the threshold; "
                    f"retaining {n_result['retained']} sequences.\n")

            f.write(f"\nTotal sequences removed: {initial_count - total}\n")
            f.write(f"\nTotal sequences retained: {total}\n")
            f.write(f"Sequence length range: {length_stats['min']}-{length_stats['max']} bp\n")
            f.write(f"Average sequence length: {length_stats['avg']:.0f} bp\n")

            if records:
                country_counts, subclade_counts, host_counts, year_counts = {}, {}, {}, {}
                for r in records:
                    country_counts[r[1] or "None"]  = country_counts.get(r[1] or "None", 0) + 1
                    subclade_counts[r[2] or "None"] = subclade_counts.get(r[2] or "None", 0) + 1
                    host_counts[r[3] or "None"]     = host_counts.get(r[3] or "None", 0) + 1
                    year_counts[r[4] or "None"]     = year_counts.get(r[4] or "None", 0) + 1

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
            f.write(f"\t✓ Saved {total} sequences to: {fasta_path}\n")
            f.write(f"\t✓ Saved metadata to: {tsv_path}\n")
            f.write(f"\t✓ Saved summary to: {summary_path}\n")

        print(f"✓ Saved summary to: {summary_path}")

    def run(self):
        start_time = datetime.datetime.now()
        command    = " ".join([os.path.basename(sys.argv[0])] + sys.argv[1:])

        if not os.path.exists(self.fasta_path):
            print(f"✗ Input file not found: {self.fasta_path}")
            sys.exit(1)

        # the source folder's own name already encodes the file prefix RetrieveSequences.py used
        # (e.g. "morocco_sequences"), so its sibling files are named "morocco_sequences_<type>"
        source_folder_name = Path(self.fasta_path).parent.name
        file_prefix = source_folder_name[:-len("_sequences")] if source_folder_name.endswith("_sequences") else source_folder_name
        source_tsv_path = Path(self.fasta_path).parent / f"{source_folder_name}_metadata.tsv"

        try:
            source_lookup = self.read_source_metadata(str(source_tsv_path))
        except (FileNotFoundError, ValueError) as error:
            print(f"✗ {error}")
            sys.exit(1)

        country_label  = self.summarize_label((v[0] or "None") for v in source_lookup.values())
        subclade_label = self.summarize_label((v[1] or "None") for v in source_lookup.values())
        country_label, subclade_label = self.resolve_filter_labels(file_prefix, country_label, subclade_label)

        # QC gets its own folder — RetrieveSequences.py's output is never modified.
        # Files are named after their containing folder, matching RetrieveSequences.py's convention.
        folder_name = f"{file_prefix}_qc"
        qc_dir = Path(self.output_dir) / folder_name
        qc_dir.mkdir(parents=True, exist_ok=True)

        qc_fasta_path = qc_dir / f"{folder_name}_sequences.fasta"
        shutil.copy(self.fasta_path, qc_fasta_path)
        self.fasta_path = str(qc_fasta_path)
        initial_count = self.count_sequences()

        tsv_path     = qc_dir / f"{folder_name}_metadata.tsv"
        summary_path = qc_dir / f"{folder_name}_summary.txt"

        print("Running QC on retrieved sequences...")
        dedup_lines, dedup_warning = self.dedup_sequences()
        dedup_retained = self.count_sequences()
        trim_warning = self.trim_alignment()
        gap_result = self.calculate_gap_content()
        n_result = self.calculate_n_content()
        n_result["excluded"], n_result["retained"] = self.remove_high_n_content_sequences(n_result)
        length_stats = self.compute_length_stats()

        accessions = self.read_accessions()

        country_set  = country_label  not in (None, "None")
        subclade_set = subclade_label not in (None, "None")
        if country_set and subclade_set:
            filter_desc = f"{country_label} (country) and {subclade_label} (subclade)"
        elif country_set:
            filter_desc = f"{country_label} (country)"
        else:
            filter_desc = f"{subclade_label} (subclade)"

        if not accessions:
            print(f"\nNo sequences were retained for {filter_desc}.")
            print("Please select a different country and/or subclade.")
            shutil.rmtree(qc_dir, ignore_errors=True)
            sys.exit(0)

        rows = self.build_metadata(accessions, source_lookup)

        print(f"✓ Saved {len(accessions)} sequences to: {qc_fasta_path}")
        self.write_metadata(rows, str(tsv_path))
        self.write_summary(start_time, command, trim_warning, dedup_lines, dedup_warning, dedup_retained,
                            gap_result, n_result, length_stats, rows, initial_count, str(qc_fasta_path), str(tsv_path), str(summary_path))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run quality control (gap trimming, duplicate removal, N-content filtering) on retrieved sequences."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="FASTA file of sequences selected from a country and/or subclade, sequences submitted "
             "by the user, or a combination of both.",
    )
    parser.add_argument(
        "-o", "--output_dir",
        default="results",
        help="Directory where the QC output folder will be saved (default: results/)",
    )
    parser.add_argument(
        "-g", "--max-gap-content",
        type=float,
        default=20.0,
        help="Maximum gap content (%%) allowed per sequence (default: 20)",
    )
    parser.add_argument(
        "-n", "--max-n-content",
        type=float,
        default=20.0,
        help="Maximum N-content (%%) allowed per sequence; sequences with N-content above this "
             "excluded (default: 20)",
    )
    args = parser.parse_args()

    QualityControl(args.input, args.output_dir, args.max_n_content, args.max_gap_content).run()
