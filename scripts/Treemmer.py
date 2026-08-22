#!/usr/bin/env python3
import argparse
import datetime
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

try:
    import matplotlib
    matplotlib.use('Agg')
    from Bio import Phylo
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

class TreemmerWorkflow:
    # RTL 0.9 default - pipeline retains 90% of phylogenetic diversity unless specified otherwise
    def __init__(self, fasta_file, output_dir="results", keep_rtl=0.9):
        self.fasta_file = fasta_file
        self.output_dir = output_dir
        self.keep_rtl = keep_rtl

    def count_sequences(self):
        with open(self.fasta_file, encoding="utf-8") as f:
            return sum(1 for line in f if line.startswith(">"))

    def check_dependencies(self):
        tools = {
            'fasttree': 'FastTree (tree building)',
        }
        missing = []
        for tool, description in tools.items():
            result = subprocess.run(['which', tool], capture_output=True)
            if result.returncode != 0:
                missing.append(description)
        if missing:
            print("\nMissing dependencies:")
            for dep in missing:
                print(f"  - {dep}")
            print("\nInstall with: brew install fasttree")
            return False
        return True

    def build_tree(self, alignment_file, tree_file):
        cmd = ['fasttree', '-nt', alignment_file]
        with open(tree_file, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)
        if result.returncode != 0:
            print(f"✗ Tree building failed: {result.stderr.decode()}")
            return False
        return True

    def root_tree(self, tree_file):
        if not BIOPYTHON_AVAILABLE:
            print("⚠️  Biopython not available — skipping midpoint rooting.")
            return False
        tree = Phylo.read(tree_file, 'newick')
        tree.root_at_midpoint()
        Phylo.write(tree, tree_file, 'newick')
        print(f"✓ Saved midpoint-rooted tree to: {tree_file}")
        return True

    def plot_tree(self, tree_file, image_file):
        if not BIOPYTHON_AVAILABLE:
            print("⚠️  Biopython not installed. Skipping tree diagram.")
            return False
        tree = Phylo.read(tree_file, "newick")
        tree.ladderize()  # flip branches so deeper clades are displayed at top

        import matplotlib.pyplot as plt
        n_leaves = tree.count_terminals()
        fig = plt.figure(figsize=(12, max(6, n_leaves * 0.3)))
        axes = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, do_show=False, show_confidence=False, axes=axes)
        axes.get_yaxis().set_visible(False)

        plt.savefig(image_file, bbox_inches="tight")
        plt.close()
        return True

    def calculate_mpd(self, tree_path):
        if not BIOPYTHON_AVAILABLE:
            return None
        tree = Phylo.read(tree_path, "newick")
        branch_lengths = [c.branch_length for c in tree.find_clades() if c.branch_length is not None]
        total_branches = len(branch_lengths)
        total_length   = sum(branch_lengths)

        # Mean Pairwise Distance (MPD): average tree-distance between every pair of tips
        terminals = tree.get_terminals()
        distances = [
            tree.distance(terminals[i], terminals[j])
            for i in range(len(terminals))
            for j in range(i + 1, len(terminals))
        ]
        mean_pairwise_distance = (sum(distances) / len(distances)) if distances else 0.0

        return {
            "branches":     total_branches,
            "total_length": total_length,
            "mpd":          mean_pairwise_distance,
        }

    def compute_length_stats(self, fasta_path):
        lengths = []
        header, parts = None, []
        with open(fasta_path, encoding="utf-8") as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    if header is not None:
                        lengths.append(len("".join(parts).replace('-', '').replace('.', '')))
                    header, parts = line, []
                elif line:
                    parts.append(line)
            if header is not None:
                lengths.append(len("".join(parts).replace('-', '').replace('.', '')))
        if not lengths:
            return {"min": 0, "max": 0, "avg": 0.0}
        return {"min": min(lengths), "max": max(lengths), "avg": sum(lengths) / len(lengths)}

    def download_treemmer(self):
        treemmer_script = 'Treemmer_v0.3.py'
        if os.path.exists(treemmer_script):
            return True
        url = "https://raw.githubusercontent.com/fmenardo/Treemmer/master/Treemmer_v0.3.py"
        result = subprocess.run(['curl', '-s', '-o', treemmer_script, url], capture_output=True)
        if result.returncode == 0:
            os.chmod(treemmer_script, 0o755)
            return True
        print(f"✗ Failed to download Treemmer")
        return False

    # returns (success, rtl_table) where rtl_table is [(rtl_value_str, n_remaining), ...]
    # parsed straight from Treemmer_v0.3.py's stdout, since it doesn't reliably write the "_LD" file
    def run_treemmer(self, tree_file, pruning_option, pruning_value):
        treemmer_script = 'Treemmer_v0.3.py'
        cmd = [sys.executable, treemmer_script, tree_file, pruning_option, str(pruning_value)]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"✗ Treemmer failed: {result.stderr}")
            return False, []

        rtl_table = []
        for line in result.stdout.splitlines():
            m = re.search(r"RTL\s*:\s*([\d.]+)\s*N_seq:\s*(\d+)", line)
            if m:
                rtl_table.append((m.group(1), int(m.group(2))))
        return True, rtl_table

    def write_representative_sequences(self, accessions, alignment_file, output_fasta):
        representatives = set(accessions)

        written = 0
        with open(alignment_file, encoding="utf-8") as infile, \
             open(output_fasta, "w", encoding="utf-8") as outfile:
            write = False
            for line in infile:
                if line.startswith(">"):
                    accession = line[1:].strip().split(",")[0].strip()
                    write = accession in representatives
                if write:
                    outfile.write(line)
                    if line.startswith(">"):
                        written += 1

        return written

    def load_metadata(self, source_folder_name):
        meta_path = Path(self.fasta_file).parent / f"{source_folder_name}_metadata.tsv"
        lookup = {}
        if not meta_path.exists():
            return lookup
        with open(meta_path, encoding="utf-8") as f:
            header = f.readline().strip().split('\t')
            header_lower = [h.lower() for h in header]
            try:
                country_idx  = header_lower.index('country')
                subclade_idx = header_lower.index('subclade')
            except ValueError:
                return lookup
            year_idx = header_lower.index('year') if 'year' in header_lower else None
            host_idx = header_lower.index('host') if 'host' in header_lower else None
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) > max(country_idx, subclade_idx):
                    lookup[parts[0]] = {
                        'country':  parts[country_idx],
                        'subclade': parts[subclade_idx],
                        'year':     parts[year_idx] if year_idx is not None and len(parts) > year_idx else '-',
                        'host':     parts[host_idx] if host_idx is not None and len(parts) > host_idx else '-',
                    }
        return lookup

    def write_metadata_list(self, accessions, metadata, path):
        with open(path, "w", encoding="utf-8") as f:
            f.write("Accession\tCountry\tSubclade\tYear\tHost\n")
            for acc in accessions:
                info     = metadata.get(acc, {})
                country  = info.get('country', '-')
                subclade = info.get('subclade', '-')
                year     = info.get('year', '-')
                host     = info.get('host', '-')
                f.write(f"{acc}\t{country}\t{subclade}\t{year}\t{host}\n")

    def enrich_trimmed_list(self, trimmed_list_path, metadata):
        with open(trimmed_list_path, encoding="utf-8") as f:
            accessions = [line.strip() for line in f if line.strip()]
        self.write_metadata_list(accessions, metadata, trimmed_list_path)
        return accessions

    def write_summary(self, start_time, command, rtl_table, not_reached_note, original_count,
                       accessions, metadata, tree_stats, length_stats, output_files, summary_path):
        total = len(accessions)

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

            f.write("\nBuilding phylogenetic tree using FastTree...\n")
            f.write("\nApplying midpoint rooting...\n")

            f.write("\nPruning using Treemmer...\n")
            if rtl_table:
                rows = [(f"{float(rtl):.4f}", n_seq) for rtl, n_seq in rtl_table]
                rtl_w = max(max(len(r) for r, _ in rows), len("RTL")) + 2
                f.write(f"{'RTL':<{rtl_w}}{'Remaining Sequences':>20}\n")
                for rtl, n_seq in rows:
                    f.write(f"{rtl:<{rtl_w}}{n_seq:>20}\n")
                if not_reached_note:
                    f.write(f"{not_reached_note}\n")
                else:
                    removed = original_count - total
                    seq_word = "sequence" if removed == 1 else "sequences"
                    f.write(f"Removed {removed} {seq_word}.\n")
            else:
                f.write(f"{not_reached_note}\n")

            f.write(f"\nTotal representative sequences: {total}\n")
            f.write(f"Sequence length range: {length_stats['min']}-{length_stats['max']} bp\n")
            f.write(f"Average sequence length: {length_stats['avg']:.0f} bp\n")

            if tree_stats:
                rtl_achieved = float(rtl_table[-1][0]) if rtl_table else 1.0
                f.write("\nPhylogenetic Diversity Summary\n")
                f.write(f"RTL achieved: {rtl_achieved:.2f}\n")
                f.write(f"Total tree length: {tree_stats['total_length']:.4f}\n")
                f.write(f"Mean pairwise distance: {tree_stats['mpd']:.4f}\n")

            if accessions:
                country_counts, subclade_counts, host_counts, year_counts = {}, {}, {}, {}
                for acc in accessions:
                    info = metadata.get(acc, {})
                    country_counts[info.get('country') or 'None']   = country_counts.get(info.get('country') or 'None', 0) + 1
                    subclade_counts[info.get('subclade') or 'None'] = subclade_counts.get(info.get('subclade') or 'None', 0) + 1
                    host_counts[info.get('host') or 'None']         = host_counts.get(info.get('host') or 'None', 0) + 1
                    year_counts[info.get('year') or 'None']         = year_counts.get(info.get('year') or 'None', 0) + 1

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

            f.write("\nResults written to: \n")
            for label, path in output_files:
                f.write(f"\t✓ Saved {label} to: {path} \n")

        print(f"✓ Saved summary to: {summary_path}")

    def write_failure_summary(self, output_dir, file_prefix, country, subclade_filter, sequence_count, remark):
        output_dir.mkdir(parents=True, exist_ok=True)
        summary_path = output_dir / f"{file_prefix}_summary.txt"
        with open(summary_path, "w", encoding="utf-8") as f:
            f.write(f"Time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Command: {' '.join([os.path.basename(sys.argv[0])] + sys.argv[1:])}\n\n")
            f.write(f"Country: {country}\n")
            f.write(f"Subclade: {subclade_filter}\n")
            f.write(f"total sequences: {sequence_count}\n")
            f.write(f"\n{remark}\n")
            f.write("\nStatus: Unsuccessful\n")
        print(f"✓ Saved summary to: {summary_path}")

    @staticmethod
    def sanitize(text):
        return "_".join(
            part for part in text.strip().lower().replace("/", "_").split() if part
        )

    # A lone folder-name token (e.g. "af1a") is ambiguous — it could be a country-only or a
    # subclade-only filter, and sanitize() can't be reversed to tell which. Peek at one metadata
    # row to disambiguate (and to recover proper casing); only guess from the bare token if no
    # metadata file is available at all.
    def summarize_label(self, values):
        distinct = set(values)
        return next(iter(distinct)) if len(distinct) == 1 else "Multiple"

    # Data content can be homogeneous (or diverse) by coincidence (e.g. every Moroccan record
    # happens to be subclade AF1a, or a user-submitted FASTA happens to span many countries)
    # regardless of whether that dimension was ever actually filtered on. The file_prefix (built
    # from sanitize(country)[_sanitize(subclade)]) is the reliable signal for which dimension(s)
    # were actually filtered on — a user-submitted FASTA (SubmitSequences.py) has a file_prefix
    # that isn't a country/subclade token at all (e.g. "ea_2024"), meaning neither was ever
    # filtered on, so report "-" (None) for both rather than surfacing what's incidentally true
    # of the data.
    def resolve_country_subclade(self, file_prefix, metadata):
        countries = [info.get("country") for info in metadata.values() if info.get("country")]
        subclades = [info.get("subclade") for info in metadata.values() if info.get("subclade")]

        country_label  = self.summarize_label(countries) if countries else None
        subclade_label = self.summarize_label(subclades) if subclades else None

        country_token  = self.sanitize(country_label) if country_label not in (None, "Multiple") else None
        subclade_token = self.sanitize(subclade_label) if subclade_label not in (None, "Multiple") else None

        if country_token and subclade_token and file_prefix == f"{country_token}_{subclade_token}":
            return country_label, subclade_label
        if country_token and file_prefix == country_token:
            return country_label, "-"
        if subclade_token and file_prefix == subclade_token:
            return "-", subclade_label

        return "-", "-"

    def run(self):
        start_time = datetime.datetime.now()
        command = " ".join([os.path.basename(sys.argv[0])] + sys.argv[1:])

        # derive the base name from the containing folder (e.g. "morocco_qc" or "morocco_sequences"),
        # not the fasta filename, since QualityControl.py/RetrieveSequences.py name their files after
        # their own folder (e.g. "morocco_qc_sequences.fasta")
        source_folder_name = Path(self.fasta_file).parent.name
        if source_folder_name.endswith("_qc"):
            file_prefix = source_folder_name[:-len("_qc")]
        elif source_folder_name.endswith("_sequences"):
            file_prefix = source_folder_name[:-len("_sequences")]
        else:
            file_prefix = source_folder_name

        folder_name = f"{file_prefix}_treemmer"
        output_dir = Path(self.output_dir) / folder_name

        print(f"Using Treemmer (v0.3): {os.path.abspath('Treemmer_v0.3.py')}\n")

        if not os.path.exists(self.fasta_file):
            print(f"✗ Input file not found: {self.fasta_file}")
            sys.exit(1)

        metadata = self.load_metadata(source_folder_name)
        m49_code, subregion_name = self.resolve_country_subclade(file_prefix, metadata)

        seq_count = self.count_sequences()
        if seq_count < 4:
            n = seq_count
            seq_word = "sequence" if n == 1 else "sequences"
            was_word = "was" if n == 1 else "were"

            country_set  = m49_code not in (None, "-")
            subclade_set = subregion_name not in (None, "-")
            if country_set and subclade_set:
                filter_desc = f"{m49_code} (country) and {subregion_name} (subclade)"
            elif country_set:
                filter_desc = f"{m49_code} (country)"
            elif subclade_set:
                filter_desc = f"{subregion_name} (subclade)"
            else:
                filter_desc = "the provided sequences"

            print(f"\nWarning: Only {n} {seq_word} {was_word} found for {filter_desc}.")
            print("Treemmer requires at least 4 sequences to identify representative sequences.")
            if country_set and subclade_set:
                print("\nTo continue, please broaden your selection by:")
                print("  • selecting additional countries and/or subclades")
                print("  • removing one or more filters")
            elif country_set or subclade_set:
                print("\nPlease broaden your selection by selecting additional countries and/or subclades.")
            else:
                print("\nPlease submit additional sequences.")

            remark = f"Only {n} {seq_word} found."
            self.write_failure_summary(output_dir, folder_name, m49_code, subregion_name, seq_count, remark)
            sys.exit(0)

        output_dir.mkdir(parents=True, exist_ok=True)
        tree_file = str(output_dir / f"{file_prefix}.nwk")

        if not self.check_dependencies():
            print("\nTo install dependencies on macOS:")
            print("  brew install fasttree")
            sys.exit(1)

        print("\nBuilding phylogenetic tree using FastTree...")
        if not self.build_tree(self.fasta_file, tree_file):
            sys.exit(1)

        print("\nApplying midpoint rooting...")
        self.root_tree(tree_file)

        if not self.download_treemmer():
            sys.exit(1)

        print("\nPruning using Treemmer...")
        treemmer_success, rtl_table = self.run_treemmer(tree_file, '-RTL', self.keep_rtl)

        if not treemmer_success:
            remark = "not enough sequences for pruning"
            self.write_failure_summary(output_dir, folder_name, m49_code, subregion_name, 0, remark)
            sys.exit(1)

        tree_base = Path(tree_file).name

        # Treemmer_v0.3.py itself writes this file as "..._trimmed_list_...";
        # rename it to "..._metadata_..." to match our naming convention.
        trimmed_list_path = output_dir / f"{tree_base}_trimmed_list_RTL_{self.keep_rtl}"
        metadata_list_path = output_dir / f"{tree_base}_metadata_RTL_{self.keep_rtl}"
        representative_fasta = str(output_dir / f"{file_prefix}_representative.fasta")
        summary_path = output_dir / f"{folder_name}_summary.txt"
        original_count = self.count_sequences()

        if trimmed_list_path.exists():
            trimmed_list_path.rename(metadata_list_path)
            accessions = self.enrich_trimmed_list(metadata_list_path, metadata)
            self.write_representative_sequences(accessions, self.fasta_file, representative_fasta)
            print(f"✓ Saved {len(accessions)} representative sequences to: {representative_fasta}")
            print(f"✓ Saved representative sequence list to: {metadata_list_path}")

            # Treemmer_v0.3.py writes this as "..._trimmed_tree_..."; rename to match our convention.
            old_trimmed_tree_path = output_dir / f"{tree_base}_trimmed_tree_RTL_{self.keep_rtl}"
            trimmed_tree_path = output_dir / f"{tree_base}_RTL_{self.keep_rtl}"
            old_trimmed_tree_path.rename(trimmed_tree_path)
            print(f"✓ Saved pruned tree to: {trimmed_tree_path}")

            plot_path = f"{trimmed_tree_path}.png"
            self.plot_tree(str(trimmed_tree_path), plot_path)
            print(f"✓ Saved pruned tree plot to: {plot_path}")

            tree_stats = self.calculate_mpd(str(trimmed_tree_path))
            length_stats = self.compute_length_stats(representative_fasta)
            output_files = [
                (f"{len(accessions)} representative sequences", representative_fasta),
                ("metadata", str(metadata_list_path)),
                ("midpoint-rooted tree", tree_file),
                ("pruned tree", str(trimmed_tree_path)),
                ("pruned tree plot", plot_path),
                ("summary", str(summary_path)),
            ]
            self.write_summary(start_time, command, rtl_table, None, original_count,
                                accessions, metadata, tree_stats, length_stats, output_files, summary_path)
        else:
            not_reached_note = f"INFO: RTL {self.keep_rtl} was not reached due to high phylogenetic diversity. All sequences were retained as representative."
            print(f"\n{not_reached_note}")
            shutil.copy(self.fasta_file, representative_fasta)
            accessions = list(metadata.keys()) if metadata else []
            print(f"\n✓ Saved {len(accessions)} representative sequences to: {representative_fasta}")

            self.write_metadata_list(accessions, metadata, metadata_list_path)
            print(f"✓ Saved representative sequence list to: {metadata_list_path}")

            # when Treemmer_v0.3.py runs through every leaf without reaching the target RTL, it
            # leaves behind its own diagnostic RTL-vs-leaves plot/data files — not part of our
            # output set (see the "Results written to" list below), so clean them up
            if rtl_table:
                for suffix in ("_res_1_LD", "_res_1_TLD.pdf"):
                    leftover = output_dir / f"{tree_base}{suffix}"
                    if leftover.exists():
                        leftover.unlink()

            plot_path = f"{tree_file}.png"
            self.plot_tree(tree_file, plot_path)
            print(f"✓ Saved midpoint-rooted tree plot to: {plot_path}")

            tree_stats = self.calculate_mpd(tree_file)
            length_stats = self.compute_length_stats(representative_fasta)
            output_files = [
                (f"{len(accessions)} representative sequences", representative_fasta),
                ("metadata", str(metadata_list_path)),
                ("midpoint-rooted tree", tree_file),
                ("tree plot", plot_path),
                ("summary", str(summary_path)),
            ]
            self.write_summary(start_time, command, rtl_table, not_reached_note, original_count,
                                accessions, metadata, tree_stats, length_stats, output_files, summary_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Generate a midpoint-rooted phylogenetic tree and prune it using Treemmer '
                     '(v0.3) to retain representative sequences.'
    )
    parser.add_argument('-i', '--input', required=True,
                        help='Input FASTA file from QualityControl.py')
    parser.add_argument('-o', '--output-dir', default='results',
                        help='Directory where the output folder will be saved (default: results/)')
    parser.add_argument('-RTL', '--keep-rtl', type=float, default=0.9,
                        help='Relative tree length to keep with Treemmer 0-1 (default: 0.9)')
    args = parser.parse_args()

    workflow = TreemmerWorkflow(args.input, args.output_dir, args.keep_rtl)

    workflow.run()
