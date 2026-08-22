# RABV primer design pipeline
This pipeline designs rabies virus (RABV) primers through:

1. Retrieving aligned RABV sequences from RABV-gDB.
2. Selecting a representative subset of sequences that retains 90% of the phylogenetic diversity.
3. Designing primers based on the selected representative sequences.

# Step 1. uses retrieve_sequences.py (/scripts/retrieve_sequences.py).

```bash
usage: retrieve_sequences.py -db DATABASE [-c COUNTRY] [-s SUBCLADE] [-o OUTPUT_DIR]

required:
  -db, --database     path to the SQLite database file

optional:
  -c, --country       country name (e.g. 'China')
  -s, --subclade      subclade name (e.g. 'CA1')
  -o, --output_dir    output directory (default: results/)
```

## Example commands:
1. retrieve aligned sequences by country only:
```bash
python scripts/retrieve_sequences.py -db rabv_V-gDB_10042026.db -c Indonesia
```

2. retrieve aligned sequences by subclade only:
```bash
python scripts/retrieve_sequences.py -db rabv_V-gDB_10042026.db -s CA1
```

3. retrieve aligned sequences by both country and subclade:
```bash
python scripts/retrieve_sequences.py -db rabv_V-gDB_10042026.db -c China -s CA1 
```

## Outputs
Step 1 generates a single output folder (_sequences) containing the following files:

| File  | Description |
| ------------- | ------------- |
| <name>_metadata.tsv  | Tab-delimited file containing the sequence ID, reference sequence, selected country, and selected subclade. |
| <name>_sequences.fasta | FASTA file containing the aligned sequences for the selected country and/or subclade. |
| <name>_statistics.txt | Summary statistics of the retrieved dataset. These include the distribution of countries within the selected subclade and/or the distribution of subclades within the selected country. |

# Step 2. uses treemmer.py (/scripts/treemmer.py) and treemmer_v0.3.py (/scripts/treemmer_v0.3.py).
`treemmer.py` depends on `treemmer_v0.3.py`, which contains the phylogenetic pruning algorithm (Menardo et al., 2018).

| Script | Description |
|--------|-------------|
| `treemmer_v0.3.py` | Performs phylogenetic pruning to retain the representative sequences that account for 90% of the total phylogenetic diversity. |
| `treemmer.py` | 1. Builds a phylogenetic tree from the aligned sequences obtained in Step 1.<br>2. Applies midpoint rooting and plots the full phylogenetic tree.<br>3. Runs `treemmer_v0.3.py` to select representative  sequences.<br>4. Exports the aligned representative sequences in FASTA format. |

If the relatatie tree length (RTL) never drops below the 90% diversity threshold — meaning even the most aggressively pruned tree still retains more diversity than the threshold — treemmer.py keeps all input sequences as the representative (100% diversity retained).

```bash
usage: python3 scripts/treemmer.py -i <fasta> [-o <output_dir>] [-RTL <0-1>] [-X <n>]

arguments:
  -i,   --input       aligned FASTA from retrieve_sequences.py (required)
  -o,   --output-dir  output directory (default: results/)
  -RTL, --keep-rtl    phylogenetic diversity to retain, 0–1 (default is 0.9)
  -X,   --keep        keep exact number of sequences instead of using RTL
```

## Example commands:
1. default (90% diversity) 
```bash
python scripts/treemmer.py -i results/indonesia_sequences/indonesia_sequences.fasta
```

2. custom diversity threshold
```bash
python scripts/treemmer.py -i results/china_ca1_sequences/china_ca1_sequences.fasta -RTL 0.85
```

3. exact number of sequences
```bash
python scripts/treemmer.py -i results/indonesia_sequences/indonesia_sequences.fasta -X 20
```

## Outputs
Step 2 generates a single output folder (_treemmer) containing the following files:

**treemmer_v0.3.py outputs**
1. Scenario 1: RTL reached the 90% diversity threshold
   
| File  | Description |
| ------------- | ------------- |
|<name>.nwk_trimmed_tree_RTL_0.9 | Newick-format tree file containing the pruned phylogenetic tree of the representative sequences. |
| <name>.nwk_trimmed_list_RTL_0.9 | Text file containing the sequence IDs of the representative sequences.  |

2. Scenario 2: RTL never reached the 90% diversity threshold

| File  | Description |
| ------------- | ------------- |
| <name>.nwk_res_1_LD  | Tab-separated log file containing the results of each Treemmer pruning iteration.  |
| <name>.nwk_res_1_TLD.pdf |  Scatter plot of RTL vs the number of remaining sequences. This shows the change in phylogenetic diversity across all pruning iterations.  |

**treemmer.py outputs**
| File  | Description |
| ------------- | ------------- |
| <name>_representative.fasta | FASTA file containing the aligned representative sequences selected by Treemmer. |
| <name>.nwk  | Newick-format tree file containing the midpoint-rooted phylogenetic tree of all aligned sequences obtained in Step 1. |
| <name>.nwk.png |  Tree diagram (mid-point rooted) showing the taxa and branch lengths of the representative  sequences selected by Treemmer.  |
| <name>_statistics.txt | Summary statistics of the aligned representative sequences selected by Treemmer. These include the distribution of countries within the selected subclade and/or the distribution of subclades within the selected country before and after Treemmer. |
