#!/usr/bin/env python3
import argparse
import base64
import csv
import html
import io
import os
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap

DEFAULT_THREE_PRIME_THRESHOLD = 1
DEFAULT_SNP_THRESHOLD = 4

# traffic-light heatmap flags: no mismatch (green), a mismatch below the flagging thresholds
# (yellow), or flagged — meets the same ≥3'/≥snp thresholds used to flag primers dataset-wide
# (red). Each cell is still annotated with its raw mismatch count regardless of color.
HEATMAP_COLORS = ["#A5D6A7", "#FFF59D", "#EF9A9A"]

# ISO alpha-3 country codes ("id" column) keyed by every name the database's own m49_country
# table recognizes (both display_name and full_name), used to label the heatmap's country bar
# with codes instead of full country names. A few aliases are added at the bottom for
# non-standard country strings actually found in meta_data.country that don't match m49's own
# spelling (typos, abbreviations, alternate names).
COUNTRY_CODES = {
    'Afghanistan': 'AFG',
    'Albania': 'ALB',
    'Algeria': 'DZA',
    'American Samoa': 'ASM',
    'Andorra': 'AND',
    'Angola': 'AGO',
    'Anguilla': 'AIA',
    'Antarctica': 'ATA',
    'Antigua and Barbuda': 'ATG',
    'Argentina': 'ARG',
    'Armenia': 'ARM',
    'Aruba': 'ABW',
    'Australia': 'AUS',
    'Austria': 'AUT',
    'Azerbaijan': 'AZE',
    'Bahamas': 'BHS',
    'Bahrain': 'BHR',
    'Bangladesh': 'BGD',
    'Barbados': 'BRB',
    'Belarus': 'BLR',
    'Belgium': 'BEL',
    'Belize': 'BLZ',
    'Benin': 'BEN',
    'Bermuda': 'BMU',
    'Bhutan': 'BTN',
    'Bolivia': 'BOL',
    'Bolivia (Plurinational State of)': 'BOL',
    'Bonaire, Sint Eustatius and Saba': 'BES',
    'Bosnia and Herzegovina': 'BIH',
    'Botswana': 'BWA',
    'Bouvet Island': 'BVT',
    'Brazil': 'BRA',
    'British Indian Ocean Territory': 'IOT',
    'British Virgin Islands': 'VGB',
    'Brunei Darussalam': 'BRN',
    'Bulgaria': 'BGR',
    'Burkina Faso': 'BFA',
    'Burundi': 'BDI',
    'Cabo Verde': 'CPV',
    'Cambodia': 'KHM',
    'Cameroon': 'CMR',
    'Canada': 'CAN',
    'Cape Verde': 'CPV',
    'Cayman Islands': 'CYM',
    'Central African Republic': 'CAF',
    'Chad': 'TCD',
    'Chile': 'CHL',
    'China': 'CHN',
    'China, Hong Kong Special Administrative Region': 'HKG',
    'China, Macao Special Administrative Region': 'MAC',
    'Christmas Island': 'CXR',
    'Cocos (Keeling) Islands': 'CCK',
    'Colombia': 'COL',
    'Comoros': 'COM',
    'Congo': 'COG',
    'Congo (Brazzaville)': 'COG',
    'Cook Islands': 'COK',
    'Costa Rica': 'CRI',
    'Croatia': 'HRV',
    'Cuba': 'CUB',
    'Curaçao': 'CUW',
    'Cyprus': 'CYP',
    'Czechia': 'CZE',
    "Côte d'Ivoire": 'CIV',
    "Democratic People's Republic of Korea": 'PRK',
    'Democratic Republic of the Congo': 'COD',
    'Denmark': 'DNK',
    'Djibouti': 'DJI',
    'Dominica': 'DMA',
    'Dominican Republic': 'DOM',
    'Ecuador': 'ECU',
    'Egypt': 'EGY',
    'El Salvador': 'SLV',
    'Equatorial Guinea': 'GNQ',
    'Eritrea': 'ERI',
    'Estonia': 'EST',
    'Ethiopia': 'ETH',
    'Falkland Islands': 'FLK',
    'Falkland Islands (Malvinas)': 'FLK',
    'Faroe Islands': 'FRO',
    'Fiji': 'FJI',
    'Finland': 'FIN',
    'France': 'FRA',
    'French Guiana': 'GUF',
    'French Polynesia': 'PYF',
    'French Southern Territories': 'ATF',
    'Gabon': 'GAB',
    'Gambia': 'GMB',
    'Georgia': 'GEO',
    'Germany': 'DEU',
    'Ghana': 'GHA',
    'Gibraltar': 'GIB',
    'Greece': 'GRC',
    'Greenland': 'GRL',
    'Grenada': 'GRD',
    'Guadeloupe': 'GLP',
    'Guam': 'GUM',
    'Guatemala': 'GTM',
    'Guernsey': 'GGY',
    'Guinea': 'GIN',
    'Guinea-Bissau': 'GNB',
    'Guyana': 'GUY',
    'Haiti': 'HTI',
    'Heard Island and McDonald Islands': 'HMD',
    'Heard Island and Mcdonald Islands': 'HMD',
    'Holy See': 'VAT',
    'Holy See (Vatican City State)': 'VAT',
    'Honduras': 'HND',
    'Hong Kong': 'HKG',
    'Hungary': 'HUN',
    'Iceland': 'ISL',
    'India': 'IND',
    'Indonesia': 'IDN',
    'Iran': 'IRN',
    'Iran (Islamic Republic of)': 'IRN',
    'Iraq': 'IRQ',
    'Ireland': 'IRL',
    'Isle of Man': 'IMN',
    'Israel': 'ISR',
    'Italy': 'ITA',
    'Jamaica': 'JAM',
    'Japan': 'JPN',
    'Jersey': 'JEY',
    'Jordan': 'JOR',
    'Kazakhstan': 'KAZ',
    'Kenya': 'KEN',
    'Kiribati': 'KIR',
    'Kosovo': 'XKX',
    'Kuwait': 'KWT',
    'Kyrgyzstan': 'KGZ',
    'Lao': 'LAO',
    "Lao People's Democratic Republic": 'LAO',
    'Latvia': 'LVA',
    'Lebanon': 'LBN',
    'Lesotho': 'LSO',
    'Liberia': 'LBR',
    'Libya': 'LBY',
    'Liechtenstein': 'LIE',
    'Lithuania': 'LTU',
    'Luxembourg': 'LUX',
    'Macao': 'MAC',
    'Madagascar': 'MDG',
    'Malawi': 'MWI',
    'Malaysia': 'MYS',
    'Maldives': 'MDV',
    'Mali': 'MLI',
    'Malta': 'MLT',
    'Marshall Islands': 'MHL',
    'Martinique': 'MTQ',
    'Mauritania': 'MRT',
    'Mauritius': 'MUS',
    'Mayotte': 'MYT',
    'Mexico': 'MEX',
    'Micronesia': 'FSM',
    'Micronesia (Federated States of)': 'FSM',
    'Moldova': 'MDA',
    'Monaco': 'MCO',
    'Mongolia': 'MNG',
    'Montenegro': 'MNE',
    'Montserrat': 'MSR',
    'Morocco': 'MAR',
    'Mozambique': 'MOZ',
    'Myanmar': 'MMR',
    'Namibia': 'NAM',
    'Nauru': 'NRU',
    'Nepal': 'NPL',
    'Netherlands': 'NLD',
    'New Caledonia': 'NCL',
    'New Zealand': 'NZL',
    'Nicaragua': 'NIC',
    'Niger': 'NER',
    'Nigeria': 'NGA',
    'Niue': 'NIU',
    'Norfolk Island': 'NFK',
    'North Korea': 'PRK',
    'Northern Mariana Islands': 'MNP',
    'Norway': 'NOR',
    'Oman': 'OMN',
    'Pakistan': 'PAK',
    'Palau': 'PLW',
    'Palestinian Territory': 'PSE',
    'Panama': 'PAN',
    'Papua New Guinea': 'PNG',
    'Paraguay': 'PRY',
    'Peru': 'PER',
    'Philippines': 'PHL',
    'Pitcairn': 'PCN',
    'Poland': 'POL',
    'Portugal': 'PRT',
    'Puerto Rico': 'PRI',
    'Qatar': 'QAT',
    'Republic of Korea': 'KOR',
    'Republic of Macedonia': 'MKD',
    'Republic of Moldova': 'MDA',
    'Romania': 'ROU',
    'Russia': 'RUS',
    'Russian Federation': 'RUS',
    'Rwanda': 'RWA',
    'Réunion': 'REU',
    'Saint Helena': 'SHN',
    'Saint Kitts and Nevis': 'KNA',
    'Saint Lucia': 'LCA',
    'Saint Martin (French Part)': 'MAF',
    'Saint Pierre and Miquelon': 'SPM',
    'Saint Vincent and Grenadines': 'VCT',
    'Saint Vincent and the Grenadines': 'VCT',
    'Saint-Barthélemy': 'BLM',
    'Saint-Martin (French part)': 'MAF',
    'Samoa': 'WSM',
    'San Marino': 'SMR',
    'Sao Tome and Principe': 'STP',
    'Saudi Arabia': 'SAU',
    'Senegal': 'SEN',
    'Serbia': 'SRB',
    'Seychelles': 'SYC',
    'Sierra Leone': 'SLE',
    'Singapore': 'SGP',
    'Sint Maarten (Dutch part)': 'SXM',
    'Slovakia': 'SVK',
    'Slovenia': 'SVN',
    'Solomon Islands': 'SLB',
    'Somalia': 'SOM',
    'South Africa': 'ZAF',
    'South Georgia and the South Sandwich Islands': 'SGS',
    'South Korea': 'KOR',
    'South Sudan': 'SSD',
    'Spain': 'ESP',
    'Sri Lanka': 'LKA',
    'State of Palestine': 'PSE',
    'Sudan': 'SDN',
    'Suriname': 'SUR',
    'Svalbard and Jan Mayen Islands': 'SJM',
    'Swaziland': 'SWZ',
    'Sweden': 'SWE',
    'Switzerland': 'CHE',
    'Syria': 'SYR',
    'Syrian Arab Republic': 'SYR',
    'Taiwan': 'TWN',
    'Tajikistan': 'TJK',
    'Tanzania': 'TZA',
    'Thailand': 'THA',
    'The former Yugoslav Republic of Macedonia': 'MKD',
    'Timor-Leste': 'TLS',
    'Togo': 'TGO',
    'Tokelau': 'TKL',
    'Tonga': 'TON',
    'Trinidad and Tobago': 'TTO',
    'Tunisia': 'TUN',
    'Turkey': 'TUR',
    'Turkmenistan': 'TKM',
    'Turks and Caicos Islands': 'TCA',
    'Tuvalu': 'TUV',
    'Uganda': 'UGA',
    'Ukraine': 'UKR',
    'United Arab Emirates': 'ARE',
    'United Kingdom': 'GBR',
    'United Kingdom of Great Britain and Northern Ireland': 'GBR',
    'United Republic of Tanzania': 'TZA',
    'United States': 'USA',
    'United States Minor Outlying Islands': 'UMI',
    'United States Virgin Islands': 'VIR',
    'United States of America': 'USA',
    'Uruguay': 'URY',
    'Uzbekistan': 'UZB',
    'Vanuatu': 'VUT',
    'Venezuela': 'VEN',
    'Venezuela (Bolivarian Republic of)': 'VEN',
    'Viet Nam': 'VNM',
    'Vietnam': 'VNM',
    'Virgin Islands, US': 'VIR',
    'Wallis and Futuna Islands': 'WLF',
    'Western Sahara': 'ESH',
    'Yemen': 'YEM',
    'Zambia': 'ZMB',
    'Zimbabwe': 'ZWE',
    'Åland Islands': 'ALA',
    # aliases: non-standard country strings seen in meta_data.country that don't match
    # m49_country's own display_name/full_name spelling (typos, abbreviations, alternate names)
    'Bosnia': 'BIH',
    'Columbia': 'COL',
    'Ivory Coast': 'CIV',
    'Korea': 'KOR',
    'Laos': 'LAO',
    'USA': 'USA',
    'USA (TX)': 'USA',
    'West Bank': 'PSE',
}

DETAIL_HEADER = ["amplicon", "primer", "position_in_primer", "ref_base", "primer_base", "types"]

# the trailing "_<variant>" is PrimalScheme3's own convention and is optional here — an
# externally-supplied primer.bed (e.g. checking someone else's scheme) may omit it entirely
# (e.g. "rabv_ea_1_LEFT" instead of "988d9a4f_24_RIGHT_1")
PRIMER_NAME_RE = re.compile(r"^(?P<prefix>.+)_(?P<num>\d+)_(?P<side>LEFT|RIGHT)(?:_(?P<variant>\d+))?$")

# 3' mismatches affect primer extension far more than a mismatch elsewhere, so they get their own,
# much lower flagging bar, checked independently of the snp bar — a primer with e.g. 1 3' mismatch
# and 6 snp mismatches still flags on the snp count alone, it doesn't need zero 3' mismatches
THREE_PRIME_LABEL = "3' mismatch"

# IUPAC nucleotide complement, "-" (gap) maps to itself
COMPLEMENT_TABLE = str.maketrans("ACGTRYSWKMBDHVN-", "TGCAYRSWMKVHDBN-")

def complement_base(base):
    return base.upper().translate(COMPLEMENT_TABLE)

# LEFT primers are always "+" strand and RIGHT always "-" strand (PrimalScheme3/primer.bed
# convention, confirmed against primer.bed's own strand column) — so for a RIGHT primer, the
# ref_base recorded by PrimerSequenceMismatch.py is BLAST's already-reverse-complemented (primer-
# oriented) sseq value, not the literal base actually present in the sequenced genome. Complementing
# it back recovers that literal base — LEFT primers need no correction, they're already plus/plus.
def is_right_primer(primer_name):
    m = PRIMER_NAME_RE.match(primer_name)
    return m is not None and m.group("side") == "RIGHT"

def literal_ref_base(base, primer_name):
    return complement_base(base) if is_right_primer(primer_name) else base

# each per-accession file (from PrimerSequenceMismatch.py) has a coverage table, a blank line, then the
# detail table — find the detail header and read everything after it as that accession's mismatches
def read_detail_rows(path):
    with open(path, encoding="utf-8") as f:
        lines = f.readlines()

    header_line = "\t".join(DETAIL_HEADER)
    try:
        start = next(i for i, line in enumerate(lines) if line.rstrip("\n") == header_line)
    except StopIteration:
        return []

    reader = csv.DictReader(lines[start + 1:], fieldnames=DETAIL_HEADER, delimiter="\t")
    return [row for row in reader if any(v.strip() for v in row.values())]

# *_mismatches.tsv also matches this script's own "<name>_mismatches.tsv" output, so glob results
# are filtered down to files that actually contain a per-accession detail header
def find_accession_files(mismatches_dir):
    header_line = "\t".join(DETAIL_HEADER)
    files = []
    for path in sorted(mismatches_dir.glob("*_mismatches.tsv")):
        with open(path, encoding="utf-8") as f:
            if any(line.rstrip("\n") == header_line for line in f):
                files.append(path)
    if not files:
        print(f"✗ No *_mismatches.tsv files found in: {mismatches_dir}")
        sys.exit(1)
    return files

def sanitize(text):
    return "_".join(part for part in text.strip().lower().replace("/", "_").split() if part)

def summarize_label(values):
    distinct = set(values)
    return next(iter(distinct)) if len(distinct) == 1 else "Multiple"

# best-effort label for the pipeline log — aggregates every data row of the metadata.tsv sibling
# of the <name>_qc/ folder next to mismatches_dir, not just the first row: a user-submitted FASTA
# (SubmitSequences.py) can have most rows blank (only accessions that matched an existing database
# entry get real metadata) or span many countries, so a single row is not representative.
#
# Data can be homogeneous (or diverse) by coincidence regardless of whether that dimension was
# ever actually filtered on, so "name" (== file_prefix) is checked against sanitize(country)/
# sanitize(subclade) — the reliable signal for what was actually filtered — before trusting the
# aggregated value; if it matches neither (e.g. a user-submitted FASTA, which never filters by
# country/subclade at all), both report "None" rather than what's incidentally true of the data.
def resolve_country_subclade(mismatches_dir):
    name = mismatches_dir.name[:-len("_mismatches")] if mismatches_dir.name.endswith("_mismatches") \
        else mismatches_dir.name
    metadata_path = mismatches_dir.parent / f"{name}_qc" / f"{name}_qc_metadata.tsv"
    if not metadata_path.exists():
        return None, None
    with open(metadata_path, encoding="utf-8") as f:
        header_lower = [h.lower() for h in f.readline().rstrip("\n").split("\t")]
        try:
            country_idx = header_lower.index("country")
            subclade_idx = header_lower.index("subclade")
        except ValueError:
            return None, None
        countries, subclades = [], []
        for line in f:
            row = line.rstrip("\n").split("\t")
            if len(row) > max(country_idx, subclade_idx):
                if row[country_idx]:
                    countries.append(row[country_idx])
                if row[subclade_idx]:
                    subclades.append(row[subclade_idx])

    country = summarize_label(countries) if countries else None
    subclade = summarize_label(subclades) if subclades else None

    country_token = sanitize(country) if country not in (None, "Multiple") else None
    subclade_token = sanitize(subclade) if subclade not in (None, "Multiple") else None

    if country_token and subclade_token and name == f"{country_token}_{subclade_token}":
        return country, subclade
    if country_token and name == country_token:
        return country, None
    if subclade_token and name == subclade_token:
        return None, subclade
    return None, None

# per-accession Country, keyed by accession — unlike resolve_country_subclade() (which collapses
# every row into one aggregate label), this keeps each accession's own value since the heatmap's
# country bar groups/labels rows individually. The mismatches folder's own name isn't a reliable
# guide to which *_qc/ folder its accessions came from — e.g. run_pipeline_verify.sh names the QC
# folder after the retrieved country ("morocco_qc") but the mismatches folder after a separate
# --output-name ("morocco_rabv_ea_mismatches") — so every *_qc_metadata.tsv under mismatches_dir's
# parent is searched for the accessions actually being looked up, rather than guessing one path
def read_accession_countries(mismatches_dir, accessions):
    needed = set(accessions)
    lookup = {}
    for metadata_path in sorted(mismatches_dir.parent.glob("*_qc/*_qc_metadata.tsv")):
        with open(metadata_path, encoding="utf-8") as f:
            header_lower = [h.lower() for h in f.readline().rstrip("\n").split("\t")]
            try:
                id_idx = header_lower.index("accession")
                country_idx = header_lower.index("country")
            except ValueError:
                continue
            for line in f:
                row = line.rstrip("\n").split("\t")
                if len(row) > max(id_idx, country_idx) and row[id_idx] in needed:
                    lookup[row[id_idx]] = row[country_idx] or None
        if len(lookup) == len(needed):
            break
    return lookup

# the country's ISO alpha-3 code, its own name if it isn't a recognized m49 country (rather than
# silently dropping it), or "unknown" when no country is on file at all (e.g. a user-submitted
# accession with no database match)
def country_bar_label(country):
    if not country:
        return "unknown"
    return COUNTRY_CODES.get(country, country)

def find_recurring_mismatches(mismatches_dir, threshold):
    accession_files = find_accession_files(mismatches_dir)

    total = len(accession_files)

    # keyed by (primer, position_in_primer) — a specific recurring site in a specific primer;
    # value is a Counter of every ref_base seen there (a site can mismatch to more than one
    # alternate base across the dataset), plus the amplicon/primer_base/type for that site
    site_base_counts = defaultdict(Counter)
    site_info = {}
    # which accessions actually carry each site — lets build_recommendations_rows compute the
    # true overlap (accessions carrying *every* site in a primer's group), not just each site's
    # own individual proportion
    site_accessions = defaultdict(set)

    for path in accession_files:
        accession = path.name[:-len("_mismatches.tsv")]
        for row in read_detail_rows(path):
            if row["position_in_primer"] == "-":
                continue  # "no alignment found" placeholder row, not a positional mismatch

            key = (row["primer"], row["position_in_primer"])
            site_base_counts[key][row["ref_base"]] += 1
            site_info[key] = (row["amplicon"], row["primer_base"], row["types"])
            site_accessions[key].add(accession)

    suggestions = []
    for (primer, position), base_counts in site_base_counts.items():
        count = sum(base_counts.values())
        mismatch_prop = count / total
        if mismatch_prop < threshold:
            continue
        amplicon, primer_base, issue_type = site_info[(primer, position)]
        # distribution of alternate bases actually found in the sequences at this site (literal,
        # not primer-oriented), most common first — generalizes to indels too, since ref_base is
        # "-" for a primer insertion
        ref_dist = ", ".join(
            f"{literal_ref_base(base, primer)}:{base_count / count:.2f}"
            for base, base_count in base_counts.most_common()
        )
        # the base the primer needs to carry to anneal to the majority literal ref base — derived
        # from the stored (primer-oriented) majority base, which already *is* that answer: BLAST's
        # strand normalization and the anneal-complement requirement cancel out for RIGHT primers,
        # and both are no-ops for LEFT primers, so no further transform is needed here
        consensus_base = base_counts.most_common(1)[0][0]
        suggestions.append({
            "amplicon": amplicon,
            "primer.issues": primer,
            "positionInPrimer": position,
            "refSeq": ref_dist,
            "primerSeq": primer_base,
            "primerConsensus": consensus_base,
            "type": issue_type.replace(THREE_PRIME_LABEL, "3'"),
            "count": count,
            "mismatch_prop": f"{mismatch_prop:.2f}",
            "accessions": site_accessions[(primer, position)],
        })

    suggestions.sort(key=lambda r: r["count"], reverse=True)
    return suggestions, total

# "a, b and c" — Oxford-comma-free join used for both recommendation bullets
def join_names(names):
    if len(names) == 1:
        return names[0]
    return ", ".join(names[:-1]) + " and " + names[-1]

# a primer is red if it has ANY recurring 3' mismatch site (a single one is enough — every row in
# suggestions already meets the ≥threshold mismatch_prop bar to be there at all, so "1 3' in ≥20%
# of the sequences" is just "any 3' row present") or ≥RED_SNP_INDEL_SITE_THRESHOLD recurring
# SNP/indel sites (pooled together); everything else in the recurring table is yellow. Within one
# primer, 3' takes priority over SNP/indel for the displayed issue text, since a single 3'
# recurring mismatch alone is already enough to flag it red regardless of how many SNP/indel sites
# it also has.
RED_SNP_INDEL_SITE_THRESHOLD = 4

# a recurring SNP/indel site below this population frequency isn't worth a redesign/alternate-
# primer call — only 3' mismatches (always kept, any one is disruptive regardless of frequency)
# skip this filter
MIN_SNP_INDEL_PROP_FOR_RECOMMENDATION = 0.40

# builds one "N <label> mismatch(es) found in [at least] X% of the sequences" clause for a
# same-type group of sites, plus the (n, prop) pair used for row sorting
def describe_group(group, label):
    n = len(group)
    if n == 1:
        prop = float(group[0]["mismatch_prop"])
        return f"1 {label} mismatch found in {prop:.0%} of the sequences", n, prop
    # each site's own mismatch_prop is only how many sequences carry *that* site — sites in the
    # group need not co-occur in the same sequences (their accession sets can be disjoint), so
    # the group's least-prevalent site is the only prevalence every one of them is guaranteed
    # to clear
    prop = min(float(r["mismatch_prop"]) for r in group)
    return f"{n} {label} mismatches found in at least {prop:.0%} of the sequences", n, prop

def build_recommendations_rows(suggestions, n_sequences):
    by_primer = defaultdict(list)
    for row in suggestions:
        by_primer[row["primer.issues"]].append(row)

    rows = []
    for primer, primer_rows in by_primer.items():
        three_prime_rows = [r for r in primer_rows if "3'" in r["type"]]
        snp_indel_rows = [r for r in primer_rows if "3'" not in r["type"]
                           and float(r["mismatch_prop"]) >= MIN_SNP_INDEL_PROP_FOR_RECOMMENDATION]

        if three_prime_rows:
            clauses = [describe_group(three_prime_rows, "3'")]
        elif snp_indel_rows:
            # SNP and indel are pooled for the red/yellow decision below (both count toward the
            # same ≥40% qualification bar and ≥RED_SNP_INDEL_SITE_THRESHOLD count), but reported
            # as separate clauses here rather than blended under one generic "SNP/indel" label
            snp_only = [r for r in snp_indel_rows if r["type"] == "snp"]
            indel_only = [r for r in snp_indel_rows if r["type"] != "snp"]
            clauses = []
            if snp_only:
                clauses.append(describe_group(snp_only, "SNP"))
            if indel_only:
                clauses.append(describe_group(indel_only, "indel"))
        else:
            continue  # only low-frequency (<40%) SNP/indel sites left — not worth a recommendation row

        issue = "; ".join(text for text, _, _ in clauses)
        n = sum(n for _, n, _ in clauses)
        prop = max(prop for _, _, prop in clauses)

        is_red = bool(three_prime_rows) or len(snp_indel_rows) >= RED_SNP_INDEL_SITE_THRESHOLD

        rows.append({
            "primer": primer,
            "issue": issue,
            "recommendation": ("Consider redesigning new primers" if is_red
                                else "Consider alternate or degenerate primers"),
            "severity": "red" if is_red else "yellow",
            "n": n,
            "prop": prop,
        })

    # red rows before yellow; within each severity, most mismatches first, then — among rows
    # with the same mismatch count — highest population proportion first (e.g. 2 SNP @ 30% before
    # 2 SNP @ 10% before 1 SNP @ 100%: mismatch count outranks prevalence), falling back to the
    # primer name only to break an exact tie
    def sort_key(r):
        return (r["severity"] != "red", -r["n"], -r["prop"], r["primer"])

    rows.sort(key=sort_key)
    return rows

# amplicon number ascending, LEFT before RIGHT, variant ascending (0 if no variant suffix at
# all) — falls back to alphabetical for any primer name that doesn't match at all
def primer_sort_key(primer_name):
    m = PRIMER_NAME_RE.match(primer_name)
    if not m:
        return (float("inf"), 0, 0, primer_name)
    return (int(m.group("num")), 0 if m.group("side") == "LEFT" else 1, int(m.group("variant") or 0), "")

# primer.bed's 4th column is the primer name (PrimerSequenceMismatch.py's read_primer_bed)
def read_primer_names(bed_path):
    names = []
    with open(bed_path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 4:
                names.append(parts[3])
    return names

def build_heatmap_matrix(mismatches_dir, three_prime_threshold, other_threshold, primer_bed=None):
    accession_files = find_accession_files(mismatches_dir)

    counts = defaultdict(lambda: defaultdict(int))
    # per-cell breakdown of the two mismatch groups the --3-prime/--snp thresholds are checked
    # against, tallied per (accession, primer) pair. snp_indel_counts covers both SNP and indel
    # mismatches together — a 3'-proximal mismatch (of either kind) is excluded here since it's
    # already independently checked via three_prime_counts, at its own much lower bar
    three_prime_counts = defaultdict(lambda: defaultdict(int))
    snp_indel_counts = defaultdict(lambda: defaultdict(int))
    accessions = []
    seen_primers = set()

    for path in accession_files:
        accession = path.name[:-len("_mismatches.tsv")]
        accessions.append(accession)
        for row in read_detail_rows(path):
            if row["position_in_primer"] == "-":
                continue  # "no alignment found" placeholder, not a positional mismatch
            primer = row["primer"]
            counts[accession][primer] += 1
            seen_primers.add(primer)
            types = row["types"].split(",")
            if THREE_PRIME_LABEL in types:
                three_prime_counts[accession][primer] += 1
            else:
                snp_indel_counts[accession][primer] += 1

    country_lookup = read_accession_countries(mismatches_dir, accessions)

    # group rows by country bar label so same-country accessions end up adjacent and the heatmap's
    # right-side bar can merge them into one labeled span instead of one row per accession — stable
    # sort keeps each country's own accessions in their original (alphabetical) order
    row_labels = [country_bar_label(country_lookup.get(acc)) for acc in accessions]
    order = sorted(range(len(accessions)), key=lambda i: row_labels[i])
    accessions = [accessions[i] for i in order]
    row_labels = [row_labels[i] for i in order]

    # with a primer.bed, show every primer that was checked (zero-mismatch ones included);
    # without one, only primers that had at least one mismatch anywhere are known
    primers = set(read_primer_names(primer_bed)) | seen_primers if primer_bed else seen_primers
    primers = sorted(primers, key=primer_sort_key)
    matrix = [[counts[acc].get(primer, 0) for primer in primers] for acc in accessions]

    # per-cell flag — 0 (no mismatch), 1 (a mismatch that doesn't meet the flagging thresholds),
    # or 2 (flagged: this primer/accession pair alone meets ≥three_prime_threshold 3' mismatches
    # or ≥other_threshold combined SNP/indel mismatches) — independent of the raw mismatch count,
    # so e.g. 5 plain SNP/indel mismatches below other_threshold still flags as 1 (yellow), not 2 (red)
    def cell_flag(acc, primer):
        if (three_prime_counts[acc].get(primer, 0) >= three_prime_threshold
                or snp_indel_counts[acc].get(primer, 0) >= other_threshold):
            return 2
        return 1 if counts[acc].get(primer, 0) else 0

    flag_matrix = [[cell_flag(acc, primer) for primer in primers] for acc in accessions]

    return accessions, primers, matrix, flag_matrix, row_labels, bool(seen_primers)

# 90°-rotated primer-name tick labels need vertical room proportional to their *length*, not just
# a fixed pad — with few accessions (a short fig_h) and long primer names (e.g. PrimalScheme3's
# "<hash>_<n>_LEFT_<variant>"), a fixed pad clips the top of every label off the canvas entirely
def primer_label_height_in(primers):
    max_len = max((len(p) for p in primers), default=0)
    return max_len * 0.11 + 0.3

# renders straight into the HTML report as a data URI instead of a standalone file — no separate
# .png sitting alongside the report duplicating what's already shown inline
def fig_to_base64(fig, **savefig_kwargs):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", **savefig_kwargs)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")

def plot_heatmap(accessions, primers, matrix, flag_matrix, row_labels):
    fig_w = max(8.0, len(primers) * 0.28)
    fig_h = max(3.0, len(accessions) * 0.5 + 1.5) + primer_label_height_in(primers)
    line_w = 1.1
    cbar_w = 0.35
    fig, (ax, line_ax, cbar_ax) = plt.subplots(
        1, 3, figsize=(fig_w + line_w + cbar_w, fig_h),
        gridspec_kw={"width_ratios": [fig_w, line_w, cbar_w], "wspace": 0.02},
    )

    # traffic light by flag, not raw mismatch magnitude: green = no mismatch, yellow = a
    # mismatch that doesn't meet the flagging thresholds (see Figure 1's caption for the exact
    # thresholds), red = this cell alone meets them
    cmap = ListedColormap(HEATMAP_COLORS)
    norm = BoundaryNorm([0, 1, 2, 3], cmap.N)
    im = ax.imshow(flag_matrix, cmap=cmap, norm=norm, aspect="auto")

    # annotate every mismatched cell with its actual count — color no longer encodes magnitude
    # (only the flagging condition does), so the number is the only place the count is shown
    for i, row in enumerate(matrix):
        for j, value in enumerate(row):
            if value:
                ax.text(j, i, str(value), ha="center", va="center", color="#212529", fontsize=7)

    ax.set_xticks(range(len(primers)))
    ax.set_xticklabels(primers, rotation=90, fontsize=7)
    ax.set_yticks(range(len(accessions)))
    ax.set_yticklabels(accessions)
    ax.tick_params(which="major", bottom=False, left=False)

    # match the country-line label and "mismatch count" colorbar label to the accession labels'
    # actual rendered font (size and family), rather than hardcoding a possibly-different value
    accession_font = ax.get_yticklabels()[0]
    accession_fontsize = accession_font.get_fontsize()
    accession_fontfamily = accession_font.get_fontfamily()[0]

    # country line — sits between the heatmap and the colorbar. One short line segment per
    # contiguous run of accessions sharing the same country (or "unknown"), each inset by a small
    # gap from its neighbors so adjacent countries read as visually distinct segments rather than
    # one unbroken line running the full height of the heatmap; same thickness as the heatmap's
    # own axis edges, labeled once at its midpoint rather than repeating on every row
    edge_lw = ax.spines["top"].get_linewidth()
    line_ax.set_xlim(0, 1)
    line_ax.set_ylim(ax.get_ylim())
    line_ax.axis("off")

    segment_gap = 0.15
    i = 0
    n = len(row_labels)
    while i < n:
        j = i
        while j + 1 < n and row_labels[j + 1] == row_labels[i]:
            j += 1
        top, bottom = i - 0.5 + segment_gap, j + 0.5 - segment_gap
        line_ax.plot([0.02, 0.02], [top, bottom], color="black", linewidth=edge_lw,
                     solid_capstyle="butt", clip_on=False)
        line_ax.text(0.09, (top + bottom) / 2, row_labels[i], ha="left", va="center",
                     fontsize=accession_fontsize, fontfamily=accession_fontfamily)
        i = j + 1

    cbar = fig.colorbar(im, cax=cbar_ax, ticks=[0.5, 1.5, 2.5])
    cbar.ax.set_yticklabels(["No mismatches", "Mismatches below thresholds", "Mismatches above thresholds"],
                             fontsize=accession_fontsize, fontfamily=accession_fontfamily)

    # fig.tight_layout() silently miscomputes margins here — an explicit cax from a 3-way subplot
    # grid isn't fully compatible with it (matplotlib warns "Axes ... not compatible with
    # tight_layout"), and the bottom margin it picks clips the top off every rotated primer-name
    # tick label. bbox_inches="tight" on save measures the actual rendered artist boxes instead
    # and isn't affected by that incompatibility.
    return fig_to_base64(fig, dpi=150, bbox_inches="tight")

def describe_dataset(country, subclade, fasta_name):
    if country and subclade:
        location = f"{country} ({subclade})"
    elif country:
        location = country
    elif subclade:
        location = subclade
    else:
        location = None

    if location and fasta_name:
        return f"{location}({fasta_name})"
    if location:
        return location
    if fasta_name:
        return fasta_name
    return "a user-provided dataset"

# columns: [(row dict key, displayed header label), ...] — lets the HTML table show a friendlier
# header than the underlying dict key without touching the TSV files that use those same keys
def rows_to_html_table(columns, rows):
    head = "".join(f"<th>{html.escape(label)}</th>" for _, label in columns)
    body = "\n".join(
        "<tr>" + "".join(f"<td>{html.escape(str(row.get(key, '')))}</td>" for key, _ in columns) + "</tr>"
        for row in rows
    )
    return f"<table>\n<thead><tr>{head}</tr></thead>\n<tbody>\n{body}\n</tbody>\n</table>"

# like rows_to_html_table, but tailored to RECOMMENDATIONS_COLUMNS: colors each <tr> by its
# "severity" (red/yellow, via REPORT_CSS) and merges the "Recommendation" cell across every row of
# a contiguous same-severity run into one rowspan'd cell, since that text is identical for every
# row in the run (build_recommendations_rows() already sorts red before yellow, so same-severity
# rows are always contiguous — no grouping/re-sort needed here)
def build_recommendations_table(rows):
    head = "<th>Primer Name</th><th>Issue</th><th>Recommendation</th>"
    body_rows = []
    i, n = 0, len(rows)
    while i < n:
        j = i
        while j + 1 < n and rows[j + 1]["severity"] == rows[i]["severity"]:
            j += 1
        group = rows[i:j + 1]
        recommendation_cell = f'<td rowspan="{len(group)}">{html.escape(group[0]["recommendation"])}</td>'
        for k, row in enumerate(group):
            cells = (f'<td>{html.escape(row["primer"])}</td>'
                     f'<td>{html.escape(row["issue"])}</td>')
            if k == 0:
                cells += recommendation_cell
            body_rows.append(f'<tr class="severity-{row["severity"]}">{cells}</tr>')
        i = j + 1
    body = "\n".join(body_rows)
    return f"<table>\n<thead><tr>{head}</tr></thead>\n<tbody>\n{body}\n</tbody>\n</table>"

REPORT_CSS = """
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
       max-width: 900px; margin: 2rem auto; padding: 0 1.5rem; color: #212529; line-height: 1.5; }
h1 { margin-bottom: 0.2rem; }
.date { color: #6c757d; margin-top: 0; margin-bottom: 2rem; }
h2 { border-bottom: 2px solid #e9ecef; padding-bottom: 0.3rem; margin-top: 2.5rem; }
.caption { color: #495057; font-size: 0.9rem; margin-top: 0.3rem; margin-bottom: 1rem; }
img { max-width: 100%; border: 1px solid #e9ecef; border-radius: 4px; }
table { border-collapse: collapse; width: 100%; font-size: 0.85rem; margin-top: 0.5rem; }
th, td { border: 1px solid #dee2e6; padding: 0.4rem 0.6rem; text-align: left; vertical-align: top; }
th { background: #f1f3f5; }
tr:nth-child(even) { background: #f8f9fa; }
tr.severity-red td { background-color: #EF9A9A; }
tr.severity-yellow td { background-color: #FFF59D; }
.empty-note { color: #6c757d; font-style: italic; }
.note { color: #6c757d; font-size: 0.85rem; margin-top: 1rem; }
"""

# self-contained (the heatmap is embedded as a base64 data URI, not a sibling file) — pulls
# together the outputs of this run() call into one summary page rather than leaving the user to
# open separate files
RECURRING_COLUMNS = [
    ("amplicon", "Amplicon"),
    ("primer.issues", "Primer Name"),
    ("positionInPrimer", "Position in Primer"),
    ("primerSeq", "PrimerRef"),
    ("PrimerAlt", "PrimerAlt"),
    ("type", "Type"),
    ("count", "Count"),
    ("mismatch_prop", "Mismatch Proportion (%)"),
]

def build_html_report(country, subclade, fasta_name, n_primers, n_sequences, threshold,
                       three_prime_threshold, other_threshold,
                       heatmap_image_b64,
                       suggestions):
    dataset = describe_dataset(country, subclade, fasta_name)
    report_date = datetime.now().strftime("%d/%m/%Y")

    if heatmap_image_b64:
        heatmap_block = f'<img src="data:image/png;base64,{heatmap_image_b64}" alt="Mismatch heatmap">'
    else:
        heatmap_block = '<p class="empty-note">No mismatches were found — no heatmap was generated.</p>'

    is_default_thresholds = (three_prime_threshold == DEFAULT_THREE_PRIME_THRESHOLD
                              and other_threshold == DEFAULT_SNP_THRESHOLD)
    # both the Overview blurb and Figure 1's caption depend on whether the run used the default
    # --3-prime/--snp thresholds or ones the user set explicitly
    overview_intro = (
        f"This report summarises primer mismatches identified across <strong>{n_primers} primers</strong> "
        f"evaluated against <strong>{n_sequences} sequences</strong> from <strong>{html.escape(dataset)}</strong>. "
        f"Potentially problematic primers shown in the mismatch heatmap (Figure 1) were flagged based on the "
    )
    if is_default_thresholds:
        overview_text = (
            overview_intro
            + f"default 3′ and SNP/indel thresholds (≥{DEFAULT_THREE_PRIME_THRESHOLD} 3′ mismatches or "
              f"≥{DEFAULT_SNP_THRESHOLD} SNP/indel mismatches) adapted from Huang et al. (2024) and "
              f"Bru et al. (2008). These primers may affect the sensitivity, specificity, and "
              f"amplification efficiency of PCR."
        )
        heatmap_caption = (
            f"Figure 1. Heatmap showing potentially problematic primers for each sequence based on "
            f"the default thresholds (≥{DEFAULT_THREE_PRIME_THRESHOLD} 3′ mismatches or "
            f"≥{DEFAULT_SNP_THRESHOLD} SNP/indel mismatches). The total number of mismatches for each "
            f"primer is indicated."
        )
    else:
        overview_text = (
            overview_intro
            + f"3′ and SNP/indel mismatch thresholds specified by the user "
              f"(≥{three_prime_threshold} 3′ or ≥{other_threshold} SNP/indel)."
        )
        heatmap_caption = (
            f"Figure 1. Heatmap showing potentially problematic primers for each sequence based on "
            f"the defined thresholds (≥{three_prime_threshold} 3′ mismatches or "
            f"≥{other_threshold} SNP/indel mismatches). The total number of mismatches for each primer is "
            f"indicated."
        )

    if suggestions:
        # "count" as count/total (n_sequences is the same accession universe find_recurring_
        # mismatches() computed mismatch_prop against), mismatch_prop as a bare percentage (e.g.
        # 80, not 0.80), and PrimerAlt as just the alternate base letters found in the sequences
        # (e.g. "AT" for a site with both A and T, most common first) — refSeq's own "base:prop"
        # breakdown dropped for this view, and "N" (ambiguous/unknown, not a real alternate base)
        # excluded. All display-only; the TSV summary keeps the raw values, N included.
        display_rows = [
            {**row,
             "count": f"{row['count']}/{n_sequences}",
             "mismatch_prop": f"{float(row['mismatch_prop']) * 100:.0f}",
             "PrimerAlt": "".join(base.split(":")[0] for base in row["refSeq"].split(", ")
                                   if base.split(":")[0] != "N")}
            for row in suggestions
        ]
        recurring_block = rows_to_html_table(RECURRING_COLUMNS, display_rows)
    else:
        recurring_block = (f'<p class="empty-note">No recurring mismatches were found (none shared by '
                            f'≥{threshold:.0%} of sequences).</p>')

    recommendations_rows = build_recommendations_rows(suggestions, n_sequences) if suggestions else []
    if recommendations_rows:
        recommendations_block = build_recommendations_table(recommendations_rows)
    else:
        recommendations_block = ('<p class="empty-note">No recommendations — no primers had a recurring '
                                  'mismatch shared by enough sequences to warrant redesign.</p>')

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Primer Mismatches Report</title>
<style>{REPORT_CSS}</style>
</head>
<body>
<h1>Primer Mismatches Report</h1>
<p class="date">{report_date}</p>

<h2>Overview</h2>
<p>{overview_text}</p>

<h2>Mismatch Heatmap</h2>
<p class="caption">{heatmap_caption}</p>
{heatmap_block}

<h2>Recurring Mismatches</h2>
<p class="caption">Table 1. Recurring mismatches observed across all sequences (mismatches at the
same position with the same alternate base occurring in ≥{threshold:.0%} of sequences).</p>
{recurring_block}

<h2>Recommendations</h2>
{recommendations_block}
<p class="note">These results have not been experimentally validated. Users should consider
relevant experimental evidence and the published literature before making any decisions.</p>

</body>
</html>
"""

def run(mismatches_dir, threshold, summary_path, report_path,
        three_prime_threshold, other_threshold, primer_bed=None,
        country_override=None, subclade_override=None, fasta_name=None):
    # an explicit --country/--subclade from the caller (the run_pipeline*.sh script, which
    # already knows exactly what the user asked for) is trusted over inferring it from the
    # mismatches folder name/metadata — that inference exists only for someone invoking this
    # script directly against a folder built outside the RABVPrime.py pipeline, where no such
    # explicit value is available. Both accept one or more names (RetrieveSequences.py's -c/-s
    # can retrieve several combined into one dataset), joined for display via join_names().
    if country_override or subclade_override:
        country = join_names(country_override) if country_override else None
        subclade = join_names(subclade_override) if subclade_override else None
    else:
        country, subclade = resolve_country_subclade(mismatches_dir)

    degeneracy_fieldnames = ["amplicon", "primer.issues", "positionInPrimer", "refSeq", "primerSeq",
                              "primerConsensus", "type", "count", "mismatch_prop"]

    print("Generating a primer mismatch heatmap…")
    accessions, primers, matrix, flag_matrix, row_labels, has_mismatches = build_heatmap_matrix(
        mismatches_dir, three_prime_threshold, other_threshold, primer_bed)
    heatmap_image_b64 = (plot_heatmap(accessions, primers, matrix, flag_matrix, row_labels)
                          if has_mismatches else None)

    # independent of the heatmap above — a scatter of one-off, non-recurring mismatches across
    # many primers/positions still gets visualized there even if nothing recurs here
    print("\nIdentifying recurring primer mismatches...")
    suggestions, _ = find_recurring_mismatches(mismatches_dir, threshold)

    if suggestions:
        with open(summary_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=degeneracy_fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows({key: row[key] for key in degeneracy_fieldnames} for row in suggestions)
        print(f"✓ Saved recurring mismatch summary to: {summary_path}")
    else:
        if summary_path.exists():
            summary_path.unlink()  # stale file from a previous run with different thresholds
        print("No recurring primer mismatches found - skipped the summary file")

    print("\nGenerating an HTML report…")
    report_html = build_html_report(
        country, subclade, fasta_name, len(primers), len(accessions), threshold,
        three_prime_threshold, other_threshold,
        heatmap_image_b64,
        suggestions,
    )
    report_path.write_text(report_html, encoding="utf-8")
    print(f"✓ Saved HTML report (with embedded heatmap) to: {report_path.resolve()}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Identify recurring mismatches across all sequences by searching for mismatches "
                     "that occur at the same position in at least 20% of sequences by default "
                     "(configurable with --recurring-threshold). Generate a mismatch heatmap, with "
                     "mismatches flagged according to the --3-prime and --snp thresholds, along with a "
                     "table of recurring mismatches and a warning summary highlighting potentially "
                     "problematic primer sites."
    )
    parser.add_argument(
        "-i", "--mismatches-dir",
        required=True,
        help="The <name>_mismatches/ folder produced by PrimerSequenceMismatch.py",
    )
    parser.add_argument(
        "--recurring-threshold",
        type=float,
        default=0.20,
        help="Minimum proportion of sequences (0-1) that must share a mismatch to be flagged "
             "(default: 0.20)",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output path prefix for the summary TSV and HTML report (default: <mismatches_dir>/<name>)",
    )
    parser.add_argument(
        "-b", "--bed",
        dest="primer_bed",
        default=None,
        help="Primer BED file containing the designed primers, either generated by PrimerDesign.py "
             "or provided by the user. This file allows primers with no mismatches to be identified "
             "and displayed on the mismatch heatmap.",
    )
    parser.add_argument(
        "--country",
        nargs="+",
        default=None,
        help="Country name(s) to show in the report's Overview section "
             "(default: inferred from the mismatches folder)",
    )
    parser.add_argument(
        "--subclade",
        nargs="+",
        default=None,
        help="Subclade name(s) to show in the report's Overview section "
             "(default: inferred from the mismatches folder)",
    )
    parser.add_argument(
        "--fasta-name",
        dest="fasta_name",
        default=None,
        help="Original FASTA filename to show in the report's Overview section",
    )
    parser.add_argument(
        "--3-prime",
        dest="three_prime_threshold",
        type=int,
        default=DEFAULT_THREE_PRIME_THRESHOLD,
        help=f"Flag a primer red in the heatmap if it has at least this many 3' mismatches against "
             f"one sequence (default: {DEFAULT_THREE_PRIME_THRESHOLD})",
    )
    parser.add_argument(
        "--snp",
        dest="other_threshold",
        type=int,
        default=DEFAULT_SNP_THRESHOLD,
        help=f"Flag a primer red in the heatmap if it has at least this many SNP/indel "
             f"mismatches against one sequence (default: {DEFAULT_SNP_THRESHOLD})",
    )
    args = parser.parse_args()

    mismatches_dir = Path(args.mismatches_dir)
    if not mismatches_dir.is_dir():
        print(f"✗ Not a directory: {mismatches_dir}")
        sys.exit(1)

    primer_bed = Path(args.primer_bed) if args.primer_bed else None
    if primer_bed and not primer_bed.is_file():
        print(f"✗ Not a file: {primer_bed}")
        sys.exit(1)

    if args.output:
        prefix = Path(args.output)
    else:
        name = mismatches_dir.name[:-len("_mismatches")] if mismatches_dir.name.endswith("_mismatches") \
            else mismatches_dir.name
        prefix = mismatches_dir / name

    summary_path = prefix.with_name(prefix.name + "_recurring_mismatches_summary.tsv")
    report_path = prefix.with_name(prefix.name + "_report.html")

    run(mismatches_dir, args.recurring_threshold, summary_path, report_path,
        args.three_prime_threshold, args.other_threshold, primer_bed,
        args.country, args.subclade, args.fasta_name)
