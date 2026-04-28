import csv
import numpy as np
from scipy import stats
from collections import defaultdict
import matplotlib.pyplot as plt
import math

# Part (a) helper functions

def read_segments(path):
    segments = []
    with open(path) as f:
        for idx, line in enumerate(f):
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, strand = line.strip().split()[:4]
            segments.append({
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "segment_order": idx + 1,
            })
    return segments

def read_genes(path):
    genes = {}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, gene = line.strip().split()[:4]
            genes[gene] = {
                "chrom": chrom,
                "start": int(start),
                "end": int(end)
            }
    return genes

def plot_ecdna_cartoon(segments, output_path=None):
    
    fig, ax = plt.subplots(figsize=(12, 2))
    y_levels = [1 if i % 2 == 0 else 0 for i in range(len(segments))]

    for i, seg in enumerate(segments):
        y = y_levels[i]
        color = "tab:blue" if seg["strand"] == "+" else "tab:red"
        ax.plot([seg["start"], seg["end"]], [y, y], color=color,
                linewidth=8, solid_capstyle="butt")

        if i > 0:
            prev = segments[i - 1]
            prev_y = y_levels[i - 1]
            ax.plot([prev["end"], seg["start"]], [prev_y, y], color="gray", lw=1)

        label = f"S{i + 1}"
        ax.text((seg["start"] + seg["end"]) / 2, y + 0.12, label,
                ha="center", va="bottom", fontsize=8)

    if len(segments) > 1:
        first = segments[0]
        last = segments[-1]
        first_y = y_levels[0]
        last_y = y_levels[-1]
        ax.plot([last["end"], first["start"]], [last_y, first_y], color="gray", lw=1)

    ax.set_xlabel("Reference genomic coordinate (hg38)")
    ax.set_yticks([])
    ax.set_title("ecDNA segments and genes")
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.ticklabel_format(style="plain", axis="x")

    min_start = min(seg["start"] for seg in segments)
    max_end = max(seg["end"] for seg in segments)
    padding = int((max_end - min_start) * 0.03) or 1
    ax.set_xlim(min_start - padding, max_end + padding)
    ax.set_ylim(-0.5, 1.3)
    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=200)
        print(f"Saved ecDNA cartoon to {output_path}")
    plt.close(fig)


# Part (b) helper functions

def load_refgene_exons(refgene_csv, ecdna_genes):
    gene_exons = defaultdict(list)
    with open(refgene_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene = row["name2"]
            if gene not in set(ecdna_genes.keys()):
                continue

            cds_start = int(row["cdsStart"])
            cds_end = int(row["cdsEnd"])
            if cds_start >= cds_end:
                continue
            chrom = row["chrom"]
            if gene not in ecdna_genes or chrom != ecdna_genes[gene]["chrom"]:
                continue

            exons = split_exons(row["exonStarts"], row["exonEnds"])
            ginfo = ecdna_genes[gene]
            for s, e in exons:
                ov = overlap_len(s, e, ginfo["start"], ginfo["end"])
                if ov > 0:
                    gene_exons[gene].append((max(s, ginfo["start"]), min(e, ginfo["end"])))
    for gene in list(gene_exons):
        gene_exons[gene] = merge_intervals(gene_exons[gene])
    return gene_exons

def split_exons(starts_str, ends_str):
    starts = [int(x) for x in starts_str.strip().rstrip(",").split(",") if x != ""]
    ends = [int(x) for x in ends_str.strip().rstrip(",").split(",") if x != ""]
    return list(zip(starts, ends))

def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [tuple(x) for x in merged]


def overlap_len(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))

def read_expression_bg(path):
    rows = []
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            if len(line.strip().split()) < 4:
                continue
            chrom, s, e, v = line.rstrip().split()[:4]
            if chrom == "chr2":
                rows.append((int(s), int(e), float(v)))
    return rows

def interval_average_signal(intervals, gene_expr):
    total_weighted = 0.0
    total_bases = 0

    for s, e in intervals:
        for bg_s, bg_e, val in gene_expr:
            ov = overlap_len(s, e, bg_s, bg_e)
            if ov > 0:
                total_weighted += val * ov
                total_bases += ov

    return total_weighted / total_bases if total_bases > 0 else 0.0

# Part (c) helper functions

def find_best_refgene_transcript(refgene_csv, gene_name, chrom, gene_start, gene_end):
    best = None
    best_score = None
    with open(refgene_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["name2"] != gene_name or row["chrom"] != chrom:
                continue
            tx_start = int(row["txStart"])
            tx_end = int(row["txEnd"])
            overlap = overlap_len(tx_start, tx_end, gene_start, gene_end)
            if overlap <= 0:
                continue
            score = (-overlap, abs(tx_start - gene_start), abs(tx_end - gene_end), -(tx_end - tx_start))
            if best_score is None or score < best_score:
                best_score = score
                best = row
    return best


def build_transcript_bins(exon_spans, bin_size=50):
    transcript_bins = []
    current_bin = []
    current_len = 0

    for start, end in exon_spans:
        while start < end:
            take = min(bin_size - current_len, end - start)
            current_bin.append((start, start + take))
            current_len += take
            start += take
            if current_len == bin_size:
                transcript_bins.append(current_bin)
                current_bin = []
                current_len = 0

    if current_bin:
        transcript_bins.append(current_bin)
    return transcript_bins


def bin_expression_count(bin_spans, gene_expr):
    total = 0.0
    for s, e in bin_spans:
        for bg_s, bg_e, val in gene_expr:
            ov = overlap_len(s, e, bg_s, bg_e)
            if ov > 0:
                total += val * (ov / (bg_e - bg_s))  # proportional share of the 50bp bedgraph bin
    return total
