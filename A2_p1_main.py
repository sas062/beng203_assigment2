from A2_p1_helpers import *
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# Part (a) 
segments = read_segments("ecDNA_segments.bed")
genes = read_genes("ecDNA_genes.bed")
# plot_ecdna_cartoon(segments, output_path="ecDNA_cartoon.png")

# Part (b)

exons = load_refgene_exons("UCSC_genome_table.csv", genes)
gene_expr = read_expression_bg("ecDNA_expression.bg")

# gene_scores = {}
# for gene, exon_intervals in exons.items():
#     gene_scores[gene] = interval_average_signal(exon_intervals, gene_expr)

# top_gene = max(gene_scores, key=gene_scores.get)
# print("Highest expression gene:", top_gene)
# print("Expression score:", gene_scores[top_gene])

# for gene in ["LPIN1", "MYCN", "TRIB2", "NTSR2", "GREB1"]:
#     print(f"{gene} expression score:", gene_scores[gene])


# Part (c)

lpin1_bounds = genes["LPIN1"]
transcript = find_best_refgene_transcript(
    "UCSC_genome_table.csv", "LPIN1", 
    lpin1_bounds["chrom"], lpin1_bounds["start"], lpin1_bounds["end"]
)

exon_spans = split_exons(transcript["exonStarts"], transcript["exonEnds"])
exon_spans = merge_intervals(exon_spans)

transcript_bins = build_transcript_bins(exon_spans, bin_size=50)
bin_results = []
higher_count = 0
lower_count = 0

all_bin_counts = [bin_expression_count(b, gene_expr) for b in transcript_bins]
lambda_est = np.mean(all_bin_counts)  # this is λ — the Poisson null rate

for i, (bin_spans, obs) in enumerate(zip(transcript_bins, all_bin_counts)):
    obs_int = int(round(obs))

    p_low  = stats.poisson.cdf(obs_int, lambda_est)          # P(X <= k)
    p_high = 1.0 - stats.poisson.cdf(obs_int - 1, lambda_est)  # P(X >= k)
    p_val  = 2 * min(p_low, p_high)

    bin_results.append((i * 50, obs, p_val))
    if p_high < 0.05:
        higher_count += 1
    elif p_low < 0.05:
        lower_count += 1
    

pos, scores, pvals = zip(*bin_results)
n_signif = sum(1 for p in pvals if p < 0.05)

plt.figure(figsize=(10, 5))
plt.plot(pos, scores, label="Mean expression", color="black", lw=1)
sig_mask = np.array(pvals) < 0.05
plt.scatter(np.array(pos)[sig_mask], np.array(scores)[sig_mask], 
            color="red", s=10, label="Significantly different")

plt.xlabel("Position along LPIN1 transcript (exons)")
plt.ylabel("Mean expression")
plt.title(f"LPIN1 Transcript Expression ({n_signif} significant bins)")
plt.legend()
plt.tight_layout()
plt.savefig("LPIN1_transcript_profile.png")

print(f"--- Poisson Analysis Results (lambda={lambda_est:.2f}) ---")
print(f"Number of significant higher-expression bins: {higher_count}")
print(f"Number of significant lower-expression bins: {lower_count}")
print(f"Total significant bins: {higher_count + lower_count}")
