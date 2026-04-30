import math
from pathlib import Path
import mmh3

READ_FASTQ = Path("BF-database/SRR873426.R1_sub.fastq.txt")
QUERY_FILE = Path("BF-Queries.txt")

READS = 10000
KMER = 25
P_FP = 0.01

# use derived lower bound on m to minimize false positive rate
def bloom_filter_size_bits(n_items, p_fp):
    return math.ceil(-(n_items * math.log(p_fp)) / (math.log(2) ** 2))

# use derived h to minimize false positive rate for given m and n
def optimal_hashes(m_bits, n_items):
    return math.ceil((m_bits / n_items) * math.log(2))

def read_queries(query_file):
    queries = []
    with query_file.open() as f:
        for line in f:
            line = line.strip()
            if line.startswith('>') or not line:
                continue
            queries.append(line)
    return queries

def read_sequences(fastq_file):
    sequences = []
    with fastq_file.open() as f:
        for i, line in enumerate(f, start=1):
            if i % 4 == 2:  # Sequence line in FASTQ format
                sequences.append(line.strip())
    return sequences

def main():
    queries = read_queries(QUERY_FILE)
    sequences = read_sequences(READ_FASTQ)
    n_kmers = READS * (len(sequences[0]) - KMER + 1)
    m_bits = bloom_filter_size_bits(n_kmers, P_FP)
    h = optimal_hashes(m_bits, n_kmers)
    print(f"Bloom filter size (bits): {m_bits:d}")
    print(f"Optimal number of hash functions: {h:d}")
    print(f"Num sequences: {len(sequences)}")
    print(f"Sequence length: {len(sequences[0])}")


if __name__ == '__main__':
    main()