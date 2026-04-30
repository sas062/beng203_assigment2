import math
from pathlib import Path
import mmh3

READ_FASTQ = Path("BF-database/SRR873426.R1_sub.fastq.txt")
QUERY_FILE = Path("BF-Queries.txt")

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

def get_kmers(sequence, k):
    if len(sequence) < k:
        return []
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def hash_bloom_filter(obj, m, h):
    for seed in range(h):
        yield mmh3.hash(obj, seed=seed, signed=False) % m

def human_query(query, bit_array, m_bits, h, k):
    for kmer in get_kmers(query, k):
        for i in hash_bloom_filter(kmer, m_bits, h):
            if bit_array[i] == 0:
                return False
    return True

def query_match_fraction(query, bit_array, m_bits, h, k):
    kmers = get_kmers(query, k)
    if not kmers:
        return 0.0

    matches = 0
    for kmer in kmers:
        present = True
        for i in hash_bloom_filter(kmer, m_bits, h):
            if bit_array[i] == 0:
                present = False
                break
        if present:
            matches += 1

    return matches / len(kmers)

def main():
    queries = read_queries(QUERY_FILE)
    sequences = read_sequences(READ_FASTQ)
    n_kmers = len(sequences) * (len(sequences[0]) - KMER + 1)
    m_bits = bloom_filter_size_bits(n_kmers, P_FP)
    h = optimal_hashes(m_bits, n_kmers)

    # Initialize Bloom filter bit array
    bit_array = [0] * m_bits
    for sequence in sequences:
        for kmer in get_kmers(sequence, KMER):
            for i in hash_bloom_filter(kmer, m_bits, h):
                bit_array[i] = 1

    for idx, query in enumerate(queries, start=1):
        frac = query_match_fraction(query, bit_array, m_bits, h, KMER)
        bool_result = human_query(query, bit_array, m_bits, h, KMER)
        print(f"query{idx} is {'human' if bool_result else 'not human'}. Match fraction = {frac:.3f}")

if __name__ == '__main__':
    main()