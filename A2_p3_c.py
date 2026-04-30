import math
import pandas as pd

READS = 10000
READ_LENGTH = 99
KMER = 25
P_FP_VALUES = [0.001, 0.01, 0.05, 0.1, 0.5]

# use derived lower bound on m to minimize false positive rate
def bloom_filter_size_bits(n_items, p_fp):
    return -(n_items * math.log(p_fp)) / (math.log(2) ** 2)

def main():
    n_kmers = READS * (READ_LENGTH - KMER + 1)
    rows = []

    for p_fp in P_FP_VALUES:
        m_bits = bloom_filter_size_bits(n_kmers, p_fp)
        rows.append({
            'P_FP': p_fp,
            'min_BF_size_kB': m_bits / 8 / 1000,
        })

    df = pd.DataFrame(rows)
    print(df.to_string(index=False, formatters={
        'min_BF_size_kB': '{:.3f}'.format,
    }))

    df.to_csv('BF_size_error.csv', index=False, float_format='%.3f')


if __name__ == '__main__':
    main()