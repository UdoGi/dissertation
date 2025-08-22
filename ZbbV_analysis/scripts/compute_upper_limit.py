import math

import pandas as pd
from scipy.stats import skellam, poisson


def compute_pval_from_rate(rate, genome_len, ml_subs_penz,
                           ml_subs_bad, observed_subs, observation_time ):
    # the expected delta if the rate would be true
    delta = rate * observation_time * genome_len

    # iterating through possible placements of d_A and d_B
    d_B = 0
    d_A = int(d_B + delta)
    l = []
    while d_B < ml_subs_penz:
        l.append([observed_subs, d_A, d_B,
                  skellam.cdf(observed_subs, d_A, d_B),
                  poisson.cdf(ml_subs_bad, d_A),
                  poisson.cdf(ml_subs_penz, d_B)])
        d_B += 1
        d_A += 1

    df = pd.DataFrame(l, columns=['observed', 'd_A', 'd_B',
                                  'skellam', 'poissonA', 'poissonB'])
    df['pval_A'] = df.poissonA.apply(lambda x: 1 - x if x > 0.5 else x)
    df['pval_B'] = df.poissonB.apply(lambda x: 1 - x if x > 0.5 else x)
    df['pvaldiff'] = (df.pval_A - df.pval_B)

    # find out when the pvals of the distr of A and B cross
    previous_row = df.iloc[0]
    for i, row in df.iloc[1:].iterrows():
        # we crossed the point where both p-values are most similar
        if row.pvaldiff < 0:
            crossing_row = row if abs(row.pvaldiff) < abs(previous_row.pvaldiff) else previous_row
            break
        previous_row = row

    return crossing_row


if __name__ == '__main__':
    genome_len = 6255
    ml_dist_penz = .013
    ml_subs_penz = int(ml_dist_penz * genome_len)
    ml_dist_bad = 0.0157
    ml_subs_bad = int(ml_dist_bad * genome_len)
    observed_subs = math.ceil((ml_dist_bad - ml_dist_penz) * genome_len)
    observation_time = 93

    for x in range(4, 9):
        rate = x * 10 ** -5
        r = compute_pval_from_rate(rate, genome_len, ml_subs_penz,
                                   ml_subs_bad, observed_subs, observation_time)
        print(f"rate {rate:.2e} d_A {r.d_A}, d_B {r.d_B}, delta {r.d_A - r.d_B} "
              f"pval_delta {r.skellam:.2e}, pval_A {r.pval_A:.2e}, pval_B {r.pval_B:.2e}")

