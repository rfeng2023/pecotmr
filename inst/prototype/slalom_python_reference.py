#!/usr/bin/env python3
"""
Reference Python implementation of SLALOM core algorithm.
Extracted from https://github.com/mkanai/slalom/blob/master/slalom.py
Used for validating the pecotmr R implementation.
"""
import sys
import numpy as np
import scipy as sp
import scipy.stats
import pandas as pd


def abf(beta, se, W=0.04):
    """Approximate Bayes Factor - exact copy from slalom.py"""
    z = beta / se
    V = se**2
    r = W / (W + V)
    lbf = 0.5 * (np.log(1 - r) + (r * z**2))
    denom = sp.special.logsumexp(lbf)
    prob = np.exp(lbf - denom)
    return lbf, prob


def get_cs(variant_idx, prob, coverage=0.95):
    """Credible set - adapted from slalom.py (uses indices instead of variant names)"""
    ordering = np.argsort(prob)[::-1]
    idx = np.where(np.cumsum(prob[ordering]) > coverage)[0][0]
    cs = ordering[: (idx + 1)]
    return cs


def slalom_python(z_scores, LD_matrix, se=None,
                  abf_prior_variance=0.04,
                  r2_threshold=0.6,
                  nlog10p_dentist_s_threshold=4.0):
    """
    Run SLALOM algorithm matching the original Python implementation.

    Parameters
    ----------
    z_scores : array-like
        Z-scores for each variant
    LD_matrix : 2D array
        LD correlation matrix (r, not r²)
    se : array-like, optional
        Standard errors. Default: all ones.

    Returns
    -------
    dict with keys: prob, pvalue, lbf, nlog10p_dentist_s, outliers,
                    lead_idx, n_r2, n_dentist_s_outlier, max_pip,
                    cs_95, cs_99
    """
    z = np.array(z_scores, dtype=float)
    R = np.array(LD_matrix, dtype=float)
    n = len(z)

    if se is None:
        se = np.ones(n)
    se = np.array(se, dtype=float)
    beta = z * se  # When se=1, beta = z

    # ABF
    lbf, prob = abf(beta, se, W=abf_prior_variance)

    # Credible sets
    cs_95 = get_cs(np.arange(n), prob, coverage=0.95)
    cs_99 = get_cs(np.arange(n), prob, coverage=0.99)

    # Lead variant selection via one-sided p-value (matching slalom.py notebook)
    # Python original: p = stats.norm.cdf(z), then idxmin
    pvalue = scipy.stats.norm.cdf(z)
    lead_idx = np.argmin(pvalue)

    # DENTIST-S statistic
    lead_z = z[lead_idx]
    r_lead = R[:, lead_idx]
    r2_lead = r_lead ** 2

    t_dentist_s = (z - r_lead * lead_z) ** 2 / (1 - r2_lead)
    t_dentist_s = np.where(t_dentist_s < 0, np.inf, t_dentist_s)
    t_dentist_s[lead_idx] = np.nan  # Explicitly set lead to NaN

    nlog10p_dentist_s = scipy.stats.chi2.logsf(t_dentist_s, df=1) / -np.log(10)

    # Outlier flags
    outliers = (r2_lead > r2_threshold) & (nlog10p_dentist_s > nlog10p_dentist_s_threshold)

    # Summary
    n_r2 = np.sum(r2_lead > r2_threshold)
    n_dentist_s_outlier = np.sum(outliers)
    max_pip = np.max(prob)

    return {
        'prob': prob,
        'pvalue': pvalue,
        'lbf': lbf,
        'nlog10p_dentist_s': nlog10p_dentist_s,
        'outliers': outliers,
        'lead_idx': lead_idx,
        'n_r2': int(n_r2),
        'n_dentist_s_outlier': int(n_dentist_s_outlier),
        'max_pip': float(max_pip),
        'cs_95': cs_95.tolist(),
        'cs_99': cs_99.tolist(),
    }


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: slalom_python_reference.py <z_scores.csv> <LD_matrix.csv> <output.csv>")
        sys.exit(1)

    z_file, ld_file, out_file = sys.argv[1], sys.argv[2], sys.argv[3]

    z_scores = np.loadtxt(z_file)
    LD_matrix = np.loadtxt(ld_file, delimiter=',')

    result = slalom_python(z_scores, LD_matrix)

    # Save per-variant results
    df = pd.DataFrame({
        'prob': result['prob'],
        'pvalue': result['pvalue'],
        'lbf': result['lbf'],
        'nlog10p_dentist_s': result['nlog10p_dentist_s'],
        'outliers': result['outliers'].astype(int),
    })
    df.to_csv(out_file, index=False)

    # Print summary
    print(f"lead_idx={result['lead_idx']}")
    print(f"n_r2={result['n_r2']}")
    print(f"n_dentist_s_outlier={result['n_dentist_s_outlier']}")
    print(f"max_pip={result['max_pip']:.15e}")
    print(f"cs_95_size={len(result['cs_95'])}")
    print(f"cs_99_size={len(result['cs_99'])}")
