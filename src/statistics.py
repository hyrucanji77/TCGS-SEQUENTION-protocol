#!/usr/bin/env python3
"""
TCGS-SEQUENTION: Statistical Functions
=======================================

Statistical functions for the three-gate protocol.

Author: Henry Arellano-Peña
Framework: TCGS-SEQUENTION
License: CC BY 4.0
"""

import numpy as np
from scipy import stats
from typing import List, Tuple, Optional


def compute_synchrony_statistic(times: np.ndarray, 
                                 delta_t: float,
                                 reference: str = 'median') -> Tuple[float, float]:
    """
    Compute the S* synchrony statistic.
    
    S* = (1/m) * sum_i I(|tau_i - tau_ref| <= delta_t)
    
    Parameters
    ----------
    times : array-like
        Emergence times for each population
    delta_t : float
        Synchrony window size
    reference : str
        Reference point: 'median' (default) or 'mean'
        
    Returns
    -------
    s_stat : float
        S* statistic (fraction synchronized)
    ref_time : float
        Reference time used
    """
    times = np.asarray(times)
    
    if reference == 'median':
        ref_time = np.median(times)
    elif reference == 'mean':
        ref_time = np.mean(times)
    else:
        raise ValueError(f"Unknown reference: {reference}")
    
    within_window = np.abs(times - ref_time) <= delta_t
    s_stat = np.mean(within_window)
    
    return float(s_stat), float(ref_time)


def permutation_test_synchrony(times: np.ndarray,
                                delta_t: float,
                                n_permutations: int = 10000,
                                seed: Optional[int] = None) -> Tuple[float, np.ndarray]:
    """
    Compute p-value for synchrony via permutation test.
    
    Parameters
    ----------
    times : array-like
        Observed emergence times
    delta_t : float
        Synchrony window size
    n_permutations : int
        Number of permutations
    seed : int, optional
        Random seed for reproducibility
        
    Returns
    -------
    p_value : float
        Permutation p-value
    null_distribution : ndarray
        Null distribution of S* values
    """
    if seed is not None:
        np.random.seed(seed)
    
    times = np.asarray(times)
    n = len(times)
    
    # Observed statistic
    observed_s, _ = compute_synchrony_statistic(times, delta_t)
    
    # Time range for null sampling
    time_min = np.min(times)
    time_max = np.max(times)
    time_range = time_max - time_min
    
    if time_range == 0:
        return 0.0, np.zeros(n_permutations)  # Perfect synchrony
    
    # Generate null distribution
    null_s = np.zeros(n_permutations)
    for i in range(n_permutations):
        # Random uniform times over expanded range
        random_times = np.random.uniform(
            time_min - time_range * 0.5,
            time_max + time_range * 0.5,
            n
        )
        null_s[i], _ = compute_synchrony_statistic(random_times, delta_t)
    
    # p-value: fraction of null >= observed
    p_value = np.mean(null_s >= observed_s)
    
    return float(p_value), null_s


def bootstrap_confidence_interval(times: np.ndarray,
                                   delta_t: float,
                                   n_bootstrap: int = 10000,
                                   ci_level: float = 0.95,
                                   seed: Optional[int] = None) -> Tuple[float, float, float]:
    """
    Compute bootstrap confidence interval for S*.
    
    Parameters
    ----------
    times : array-like
        Emergence times
    delta_t : float
        Synchrony window size
    n_bootstrap : int
        Number of bootstrap samples
    ci_level : float
        Confidence level (e.g., 0.95 for 95% CI)
    seed : int, optional
        Random seed
        
    Returns
    -------
    s_stat : float
        Point estimate of S*
    ci_lower : float
        Lower confidence bound
    ci_upper : float
        Upper confidence bound
    """
    if seed is not None:
        np.random.seed(seed)
    
    times = np.asarray(times)
    n = len(times)
    
    # Point estimate
    s_stat, _ = compute_synchrony_statistic(times, delta_t)
    
    # Bootstrap
    boot_s = np.zeros(n_bootstrap)
    for i in range(n_bootstrap):
        idx = np.random.choice(n, size=n, replace=True)
        boot_s[i], _ = compute_synchrony_statistic(times[idx], delta_t)
    
    # Confidence interval
    alpha = 1 - ci_level
    ci_lower = np.percentile(boot_s, 100 * alpha / 2)
    ci_upper = np.percentile(boot_s, 100 * (1 - alpha / 2))
    
    return float(s_stat), float(ci_lower), float(ci_upper)


def transfer_entropy(x: np.ndarray, 
                     y: np.ndarray,
                     lag: int = 1,
                     bins: int = 10) -> float:
    """
    Estimate transfer entropy TE(X -> Y).
    
    TE measures the reduction in uncertainty about Y's future
    when knowing X's past, beyond knowing Y's own past.
    
    Parameters
    ----------
    x : array-like
        Source time series
    y : array-like
        Target time series
    lag : int
        Time lag for conditioning
    bins : int
        Number of bins for discretization
        
    Returns
    -------
    te : float
        Transfer entropy estimate
    """
    x = np.asarray(x)
    y = np.asarray(y)
    
    if len(x) != len(y):
        raise ValueError("Time series must have equal length")
    
    n = len(x) - lag
    if n < 10:
        return 0.0  # Insufficient data
    
    # Discretize
    x_bins = np.digitize(x, np.linspace(x.min(), x.max(), bins))
    y_bins = np.digitize(y, np.linspace(y.min(), y.max(), bins))
    
    # Build joint distributions
    y_future = y_bins[lag:]
    y_past = y_bins[:-lag]
    x_past = x_bins[:-lag]
    
    # Estimate entropies using histograms
    def entropy(data):
        counts = np.bincount(data.flatten())
        probs = counts[counts > 0] / counts.sum()
        return -np.sum(probs * np.log2(probs + 1e-10))
    
    # Joint entropy estimates
    h_y_future_y_past = entropy(y_future * bins + y_past)
    h_y_future_y_past_x_past = entropy(y_future * bins**2 + y_past * bins + x_past)
    h_y_past = entropy(y_past)
    h_y_past_x_past = entropy(y_past * bins + x_past)
    
    # Transfer entropy
    te = h_y_future_y_past - h_y_future_y_past_x_past - h_y_past + h_y_past_x_past
    
    return max(0.0, float(te))  # TE is non-negative


def compute_fst(allele_freqs: np.ndarray) -> float:
    """
    Compute Weir-Cockerham FST estimator.
    
    Parameters
    ----------
    allele_freqs : ndarray
        Allele frequencies, shape (n_populations, n_loci)
        
    Returns
    -------
    fst : float
        Global FST estimate
    """
    allele_freqs = np.asarray(allele_freqs)
    
    if allele_freqs.ndim == 1:
        allele_freqs = allele_freqs.reshape(-1, 1)
    
    n_pops, n_loci = allele_freqs.shape
    
    # Mean frequency
    p_mean = np.mean(allele_freqs, axis=0)
    
    # Variance among populations
    var_p = np.var(allele_freqs, axis=0, ddof=1)
    
    # Expected heterozygosity
    h_t = 2 * p_mean * (1 - p_mean)
    h_s = 2 * allele_freqs * (1 - allele_freqs)
    h_s_mean = np.mean(h_s, axis=0)
    
    # FST per locus
    with np.errstate(divide='ignore', invalid='ignore'):
        fst_loci = np.where(h_t > 0, (h_t - h_s_mean) / h_t, 0)
    
    # Global FST
    fst = np.nanmean(fst_loci)
    
    return float(fst)


def compute_ibd(genotypes1: np.ndarray, genotypes2: np.ndarray) -> float:
    """
    Estimate identity-by-descent (IBD) between two individuals.
    
    Simple estimator based on allele sharing.
    
    Parameters
    ----------
    genotypes1, genotypes2 : ndarray
        Genotype vectors (0, 1, 2 encoding)
        
    Returns
    -------
    ibd : float
        IBD estimate (0 to 1)
    """
    g1 = np.asarray(genotypes1)
    g2 = np.asarray(genotypes2)
    
    if len(g1) != len(g2):
        raise ValueError("Genotype vectors must have equal length")
    
    # Allele sharing score
    sharing = 1 - np.abs(g1 - g2) / 2
    ibd = np.mean(sharing)
    
    return float(ibd)


def fisher_combine_pvalues(pvalues: List[float]) -> float:
    """
    Combine p-values using Fisher's method.
    
    Parameters
    ----------
    pvalues : list of float
        Individual p-values
        
    Returns
    -------
    combined_p : float
        Combined p-value
    """
    pvalues = np.asarray(pvalues)
    pvalues = pvalues[pvalues > 0]  # Remove zeros
    
    if len(pvalues) == 0:
        return 1.0
    
    # Fisher's statistic
    chi2_stat = -2 * np.sum(np.log(pvalues))
    df = 2 * len(pvalues)
    
    # Combined p-value
    combined_p = 1 - stats.chi2.cdf(chi2_stat, df)
    
    return float(combined_p)


def compute_aic(log_likelihood: float, n_parameters: int) -> float:
    """
    Compute Akaike Information Criterion.
    
    AIC = 2k - 2ln(L)
    
    Parameters
    ----------
    log_likelihood : float
        Log-likelihood of the model
    n_parameters : int
        Number of estimated parameters
        
    Returns
    -------
    aic : float
        AIC value
    """
    return 2 * n_parameters - 2 * log_likelihood


def compute_bic(log_likelihood: float, n_parameters: int, n_observations: int) -> float:
    """
    Compute Bayesian Information Criterion.
    
    BIC = k*ln(n) - 2ln(L)
    
    Parameters
    ----------
    log_likelihood : float
        Log-likelihood of the model
    n_parameters : int
        Number of estimated parameters
    n_observations : int
        Number of observations
        
    Returns
    -------
    bic : float
        BIC value
    """
    return n_parameters * np.log(n_observations) - 2 * log_likelihood
