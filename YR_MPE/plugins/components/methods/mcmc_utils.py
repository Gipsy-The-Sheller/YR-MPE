#!/usr/bin/env python3
"""
MCMC utility functions for calculating Effective Sample Size (ESS) and 
Highest Posterior Density (HPD) intervals.
"""

import numpy as np


def calculate_ESS(samples):
    """
    Calculate Effective Sample Size (ESS) using numpy for acceleration.
    Implements Geyer's initial monotone sequence criterion to match Tracer's implementation.
    
    The basic ess (N_eff) diagnostic is computed by:
    N_eff = N / tau_hat
    tau_hat = -1 + 2 * sum_{k=0}^K P_k
    where P_k = rho_{2k} + rho_{2k+1}
    and K is the last integer for which P_k is still positive and monotone decreasing.
    
    Parameters:
    samples: array-like, MCMC samples
    
    Returns:
    float: Effective Sample Size
    """
    if len(samples) < 4:  # Need at least 4 samples for valid ESS calculation
        return len(samples)
    
    samples = np.asarray(samples)
    n = len(samples)
    
    # Calculate mean and variance
    mean_val = np.mean(samples)
    centered = samples - mean_val
    var = np.var(centered, ddof=1)
    
    if var == 0:
        return n
    
    # Calculate autocorrelation directly (more accurate than FFT for this purpose)
    # We only need autocorrelations up to a reasonable lag
    max_lag = min(n // 2, 1000)  # Limit to avoid excessive computation
    acf = np.zeros(max_lag + 1)
    acf[0] = 1.0  # Lag 0 autocorrelation is always 1
    
    # Calculate autocorrelations for lags 1 to max_lag
    for lag in range(1, max_lag + 1):
        if lag >= n:
            break
        # Calculate autocorrelation at this lag
        numerator = np.sum(centered[:n-lag] * centered[lag:])
        denominator = var * (n - lag)
        if denominator != 0:
            acf[lag] = numerator / denominator
        else:
            acf[lag] = 0.0
    
    # Implement Geyer's initial monotone sequence criterion
    # Sum pairs: P_k = rho_{2k} + rho_{2k+1}
    sum_pairs = 0.0
    prev_pair_sum = np.inf  # Start with infinity to ensure first pair is accepted
    
    k = 0
    while True:
        # Get indices for the current pair
        idx_2k = 2 * k
        idx_2k1 = 2 * k + 1
        
        # Check if we have enough lags
        if idx_2k >= len(acf):
            break
            
        # Get rho_{2k}
        rho_2k = acf[idx_2k]
        
        # Get rho_{2k+1} if available, otherwise 0
        rho_2k1 = acf[idx_2k1] if idx_2k1 < len(acf) else 0.0
        
        current_pair_sum = rho_2k + rho_2k1
        
        # Stop if the pair sum becomes negative
        if current_pair_sum < 0:
            break
            
        # Stop if the sequence is no longer monotone decreasing
        # (current pair sum should be <= previous pair sum)
        if current_pair_sum > prev_pair_sum + 1e-10:
            break
            
        sum_pairs += current_pair_sum
        prev_pair_sum = current_pair_sum
        k += 1
        
        # Safety check to avoid infinite loop
        if k > len(acf) // 2:
            break
    
    # tau_hat = -1 + 2 * sum_pairs
    tau_hat = -1.0 + 2.0 * sum_pairs
    
    # Avoid division by zero or negative values
    if tau_hat <= 0:
        return float(n)
    
    # ESS = n / tau_hat
    ess = n / tau_hat
    
    # Ensure ESS is within reasonable bounds
    return max(1.0, min(ess, float(n)))


def calculate_95HPD(samples):
    """
    Calculate 95% Highest Posterior Density (HPD) interval using numpy.
    
    Parameters:
    samples: array-like, MCMC samples
    
    Returns:
    tuple: (lower_bound, upper_bound) of 95% HPD interval
    """
    if len(samples) < 2:
        return (samples[0], samples[0]) if len(samples) == 1 else (0, 0)
    
    samples = np.asarray(samples)
    n = len(samples)
    
    # Sort the samples
    sorted_samples = np.sort(samples)
    
    # For 95% HPD, we want to exclude 5% of the samples
    # This means we need to find the shortest interval containing 95% of samples
    hpd_size = int(np.ceil(0.95 * n))
    
    if hpd_size >= n:
        return (sorted_samples[0], sorted_samples[-1])
    
    # Find the shortest interval containing hpd_size samples
    intervals = sorted_samples[hpd_size-1:] - sorted_samples[:n-hpd_size+1]
    min_idx = np.argmin(intervals)
    
    lower_bound = sorted_samples[min_idx]
    upper_bound = sorted_samples[min_idx + hpd_size - 1]
    
    return (lower_bound, upper_bound)


if __name__ == "__main__":
    """Simple test of the functions."""
    np.random.seed(42)
    
    # Test ESS
    independent_samples = np.random.normal(-825078, 100, 1000)
    ess = calculate_ESS(independent_samples)
    print(f"Independent samples ESS: {ess:.2f}")
    
    # Test 95% HPD
    hpd = calculate_95HPD(independent_samples)
    print(f"95% HPD: [{hpd[0]:.2f}, {hpd[1]:.2f}]")

def calculate_psrf(chain1, chain2):
    """Gelman--Rubin potential scale reduction factor (PSRF, R-hat).

    The two chains are truncated to the shorter post-burn-in length, which is
    the conservative choice when the runs were sampled at different
    frequencies.  A value near 1 indicates that the between-chain variance is
    small relative to the within-chain variance.  The function returns
    ``None`` when fewer than two paired samples are available.

    Parameters
    ----------
    chain1, chain2 : array-like
        Post-burn-in samples for one parameter from two independent chains.

    Returns
    -------
    float or None
        PSRF, ``inf`` for two constant but different chains, or ``None`` when
        the statistic is undefined.
    """
    first = np.asarray(chain1, dtype=float).ravel()
    second = np.asarray(chain2, dtype=float).ravel()
    first = first[np.isfinite(first)]
    second = second[np.isfinite(second)]

    length = min(len(first), len(second))
    if length < 2:
        return None

    first = first[:length]
    second = second[:length]

    means = np.array([first.mean(), second.mean()])
    variances = np.array([first.var(ddof=1), second.var(ddof=1)])
    within = float(variances.mean())

    if within <= 0.0:
        if np.allclose(first, second):
            return 1.0
        return float("inf")

    between = float(length * means.var(ddof=1))
    variance_hat = ((length - 1.0) / length) * within + between / length
    return float(np.sqrt(max(variance_hat / within, 0.0)))


def is_log_statistic(name):
    """Best-effort detection of log-scale trace columns (lnL, logPopSize...)."""
    if not name:
        return False
    normalized = "".join(ch for ch in name.lower() if ch.isalnum())
    if normalized.startswith(("ln", "log")):
        return True
    return any(marker in normalized
               for marker in ("lnl", "logl", "loglik", "lnlik", "loglike", "lnlike"))


def calculate_harmonic_mean(samples, effective_n=None):
    """Log-space harmonic mean and its delta-method standard error.

    Every likelihood trace is reported in log space, so the statistic is
    computed as log(HM(exp(x))): exponentiate, take the harmonic mean in
    linear space, then log the result. The standard error is transformed with
    the delta method and stays on the log scale.

    This is a descriptive statistic only; the harmonic-mean marginal-
    likelihood estimator is numerically unstable and not recommended for
    model comparison.
    """
    values = np.asarray(samples, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return None, None

    # Stable log-sum-exp form of HM(exp(x)).
    scaled = -values
    maximum = float(np.max(scaled))
    weights = np.exp(scaled - maximum)
    mean_weight = float(np.mean(weights))
    harmonic = -(maximum + float(np.log(mean_weight)))
    if len(values) < 2:
        return harmonic, 0.0
    count = len(values)
    if effective_n is not None:
        count = max(1.0, min(float(effective_n), float(count)))
    variance = float(np.var(weights, ddof=1))
    stderr = float(np.sqrt(variance / count) / mean_weight)
    return harmonic, stderr
