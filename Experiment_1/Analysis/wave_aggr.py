# Requirements:
#   pip install numpy pandas pywavelets matplotlib
# Notes:
# - Aggregates across all participants by averaging CWT power over participants.
# - Monte-Carlo significance (n_sim=200) is computed for the group-mean power
#   by simulating AR(1) noise per participant using their own AR(1) params.
# - Detrending per participant via polynomial (deg=2 by default).
# - Reconstruction from significant regions is done from mean CWT coefficients.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pywt
from glob import glob

# ----------------------------
# Data loading (same as before)
# ----------------------------
def load_data():
    data_names = glob('../Data/*.csv')
    data_files = pd.concat([pd.read_csv(f) for f in data_names])

    data_aggr = (data_files[data_files['Trial_type'] == 'experiment']
                 .query('Rt > 0.0')
                 [['PART_ID', 'CSI', 'Corr']]
                 .assign(t=lambda x: x['CSI'] * 16.6)
                 .groupby(['PART_ID', 't'])
                 .agg({'Corr': lambda x: 100.0 * x.mean()})
                 .rename(columns={'Corr': 'mean_corr'})
                 .reset_index())
    return data_aggr

# ----------------------------
# Utilities (as before + robust)
# ----------------------------
def detrend_poly(t, y, deg=2):
    p = np.polyfit(t, y, deg)
    return y - np.polyval(p, t)

def ar1_params(x):
    x = np.asarray(x) - np.mean(x)
    if len(x) < 3 or np.allclose(np.var(x), 0.0):
        return 0.0, np.std(x)  # fallback
    phi = np.corrcoef(x[:-1], x[1:])[0, 1]
    if not np.isfinite(phi):
        phi = 0.0
    phi = float(np.clip(phi, -0.99, 0.99))
    var = np.var(x, ddof=1)
    sigma_e = np.sqrt(max(1e-12, var * (1 - phi**2)))
    return phi, sigma_e

def simulate_ar1(phi, sigma_e, n, burnin=200, rng=None):
    # Accept either a numpy Generator/BitGenerator or a seed/None
    if rng is None or not hasattr(rng, "normal"):
        rng = np.random.default_rng(rng)
    e = rng.normal(0.0, sigma_e, size=n + burnin)
    y = np.zeros(n + burnin)
    for i in range(1, n + burnin):
        y[i] = phi * y[i - 1] + e[i]
    return y[burnin:]

def wavelet_scales_for_periods(period_min, period_max, dj, dt, wavelet_name='morl'):
    # For PyWavelets: period ~= (scale * dt) / fc
    fc = pywt.central_frequency(wavelet_name)
    s0 = (period_min * fc) / dt
    smax = (period_max * fc) / dt
    J = int(np.floor(np.log2(smax / s0) / dj))
    scales = s0 * 2.0 ** (np.arange(J + 1) * dj)
    periods = (scales * dt) / fc
    return scales, periods

def cwt_power(x, scales, wavelet_name='morl', dt=1.0):
    coeffs, _ = pywt.cwt(x, scales, wavelet_name, sampling_period=dt)  # (n_scales, n_time)
    power = np.abs(coeffs) ** 2
    return coeffs, power

def monte_carlo_sig_group(n_sim, phis, sigmas, n, scales, wavelet_name, dt, alpha=0.05, seed=0):
    # Simulate AR(1) per participant; average power across participants per sim.
    rng = np.random.default_rng(seed)
    n_part = len(phis)
    sim_q = []
    for k in range(n_sim):
        mean_pwr = None
        for i in range(n_part):
            y = simulate_ar1(phis[i], sigmas[i], n, rng=rng)
            _, pwr = cwt_power(y, scales, wavelet_name, dt)
            if mean_pwr is None:
                mean_pwr = pwr
            else:
                mean_pwr += pwr
        mean_pwr /= float(n_part)
        sim_q.append(mean_pwr)
    sim_q = np.stack(sim_q, axis=0)  # (n_sim, n_scales, n_time)
    thr = np.quantile(sim_q, 1 - alpha, axis=0)
    return thr

def inverse_cwt_from_masked(coeffs_masked, scales, wavelet_name='morl', dt=1.0, x_ref=None):
    if hasattr(pywt, 'icwt'):
        return pywt.icwt(coeffs_masked, scales, wavelet_name, sampling_period=dt)

    import warnings
    warnings.warn(
        "pywt.icwt not found. Using approximate inverse CWT; amplitudes may be scaled.",
        RuntimeWarning,
    )
    log2s = np.log2(scales)
    w = np.gradient(log2s)
    x_hat = np.sum(np.real(coeffs_masked) / np.sqrt(scales)[:, None] * w[:, None], axis=0)
    x_hat = x_hat - np.mean(x_hat)
    if x_ref is not None:
        s_hat = float(np.std(x_hat))
        s_ref = float(np.std(x_ref))
        if s_hat > 0:
            x_hat = x_hat * (s_ref / s_hat)
    return x_hat

# ----------------------------
# Build common time grid across participants
# ----------------------------
def build_common_panel(data_aggr, min_coverage=1.0):
    """
    Returns:
      t (n_time,), X (n_time, n_part), part_ids (n_part,)
    Keeps only time points (rows) where at least min_coverage fraction of participants have data.
    If min_coverage=1.0, this is strict intersection across all participants.
    """
    pivot = (data_aggr
             .pivot_table(index='t', columns='PART_ID', values='mean_corr', aggfunc='mean')
             .sort_index())
    # Keep only rows with enough coverage
    coverage = pivot.notna().mean(axis=1)
    pivot = pivot.loc[coverage >= float(min_coverage)]
    # After enforcing coverage, interpolate any remaining internal gaps and drop residual NaNs (should be none if min_coverage=1.0)
    pivot = pivot.interpolate(method='index', axis=0, limit_direction='both')
    pivot = pivot.dropna(axis=0, how='any')  # final safety
    t = pivot.index.to_numpy()
    X = pivot.to_numpy()  # (n_time, n_part)
    part_ids = pivot.columns.to_numpy()
    return t, X, part_ids

# ----------------------------
# Group analysis over all participants
# ----------------------------
def analyze_all_parts(
    data_aggr,
    deg=2,
    dt=1.0,
    dj=1/60,
    lowerPeriod=1.0,
    upperPeriod=16.0,
    n_sim=200,
    alpha=0.05,
    wavelet_name='morl',
    seed=0,
    min_coverage=1.0
):
    # Build common time grid and matrix
    t, X, part_ids = build_common_panel(data_aggr, min_coverage=min_coverage)  # X: (n_time, n_part)
    n_time, n_part = X.shape

    # Detrend each participant
    X_detr = np.empty_like(X, dtype=float)
    for j in range(n_part):
        X_detr[:, j] = detrend_poly(t, X[:, j], deg=deg)

    # Scales/periods
    scales, periods = wavelet_scales_for_periods(lowerPeriod, upperPeriod, dj, dt, wavelet_name)

    # CWT per participant
    coeffs_list = []
    power_list = []
    for j in range(n_part):
        c, p = cwt_power(X_detr[:, j], scales, wavelet_name, dt)
        coeffs_list.append(c)  # (n_scales, n_time)
        power_list.append(p)
    coeffs_stack = np.stack(coeffs_list, axis=0)     # (n_part, n_scales, n_time)
    power_stack = np.stack(power_list, axis=0)       # (n_part, n_scales, n_time)

    # Group means
    coeffs_mean = np.mean(coeffs_stack, axis=0)      # (n_scales, n_time)
    power_mean = np.mean(power_stack, axis=0)        # (n_scales, n_time)
    x_mean = np.mean(X_detr, axis=1)                 # (n_time,)

    # AR(1) params per participant (for MC)
    phis, sigmas = zip(*[ar1_params(X_detr[:, j]) for j in range(n_part)])

    # Monte Carlo significance for group-mean power
    thr = monte_carlo_sig_group(
        n_sim=n_sim,
        phis=phis,
        sigmas=sigmas,
        n=n_time,
        scales=scales,
        wavelet_name=wavelet_name,
        dt=dt,
        alpha=alpha,
        seed=seed
    )

    sig_mask = power_mean >= thr

    # Reconstruction from significant regions using mean coefficients
    coeffs_sig = np.where(sig_mask, coeffs_mean, 0.0)
    x_rec = inverse_cwt_from_masked(coeffs_sig, scales, wavelet_name, dt, x_ref=x_mean)

    out = {
        't': t,
        'periods': periods,
        'scales': scales,
        'power_mean': power_mean,
        'thr': thr,
        'sig_mask': sig_mask,
        'coeffs_mean': coeffs_mean,
        'x_detrended_mean': x_mean,
        'x_rec_sig_mean': x_rec,
        'part_ids': part_ids
    }
    return out

# ----------------------------
# Plotting: single plot for all participants (mean power)
# ----------------------------
def plot_wavelet_group(out,
                       timelab="PDI in [ms]",
                       periodlab="Period [data points per cycle]",
                       cmap='viridis'):
    t = out['t']
    P = out['power_mean']
    thr = out['thr']
    periods = out['periods']
    sig_mask = out['sig_mask']

    # Robust color scaling
    vmin = np.quantile(P, 0.01)
    vmax = np.quantile(P, 0.99)

    fig, ax = plt.subplots(figsize=(8, 4))
    im = ax.imshow(P,
                   extent=[t.min(), t.max(), periods.max(), periods.min()],  # small period at top
                   aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xlabel(timelab)
    ax.set_ylabel(periodlab)

    # Proper significance contour
    T, Per = np.meshgrid(t, periods)  # shapes match P (n_scales, n_time)
    if np.any(sig_mask):
        ax.contour(T, Per, sig_mask.astype(float),
                   levels=[0.5], colors='white', linewidths=0.8)
        # Legend proxy (avoid accessing cs.collections for compatibility)
        ax.plot([], [], color='white', lw=0.8, label='significant (group, p<0.05)')

    cbar = plt.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label('Mean wavelet power across participants')
    ax.set_title('Group wavelet power spectrum (Monte Carlo significance)')
    plt.tight_layout()
    plt.show()

def plot_reconstruction_group(out, timelab="PDI in [ms]"):
    t = out['t']
    x_mean = out['x_detrended_mean']
    x_rec = out['x_rec_sig_mean']

    fig, ax = plt.subplots(figsize=(8, 3))
    ax.plot(t, x_mean, label='Mean detrended x', lw=1.5, color='0.4')
    ax.plot(t, x_rec, label='Reconstruction (significant only)', lw=2.0, color='tab:blue')
    ax.set_xlabel(timelab)
    ax.legend(loc='lower left', frameon=False)
    ax.set_title('Group reconstruction from significant wavelet power')
    plt.tight_layout()
    plt.show()

# ----------------------------
# Example usage
# ----------------------------
if __name__ == "__main__":
    data_aggr = load_data()

    out_grp = analyze_all_parts(
        data_aggr=data_aggr,
        deg=2,
        dt=1.0,
        dj=1/60,
        lowerPeriod=1,
        upperPeriod=16,
        n_sim=200,         # reduce if slow
        alpha=0.05,
        wavelet_name='morl',
        seed=0,
        min_coverage=1.0   # if too strict, try 0.9 or 0.8
    )

    plot_wavelet_group(out_grp,
                       timelab="PDI in [ms]",
                       periodlab="Period [data points per cycle]")

    # Optional: show reconstruction of mean signal
    plot_reconstruction_group(out_grp, timelab="PDI in [ms]")