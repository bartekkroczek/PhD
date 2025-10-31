# Requirements:
#   pip install numpy pandas pywavelets matplotlib
# Notes:
# - Uses PyWavelets for CWT and Monte‑Carlo significance (n_sim=200) against AR(1) noise.
# - Detrending by subtracting a polynomial of chosen degree (deg=2 by default).
# - Reconstructs signal from significant regions only via inverse CWT.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pywt

# ----------------------------
# Inputs (replace with your data/vars)
# ----------------------------
from glob import glob


def load_data():
    # Read all CSV files from ../Data directory
    data_names = glob('../Data/*.csv')
    data_files = pd.concat([pd.read_csv(f) for f in data_names])

    # Process data similar to R code
    data_aggr = (data_files[data_files['Trial_type'] == 'experiment']
                 .query('Rt > 0.0')
                 [['PART_ID', 'CSI', 'Corr']]
                 .assign(t=lambda x: x['CSI'] * 16.6)
                 .groupby(['PART_ID', 't'])
                 .agg({'Corr': lambda x: 100.0 * x.mean()})
                 .rename(columns={'Corr': 'mean_corr'})
                 .reset_index())
    return data_aggr


# Example usage:
data_aggr = load_data()
random_part_id = data_aggr['PART_ID'].iloc[0]  # Get first participant ID
random_part_id = "145K21"

def detrend_poly(t, y, deg=2):
    p = np.polyfit(t, y, deg)
    return y - np.polyval(p, t)

def ar1_params(x):
    x = np.asarray(x) - np.mean(x)
    # lag-1 autocorrelation
    phi = np.corrcoef(x[:-1], x[1:])[0, 1]
    phi = np.clip(phi, -0.99, 0.99)
    var = np.var(x, ddof=1)
    sigma_e = np.sqrt(max(1e-12, var * (1 - phi**2)))
    return phi, sigma_e

def simulate_ar1(phi, sigma_e, n, burnin=200):
    e = np.random.normal(0.0, sigma_e, size=n + burnin)
    y = np.zeros(n + burnin)
    for i in range(1, n + burnin):
        y[i] = phi * y[i - 1] + e[i]
    return y[burnin:]

def wavelet_scales_for_periods(period_min, period_max, dj, dt, wavelet_name='morl'):
    # For PyWavelets: period ~= (scale * dt) / fc, where fc = central_frequency(wavelet)
    fc = pywt.central_frequency(wavelet_name)
    s0 = (period_min * fc) / dt
    smax = (period_max * fc) / dt
    J = int(np.floor(np.log2(smax / s0) / dj))
    scales = s0 * 2.0 ** (np.arange(J + 1) * dj)
    periods = (scales * dt) / fc
    return scales, periods

def cwt_power(x, scales, wavelet_name='morl', dt=1.0):
    coeffs, _ = pywt.cwt(x, scales, wavelet_name, sampling_period=dt)  # shape: (n_scales, n_time)
    power = np.abs(coeffs) ** 2
    return coeffs, power

def monte_carlo_sig(power_fun, n_sim, phi, sigma_e, n, scales, wavelet_name, dt, alpha=0.05, rng=None):
    rng = np.random.default_rng(rng)
    # Estimate threshold at each (scale, time) as (1-alpha) quantile over simulations
    # To avoid huge memory, accumulate online via order statistics (simple approach: stack).
    sim_q = []
    for k in range(n_sim):
        y = simulate_ar1(phi, sigma_e, n)
        _, pwr = power_fun(y, scales, wavelet_name, dt)
        sim_q.append(pwr)
    sim_q = np.stack(sim_q, axis=0)  # (n_sim, n_scales, n_time)
    thr = np.quantile(sim_q, 1 - alpha, axis=0)  # (n_scales, n_time)
    return thr

def inverse_cwt_from_masked(coeffs_masked, scales, wavelet_name='morl', dt=1.0, x_ref=None):
    """Inverse CWT with compatibility fallback.

    If pywt.icwt is available, use it. Otherwise, approximate the inverse
    using a log-scale weighted sum of real(CWT)/sqrt(scale) and rescale to
    match the reference series variance (if provided).
    """
    # Native inverse if available (PyWavelets >= 1.1/1.2 depending on distro)
    if hasattr(pywt, 'icwt'):
        return pywt.icwt(coeffs_masked, scales, wavelet_name, sampling_period=dt)

    # Fallback: approximate inverse following Torrence & Compo-like weighting.
    # Note: amplitude may differ by a constant factor; we variance-match to x_ref.
    try:
        import warnings
        warnings.warn(
            "pywt.icwt not found. Using approximate inverse CWT; amplitudes may be scaled.",
            RuntimeWarning,
        )
    except Exception:
        pass

    # We integrate over log2(scales): weights approximate dj
    log2s = np.log2(scales)
    w = np.gradient(log2s)  # shape (n_scales,)

    # Sum over scales: real part and 1/sqrt(scale) normalization
    x_hat = np.sum(np.real(coeffs_masked) / np.sqrt(scales)[:, None] * w[:, None], axis=0)

    # Mean-center
    x_hat = x_hat - np.mean(x_hat)

    # Match variance to reference detrended series if provided
    if x_ref is not None:
        s_hat = float(np.std(x_hat))
        s_ref = float(np.std(x_ref))
        if s_hat > 0:
            x_hat = x_hat * (s_ref / s_hat)
    return x_hat

# ----------------------------
# Pipeline
# ----------------------------
def analyze_one_part(data_aggr, random_part_id,
                     deg=2,
                     dt=1.0,
                     dj=1/60,
                     lowerPeriod=1.0,
                     upperPeriod=16.0,
                     n_sim=200,
                     alpha=0.05,
                     wavelet_name='morl',
                     seed=0):
    df = (data_aggr.query("PART_ID == @random_part_id")
                    .sort_values('t')
                    .reset_index(drop=True))
    t = df['t'].to_numpy()
    x_raw = df['mean_corr'].to_numpy()

    # Detrend (poly degree = deg)
    x = detrend_poly(t, x_raw, deg=deg)

    # Scales/periods
    scales, periods = wavelet_scales_for_periods(lowerPeriod, upperPeriod, dj, dt, wavelet_name)

    # CWT
    coeffs, power = cwt_power(x, scales, wavelet_name, dt)

    # Monte Carlo significance (AR1 null)
    phi, sigma_e = ar1_params(x)
    thr = monte_carlo_sig(cwt_power, n_sim, phi, sigma_e, len(x), scales, wavelet_name, dt, alpha, rng=seed)

    sig_mask = power >= thr
    coeffs_sig = np.where(sig_mask, coeffs, 0.0)
    x_rec = inverse_cwt_from_masked(coeffs_sig, scales, wavelet_name, dt, x_ref=x)

    out = {
        't': t,
        'x_detrended': x,
        'coeffs': coeffs,
        'power': power,
        'periods': periods,
        'scales': scales,
        'thr': thr,
        'sig_mask': sig_mask,
        'x_rec_sig': x_rec
    }
    return out

# ----------------------------
# Plotting similar to wt.image and reconstruct
# ----------------------------
def plot_wavelet(out,
                 timelab="PDI in [ms]",
                 periodlab="Period [data points per cycle]",
                 time_tick_pos=(0, 12, 24, 36, 48),
                 time_tick_labels=(200, 400, 600, 800, 1000),
                 cmap='viridis'):
    t = out['t']
    power = out['power']
    thr = out['thr']
    periods = out['periods']
    sig_mask = out['sig_mask']

    # Robust color scaling similar to "quantile" with many levels
    vmin = np.quantile(power, 0.01)
    vmax = np.quantile(power, 0.99)

    fig, ax = plt.subplots(figsize=(8, 4))
    im = ax.imshow(power,
                   extent=[0, len(t)-1, periods.max(), periods.min()],
                   aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xlabel(timelab)
    ax.set_ylabel(periodlab)

    # Custom time ticks (optional; adapt as needed)
    valid_ticks = [p for p in time_tick_pos if 0 <= p < len(t)]
    ax.set_xticks(valid_ticks)
    if len(valid_ticks) == len(time_tick_labels):
        ax.set_xticklabels(time_tick_labels)

    # Significance contour at alpha (power >= thr)
    # Draw a thin contour where mask transitions true/false
    try:
        import matplotlib
        cs = ax.contour(sig_mask.astype(float),
                        levels=[0.5],
                        colors='white',
                        linewidths=0.8,
                        extent=[0, len(t)-1, periods.max(), periods.min()])
        cs.collections[0].set_label(f'significant (p<{0.05})')
    except Exception:
        pass

    cbar = plt.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label('Wavelet power')
    ax.set_title('Wavelet power spectrum (Monte Carlo significance)')
    plt.tight_layout()
    plt.show()

def plot_reconstruction(out, timelab="PDI in [ms]"):
    t = out['t']
    x = out['x_detrended']
    x_rec = out['x_rec_sig']

    fig, ax = plt.subplots(figsize=(8, 3))
    ax.plot(t, x, label='Detrended x', lw=1.5, color='0.4')
    ax.plot(t, x_rec, label='Reconstruction (significant only)', lw=2.0, color='tab:blue')
    ax.set_xlabel(timelab)
    ax.legend(loc='lower left', frameon=False)
    ax.set_title('Reconstruction from significant wavelet power')
    plt.tight_layout()
    plt.show()

# ----------------------------
# Example usage
# ----------------------------
out = analyze_one_part(
    data_aggr=data_aggr,
    random_part_id=random_part_id,
    deg=2,              # polynomial degree for detrending
    dt=1.0,
    dj=1/60,
    lowerPeriod=1,
    upperPeriod=16,
    n_sim=200,
    alpha=0.05,         # matches siglvl=0.05
    wavelet_name='morl',
    seed=0
)
plot_wavelet(out,
             timelab="PDI in [ms]",
             periodlab="Period [data points per cycle]")
plot_reconstruction(out, timelab="PDI in [ms]")