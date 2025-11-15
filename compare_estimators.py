import numpy as np
import matplotlib.pyplot as plt

def compare_estimators_over_theta(scales, xi_ls, err_ls, xi_dp, err_dp, 
                       xlabel=r'$\theta$ [deg]', title="LS vs DP"):
    """
    Compare Landy–Szalay and Davis–Peebles correlation results with error propagation.

    Parameters
    ----------
    scales : array
        Separation bins (e.g. r or theta).
    xi_ls : array
        Correlation function from Landy–Szalay estimator.
    err_ls : array
        1σ errors on xi_ls.
    xi_dp : array
        Correlation function from Davis–Peebles estimator.
    err_dp : array
        1σ errors on xi_dp.
    xlabel : str
        Label for the x-axis (default: r or θ).
    title : str
        Figure title.

    Returns
    -------
    ratio, ratio_err, residual, residual_err
    """
    xi_ls = np.asarray(xi_ls)
    xi_dp = np.asarray(xi_dp)
    err_ls = np.asarray(err_ls)
    err_dp = np.asarray(err_dp)

    # Residual and its error
    residual = xi_dp - xi_ls
    residual_err = np.sqrt(err_dp**2 + err_ls**2)

    # Ratio and its error (avoid division by zero)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = xi_dp / xi_ls
        ratio_err = np.abs(ratio) * np.sqrt(
            (err_dp/xi_dp)**2 + (err_ls/xi_ls)**2
        )

        # Mask invalid entries (NaN/inf when xi=0)
        ratio[~np.isfinite(ratio)] = np.nan
        ratio_err[~np.isfinite(ratio_err)] = np.nan

    # --- Plotting ---
    fig, axes = plt.subplots(3, 1, figsize=(6, 11), sharex=True,
                             gridspec_kw={'height_ratios': [2, 1, 1]})

    # Correlation functions
    axes[0].errorbar(scales, xi_ls, yerr=err_ls, fmt='o-', label='Landy–Szalay')
    axes[0].errorbar(scales, xi_dp, yerr=err_dp, fmt='s--', label='Davis–Peebles')
    axes[0].set_ylabel(r'$\omega_{\rm fg-bg}(\theta)$')
    axes[0].legend()
    #axes[0].set_xscale('log')
    axes[0].set_title(title)

    # Ratio
    axes[1].errorbar(scales, ratio, yerr=ratio_err, fmt='o-', color='purple')
    axes[1].axhline(1.0, color='k', ls='--')
    axes[1].set_ylabel(r"$\omega_{DP}/\omega_{LS}$")
    #axes[1].set_ylim(-2, 2)
    #axes[1].set_xscale('log')
    # Residual
    axes[2].errorbar(scales, residual, yerr=residual_err, fmt='o-', color='darkgreen')
    axes[2].axhline(0.0, color='k', ls='--')
    axes[2].set_xlabel(xlabel)
    axes[2].set_ylabel(r"$\omega_{DP}-\omega_{LS}$")
    #axes[2].set_xscale('log')
    plt.xscale('log')
    plt.tight_layout()
    plt.show()

    return ratio, ratio_err, residual, residual_err

def compare_estimators_over_redshift(z_bins, dp_vals, dp_err, ls_vals, ls_err, title="DP vs LS Comparison"):
    """
    Compare Davis-Peebles (DP) and Landy-Szalay (LS) cross-correlation results,
    and plot the ratios and residuals as a function of redshift.

    Parameters
    ----------
    z_bins : array-like
        Central redshift values (or bin centers).
    dp_vals : array-like
        Cross-correlation values from Davis-Peebles estimator per redshift bin.
    dp_err : array-like
        Errors on DP values.
    ls_vals : array-like
        Cross-correlation values from Landy-Szalay estimator per redshift bin.
    ls_err : array-like
        Errors on LS values.
    title : str
        Title for the plot.

    Returns
    -------
    ratio, ratio_err, residual, residual_err : np.ndarray
        Arrays of ratios and residuals with propagated errors.
    """

    dp_vals = np.array(dp_vals, dtype=float)
    dp_err = np.array(dp_err, dtype=float)
    ls_vals = np.array(ls_vals, dtype=float)
    ls_err = np.array(ls_err, dtype=float)
    z_bins = np.array(z_bins, dtype=float)

    # --- Residuals ---
    residual = dp_vals - ls_vals
    residual_err = np.sqrt(dp_err**2 + ls_err**2)

    # --- Ratios ---
    ratio = np.full_like(dp_vals, np.nan)
    ratio_err = np.full_like(dp_vals, np.nan)
    valid = ls_vals != 0
    ratio[valid] = dp_vals[valid] / ls_vals[valid]
    ratio_err[valid] = np.sqrt(
        (dp_err[valid] / ls_vals[valid])**2 +
        (dp_vals[valid] * ls_err[valid] / ls_vals[valid]**2)**2
    )

    # --- Plot ---
    fig, ax = plt.subplots(1, 2, figsize=(10, 5), sharex=True)
    plt.suptitle(title)
    # Residuals
    ax[0].errorbar(z_bins, residual, yerr=residual_err, fmt='o', capsize=4, label="DP - LS")
    ax[0].axhline(0, color='k', ls='--', lw=1)
    ax[0].set_ylabel("Residual (DP - LS)")
    ax[0].set_xlabel("Redshift bin [z]")
    #ax[0].set_title(title)
    ax[0].legend()

    # Ratios
    ax[1].errorbar(z_bins, ratio, yerr=ratio_err, fmt='o', capsize=4, label="DP / LS")
    ax[1].axhline(1, color='k', ls='--', lw=1)
    ax[1].set_ylabel("Ratio (DP/LS)")
    ax[1].set_xlabel("Redshift bin [z]")
    ax[1].legend()

    plt.tight_layout()
    plt.show()

    return ratio, ratio_err, residual, residual_err