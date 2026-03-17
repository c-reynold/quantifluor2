from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def weighted_curve_fit(
    model_func,
    x,
    y,
    p0=None,
    weighting: str | None = None,
    bounds=(-np.inf, np.inf),
    absolute_sigma: bool = False,
    make_plot: bool = True,
    n_curve_points: int = 500,
    title: str = "Curve Fit",
):
    """
    Fit a custom model to x/y data with optional weighting and residual plots.

    Parameters
    ----------
    model_func : callable
        Function of the form f(x, *params).
    x : array-like
        X values for fit.
    y : array-like
        Y values for fit.
    p0 : list/tuple/array, optional
        Initial parameter guesses. Should leave empty unill you have 
    weighting : str or None
        Weighting scheme. Options:
            None       -> unweighted
            "1/y^2"    -> weights proportional to 1 / y^2
            "1/y"      -> weights proportional to 1 / y
    bounds : 2-tuple
        Lower and upper parameter bounds for curve_fit.
    absolute_sigma : bool
        Passed to scipy.optimize.curve_fit.
    make_plot : bool
        Whether to generate fit and residual plots.
    n_curve_points : int
        Number of points used to draw the smooth fitted curve.
    title : str
        Plot title.

    Returns
    -------
    results : dict
        Dictionary containing:
            - popt: fitted parameter values
            - pcov: covariance matrix
            - perr: parameter standard errors
            - y_fit: fitted y values at original x points
            - residuals: y - y_fit
            - sigma: sigma array used by curve_fit
    """

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if x.shape != y.shape:
        raise ValueError("Number of x and y fit data coordinates don't match.")

    # Remove NaN/inf rows
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]

    if len(x) < 3:
        raise ValueError("Need at least 3 non-empty and valid data point pairs to fit.")

    sigma = _build_sigma(y=y, weighting=weighting)

    popt, pcov = curve_fit(
        f=model_func,
        xdata=x,
        ydata=y,
        p0=p0,
        sigma=sigma,
        absolute_sigma=absolute_sigma,
        bounds=bounds,
        maxfev=10000,
    )

    y_fit = model_func(x, *popt)
    residuals = y - y_fit
    perr = np.sqrt(np.diag(pcov))

    if make_plot:
        _plot_fit_and_residuals(
            model_func=model_func,
            x=x,
            y=y,
            popt=popt,
            residuals=residuals,
            sigma=sigma,
            n_curve_points=n_curve_points,
            title=title,
        )

    return {
        "popt": popt,
        "pcov": pcov,
        "perr": perr,
        "y_fit": y_fit,
        "residuals": residuals,
        "sigma": sigma,
    }


def _build_sigma(y: np.ndarray, weighting: str | None):
    """
    Build sigma for scipy curve_fit.

    Important:
    curve_fit uses sigma as standard deviations, and minimizes:
        sum(((y - f(x)) / sigma) ** 2)

    So:
    - weighting = 1 / y^2  means sigma should be proportional to y
    - weighting = 1 / y    means sigma should be proportional to sqrt(y), but
      for simplicity many people loosely use sigma ~ sqrt(|y|) when variance
      scales with signal. Here, for "1/y", we use
      sigma = sqrt(|y|).

    For 1 / y^2 weighting, we use:
        sigma = abs(y)

    To avoid division-by-zero or tiny unstable sigma values, a floor is applied.
    """
    if weighting is None:
        return None

    y_abs = np.abs(y)

    # Small floor to avoid zero sigma values
    floor = max(np.nanmax(y_abs) * 1e-12, 1e-12)
    y_safe = np.maximum(y_abs, floor)

    if weighting == "1/y^2":
        return y_safe

    if weighting == "1/y":
        return np.sqrt(y_safe)

    raise ValueError(
        "Unsupported weighting. Use None, '1/y^2', or '1/y'."
    )


def _plot_fit_and_residuals(
    model_func,
    x,
    y,
    popt,
    residuals,
    sigma,
    n_curve_points,
    title,
):
    """
    Plot data + fit and residuals.
    """
    x_curve = np.linspace(np.min(x), np.max(x), n_curve_points)
    y_curve = model_func(x_curve, *popt)

    fig, axes = plt.subplots(
        2, 1, figsize=(8, 8), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]}
    )

    ax_fit = axes[0]
    ax_res = axes[1]

    if sigma is not None:
        ax_fit.errorbar(
            x, y, yerr=sigma, fmt="o", capsize=3, label="Data"
        )
    else:
        ax_fit.plot(x, y, "o", label="Data")

    ax_fit.plot(x_curve, y_curve, "-", label="Fit")
    ax_fit.set_ylabel("y")
    ax_fit.set_title(title)
    ax_fit.legend()
    ax_fit.grid(True, alpha=0.3)

    ax_res.axhline(0, linestyle="--")
    ax_res.plot(x, residuals, "o")
    ax_res.set_xlabel("x")
    ax_res.set_ylabel("Residual")
    ax_res.grid(True, alpha=0.3)

    
    ax_fit.set_xscale("log")
    ax_fit.set_yscale("log")

    ax_res.set_xscale("log")

    #right now, the plot is just shown. It can be saved as a PNG, but not sure exactly how to do that in the moment with the way that this function is wrapped in another function in the fitting.py script.
    plt.tight_layout()
    plt.show()



def four_pl_interpolate_concentration(y,bottom,top,ec50,hill_slope):
    """
    Given a flouresence value, and fit parameters, calculate the interpolated concentration for a 4PL fit
    """

    conc=ec50*((y-bottom)/(top-y)) ** (1/hill_slope)

    return conc
    