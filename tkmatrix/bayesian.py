import os.path
import sys
import numpy as np
import matplotlib.colors as mcolors
np_mod = sys.modules.get("numpy")
import scipy.integrate as si
# --- Patch numpy.int (for old code that imports it) ---
if not hasattr(np, "int"):
    np.int = int
# --- Patch scipy.integrate.trapz (for old code that imports it) ---
if not hasattr(si, "trapz"):
    # just delegate to numpy.trapz – signature-compatible for typical use
    def _trapz(y, x=None, dx=1.0, axis=-1):
        return np.trapz(y, x=x, dx=dx, axis=axis)
    si.trapz = _trapz

from spright import RMRelation
import pandas as pd
import astropy.units as u

import matplotlib.pyplot as plt
from lcbuilder.helper import LcbuilderHelper
from lcbuilder.star.HabitabilityCalculator import HabitabilityCalculator
from lcbuilder.star.starinfo import StarInfo

# Earth radius in units of solar radius
R_EARTH_OVER_R_SUN = 0.0091577


def make_completeness_function(radius_grid_Re, completeness_grid):
    """
    Build a simple 1D linear completeness interpolator C(R_p) using numpy.

    Parameters
    ----------
    radius_grid_Re : array-like
        Radii in Earth radii where MATRIX has been evaluated.
    completeness_grid : array-like
        Corresponding detection probabilities in [0, 1].

    Returns
    -------
    callable
        C(R_p) that returns completeness for any R_p (Earth radii),
        extrapolating flat outside the provided range.
    """
    radius_grid_Re = np.asarray(radius_grid_Re)
    completeness_grid = np.asarray(completeness_grid)

    def C(Rp_Re):
        """1D linear interpolator for detection completeness at a given planet radius.

        Parameters
        ----------
        Rp_Re : float or array-like
            Planet radius in Earth radii.

        Returns
        -------
        float or ndarray
            Completeness probability in [0, 1].
        """
        Rp_Re = np.asarray(Rp_Re)
        # np.interp extrapolates flat at the extremes
        return np.interp(Rp_Re, radius_grid_Re, completeness_grid)

    return C


def sample_radius_prior_parviainen(m_sin_i, m_sin_i_std, size=1, rng=None):
    """
    EXAMPLE radius prior sampler.

    Replace this with your actual implementation that uses the
    Parviainen et al. (2023) mass–radius relation.

    For now, we just assume a log-normal in R (Earth radii) centered on some
    R0(m) with sigma in ln R. You must calibrate R0 and sigma to your case.
    """
    rmr = RMRelation()
    mds = rmr.predict_radius(mass=(m_sin_i, m_sin_i / 10), nsamples=size)
    return mds.samples


def sample_posterior_no_transit(
        N_samples,
        m_sin_i,
        a_over_Rstar,
        Rstar_Rsun,
        completeness_func,
        radius_sampler,
        batch_size=100000,
        rng=None,
):
    """
    Draw samples from the posterior p(i, R_p, M | 'no detected transit'),
    combining:
      - isotropic inclination prior p(i) ∝ sin i,
      - radius prior given by `radius_sampler`,
      - MATRIX completeness C(R_p) at the known RV period,
      - condition 'no transit detected'.

    Parameters
    ----------
    N_samples : int
        Number of posterior samples to return.
    m_sin_i : float
        Minimum mass from RV in arbitrary units (e.g. Earth masses).
        All masses will be returned in the same units.
    a_over_Rstar : float
        Semi-major axis in units of stellar radius (a / R_star).
    Rstar_Rsun : float
        Stellar radius in solar radii.
    completeness_func : callable
        Function C(R_p) that takes planet radius in Earth radii and returns
        detection probability in [0, 1] at that period for your TESS+MATRIX setup.
    radius_sampler : callable
        Function R_p = radius_sampler(size, rng) that returns radii (Earth radii)
        drawn from your radius prior (e.g. Parviainen+2023 posterior).
    batch_size : int, optional
        How many prior samples to draw per loop when thinning by non-detection.
    rng : np.random.Generator, optional
        Random number generator. If None, a new default_rng() is created.

    Returns
    -------
    dict
        {
          "inclinations_rad": array,
          "radii_Re": array,
          "masses": array  # same mass units as m_sin_i
        }
    """
    if rng is None:
        rng = np.random.default_rng()

    accepted_i = []
    accepted_Rp = []

    while len(accepted_i) < N_samples:
        # 1) Sample inclinations from isotropic prior: cos(i) ~ U(0, 1)
        u = rng.uniform(0.0, 1.0, size=batch_size)
        i = np.arccos(u)  # radians, in [0, pi/2]

        # 2) Sample radii from prior (Earth radii)
        Rp_Re = radius_sampler(size=batch_size, rng=rng)

        # 3) Compute transit threshold i_tr (radians)
        # R_p / R_star:
        Rp_over_Rstar = (Rp_Re * R_EARTH_OVER_R_SUN) / Rstar_Rsun

        # Argument of arccos can be slightly >1 due to numerical issues, clip it:
        arg = (1.0 + Rp_over_Rstar) / a_over_Rstar
        arg = np.clip(arg, 0.0, 1.0)
        i_tr = np.arccos(arg)  # radians

        # 4) Determine whether each configuration is geometrically transiting
        transiting = i >= i_tr

        # 5) Compute probability of 'no detection' for each configuration
        C = completeness_func(Rp_Re)  # detection prob if transiting
        # P(no detection | i, R) = 1 if non-transiting; = 1 - C(R) if transiting
        p_no_det = np.where(transiting, 1.0 - C, 1.0)

        # 6) Rejection sampling: keep configurations consistent with "no detection"
        # Draw uniforms and accept if u < p_no_det
        u2 = rng.uniform(0.0, 1.0, size=batch_size)
        accept_mask = (u2 < p_no_det)

        accepted_i.extend(i[accept_mask])
        accepted_Rp.extend(Rp_Re[accept_mask])

    # Trim to requested N_samples
    accepted_i = np.array(accepted_i[:N_samples])
    accepted_Rp = np.array(accepted_Rp[:N_samples])

    # 7) Compute true masses: M = (m sin i) / sin(i)
    # (Mass units are whatever units m_sin_i is in)
    sin_i = np.sin(accepted_i)
    # Guard against sin(i)~0; geometrically that almost never happens for exoplanets
    sin_i = np.clip(sin_i, 1e-6, None)
    M_true = m_sin_i / sin_i

    return {
        "inclinations_rad": accepted_i,
        "radii_Re": accepted_Rp,
        "masses": M_true,
    }

def geometric_mass_lower_limit(m_sin_i, a_over_Rstar, Rp_Re, Rstar_Rsun):
    """
    Closed-form lower limit on true mass assuming:
      - completeness ~1 for the relevant R_p,
      - any transiting configuration would be detected and is excluded.

    Parameters
    ----------
    m_sin_i : float
        RV minimum mass (same units as output).
    a_over_Rstar : float
        a / R_star.
    Rp_Re : float
        Representative planet radius (Earth radii), e.g. median of radius prior.
    Rstar_Rsun : float
        Stellar radius (solar radii).

    Returns
    -------
    M_min : float
        Lower limit on true mass, same units as m_sin_i.
    i_tr_rad : float
        Transit threshold inclination in radians.
    """
    Rp_over_Rstar = (Rp_Re * R_EARTH_OVER_R_SUN) / Rstar_Rsun
    arg = (1.0 + Rp_over_Rstar) / a_over_Rstar
    arg = np.clip(arg, 0.0, 1.0)
    i_tr_rad = np.arccos(arg)
    M_min = m_sin_i / np.sin(i_tr_rad)
    return M_min, i_tr_rad

def prior_geometric_transit_probability(
    N_samples, a_over_Rstar, Rstar_Rsun, radius_sampler, rng=None
):
    """
    Estima la probabilidad geométrica de tránsito P_geom
    muestreando del prior p(i,R) = sin(i) * p(R).
    """
    if rng is None:
        rng = np.random.default_rng()

    # prior en i: cos(i) ~ U(0,1)
    u = rng.uniform(0.0, 1.0, size=N_samples)
    i = np.arccos(u)

    # prior en radio
    Rp_Re = radius_sampler(size=N_samples, rng=rng)

    # condición de tránsito
    Rp_over_Rstar = (Rp_Re * R_EARTH_OVER_R_SUN) / Rstar_Rsun
    arg = (1.0 + Rp_over_Rstar) / a_over_Rstar
    arg = np.clip(arg, 0.0, 1.0)
    i_tr = np.arccos(arg)

    transiting = i >= i_tr
    return np.mean(transiting.astype(float))

def posterior_transit_probability(i_post, Rp_post, a_over_Rstar, Rstar_Rsun):
    """
    Compute P(transit | no detection) from posterior samples.

    Parameters
    ----------
    i_post : array
        Posterior inclinations in radians (samples from p(i,R | O)).
    Rp_post : array
        Posterior planet radii in Earth radii (same length as i_post).
    a_over_Rstar : float
        Semi-major axis in units of stellar radius (a / R_star).
    Rstar_Rsun : float
        Stellar radius in solar radii.

    Returns
    -------
    float
        Posterior probability that the planet actually transits,
        given the non-detection with TESS+MATRIX.
    """
    Rp_over_Rstar = (Rp_post * R_EARTH_OVER_R_SUN) / Rstar_Rsun
    arg = (1.0 + Rp_over_Rstar) / a_over_Rstar
    arg = np.clip(arg, 0.0, 1.0)
    i_tr = np.arccos(arg)  # radians

    transiting_mask = i_post >= i_tr
    P_transit_given_no_det = np.mean(transiting_mask.astype(float))

    return P_transit_given_no_det

def plot_radius_completeness(df_dirs, star_mass, star_radius, planet_masses, planet_masses_err, sigma=95, colors=['blue']):
    """Plot recovery rate vs planet radius for one or more planet mass scenarios.

    For each planet mass, samples the radius prior (Parviainen et al. 2023),
    computes median and confidence intervals, and plots the completeness curve
    from the MATRIX results. Saves the figure as ``completeness.png``.

    Parameters
    ----------
    df_dirs : list of str
        Paths to the MATRIX result CSV files (``a_tls_report.csv``).
    star_mass : float
        Stellar mass in solar masses.
    star_radius : float
        Stellar radius in solar radii.
    planet_masses : list of float
        Planet minimum masses (m sin i) in Earth masses.
    planet_masses_err : list of float
        Uncertainties on the planet masses.
    sigma : float, optional
        Confidence level for the radius range (e.g. 95 for 95%). Default 95.
    colors : list of str, optional
        Matplotlib color codes for each planet. Default ``['blue']``.
    """
    mp_errs = planet_masses_err
    plt.figure(figsize=(8,5))
    for index, mp in enumerate(mps):
        df = pd.read_csv(df_dirs[index])
        rp_samples = sample_radius_prior_parviainen(mp, mp_errs[index], size=100000)
        rp = np.median(rp_samples)
        rp_low_err = rp - np.percentile(rp_samples, (100 - sigma) / 2)
        rp_up_err = np.percentile(rp_samples, 100 - ((100 - sigma) / 2)) - rp
        vmin = rp - rp_low_err
        vmax = rp + rp_up_err
        vmax = np.min([df.loc[:, 'radius'].max(), vmax])
        df['found'] = df['found'].astype(bool)
        grouped = df.groupby('radius')['found'].mean().reset_index()
        #plt.fill_between([vmin, vmax], [1], color=radius_region_color[index], alpha=0.3, label="Radius range")
        if index == 0:
            plt.axvline(vmin, color='black', linestyle="--", linewidth=1.5)
            plt.axvline(vmax, color='black', linestyle="--", linewidth=1.5)
        plt.plot(grouped['radius'], grouped['found'], marker='o', color=colors[index])
    plt.xlabel("Radius ($R_\oplus$)", fontsize=18)
    plt.ylabel("Recovery rate", fontsize=18)
    plt.tick_params(axis="both", labelsize=16)
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(os.path.dirname(df_dirs[0]) + "/" + "completeness.png")
    plt.close()

def compute_bayesian_transit_probabilities(df_dir, star_mass, star_radius, planet_mass, planet_mass_err, sigma=95):
    """Compute Bayesian posterior transit probability and corrected minimum mass.

    Combines an RV mass constraint with MATRIX completeness and a mass-radius
    relation prior to derive posterior distributions for inclination, true mass,
    and radius, plus transit probabilities (prior and posterior).

    Parameters
    ----------
    df_dir : str
        Path to the MATRIX result CSV (``a_tls_report.csv``).
    star_mass : float
        Stellar mass in solar masses.
    star_radius : float
        Stellar radius in solar radii.
    planet_mass : float
        Planet minimum mass (m sin i) in Earth masses.
    planet_mass_err : float
        Uncertainty on the planet mass.
    sigma : float, optional
        Confidence level used in summary stats. Default 95.

    Returns
    -------
    P_geom : float
        Prior geometric transit probability.
    P_tr_posterior : float
        Posterior transit probability given non-detection.
    M_post : ndarray
        Posterior true mass samples in Earth masses.
    i_deg : ndarray
        Posterior inclination samples in degrees.
    i_tr_limit : float
        Transit threshold inclination in degrees.
    corrected_min_mass : float
        Corrected minimum mass (2.5th percentile of posterior).
    C_samples : ndarray
        Completeness samples at the radius posterior.
    """
    period = df['period'].unique().mean()
    a, _, _ = HabitabilityCalculator().calculate_semi_major_axis(period, 0.01, 0.01, star_mass,
                                                           star_mass - star_mass / 10,
                                                           star_mass / 10)
    rs = star_radius
    mp = planet_mass
    mp_err = planet_mass_err
    df['found'] = df['found'].astype(bool)
    grouped = df.groupby('radius')['found'].mean().reset_index()
    C = make_completeness_function(grouped['radius'], grouped['found'])
    def my_radius_sampler(size, rng=None):
        """Wrap :func:`sample_radius_prior_parviainen` with fixed mass and uncertainty.

        Parameters
        ----------
        size : int
            Number of samples to draw.
        rng : numpy.random.Generator, optional
            Random number generator.

        Returns
        -------
        ndarray
            Sampled planet radii in Earth radii.
        """
        return sample_radius_prior_parviainen(m_sin_i, m_sin_i_std, size=size, rng=rng)
    m_sin_i = mp
    m_sin_i_std = mp_err
    a_Rs = a / LcbuilderHelper.convert_from_to(rs, u.R_sun, u.AU)
    post = sample_posterior_no_transit(
        N_samples=50000,
        m_sin_i=m_sin_i,
        a_over_Rstar=a_Rs,
        Rstar_Rsun=rs,
        completeness_func=C,
        radius_sampler=my_radius_sampler,
    )
    i_post = post["inclinations_rad"]
    M_post = post["masses"]
    Rp_post = post["radii_Re"]
    def summary_stats(x):
        """Compute 0.15th, 50th, and 99.85th percentiles.

        Parameters
        ----------
        x : array-like
            Input data array.

        Returns
        -------
        ndarray
            Array of [0.15th, 50th, 99.85th] percentiles.
        """
        return np.percentile(x, [0.15, 50, 99.85])
    i_deg = np.degrees(i_post)
    print("i [deg] 0.15/50/99.85:", summary_stats(i_deg))
    print("M_post 0.15/50/99.85:", summary_stats(M_post))
    P_geom = prior_geometric_transit_probability(
        N_samples=200000,
        a_over_Rstar=a_Rs,
        Rstar_Rsun=rs,
        radius_sampler=my_radius_sampler,
    )
    P_tr_posterior = posterior_transit_probability(
        i_post=i_post,
        Rp_post=Rp_post,
        a_over_Rstar=a_Rs,
        Rstar_Rsun=rs,
    )
    rp_samples = sample_radius_prior_parviainen(planet_mass, planet_mass_err, size=100000)
    rp = np.median(rp_samples)
    i_tr_limit = np.rad2deg(np.arccos((rs + LcbuilderHelper.convert_from_to(rp, u.R_earth, u.R_sun)) / LcbuilderHelper.convert_from_to(a, u.AU, u.R_sun)))
    print("Max i [deg]", i_tr_limit)
    print("Prior transit probability = ", P_geom)
    print("Posterior transit probability P(transit | no detection) = ", P_tr_posterior)
    corrected_min_mass = np.percentile(M_post, 2.5)
    print("m*sin(i) = ", m_sin_i)
    print("Corrected minimum mass = ", np.percentile(M_post, 2.5))
    C_samples = C(rp_samples)
    print("Completeness stats 0.15/50/99.85:", summary_stats(C_samples))
    return P_geom, P_tr_posterior, M_post, i_deg, i_tr_limit, corrected_min_mass, C_samples
