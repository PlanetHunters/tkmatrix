"""Unit tests for the RV engine (issue #34). Run offline, using independent oracles."""
import numpy as np
import pandas as pd
import pytest
import astropy.units as u
import lightkurve

from tkmatrix.rv import RecoverPeriodInput, RvFitter

MSTAR = 0.5
SEARCH_MAX_PERIOD = 25.0


def published_semiamplitude(mass_earth, period_days, star_mass_sun):
    """K for a circular edge-on orbit, from the 28.4329 m/s relation of Lovis & Fischer (2010)."""
    mass_jupiter = mass_earth / 317.8284
    period_years = period_days / 365.25
    return 28.4329 * mass_jupiter * star_mass_sun ** (-2 / 3) * period_years ** (-1 / 3)


def series(seed=1, points=80, span=300.0, noise=1.0):
    rng = np.random.default_rng(seed)
    time = np.sort(rng.uniform(2458000.0, 2458000.0 + span, points))
    return time, rng.normal(0, noise, points), np.full(points, noise)


def frame(time, rv, rv_err):
    return pd.DataFrame({"bjd": time, "rv": rv, "rv_err": rv_err})


def test_ravex_is_available_from_the_submodule():
    from tkmatrix.ravex.src.ravex.core import MultiPlanetSystem
    assert callable(MultiPlanetSystem.inject_planet_in_series)
    assert callable(MultiPlanetSystem.recover_periodic_signal)


@pytest.mark.parametrize("period,mass", [(3.5, 2.0), (11.3, 8.0), (40.0, 50.0)])
def test_injection_matches_published_circular_orbit(period, mass):
    time = 2458000.0 + np.linspace(0, 60, 400)
    t0 = 2458001.3
    injected = RvFitter.inject_rv(time, MSTAR * u.M_sun, 0.5 * u.R_sun, mass * u.M_earth, period, t0)
    k = published_semiamplitude(mass, period, MSTAR)
    expected = -k * np.sin(2 * np.pi * (time - t0) / period)
    assert np.sqrt(np.mean((injected - expected) ** 2)) / k < 1e-3


def test_injection_crosses_zero_descending_at_the_transit_epoch():
    t0, period = 2458001.3, 7.0
    around = np.array([t0 - 0.01, t0, t0 + 0.01])
    before, at, after = RvFitter.inject_rv(around, MSTAR * u.M_sun, 0.5 * u.R_sun, 10 * u.M_earth, period, t0)
    assert before > 0 > after
    assert abs(at) < 1e-3 * before / 0.01


def test_injection_accepts_reduced_time_stamps():
    time = np.linspace(58400.0, 58460.0, 50)
    injected = RvFitter.inject_rv(time, MSTAR * u.M_sun, 0.5 * u.R_sun, 5 * u.M_earth, 5.0, 58401.0)
    assert np.all(np.isfinite(injected))
    assert np.abs(injected).max() > 0


@pytest.mark.parametrize("offset", [0.0, 100.0, -5678.0])
def test_recover_period_handles_rv_offsets(offset):
    time, _, rv_err = series()
    k, omega, period = 3.0, 1.2, 6.3
    rv = offset + RvFitter.sinfunc(time, k, omega, period)
    found_period, found_k, found_omega, k_err, omega_err, chi2 = RvFitter.recover_period(
        RecoverPeriodInput(time, rv, rv_err, period))
    assert found_period == pytest.approx(period)
    assert found_k == pytest.approx(k, rel=1e-6)
    assert found_omega == pytest.approx(omega, rel=1e-6)
    assert chi2 == pytest.approx(0, abs=1e-12)


def test_recover_period_keeps_the_sinfunc_convention():
    time, noise, rv_err = series()
    _, k, omega, k_err, omega_err, _ = RvFitter.recover_period(
        RecoverPeriodInput(time, RvFitter.sinfunc(time, 4.0, 5.5, 9.0) + noise, rv_err, 9.0))
    assert k > 0
    assert 0 <= omega < 2 * np.pi
    assert 0 < k_err < 1
    assert 0 < omega_err < 1


def test_mask_with_mass_fixes_the_amplitude_from_the_published_relation():
    time, noise, rv_err = series()
    period, mass, omega = 6.3, 10.0, 2.0
    k_mask = published_semiamplitude(mass, period, MSTAR)
    k_true = 2.5 * k_mask
    rv = RvFitter.sinfunc(time, k_true, omega, period) + noise
    _, k, found_omega, k_err, omega_err, _ = RvFitter.recover_period(
        RecoverPeriodInput(time, rv, rv_err, period, mass, MSTAR))
    assert k == pytest.approx(k_mask, rel=2e-3)
    assert k_err == 0
    assert found_omega == pytest.approx(omega, abs=0.1)
    masked = RvFitter.mask_signals(frame(time, rv, rv_err), [{"P": period, "M": mass}], MSTAR)
    remaining = RvFitter.recover_period(RecoverPeriodInput(time, masked, rv_err, period))[1]
    assert remaining == pytest.approx(k_true - k_mask, abs=0.5)


def test_fixed_amplitude_phase_is_the_constrained_optimum_near_a_window_alias():
    rng = np.random.default_rng(9)
    time = 2458000.0 + np.arange(40) + rng.normal(0, 0.02, 40)
    rv_err = rng.uniform(0.5, 3.0, 40)
    period, mass = 1.002, 30.0
    k = RvFitter.compute_semiamplitude_from_mmin(mass, period, MSTAR)
    rv = RvFitter.sinfunc(time, k, 1.0, period) + rng.normal(0, 1, 40) * rv_err
    chi2 = RvFitter.recover_period(RecoverPeriodInput(time, rv, rv_err, period, mass, MSTAR))[5]
    scan = []
    for omega in np.linspace(0, 2 * np.pi, 7200, endpoint=False):
        residuals = rv - RvFitter.sinfunc(time, k, omega, period)
        residuals = residuals - np.average(residuals, weights=rv_err ** -2)
        scan.append(np.sum((residuals / rv_err) ** 2))
    assert chi2 <= np.min(scan) + 1e-6


def test_fit_errors_are_the_formal_ones():
    time, _, rv_err = series(points=100)
    rng = np.random.default_rng(12)
    fits = [RvFitter.recover_period(RecoverPeriodInput(time, RvFitter.sinfunc(time, 4.0, 2.0, 7.7) +
                                                       rng.normal(0, 1.0, 100), rv_err, 7.7)) for _ in range(400)]
    amplitudes, phases = np.array([fit[1] for fit in fits]), np.array([fit[2] for fit in fits])
    assert np.mean([fit[3] for fit in fits]) == pytest.approx(np.std(amplitudes), rel=0.15)
    assert np.mean([fit[4] for fit in fits]) == pytest.approx(np.std(phases), rel=0.15)


def test_fit_errors_are_not_scaled_by_the_fit_quality():
    time, _, rv_err = series(points=100)
    rng = np.random.default_rng(13)
    fits = [RvFitter.recover_period(RecoverPeriodInput(time, RvFitter.sinfunc(time, 12.0, 2.0, 7.7) +
                                                       rng.normal(0, 3.0, 100), rv_err, 7.7)) for _ in range(400)]
    assert np.mean([fit[3] for fit in fits]) == pytest.approx(np.std([fit[1] for fit in fits]) / 3, rel=0.15)


def test_fixed_amplitude_phase_error_matches_the_scatter():
    time, _, rv_err = series(points=100)
    rng = np.random.default_rng(14)
    period, mass, omega = 7.7, 30.0, 2.0
    k = RvFitter.compute_semiamplitude_from_mmin(mass, period, MSTAR)
    fits = [RvFitter.recover_period(RecoverPeriodInput(time, RvFitter.sinfunc(time, k, omega, period) +
                                                       rng.normal(0, 1.0, 100), rv_err, period, mass, MSTAR))
            for _ in range(200)]
    assert np.mean([fit[4] for fit in fits]) == pytest.approx(np.std([fit[2] for fit in fits]), rel=0.2)


def test_mask_without_mass_removes_the_signal():
    time, noise, rv_err = series()
    rv = RvFitter.sinfunc(time, 6.0, 2.0, 6.3) + noise
    masked = RvFitter.mask_signals(frame(time, rv, rv_err), [{"P": 6.3}], MSTAR)
    assert np.std(masked) == pytest.approx(np.std(noise), rel=0.2)


def test_strongest_signal_keeps_tuple_shape_and_finds_the_period():
    time, noise, rv_err = series()
    k, period = 4.0, 6.3
    result = RvFitter.recover_strongest_signal(frame(time, RvFitter.sinfunc(time, k, 1.0, period) + noise, rv_err),
                                               max_period=SEARCH_MAX_PERIOD, star_mass=MSTAR)
    assert len(result) == 10
    rv_data, found_period, found_k, omega, msin, chi2, argmax, power, snr, sde = result
    assert found_period == pytest.approx(period, rel=0.01)
    assert found_k == pytest.approx(k, rel=0.2)
    assert msin == pytest.approx(RvFitter.compute_mmin_from_semiamplitude([found_period], [found_k], MSTAR)[0])
    assert power[argmax] == np.max(power)


@pytest.mark.parametrize("seed", [2, 3, 4])
def test_snr_and_sde_keep_the_amplitude_periodogram_scale(seed):
    """With equal errors the weighted search must match lightkurve's amplitude statistics."""
    time, noise, rv_err = series(seed=seed)
    rv = noise + RvFitter.sinfunc(time, 1.0, 1.0, 6.3)
    *_, snr, sde = RvFitter.recover_strongest_signal(frame(time, rv, rv_err), max_period=SEARCH_MAX_PERIOD,
                                                     star_mass=MSTAR)
    periodogram = lightkurve.LightCurve(time=time, flux=rv - rv.mean()).to_periodogram(
        oversample_factor=10, minimum_period=0.4, maximum_period=SEARCH_MAX_PERIOD)
    periods, power = periodogram.period.value, periodogram.power.value
    best = periods[np.argmax(power)]
    window = (SEARCH_MAX_PERIOD - 0.5) / 20
    around = (periods > best - window / 2) & (periods < best + window / 2)
    assert sde == pytest.approx(np.max(power / np.std(power)), rel=0.05)
    assert snr == pytest.approx(np.max(power) / np.median(power[around]), rel=0.1)


def test_recover_signal_prewhitens_a_stronger_unrelated_signal():
    time, noise, rv_err = series(points=120)
    rv = RvFitter.sinfunc(time, 12.0, 0.4, 2.77) + RvFitter.sinfunc(time, 3.0, 1.0, 9.1) + noise
    found, run, snr, sde, period, omega, msin = RvFitter.recover_signal(
        frame(time, rv, rv_err), 9.1, 3, 0.5, None, MSTAR, 'log', 100, 0.5, SEARCH_MAX_PERIOD, 5, 5, 3)
    assert found
    assert run == 1
    assert period == pytest.approx(9.1, rel=0.01)


def test_recover_signal_applies_masks_before_searching():
    time, noise, rv_err = series(points=120)
    rv = RvFitter.sinfunc(time, 12.0, 0.4, 2.77) + RvFitter.sinfunc(time, 3.0, 1.0, 9.1) + noise
    found, run, *_ = RvFitter.recover_signal(frame(time, rv, rv_err), 9.1, 3, 0.5, [{"P": 2.77}], MSTAR, 'log', 100,
                                             0.5, SEARCH_MAX_PERIOD, 5, 5, 3)
    assert found
    assert run == 0


def test_found_does_not_depend_on_the_thresholds():
    time, noise, rv_err = series()
    rv = RvFitter.sinfunc(time, 6.0, 1.0, 6.3) + noise
    found, run, *_ = RvFitter.recover_signal(frame(time, rv, rv_err), 6.3, 3, 0.5, None, MSTAR, 'log', 100, 0.5,
                                             SEARCH_MAX_PERIOD, 1000, 1000, 3)
    assert found
    assert run == 0


@pytest.mark.parametrize("period_min, expected", [(0.4, 0.45), (0.5, None)])
def test_search_starts_at_the_given_minimum_period(period_min, expected):
    time, noise, rv_err = series(points=200, noise=0.2)
    rv = noise + 5.0 * np.sin(2 * np.pi * time / 0.45)
    strongest = RvFitter.recover_strongest_signal(frame(time, rv, rv_err), 'log', 1000, period_min, SEARCH_MAX_PERIOD,
                                                  [], MSTAR)[1]
    assert strongest >= period_min
    if expected is not None:
        assert strongest == pytest.approx(expected, rel=0.001)


@pytest.mark.parametrize("tolerance, found", [(0.01, False), (0.05, True)])
def test_found_uses_the_period_match_tolerance(tolerance, found):
    time, noise, rv_err = series(noise=0.2)
    rv = noise + 5.0 * np.sin(2 * np.pi * time / 5.15)
    result = RvFitter.recover_signal(frame(time, rv, rv_err), 5.0, 3, 0.5, [], MSTAR, 'log', 1000, 0.5,
                                     SEARCH_MAX_PERIOD, 5, 5, 1, period_match_tolerance=tolerance)
    assert result[0] == found
    assert result[4] == pytest.approx(5.15, rel=0.005)
