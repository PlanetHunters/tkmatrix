import os
from multiprocessing import Pool

from astropy import units as u
import numpy
from scipy.optimize import minimize_scalar

try:
    from tkmatrix.ravex.src.ravex.core import MultiPlanetSystem
except ModuleNotFoundError as error:
    if error.name is None or not error.name.startswith("tkmatrix.ravex"):
        raise
    raise ModuleNotFoundError("The ravex submodule is missing. Run `git submodule update --init --recursive`.",
                              name=error.name) from error

SEARCH_SAMPLES_PER_PEAK = 10
# The GLS search does not depend on the stellar mass, so a single system serves every recovery.
SEARCH_SYSTEM = MultiPlanetSystem(1 * u.M_sun, [])


class RecoverPeriodInput:
    """
    Used as input for the RV recovery stage
    """
    def __init__(self, time, rv_data, rv_err, period, msin=None, star_mass=None):
        """Initialise recovery period input.

        Parameters
        ----------
        time : numpy.ndarray
            Time series array.
        rv_data : numpy.ndarray
            Radial velocity data.
        rv_err : numpy.ndarray
            Radial velocity error values.
        period : float
            Orbital period to test, in days.
        msin : float or None, optional
            Minimum mass times sin(i). Defaults to None.
        star_mass : float or None, optional
            Star mass. Defaults to None.
        """
        self.time = time
        self.rv_data = rv_data
        self.rv_err = rv_err
        self.period = period
        self.msin = msin
        self.star_mass = star_mass


class RvFitter:
    """Main class used for RV injection and recovery"""
    def __init__(self) -> None:
        """Initialise the RV fitter."""
        super().__init__()

    @staticmethod
    def inject_rv(times, star_mass, rstar, planet_mass, period, t0):
        """
        Creates a synthetic signal of a planet with a given mass, period and t0

        :param times: the time series
        :param star_mass: the star mass
        :param rstar: the star radius. Unused, as the RV of a circular orbit does not depend on it
        :param planet_mass: the synthetic planet mass
        :param period: the synthetic planet period
        :param t0: the synthetic planet epoch
        :return: the injected RV data
        """
        # ravex places conjunction at a true anomaly of 90 degrees, which is the transit only for omega = 0.
        return MultiPlanetSystem(star_mass, []).inject_planet_in_series(
            numpy.asarray(times, dtype=float), planet_mass, period, time_conjunction=float(t0),
            argument_periapse=0 * u.deg)

    @staticmethod
    def sinfunc(t, k, omega, period):
        """
        Computes the RV measurement for the given times, semi-amplitude omega and period.

        :param t: the time series
        :param k: the semi-amplitude
        :param omega: the omega value
        :param period: the period of the planet
        :return: the RV expected measurements
        """
        return k * numpy.cos(omega + t * 2 * numpy.pi / period)

    @staticmethod
    def compute_semiamplitude_from_mmin(mmin, period, star_mass):
        """
        Calculates the semi-amplitude from the minimum mass, period and star mass

        :param mmin: the mass min value
        :param period: the period
        :param star_mass: the star mass
        :return: the semi-amplitude value
        """
        return mmin / ((star_mass ** (2. / 3.)) * (period ** (1. / 3.)) / 0.64)

    @staticmethod
    def compute_mmin_from_semiamplitude(period_grid, k_grid, star_mass):
        """
        Calculates the minimum mass from a grid of periods and semiamplitudes

        :param period_grid: the period grid
        :param k_grid: the semi-amplitude grid
        :param star_mass: the star mass
        :return: the minimum masses grid
        """
        m_min = []
        for i in numpy.arange(0, len(period_grid), 1):
            m_min_value = k_grid[i] * (star_mass ** (2./3.)) * (period_grid[i] ** (1./3.)) / 0.64
            m_min.append(m_min_value)
        return m_min

    @staticmethod
    def running_median(data, kernel):
        """Returns sliding median of width 'kernel' and same length as data

        :param data: the data for the running_median to be applied to
        :param kernel: the kernel size for the running_median
        :return: the data after the running median
        """
        idx = numpy.arange(kernel) + numpy.arange(len(data) - kernel + 1)[:, None]
        med = numpy.median(data[idx], axis=1)

        # Append the first/last value at the beginning/end to match the length of
        # data and returned median
        first_values = med[0]
        last_values = med[-1]
        missing_values = len(data) - len(med)
        values_front = int(missing_values * 0.5)
        values_end = missing_values - values_front
        med = numpy.append(numpy.full(values_front, first_values), med)
        med = numpy.append(med, numpy.full(values_end, last_values))
        return med

    @staticmethod
    def spectra(chi2, oversampling_factor=1, kernel_size=30):
        """
        Computes the SDE for a residuals series.

        :param chi2: the residuals list
        :param oversampling_factor: the oversampling factor
        :param kernel_size: the kernel size for the SDE computation
        :return: all the computed data
        """
        SR = numpy.min(chi2) / chi2
        SDE_raw = (1 - numpy.mean(SR)) / numpy.std(SR)

        # Scale SDE_power from 0 to SDE_raw
        power_raw = SR - numpy.mean(SR)  # shift down to the mean being zero
        scale = SDE_raw / numpy.max(power_raw)  # scale factor to touch max=SDE_raw
        power_raw = power_raw * scale

        # Detrended SDE, named "power"
        kernel = oversampling_factor * kernel_size
        if kernel % 2 == 0:
            kernel = kernel + 1
        if len(power_raw) > 2 * kernel:
            my_median = RvFitter.running_median(power_raw, kernel)
            power = power_raw - my_median
            # Re-normalize to range between median = 0 and peak = SDE
            # shift down to the mean being zero
            power = power - numpy.mean(power)
            SDE = numpy.max(power / numpy.std(power))
            # scale factor to touch max=SDE
            scale = SDE / numpy.max(power)
            power = power * scale
        else:
            power = power_raw
            SDE = SDE_raw

        return SR, power_raw, power, SDE_raw, SDE

    @staticmethod
    def recover_signal(rv_df, period, msin, omega, rv_masks, star_mass, period_grid_geom='lin', steps_period=None,
                       period_min=0.5, max_period=None, min_snr=5, min_sde=5, max_runs=5, cpus=os.cpu_count() - 1,
                       period_match_tolerance=0.01):
        """Iteratively recover a known injected RV signal.

        Parameters
        ----------
        rv_df : pandas.DataFrame
            DataFrame with columns 'bjd', 'rv', 'rv_err'.
        period : float
            Known injected period.
        msin : float
            Known minimum mass times sin(i).
        omega : float
            Known omega value.
        rv_masks : list
            List of RV masks to apply before recovery.
        star_mass : float
            Star mass.
        period_grid_geom : str, optional
            Period grid geometry, 'lin' or 'log'. Defaults to 'lin'.
        steps_period : int or None, optional
            Number of period steps. Defaults to None.
        period_min : float, optional
            Minimum period to search. Defaults to 0.5.
        max_period : float or None, optional
            Maximum period to search. Defaults to None.
        min_snr : float, optional
            Minimum SNR threshold to stop. Defaults to 5.
        min_sde : float, optional
            Minimum SDE threshold to stop. Defaults to 5.
        max_runs : int, optional
            Maximum number of recovery runs. Defaults to 5.
        cpus : int, optional
            Number of CPUs. Defaults to os.cpu_count() - 1.
        period_match_tolerance : float, optional
            Relative period difference below which the signal counts as found. Defaults to 0.01.

        Returns
        -------
        tuple
            (found, run, snr, SDE, signal_period, signal_omega, signal_msin)
        """
        df = rv_df.copy()
        sde = numpy.inf
        snr = numpy.inf
        found = False
        signal_period = 0
        signal_k = 0
        signal_omega = 0
        signal_msin = 0
        for run in numpy.arange(0, max_runs):
            rv_data, signal_period, signal_k, signal_omega, signal_msin, least_squares_grid, argmax, power, snr, SDE = \
                RvFitter.recover_strongest_signal(df, period_grid_geom, steps_period, period_min, max_period, rv_masks,
                                         star_mass, cpus)
            period_diff = numpy.abs(period - signal_period)
            found = period_diff / period < period_match_tolerance
            if found or snr < min_snr or SDE < min_sde:
                break
            else:
                rv_fit = RvFitter.sinfunc(df["bjd"], signal_k, signal_omega, signal_period)
                df["rv"] = df["rv"] - rv_fit
        return found, run, snr, SDE, signal_period, signal_omega, signal_msin

    @staticmethod
    def fit_phase_for_amplitude(angle, rv_data, rv_err, k):
        """Find the phase of a sinusoid of known amplitude, with the offset profiled out.

        Parameters
        ----------
        angle : numpy.ndarray
            Orbital angle of every measurement, in radians.
        rv_data : numpy.ndarray
            Radial velocity data.
        rv_err : numpy.ndarray
            Radial velocity error values.
        k : float
            The known semi-amplitude.

        Returns
        -------
        tuple
            (omega, omega_err)
        """
        weights = rv_err ** -2

        def centered(values):
            return values - numpy.average(values, axis=-1, weights=weights)[..., None]

        def chi2(omega):
            residuals = centered(rv_data - k * numpy.cos(angle + numpy.atleast_1d(omega)[:, None]))
            return numpy.sum(weights * residuals ** 2, axis=1)

        omega_grid = numpy.linspace(0, 2 * numpy.pi, 360, endpoint=False)
        best_omega = omega_grid[numpy.argmin(chi2(omega_grid))]
        step = omega_grid[1]
        omega = minimize_scalar(lambda value: chi2(value)[0], bounds=(best_omega - step, best_omega + step),
                                method='bounded').x % (2 * numpy.pi)
        curvature = numpy.sum(weights * centered(k * numpy.sin(angle + omega)) ** 2)
        return omega, 1 / numpy.sqrt(curvature) if curvature > 0 else numpy.inf

    @staticmethod
    def recover_period(input: RecoverPeriodInput):
        """Fit a sinusoidal RV model with a free offset for a single period.

        When the input gives a mass and a star mass, the amplitude is fixed from them and only the phase and the
        offset are fitted. The returned errors are formal: they trust `rv_err` and are not scaled by the fit quality.

        Parameters
        ----------
        input : RecoverPeriodInput
            Input data containing time, RV data, errors, and period.

        Returns
        -------
        tuple
            (period, k, omega, k_err, omega_err, least_squares)
        """
        period = input.period
        angle = 2 * numpy.pi * ((input.time / period) % 1)
        design = numpy.column_stack([numpy.ones_like(angle), numpy.cos(angle), numpy.sin(angle)])
        weighted_design = design / input.rv_err[:, None]
        try:
            covariance = numpy.linalg.inv(weighted_design.T @ weighted_design)
        except numpy.linalg.LinAlgError:
            return period, 0, 0, 0, 0, numpy.inf
        gamma, cos_amplitude, sin_amplitude = covariance @ (weighted_design.T @ (input.rv_data / input.rv_err))
        k = numpy.hypot(cos_amplitude, sin_amplitude)
        omega = numpy.arctan2(-sin_amplitude, cos_amplitude) % (2 * numpy.pi)
        if k == 0:
            k_err = numpy.sqrt(covariance[1, 1])
            omega_err = numpy.inf
        else:
            gradient_k = numpy.array([cos_amplitude, sin_amplitude]) / k
            gradient_omega = numpy.array([sin_amplitude, -cos_amplitude]) / k ** 2
            k_err = numpy.sqrt(gradient_k @ covariance[1:, 1:] @ gradient_k)
            omega_err = numpy.sqrt(gradient_omega @ covariance[1:, 1:] @ gradient_omega)
        if input.msin is not None and input.star_mass is not None:
            k = RvFitter.compute_semiamplitude_from_mmin(input.msin, period, input.star_mass)
            k_err = 0
            omega, omega_err = RvFitter.fit_phase_for_amplitude(angle, input.rv_data, input.rv_err, k)
            residuals = input.rv_data - RvFitter.sinfunc(input.time, k, omega, period)
            gamma = numpy.average(residuals, weights=input.rv_err ** -2)
        model = gamma + RvFitter.sinfunc(input.time, k, omega, period)
        least_squares = numpy.sum((input.rv_data - model) ** 2 / input.rv_err ** 2)
        return period, k, omega, k_err, omega_err, least_squares

    @staticmethod
    def recover_periods(rv_df, period_grid_geom='lin', steps_period=None, period_min=0.5, max_period=None,
                        rv_masks=None, star_mass=None, cpus=os.cpu_count() - 1):
        """Recover RV signals across a grid of periods.

        Parameters
        ----------
        rv_df : pandas.DataFrame
            DataFrame with columns 'bjd', 'rv', 'rv_err'.
        period_grid_geom : str, optional
            Period grid geometry, 'lin' or 'log'. Defaults to 'lin'.
        steps_period : int or None, optional
            Number of period steps. Defaults to None.
        period_min : float, optional
            Minimum period to search. Defaults to 0.5.
        max_period : float or None, optional
            Maximum period to search. Defaults to None.
        rv_masks : list or None, optional
            List of RV masks to apply. Defaults to None.
        star_mass : float or None, optional
            Star mass. Defaults to None.
        cpus : int, optional
            Number of CPUs. Defaults to os.cpu_count() - 1.

        Returns
        -------
        tuple
            (rv_data, period_grid, k_grid, omega_grid, msin_grid,
             least_squares_grid, argmax_sde, power, snr, SDE)
        """
        if rv_masks is None:
            rv_masks = []
        if steps_period is None:
            steps_period = int(rv_df["bjd"].max() - rv_df["bjd"].min())
        if max_period is None:
            max_period = (rv_df["bjd"].max() - rv_df["bjd"].min()) * 2
        time = rv_df["bjd"].to_numpy()
        rv_data = rv_df["rv"].to_numpy()
        rv_err = rv_df["rv_err"].to_numpy()
        rv_data = RvFitter.mask_signals(rv_df, rv_masks, star_mass)
        period_grid = numpy.linspace(period_min, max_period, steps_period) if period_grid_geom == "lin" \
            else numpy.logspace(numpy.log10(period_min), numpy.log10(max_period), steps_period)
        k_grid = []
        k_err_grid = []
        omega_grid = []
        omega_err_grid = []
        least_squares_grid = []
        inputs = []
        for period in period_grid:
            inputs.append(RecoverPeriodInput(time, rv_data, rv_err, period))
        with Pool(processes=cpus) as pool:
            recovered_outputs = pool.map(RvFitter.recover_period, inputs)
        period_grid = []
        for output in recovered_outputs:
            period_grid.append(output[0])
            k_grid.append(output[1])
            omega_grid.append(output[2])
            k_err_grid.append(output[3])
            omega_err_grid.append(output[4])
            least_squares_grid.append(output[5])
        least_squares_grid = numpy.array(least_squares_grid)
        least_squares_grid[least_squares_grid == 0] = numpy.max(least_squares_grid)
        msin_grid = RvFitter.compute_mmin_from_semiamplitude(period_grid, k_grid, star_mass)
        SR, power_raw, power, SDE_raw, SDE = RvFitter.spectra(numpy.array(least_squares_grid))
        argmax_sde = numpy.nanargmax(power)
        max_msin = msin_grid[argmax_sde]
        period_max_msin = period_grid[argmax_sde]
        snr_window_size = (max_period - period_min) / 20
        indexes_around = numpy.argwhere((period_grid > period_max_msin - snr_window_size / 2) &
                                        (period_grid < period_max_msin + snr_window_size / 2)).flatten()
        median_power_around = numpy.nanmedian(numpy.array(msin_grid)[indexes_around])
        snr = msin_grid[argmax_sde] / median_power_around
        return rv_data, period_grid, k_grid, omega_grid, msin_grid, least_squares_grid, argmax_sde, power, snr, SDE

    @staticmethod
    def recover_strongest_signal(rv_df, period_grid_geom='lin', steps_period=None, period_min=0.5, max_period=None,
                        rv_masks=None, star_mass=None, cpus=os.cpu_count() - 1):
        """Recover the strongest RV signal using a periodogram.

        Parameters
        ----------
        rv_df : pandas.DataFrame
            DataFrame with columns 'bjd', 'rv', 'rv_err'.
        period_grid_geom : str, optional
            Period grid geometry, 'lin' or 'log'. Defaults to 'lin'.
        steps_period : int or None, optional
            Number of period steps. Defaults to None.
        period_min : float, optional
            Minimum period to search. Defaults to 0.5.
        max_period : float or None, optional
            Maximum period to search. Defaults to None.
        rv_masks : list or None, optional
            List of RV masks to apply. Defaults to None.
        star_mass : float or None, optional
            Star mass. Defaults to None.
        cpus : int, optional
            Number of CPUs. Defaults to os.cpu_count() - 1.

        Returns
        -------
        tuple
            (rv_data, strongest_period, k, omega, msin, least_squares,
             argmax_sde, power, snr, sde)
        """
        if rv_masks is None:
            rv_masks = []
        if steps_period is None:
            steps_period = int(rv_df["bjd"].max() - rv_df["bjd"].min())
        if max_period is None:
            max_period = (rv_df["bjd"].max() - rv_df["bjd"].min()) * 2
        time = rv_df["bjd"].to_numpy()
        rv_err = rv_df["rv_err"].to_numpy()
        rv_data = RvFitter.mask_signals(rv_df, rv_masks, star_mass)
        recovered = SEARCH_SYSTEM.recover_periodic_signal(time, rv_data, rv_err, min_period=period_min,
                                                          max_period=max_period,
                                                          samples_per_peak=SEARCH_SAMPLES_PER_PEAK)
        period_grid = recovered['period_grid_days']
        # GLS power is quadratic in the amplitude. Its root keeps SNR and SDE on the scale of their thresholds.
        power = numpy.sqrt(numpy.clip(recovered['power_grid'], 0, None))
        sde_power = power / numpy.std(power)
        argmax_sde = numpy.nanargmax(sde_power)
        sde = sde_power[argmax_sde]
        strongest_period, k, omega, k_err, omega_err, least_squares = \
            RvFitter.recover_period(RecoverPeriodInput(time, rv_data, rv_err, period_grid[argmax_sde]))
        msin = RvFitter.compute_mmin_from_semiamplitude([strongest_period], [k], star_mass)[0]
        snr_window_size = (max_period - period_min) / 20
        indexes_around = numpy.argwhere((period_grid > strongest_period - snr_window_size / 2) &
                                        (period_grid < strongest_period + snr_window_size / 2)).flatten()
        snr = power[argmax_sde] / numpy.nanmedian(power[indexes_around])
        return rv_data, strongest_period, k, omega, msin, least_squares, argmax_sde, power, snr, sde

    @staticmethod
    def mask_signals(rv_df, rv_masks, star_mass):
        """Subtract RV signals specified by masks from the data.

        Parameters
        ----------
        rv_df : pandas.DataFrame
            DataFrame with columns 'bjd', 'rv', 'rv_err'.
        rv_masks : list
            List of dicts with keys 'P' (period) and optionally 'M' (mass).
        star_mass : float
            Star mass.

        Returns
        -------
        numpy.ndarray
            RV data with masked signals subtracted.
        """
        time = rv_df["bjd"].to_numpy()
        rv_data = rv_df["rv"].to_numpy()
        rv_err = rv_df["rv_err"].to_numpy()
        for rv_mask in rv_masks:
            period, k, omega, k_err, omega_err, least_squares = \
                RvFitter.recover_period(RecoverPeriodInput(time, rv_data, rv_err, rv_mask['P'],
                                                           rv_mask['M'] if 'M' in rv_mask else None, star_mass))
            rv_subtract = RvFitter.sinfunc(time, k, omega, period)
            rv_data = rv_data - rv_subtract
        return rv_data
