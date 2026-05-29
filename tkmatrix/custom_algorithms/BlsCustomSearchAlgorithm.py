from abc import ABC, abstractmethod

import lightkurve
import wotan
import numpy as np
from foldedleastsquares import transit_mask, cleaned_array
from lcbuilder.HarmonicSelector import HarmonicSelector
from lcbuilder.helper import LcbuilderHelper
from lcbuilder.star.starinfo import StarInfo

from tkmatrix.custom_algorithms.custom_search import CustomSearchAlgorithm


class BlsCustomSearchAlgorithm(CustomSearchAlgorithm):
    """
    Custom search algorithm to be implemented if injecting a custom search class in the properties file
    """
    def __init__(self):
        """Initialise the BLS custom search algorithm."""
        super().__init__()

    def search(self, time, flux, rstar, rstar_min, rstar_max, mass, mstar_min, mstar_max,
               ab, epoch, period, min_period, max_period, min_snr, cores, transit_template, detrend_method, ws,
               transits_min_count, signal_selection_mode, run_limit, oversampling, period_match_tolerance,
               epoch_match_tolerance):
        """Search for a transit signal using the BLS algorithm.

        Parameters
        ----------
        time : numpy.ndarray
            Time series array.
        flux : numpy.ndarray
            Flux values corresponding to the time series.
        rstar : float
            Star radius.
        rstar_min : float
            Minimum star radius value.
        rstar_max : float
            Maximum star radius value.
        mass : float
            Star mass.
        mstar_min : float
            Minimum star mass value.
        mstar_max : float
            Maximum star mass value.
        ab : tuple
            Quadratic limb darkening coefficients.
        epoch : float
            Epoch (t0) of the signal to be spotted.
        period : float
            Period of the signal to be spotted.
        min_period : float
            Minimum period for the period grid.
        max_period : float
            Maximum period for the period grid.
        min_snr : float
            SNR threshold to stop searching.
        cores : int
            Number of processes for parallel computation.
        transit_template : str
            Transit template to use for the search.
        detrend_method : str
            Detrending method.
        ws : float
            Window size for detrending.
        transits_min_count : int
            Minimum number of transits for a valid signal.
        signal_selection_mode : str
            Signal selection mode (e.g. 'period-epoch').
        run_limit : int
            Maximum number of runs.
        oversampling : float
            Period grid oversampling factor.
        period_match_tolerance : float
            Tolerance for period matching.
        epoch_match_tolerance : float
            Tolerance for epoch matching.

        Returns
        -------
        tuple
            (found_signals, snrs, sdes, runs, durations, periods, t0s)
        """

        snr = 1e12
        found_signal = False
        time, flux = cleaned_array(time, flux)
        run = 0
        if ws > 0:
            flux = wotan.flatten(time, flux, window_length=ws, return_trend=False, method=detrend_method,
                                 break_tolerance=0.5)
        found_signals = []
        snrs = []
        sdes = []
        runs = []
        durations = []
        periods = []
        t0s = []
        while snr >= min_snr and not found_signal and (run_limit > 0 and run < run_limit):
            star_info = StarInfo(mass=mass, radius=rstar)
            tls_period_grid, oversampling = LcbuilderHelper.calculate_period_grid(time, min_period, max_period,
                                                                                  oversampling, star_info,
                                                                                  transits_min_count)
            lc = lightkurve.LightCurve(time=time, flux=flux)
            results = lc.to_periodogram(method='bls', period=tls_period_grid, frequency_factor=oversampling * 100)
            max_power_index = np.argmax(results.power)
            sde = results.power[max_power_index].value / np.nanmedian(results.power).value
            t0 = results.transit_time_at_max_power.value
            duration = results.duration_at_max_power.value
            found_period = results.period_at_max_power.value
            intransit_result = transit_mask(time, found_period, 2 * duration, t0)
            real_intransit_result = transit_mask(time, found_period, duration, t0)
            oot_flux = flux[~real_intransit_result]
            it_flux = flux[real_intransit_result]
            snr = np.abs(np.nanmean(1 - it_flux)) / np.nanstd(oot_flux) * (len(it_flux) ** 0.5)
            if snr >= min_snr:
                time = time[~intransit_result]
                flux = flux[~intransit_result]
                time, flux = cleaned_array(time, flux)
                if results.transit_time is not None and len(results.transit_time) > 0:
                    print(f"Selecting signal with mode {signal_selection_mode}")
                    if signal_selection_mode == 'period-epoch':
                        found_signal = HarmonicSelector.is_harmonic(t0, epoch, found_period, period,
                                                                    epoch_match_tolerance, period_match_tolerance)
                    else:
                        found_signal = HarmonicSelector.multiple_of(found_period, period, period_match_tolerance) != 0
                    if found_signal:
                        found_signals = found_signals + [found_signal]
                        snrs = snrs + [snr]
                        sdes = sdes + [sde]
                        runs = runs + [run]
                        durations = durations + [duration]
                        periods = periods + [found_period]
                        t0s = t0s + [t0]
                        break
            found_signals = found_signals + [found_signal]
            snrs = snrs + [snr]
            sdes = sdes + [sde]
            runs = runs + [run]
            durations = durations + [duration]
            periods = periods + [found_period]
            t0s = t0s + [t0]
            run = run + 1
        return found_signals, snrs, sdes, runs, durations, periods, t0s
