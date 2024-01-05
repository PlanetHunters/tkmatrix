from abc import ABC, abstractmethod


class CustomSearchAlgorithm(ABC):
    """
    Custom search algorithm to be implemented if injecting a custom search class in the properties file
    """
    def __init__(self):
        pass

    @abstractmethod
    def search(self, time, flux, rstar, rstar_min, rstar_max, mass, mstar_min, mstar_max,
               ab, epoch, period, min_period, max_period, min_snr, cores,
               transit_template, detrend_method, ws, transits_min_count,
               signal_selection_mode, run_limit, oversampling):

        """
        Searches for a signal with a given epoch and period in the signal given by time and flux.

        :param time: the time series
        :param flux: the flux associated to the time series
        :param rstar: the star radius
        :param rstar_min: minimum value of the star radius
        :param rstar_max: maximum value of the star radius
        :param mass: the star mass
        :param mstar_min: the minimum value of the star mass
        :param mstar_max: the maximum value of the star mass
        :param ab: quadratic limb darkening parameters of the target
        :param epoch: t0 of the signal to be spotted
        :param period: the period of the signal to be spotted
        :param min_period: the minimum period for the period grid
        :param max_period: the maximum period for the period grid
        :param min_snr: the SNR threshold to stop searching
        :param cores: the number of processes to split the computation
        :param transit_template: the transit template to use for the search
        :param detrend_method: the strategy for detrending
        :param ws: the window size for the detrend to be applied
        :param transits_min_count: the minimum number of transits for a signal to be valid
        :param signal_selection_mode: the way for retrieving the signal
        :param run_limit: the number of runs to limit the search
        :param oversampling: the density of the period grid
        """
        pass
