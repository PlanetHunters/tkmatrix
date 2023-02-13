import os
from multiprocessing import Pool

import lightkurve
from astropy import units as u
import astropy.constants as ac
import matplotlib.pyplot as plt
import numpy
import ellc
import pandas
import scipy
from scipy.optimize import differential_evolution


class RecoverPeriodInput:
    """
    Used as input for the RV recovery stage
    """
    def __init__(self, time, rv_data, rv_err, period, msin=None, star_mass=None):
        self.time = time
        self.rv_data = rv_data
        self.rv_err = rv_err
        self.period = period
        self.msin = msin
        self.star_mass = star_mass


class RvFitter:
    """Main class used for RV injection and recovery"""
    def __init__(self) -> None:
        super().__init__()

    @staticmethod
    def inject_rv(times, star_mass, rstar, planet_mass, period, t0):
        """
        Creates a synthetic signal of a planet with a given mass, period and t0

        :param times: the time series
        :param star_mass: the star mass
        :param rstar: the star radius
        :param planet_mass: the synthetic planet mass
        :param period: the synthetic planet period
        :param t0: the synthetic planet epoch
        :return: the injected RV data
        """
        rplanet = rstar / 50
        planet_mass = planet_mass.to(u.M_sun)
        period_days = period * u.day
        a = numpy.cbrt((ac.G * star_mass * period_days ** 2) / (4 * numpy.pi ** 2)).to(u.au)
        model_kms = ellc.rv(
            t_obs=times,
            radius_1=rstar.to(u.au) / a,  # star radius convert from AU to in units of a
            radius_2=rplanet / a,  # convert from Rearth (equatorial) into AU and then into units of a
            sbratio=0,
            incl=90,
            t_zero=t0,
            period=period_days.value,
            a=a.to(u.R_sun).value,
            q=planet_mass / star_mass,
            f_c=None, f_s=None,
            gdc_1=None, gdc_2=None,
            didt=None,
            domdt=None,
            rotfac_1=1, rotfac_2=1,
            hf_1=1.5, hf_2=1.5,
            bfac_1=None, bfac_2=None,
            heat_1=None, heat_2=None,
            lambda_1=None, lambda_2=None,
            vsini_1=None, vsini_2=None, n_int=None,
            grid_1='default', grid_2='default',
            ld_1='quad', ld_2=None,
            shape_1='sphere', shape_2='sphere',
            spots_1=None, spots_2=None, flux_weighted=False, verbose=0)
        injected_rv = numpy.full(len(times), 0) + model_kms[0] * 1000 # mult by 1000 because the value comes in km/s
        return injected_rv

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
                       period_min=0.5, max_period=None, min_snr=5, min_sde=5, max_runs=5, cpus=os.cpu_count() - 1):
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
            msin_diff = numpy.abs(msin - signal_msin)
            omega_diff = numpy.abs(omega - signal_omega)
            omega_diff = omega_diff if omega_diff < numpy.pi else numpy.pi - numpy.abs(omega - signal_omega)
            found = period_diff / period < 0.01 #and msin_diff / msin < 0.5

            # rv_fit = RvFitter.sinfunc(df["bjd"], signal_k, signal_omega, signal_period)
            # import lightkurve
            # df1 = df.copy()
            # lc = lightkurve.LightCurve(time=df1['bjd'], flux=df1['rv'])
            # periodogram = lc.to_periodogram(oversample_factor=10, maximum_period=25, minimum_period=0.4)
            # periodogram.plot(scale='log')
            # plt.show()
            #plt.plot(periodogram.period.value, periodogram.power.value / numpy.std(periodogram.power))
            #plt.show()
            # df1["rv_fit"] = RvFitter.sinfunc(df1["bjd"], signal_k, signal_omega, signal_period)
            # df1["bjd_fold"] = (df1["bjd"] - df1["bjd"].min()) % signal_period
            # df1 = df1.sort_values(by=['bjd_fold'], ascending=True)
            # plt.errorbar(x=df1["bjd_fold"], y=df1["rv"], yerr=df1["rv_err"], fmt="o", color="blue", ecolor="cyan",
            #              label="RV (m/s)")
            # plt.plot(df1["bjd_fold"], df1["rv_fit"], color="black", label="X-axis motion")
            # plt.title("Fit for P=" + str(signal_period) + "d")
            # plt.xlabel("Phase")
            # plt.ylabel("RV (m/s)")
            # plt.show()

            if found or snr < min_snr or SDE < min_sde:
                break
            else:
                rv_fit = RvFitter.sinfunc(df["bjd"], signal_k, signal_omega, signal_period)
                df["rv"] = df["rv"] - rv_fit
        return found, run, snr, SDE, signal_period, signal_omega, signal_msin

    @staticmethod
    def recover_period(input: RecoverPeriodInput):
        k = 0
        k_err = 0
        omega = 0
        omega_err = 0
        least_squares = numpy.inf
        period = input.period
        try:
            if input.msin is not None and input.star_mass is not None:
                k_p0 = RvFitter.compute_semiamplitude_from_mmin(input.msin, input.period, input.star_mass)
                k_p0_max = k + 0.01
            else:
                k_p0_max = numpy.max(input.rv_data)
                k_p0 = 0
            for omega_p0 in numpy.linspace(0, numpy.pi * 2 * 0.75, 4):
                popt, pcov = scipy.optimize.curve_fit(RvFitter.sinfunc, input.time, input.rv_data, sigma=input.rv_err,
                                                      p0=[k_p0_max, omega_p0, input.period],
                                                      bounds=([k_p0, 0., input.period],
                                                              [k_p0_max, 2 * numpy.pi, input.period + 0.00001]))
                perr = numpy.sqrt(numpy.diag(pcov))
                k_err_tmp, omega_err_tmp, period_err_tmp = perr
                k_tmp, omega_tmp, period_tmp = popt
                least_squares_tmp = numpy.sum((input.rv_data -
                                           RvFitter.sinfunc(input.time, k_tmp, omega_tmp, period_tmp)) ** 2 /
                                          (input.rv_err ** 2))

                # rv_fit = RvFitter.sinfunc(input.time, k_tmp, omega_tmp, period_tmp)
                # import pandas as pd
                # df1 = pd.DataFrame(columns=['bjd', 'rv', 'rv_err'])
                # df1['bjd'] = input.time
                # df1['rv'] = input.rv_data
                # df1['rv_err'] = input.rv_err
                # df1["rv_fit"] = RvFitter.sinfunc(df1["bjd"], k_tmp, omega_tmp, period_tmp)
                # df1["bjd_fold"] = (df1["bjd"] - df1["bjd"].min()) % period_tmp
                # df1 = df1.sort_values(by=['bjd_fold'], ascending=True)
                # plt.errorbar(x=df1["bjd_fold"], y=df1["rv"], yerr=df1["rv_err"], fmt="o", color="blue", ecolor="cyan",
                #              label="RV (m/s)")
                # plt.plot(df1["bjd_fold"], df1["rv_fit"], color="black", label="X-axis motion")
                # plt.title("Fit for P=" + str(period_tmp) + "d")
                # plt.xlabel("Phase")
                # plt.ylabel("RV (m/s)")
                # plt.show()

                if least_squares_tmp < least_squares:
                    k = k_tmp
                    omega = omega_tmp
                    k_err = k_err_tmp
                    omega_err = omega_err_tmp
                    least_squares = least_squares_tmp
                    period = period_tmp
        except RuntimeError as e:
            pass
        return (period, k, omega, k_err, omega_err, least_squares)

    @staticmethod
    def recover_periods(rv_df, period_grid_geom='lin', steps_period=None, period_min=0.5, max_period=None,
                        rv_masks=None, star_mass=None, cpus=os.cpu_count() - 1):
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
        k = 0
        k_err = 0
        omega = 0
        omega_err = 0
        least_squares = 0
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
        lc = lightkurve.LightCurve(time=time, flux=rv_data)
        periodogram = lc.to_periodogram(oversample_factor=10, minimum_period=0.4, maximum_period=max_period)
        period_grid = periodogram.period.value
        power = periodogram.power.value
        sde_power = power / numpy.std(power)
        argmax_sde = numpy.nanargmax(sde_power)
        max_sde = sde_power[argmax_sde]
        strongest_period = period_grid[argmax_sde] #periodogram.periods[argmax_sde]
        k = 0
        k_err = 0
        omega = 0
        omega_err = 0
        least_squares = 0
        inputs = []
        strongest_period, k, omega, k_err, omega_err, least_squares = \
            RvFitter.recover_period(RecoverPeriodInput(time, rv_data, rv_err, strongest_period))
        msin = RvFitter.compute_mmin_from_semiamplitude([strongest_period], [k], star_mass)[0]
        snr_window_size = (max_period - period_min) / 20
        indexes_around = numpy.argwhere((period_grid > strongest_period - snr_window_size / 2) &
                                        (period_grid < strongest_period + snr_window_size / 2)).flatten()

        median_power_around = numpy.nanmedian(numpy.array(power)[indexes_around])
        snr = power[argmax_sde] / median_power_around
        sde = max_sde
        return rv_data, strongest_period, k, omega, msin, least_squares, argmax_sde, power, snr, sde

    @staticmethod
    def mask_signals(rv_df, rv_masks, star_mass):
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
