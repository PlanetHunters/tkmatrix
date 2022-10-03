import multiprocessing
from multiprocessing import Pool

import ellc
import numpy
import pandas
import scipy
import matplotlib.pyplot as plt
import astropy.constants as ac
import astropy.units as u
from scipy.optimize import differential_evolution
import lightkurve
from tkmatrix.tkmatrix_class import MATRIX


#MATRIX.plot_results("TIC 142748283", "/home/martin/Desktop/run_tests/tests/TIC142748283_ir_3", 1)
#MATRIX.plot_diff("TIC 169285097", "/home/martin/Desktop/run_tests/tests/TIC169285097_ir_sasg",
#                 "/home/martin/Desktop/run_tests/tests/TIC169285097_ir_none",
#                 "./")


class RvFitter:
    def __init__(self) -> None:
        super().__init__()

    @staticmethod
    def inject_rv(times, star_mass, rstar, planet_mass, rplanet, period, t0):
        star_mass = star_mass * u.M_sun
        rstar = rstar * u.R_sun
        rplanet = rplanet * u.R_sun
        planet_mass = planet_mass * u.M_sun
        period_days = period * u.day
        a = numpy.cbrt((ac.G * star_mass * period_days ** 2) / (4 * numpy.pi ** 2)).to(u.au)
        model = ellc.rv(
            t_obs=times,
            radius_1=rstar.to(u.au) / a,  # star radius convert from AU to in units of a
            radius_2=rplanet.to(u.au) / a,  # convert from Rearth (equatorial) into AU and then into units of a
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
            spots_1=None, spots_2=None, flux_weighted=False, verbose=1)
        injected_rv = numpy.full(len(rv_df["bjd"]), 0) + model[1]
        return injected_rv

    @staticmethod
    def sinfunc(t, k, omega, period):
        return k * numpy.cos(omega + t * 2 * numpy.pi / period)

    @staticmethod
    def compute_mmin_from_semiamplitude(period_grid, k_grid, star_mass):
        m_min = []
        for i in numpy.arange(0, len(period_grid), 1):
            m_min_value = k_grid[i] * (star_mass ** (2./3.)) * (period_grid[i] ** (1./3.)) / 0.64
            m_min.append(m_min_value)
        return m_min

    @staticmethod
    def running_median(data, kernel):
        """Returns sliding median of width 'kernel' and same length as data """
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
    def recover_signal(rv_df, period, msin, omega, period_grid_geom='lin', steps_period=None, period_min=0.5, max_period=None,
                       min_snr=5, min_sde=5, max_runs=5):
        sde = numpy.inf
        snr = numpy.inf
        found = False
        for run in numpy.arange(0, max_runs):
            period_grid, k_grid, omega_grid, msin_grid, least_squares_grid, argmax, power, snr, SDE = \
                RvFitter.recover_periods(rv_df, period_grid_geom, steps_period, period_min, max_period)
            signal_omega = omega_grid[argmax]
            signal_k = k_grid[argmax]
            signal_msin = msin_grid[argmax]
            signal_period = period_grid[argmax]
            period_diff = numpy.abs(period - signal_period)
            msin_diff = numpy.abs(msin - signal_msin)
            omega_diff = numpy.abs(omega - signal_omega)
            omega_diff = omega_diff if omega_diff < numpy.pi else numpy.pi - numpy.abs(omega - signal_omega)
            found = period_diff < period * 0.01
            if found or snr < min_snr or SDE < min_sde:
                break
            else:
                rv_fit = RvFitter.sinfunc(rv_df["bjd"], signal_k, signal_omega, signal_period)
                rv_df["rv"] = rv_df["rv"] - rv_fit
        return found, snr, sde, signal_period, signal_omega, signal_msin

    @staticmethod
    def recover_period(input):
        k = 0
        k_err = 0
        omega = 0
        omega_err = 0
        least_squares = 0
        try:
            popt, pcov = scipy.optimize.curve_fit(RvFitter.sinfunc, rv_df["bjd"], input.rv_data, sigma=rv_df["rv_err"],
                                                  bounds=([0., 0., input.period], [numpy.inf, numpy.pi, input.period + 0.00001]))
            perr = numpy.sqrt(numpy.diag(pcov))
            k_err, omega_err, period_err = perr
            k, omega, period = popt
            least_squares = numpy.sum((input.rv_data - RvFitter.sinfunc(rv_df["bjd"], k, omega, period)) ** 2 /
                                      (input.rv_err ** 2))
        except RuntimeError as e:
            print("Runtime error for period " + str(input.period))
        return (k, omega, k_err, omega_err, least_squares)

    @staticmethod
    def recover_periods(rv_df, period_grid_geom='lin', steps_period=None, period_min=0.5, max_period=None,
                        cpus=multiprocessing.cpu_count() - 1):
        if steps_period is None:
            steps_period = int(rv_df["bjd"].max() - rv_df["bjd"].min())
        if max_period is None:
            max_period = (rv_df["bjd"].max() - rv_df["bjd"].min()) * 2
        period_grid = numpy.linspace(min_period, max_period, steps_period) if period_grid_geom == "lin" \
            else numpy.logspace(numpy.log10(min_period), numpy.log10(max_period), steps_period)
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
        time = rv_df["bjd"].to_numpy()
        rv_data = rv_df["rv"].to_numpy()
        rv_err = rv_df["rv_err"].to_numpy()
        inputs = []
        for period in period_grid:
            inputs.append(RecoverPeriodInput(time, rv_data, rv_err, period))
        with Pool(processes=cpus) as pool:
            recovered_outputs = pool.map(RvFitter.recover_period, inputs)
        for output in recovered_outputs:
            k_grid.append(output[0])
            omega_grid.append(output[1])
            k_err_grid.append(output[2])
            omega_err_grid.append(output[3])
            least_squares_grid.append(output[4])
        least_squares_grid = numpy.array(least_squares_grid)
        least_squares_grid[least_squares_grid == 0] = numpy.max(least_squares_grid)
        msin_grid = RvFitter.compute_mmin_from_semiamplitude(period_grid, k_grid, star_mass)
        SR, power_raw, power, SDE_raw, SDE = RvFitter.spectra(numpy.array(least_squares_grid))
        sde_power = power / numpy.std(power)
        argmax_sde = numpy.nanargmax(sde_power)
        max_msin = msin_grid[argmax_sde]
        period_max_msin = period_grid[argmax_sde]
        snr_window_size = (max_period - period_min) / 20
        indexes_around = numpy.argwhere((period_grid > period_max_msin - snr_window_size / 2) &
                                        (period_grid < period_max_msin + snr_window_size / 2)).flatten()
        median_power_around = numpy.nanmedian(numpy.array(msin_grid)[indexes_around])
        snr = msin_grid[argmax_sde] / median_power_around
        return period_grid, k_grid, omega_grid, msin_grid, least_squares_grid, argmax_sde, power, snr, SDE

class RecoverPeriodInput:
    def __init__(self, time, rv_data, rv_err, period):
        self.time = time
        self.rv_data = rv_data
        self.rv_err = rv_err
        self.period = period

import pandas as pd
import numpy as np
harps = pd.read_csv("examples/toi175harps.csv")
harps['rv'] = harps['rv'] - np.nanmean(harps['rv'])
espresso = pd.read_csv("examples/toi175espresso.csv")
espresso['rv'] = espresso['rv'] - np.nanmean(espresso['rv'])
rv_df = pd.DataFrame(["bjd", "rv", "rv_err"])
rv_df = rv_df.append(harps.loc[:, ['bjd', 'rv', 'rv_err']])
rv_df = rv_df.append(espresso.loc[:, ['bjd', 'rv', 'rv_err']])
star_mass = 0.31
star_radius = 0.31
#rv_df = pandas.read_csv("examples/rv.csv")
# star_mass = 0.463011
# star_radius = 0.465559
planet_mass = 0.01
planet_radius = 0.05
min_period = 0.5
max_period = 20
steps_period = 2000
period_grid_geom = 'lin'
period_grid = numpy.linspace(min_period, max_period, steps_period) if period_grid_geom == "lin" \
            else numpy.logspace(numpy.log10(min_period), numpy.log10(max_period), steps_period)
rv_df.drop(0, axis=1, inplace=True)
rv_df = rv_df.dropna()
rv_df = rv_df.sort_values(by=['bjd'], ascending=True)
period_grid, k_grid, omega_grid, msin_grid, least_squares_grid, argmax_sde, power, snr, SDE = \
   RvFitter.recover_periods(rv_df, steps_period=steps_period, max_period=max_period)
inject_period = 54
#injected_rv = RvFitter.inject_rv(rv_df["bjd"], star_mass, star_radius, planet_mass, planet_mass, inject_period, rv_df["bjd"][0])
#rv_df["rv"] = injected_rv
#RvFitter.recover_signal(rv_df, inject_period, 3, 0.5, 'lin', steps_period, min_period, max_period, 3, 5, 5)

from scipy.stats import binned_statistic
bin_means_msin = binned_statistic(period_grid, msin_grid, bins=50)[0]
bin_means_period = binned_statistic(period_grid, period_grid, bins=50)[0]
plt.bar(bin_means_period, bin_means_msin, width=1.0, color="black")
plt.xscale('log', base=10)
plt.title("Mmin detections")
plt.xlabel("Period (d)")
plt.ylabel("Mmin")

plt.show()
plt.plot(period_grid, msin_grid, color="black")
plt.xscale('log', base=10)
plt.title("Mmin detections")
plt.xlabel("Period (d)")
plt.ylabel("Mmin")
plt.show()

plt.plot(period_grid, power / numpy.std(power), color="black", label="X-axis motion")
plt.xscale('log', base=10)
plt.title("SDE detections")
plt.xlabel("Period (d)")
plt.ylabel("SDE")
plt.show()

top_period = period_grid[argmax_sde]
rv_df["rv_fit"] = RvFitter.sinfunc(rv_df["bjd"], k_grid[argmax_sde], omega_grid[argmax_sde], top_period)
rv_df["bjd_fold"] = (rv_df["bjd"] - rv_df["bjd"].min()) % top_period
rv_df = rv_df.sort_values(by=['bjd_fold'], ascending=True)
plt.errorbar(x=rv_df["bjd_fold"], y=rv_df["rv"], yerr=rv_df["rv_err"], fmt="o", color="blue", ecolor="cyan",
             label="RV (m/s)")
plt.plot(rv_df["bjd_fold"], rv_df["rv_fit"], color="black", label="X-axis motion")
plt.title("SDE detections")
plt.xlabel("Phase")
plt.ylabel("RV (m/s)")
plt.show()

lc = lightkurve.LightCurve(time=rv_df['bjd'], flux=rv_df['rv'])
periodogram = lc.to_periodogram(oversample_factor=10, maximum_period=max_period, minimum_period=min_period)
periodogram.plot()
plt.show()

plt.scatter(rv_df["bjd"], injected_rv, color="red", label="X-axis motion")
plt.scatter(rv_df["bjd"], rv_df["rv"], color="gray", label="X-axis motion")
plt.errorbar(x=rv_df["bjd"], y=injected_rv, yerr=rv_df["rv_err"], fmt="o", color="blue", ecolor="cyan",
             label="RV (m/s)")
plt.plot(plot_time, fit_rv, color="orange", label="RV Fit")
plt.title("RV best fit")
plt.xlabel("Time (d)")
plt.ylabel("RV (m/s)")
plt.show()
a = 1
