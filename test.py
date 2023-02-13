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

from tkmatrix.rv import RvFitter
from tkmatrix.tkmatrix_class import MATRIX


#MATRIX.plot_results("TIC 142748283", "/home/martin/Desktop/run_tests/tests/TIC142748283_ir_3", 1)
#MATRIX.plot_diff("TIC 169285097", "/home/martin/Desktop/run_tests/tests/TIC169285097_ir_sasg",
#                 "/home/martin/Desktop/run_tests/tests/TIC169285097_ir_none",
#                 "./")

#
# MATRIX.plot_results("TOI 175", "TIC307210830_ir_3/", binning=1, xticks=None, yticks=None, period_grid_geom="lin",
#                      radius_grid_geom="lin", is_rv=True)

import pandas as pd
import numpy as np
star_mass = 0.31
star_radius = 0.31
rv_masks = [{"P": 2.2531136}, {"P": 3.6906777}, {"P": 7.4507245}, {"P": 12.796}]
v0 = -5578.51
# harps = pd.read_csv("examples/toi175harps.csv")
# harps['rv'] = harps['rv'] - np.nanmean(harps['rv'])
# espresso = pd.read_csv("examples/toi175espresso.csv")
# espresso['rv'] = espresso['rv'] - np.nanmean(espresso['rv'])
# rv_df = pd.DataFrame(["bjd", "rv", "rv_err"])
# rv_df = rv_df.append(harps.loc[:, ['bjd', 'rv', 'rv_err']])
# rv_df = rv_df.append(espresso.loc[:, ['bjd', 'rv', 'rv_err']])
rv_df = pandas.read_csv("examples/toi175rv.csv")
# star_mass = 0.463011
# star_radius = 0.465559
planet_mass = 0.00003
planet_radius = 0.1
min_period = 0.5
max_period = 20
steps_period = 3000
period_grid_geom = 'log'
period_grid = numpy.linspace(min_period, max_period, steps_period) if period_grid_geom == "lin" \
            else numpy.logspace(numpy.log10(min_period), numpy.log10(max_period), steps_period)
#rv_df.drop(0, axis=1, inplace=True)
rv_df = rv_df.dropna()
rv_df = rv_df.sort_values(by=['bjd'], ascending=True)
#rv_df.to_csv("examples/toi175rv.csv")

lc = lightkurve.LightCurve(time=rv_df['bjd'], flux=rv_df['rv'])
periodogram = lc.to_periodogram(oversample_factor=10, maximum_period=max_period, minimum_period=min_period)
periodogram.plot(scale='log')
plt.show()

rv_data, period_grid, k_grid, omega_grid, msin_grid, least_squares_grid, argmax_sde, power, snr, SDE = \
   RvFitter.recover_periods(rv_df, period_grid_geom=period_grid_geom, steps_period=steps_period, max_period=max_period,
                            rv_masks=rv_masks, star_mass=star_mass)
rv_df['rv'] = rv_data
inject_period = 12.75
injected_rv = RvFitter.inject_rv(rv_df["bjd"].to_numpy(), star_mass * u.M_sun, star_radius * u.R_sun, planet_mass / 3.00273e-6 * u.M_earth, inject_period, rv_df["bjd"].to_numpy()[0])
rv_df["rv"] = injected_rv
found, run, snr, SDE, signal_period, signal_omega, signal_msin =\
    RvFitter.recover_signal(rv_df, inject_period, planet_mass / 3.00273e-6, 0.5, rv_masks, star_mass,
                            'lin', steps_period, min_period, max_period, 3, 5, 5)

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
plt.title("Fit for P=" + str(top_period) + "d, M=" + str(msin_grid[argmax_sde]) + "Me")
plt.xlabel("Phase")
plt.ylabel("RV (m/s)")
plt.show()

lc = lightkurve.LightCurve(time=rv_df['bjd'], flux=rv_df['rv'])
periodogram = lc.to_periodogram(oversample_factor=10, maximum_period=max_period, minimum_period=min_period)
periodogram.plot(scale='log')
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
