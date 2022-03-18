import os

import math
import glob
import numpy as np
import pylab as pl
import pandas as pd
from iminuit import Minuit
import iminuit
#from astropy import units as u
#from astropy import constants as const
#PJA uncomments these two previous imports
from astropy import units as u
from astropy import constants as const
import time

#from scipy import stats
#from scipy import optimize
#from scipy.optimize import minimize

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as colors


# ## Classes
class nf(float):
    def __repr__(self):
        s = f'{self:.1f}'
        return f'{self:.0f}' if s[-1] == '0' else s


class RvFitter:

    def __init__(self) -> None:
        super().__init__()

    @staticmethod
    def axisToPeriod(a, M_star):
        """
        Converts from semimajor axis to period in days
        :param M_star: the star mass
        :return: the period in days
        """
        return np.sqrt(4. * (np.pi**2.) * (a ** 3.) / ((const.G.to(u.AU ** 3 * u.solMass ** (-1) * u.d ** (-2))).value * M_star))

    @staticmethod
    def grid_logSpaced(num_measures, order_ini, order_fin_in):
        """
        Creates a log-10 period grid
        :param num_measures: number of measurements (density of the grid)
        :param order_ini: the initial power of the logarithm
        :param order_fin_in: the final power of the logarithm
        :return: the period grid
        """
        grid_list = []
        for order_fin in np.arange(order_ini + 1, order_fin_in, 1):
            if order_fin < order_fin_in - 1:
                grid_list.extend(np.logspace(order_ini, order_fin, num_measures, endpoint=False))
            else:
                grid_list.extend(np.logspace(order_ini, order_fin, num_measures, endpoint=True))
            order_ini = order_fin
        return grid_list

    @staticmethod
    def plot_rv(file_output, rv, rv_error, bjd):
        """
        Plots the RV data against the BJD
        :param file_output: Output file
        :param rv: the RV values
        :param rv_error: the RV error values
        :param bjd: the BJD data
        """
        fig = plt.figure(figsize=(10, 6), dpi=200) # Create a figure of size 8x6 inches, 200 dots per inch
        ax = fig.add_subplot(1, 1, 1)
        ax.errorbar(bjd, rv, yerr=rv_error, fmt="o")
        ax.set_title("RVs with RMS (m/s) = {}".format(round(np.sqrt(np.mean(rv ** 2)), 2)))
        ax.set_xlabel('BJD (days)', fontsize=15)
        ax.set_ylabel('RV (m/s)', fontsize=15)
        plt.savefig(file_output)
        plt.show()

    @staticmethod
    def plot_fit(tau, m_min, fit_file):
        """
        Plot the fits results
        :param tau: period grid
        :param m_min: the Mass_min values
        :param fit_file: the file of the output
        """
        fig = plt.figure(figsize=(10, 6), dpi=200)
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(tau, m_min, "-", color='Black', linewidth=0.75, label="Thin fit, M$_m$(P)")
        plt.plot([0, 50], [3, 3], color="red", linewidth=1.0, linestyle="-.")
        plt.plot([0, 50], [10, 10], color="red", linewidth=1.0, linestyle="-.")
        plt.plot([0, 50], [30, 30], color="red", linewidth=1.0, linestyle="-.")
        plt.plot([0, 50], [100, 100], color="red", linewidth=1.0, linestyle="-.")
        plt.plot([0, 50], [300, 300], color="red", linewidth=1.0, linestyle="-.")
        plt.plot([0, 50], [1000, 1000], color="red", linewidth=1.0, linestyle="-.")
        plt.plot([50, 50], [3, 1000], color="red", linewidth=1.0, linestyle="-.")
        plt.text(60, 300, '50 days', fontsize=10, rotation=90, rotation_mode='anchor')
        ax.set_title("M_min fits")
        ax.set_xlabel('$P_{planeta} (days)$', fontsize=15)
        ax.set_ylabel('$Max(Msin i/M_{\oplus})$', fontsize=15)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(1, 10 ** 4)
        ax.set_ylim(0.1, 13 ** 3)
        ax.legend(loc='upper right', fontsize=10)
        plt.savefig(fit_file)
        plt.show()

    @staticmethod
    def compute_mmin_from_semiamplitude(k_fit, m_s, tau):
        """
        Computes the M_min values for each k_fit-tau pairs.
        :param k_fit: the fitted semi-amplitude
        :param m_s: the star mass
        :param tau: the period grid
        :return: the M_min values
        """
        m_min = []
        for i in np.arange(0, len(tau), 1):
            m_min_value = k_fit[i] * (m_s ** (2./3.)) * (tau[i] ** (1./3.)) / 0.64
            m_min.append(m_min_value)
        return m_min

    @staticmethod
    def least_squares_search_comp(period, k, omega, bjd, rv, rv_error):
        """
        Function used to fit the semi-amplitude k.
        :param period: the period for the fit
        :param k: the semi-amplitude
        :param omega: the initial phase
        :param bjd: the BJD time
        :param rv: the RV data
        :param rv_error: the RV error data
        :return: the least-squares sum to be minimized
        """
        return sum((rv - k * np.cos(omega + bjd * 2 * np.pi / period)) ** 2 / rv_error ** 2)

    @staticmethod
    def k_fit_funct(path_output, massStarList, rv, rv_error, bjd, tau):
        """
        Function that calculates a fit for RV data using iMINUIT minimizator (https://iminuit.readthedocs.io). There are also
        some other functions defined inside for plotting and other purpose. We decided to define them here because this allows
        better executions times than define them outside

        Parameters:
        -----------
        path_output (str)
        pathOutFit (str)
        massStarList (list)
        len_file (int)
        rv(list of list)
        rv_error(list of list)
        bjd(listof list)
        tau(list)
        tau_coarse (list)
        tau_coarse_plot(list)
        names_plot (list)

        Returns:
        -----------

        k_fit_list, mMinList, mEnvelopeList, mEnvelopePlotList
        and some plots in the indicated paths
        """
        k_fit = []
        rv_meas = np.array(rv)
        rv_meas_error = np.array(rv_error)
        bjd = np.array(bjd)
        plot_file = path_output + "/RV.png"
        fit_file = path_output + "/M_min_fit.png"
        RvFitter.plot_rv(plot_file, rv_meas, rv_meas_error, bjd)
        #TODO parallelize
        for index, orb_period in enumerate(tau):
            min_period = tau[index]
            max_period = tau[index]
            m = Minuit(RvFitter.least_squares_search_comp,
                       period=orb_period, k=1, omega=1, limit_period=(min_period, max_period), rv=rv_meas,
                       rv_error=rv_meas_error, pedantic=False)
            m.migrad()  # finds minimum of least_squares function
            m.hesse()  # computes errors
            if not m.get_fmin()[7]:  # if enters here it means that minimization do not converge
                # CAREFULL if it does not enter does not mean the minimization is good...
                print("Warning do not converge. Revise the program or the star")
                print("Period = ", orb_period)
                #                pprint(m.get_fmin())
                break
            k_fit.append(np.abs(m.values["k"]))
        m_min = RvFitter.compute_mmin_from_semiamplitude(k_fit, massStarList, tau)
        #TODO store M_min together to K
        RvFitter.plot_fit(tau, m_min, fit_file)
        return k_fit, m_min
