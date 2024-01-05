import logging
import os

import batman
from ellc import lc
import numpy as np
import pandas as pd
import astropy.constants as ac
import astropy.units as u


class InjectModel:
    """
    Class used to create a synthetic light curve.
    """
    def __init__(self, inject_dir, time, flux, flux_err, rstar, mstar, t0, period, rplanet, exposure_time, ab):
        self.inject_dir = inject_dir
        self.time = time
        self.flux = flux
        self.flux_err = flux_err
        self.rstar = rstar
        self.mstar = mstar
        self.t0 = t0
        self.period = period
        self.rplanet = rplanet
        self.exposure_time = exposure_time
        self.ab = ab

    @staticmethod
    def _map_planet_radius_to_mass(radius):
        max_mass = 40 * u.M_jup
        max_rad = 3 * u.R_jup
        min_mass = 1 * u.M_jup
        if radius < 1 * u.R_jup:
            mass = InjectModel.mass_from_radius(radius.value) * u.M_earth
        elif radius < max_rad:
            mass = radius.to(u.R_jup) * max_mass / max_rad
        else:
            mass = max_mass
        if mass > max_mass:
            mass = max_mass
        return mass.to(u.M_earth).value
    @staticmethod
    def make_model(inject_model):
        """
        Creates an injected transit curve given the provided InjectModel parameters instance.

        :param inject_model: the InjectModel object to be used as transit source. The result is written into a csv file.
        """
        logging.info('P = ' + str(inject_model.period) + ' days, Rp = ' +
                     str(inject_model.rplanet) + ", T0 = " +
                     str(inject_model.t0))
        P1 = inject_model.period * u.day
        a = np.cbrt((ac.G * inject_model.mstar * P1 ** 2) / (4 * np.pi ** 2)).to(u.au)
        texpo = inject_model.exposure_time / 60. / 60. / 24.
        planet_mass = InjectModel._map_planet_radius_to_mass(inject_model.rplanet) * u.M_earth
        mass_ratio = (planet_mass.to(u.M_sun) / inject_model.mstar).value
        model = lc(
            t_obs=inject_model.time,
            radius_1=inject_model.rstar.to(u.au) / a,  # star radius convert from AU to in units of a
            radius_2=inject_model.rplanet.to(u.au) / a,  # convert from Rearth (equatorial) into AU and then into units of a
            sbratio=0,
            incl=90,
            light_3=0,
            t_zero=inject_model.t0,
            period=inject_model.period,
            a=None,
            q=mass_ratio,
            f_c=None, f_s=None,
            ldc_1=inject_model.ab, ldc_2=None,
            gdc_1=None, gdc_2=None,
            didt=None,
            domdt=None,
            rotfac_1=1, rotfac_2=1,
            hf_1=1.5, hf_2=1.5,
            bfac_1=None, bfac_2=None,
            heat_1=None, heat_2=None,
            lambda_1=None, lambda_2=None,
            vsini_1=None, vsini_2=None,
            t_exp=texpo, n_int=None,
            grid_1='default', grid_2='default',
            ld_1='quad', ld_2=None,
            shape_1='sphere', shape_2='sphere',
            spots_1=None, spots_2=None,
            exact_grav=False, verbose=0)
        if not(model[0] > 0):
            params = batman.TransitParams()  # object to store transit parameters
            params.t0 = inject_model.t0  # time of inferior conjunction
            params.per = inject_model.period  # orbital period
            params.rp = inject_model.rplanet.to(u.R_sun).value / inject_model.rstar.value
            params.a = (a.to(u.R_sun) / inject_model.rstar).value
            params.inc = 90.  # orbital inclination (in degrees)
            params.ecc = 0.  # eccentricity
            params.w = 90.  # longitude of periastron (in degrees)
            params.limb_dark = "quadratic"  # limb darkening model
            params.u = inject_model.ab  # limb darkening coefficients [u1, u2, u3, u4]
            model = batman.TransitModel(params, inject_model.time)  # initializes model
            model = model.light_curve(params)
        if not (model[0] > 0):
            result_flux = []
            result_time = []
            result_flux_err = []
        else:
            flux_t = np.array(inject_model.flux) + model - 1.
            result_flux = flux_t
            result_flux_err = np.array(inject_model.flux_err)
            result_time = np.array(inject_model.time)
        file_name = os.path.join(inject_model.inject_dir + '/P' + f'{inject_model.period:06.2f}' + '_R' +
                                 f'{inject_model.rplanet.value:05.2f}' + '_T' + str(inject_model.t0) + '.csv')
        lc_df = pd.DataFrame(columns=['#time', 'flux', 'flux_err'])
        lc_df['#time'] = result_time
        lc_df['flux'] = result_flux
        lc_df['flux_err'] = result_flux_err
        lc_df.to_csv(file_name, index=False)

    @staticmethod
    def mass_from_radius(radius):
        """
        Computation of mass-radius relationship from
        Bashi D., Helled R., Zucker S., Mordasini C., 2017, A&A, 604, A83. doi:10.1051/0004-6361/201629922

        :param radius: the radius value in earth radius
        :return: the mass in earth masses
        """
        return radius ** (1 / 0.55) if radius <= 12.1 else radius ** (1 / 0.01)
