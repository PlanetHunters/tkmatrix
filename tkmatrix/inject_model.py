import logging
import os

import ellc
import numpy as np
import pandas as pd
import astropy.constants as ac
import astropy.units as u


class InjectModel:
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
    def make_model(inject_model):
        logging.info('P = ' + str(inject_model.period) + ' days, Rp = ' + str(inject_model.rplanet) + ", T0 = " +
                     str(inject_model.t0))
        P1 = inject_model.period * u.day
        a = np.cbrt((ac.G * inject_model.mstar * P1 ** 2) / (4 * np.pi ** 2)).to(u.au)
        texpo = inject_model.exposure_time / 60. / 60. / 24.
        model = ellc.lc(
            t_obs=inject_model.time,
            radius_1=inject_model.rstar.to(u.au) / a,  # star radius convert from AU to in units of a
            radius_2=inject_model.rplanet.to(u.au) / a,  # convert from Rearth (equatorial) into AU and then into units of a
            sbratio=0,
            incl=90,
            light_3=0,
            t_zero=inject_model.t0,
            period=inject_model.period,
            a=None,
            q=1e-6,
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
            exact_grav=False, verbose=1)
        if model[0] > 0:
            flux_t = inject_model.flux + model - 1.
            result_flux = flux_t
            result_flux_err = inject_model.flux_err
            result_time = inject_model.time
        else:
            result_flux = []
            result_time = []
            result_flux_err = []
        file_name = os.path.join(inject_model.inject_dir + '/P' + str(inject_model.period) + '_R' +
                                 str(inject_model.rplanet.value) + '_' + str(inject_model.t0) + '.csv')
        lc_df = pd.DataFrame(columns=['#time', 'flux', 'flux_err'])
        lc_df['#time'] = result_time
        lc_df['flux'] = result_flux
        lc_df['flux_err'] = result_flux_err
        lc_df.to_csv(file_name, index=False)
