import logging
import os
import pandas as pd

from tkmatrix.rv import RvFitter


class InjectRvModel:
    """
    Used to inject a Radial Velocity signal
    """
    def __init__(self, inject_dir, time, rv, rv_err, rstar, mstar, t0, period, mplanet):
        self.inject_dir = inject_dir
        self.time = time
        self.rv = rv
        self.rv_err = rv_err
        self.rstar = rstar
        self.mstar = mstar
        self.t0 = t0
        self.period = period
        self.mplanet = mplanet

    @staticmethod
    def make_model(inject_model):
        """
        Creates an injected RV curve given the provided InjectRvModel parameters instance.

        :param inject_model: the InjectRvModel object to be used as RV source. The result is written into a csv file.
        """
        logging.info('RV P = ' + str(inject_model.period) + ' days, Mp = ' + str(inject_model.mplanet) + ", T0 = " +
                     str(inject_model.t0))
        injected_rv = RvFitter.inject_rv(inject_model.time, inject_model.mstar, inject_model.rstar,
                                         inject_model.mplanet, inject_model.period, inject_model.t0)
        file_name = os.path.join(inject_model.inject_dir + '/RV_P' + str(inject_model.period) + '_M' +
                                 str(inject_model.mplanet.value) + '_' + str(inject_model.t0) + '.csv')
        lc_df = pd.DataFrame(columns=['bjd', 'rv', 'rv_err'])
        lc_df['bjd'] = inject_model.time
        lc_df['rv'] = inject_model.rv + injected_rv
        lc_df['rv_err'] = inject_model.rv_err
        lc_df.to_csv(file_name, index=False)
