import traceback
import numpy as np
import ellc
import matplotlib.pyplot as plt
from tkmatrix.objectinfo.MissionInputObjectInfo import MissionInputObjectInfo
from tkmatrix.objectinfo.preparer.MissionLightcurveBuilder import MissionLightcurveBuilder
from tkmatrix.objectinfo.MissionObjectInfo import MissionObjectInfo
import wotan
from transitleastsquares import transitleastsquares
from transitleastsquares import transit_mask, cleaned_array
import astropy.constants as ac
import astropy.units as u
import lightkurve as lk
import os
import re
import pandas as pd
import shutil


class TIRMA:
    """
    TIRMA: Transit Inject and Recovery with Multi-phase Analysis
    """
    object_info = None

    def __init__(self, target, sectors, dir):
        assert target is not None and isinstance(target, str)
        assert sectors is not None and (sectors == 'all' or isinstance(sectors, list))
        self.id = target
        self.dir = dir
        self.sectors = sectors

    def retrieve_object_data(self):
        self.object_info = MissionObjectInfo(self.id, self.sectors)
        lightcurve_builder = MissionLightcurveBuilder()
        self.lc, self.lc_data, self.star_info, self.transits_min_count, self.sectors, self.quarters = \
            lightcurve_builder.build(self.object_info, None)
        self.ab = self.star_info.ld_coefficients
        self.mass = self.star_info.mass
        self.massmin = self.star_info.mass_min
        self.massmax = self.star_info.mass_max
        self.radius = self.star_info.radius
        self.radiusmin = self.star_info.radius_min
        self.radiusmax = self.star_info.radius_max
        # units for ellc
        self.rstar = self.star_info.radius * u.R_sun
        self.mstar = self.star_info.mass * u.M_sun
        self.mstar_min = self.star_info.mass_min * u.M_sun
        self.mstar_max = self.star_info.mass_max * u.M_sun
        self.rstar_min = self.star_info.radius_min * u.R_sun
        self.rstar_max = self.star_info.radius_max * u.R_sun
        print('\n STELLAR PROPERTIES FOR THE SIGNAL SEARCH')
        print('================================================\n')
        print('limb-darkening estimates using quadratic LD (a,b)=', self.ab)
        print('mass =', format(self.mstar,'0.5f'))
        print('mass_min =', format(self.mstar_min,'0.5f'))
        print('mass_max =', format(self.mstar_max,'0.5f'))
        print('radius =', format(self.rstar,'0.5f'))
        print('radius_min =', format(self.rstar_min,'0.5f'))
        print('radius_max =', format(self.rstar_max,'0.5f'))

    def build_inject_dir(self):
        return self.dir + "/" + self.object_info.mission_id().replace(" ", "") + "_ir/"

    def inject(self, phases, min_period, max_period, step_period, min_radius, max_radius, step_radius):
        assert phases is not None and isinstance(phases, int)
        assert min_period is not None and isinstance(min_period, (int, float))
        assert max_period is not None and isinstance(max_period, (int, float))
        assert step_period is not None and isinstance(step_period, (int, float))
        assert min_radius is not None and isinstance(min_radius, (int, float))
        assert max_radius is not None and isinstance(max_radius, (int, float))
        assert step_radius is not None and isinstance(step_radius, (int, float))
        self.retrieve_object_data()
        lc_new = lk.LightCurve(time=self.lc.time, flux=self.lc.flux, flux_err=self.lc.flux_err)
        clean = lc_new.remove_outliers(sigma_lower=float('inf'), sigma_upper=3)
        flux0 = clean.flux.value
        time = clean.time.value
        flux_err = clean.flux_err.value
        inject_dir = self.build_inject_dir()
        if os.path.isdir(inject_dir):
            shutil.rmtree(inject_dir)
        os.mkdir(inject_dir)
        for period in np.arange(min_period, max_period + 0.01, step_period):
            for t0 in np.arange(time[60], time[60] + period - 0.1, period / phases):
                for rplanet in np.arange(min_radius, max_radius + 0.01, step_radius):
                    rplanet = np.around(rplanet, decimals=2) * u.R_earth
                    print('\n')
                    print('P = ' + str(period) + ' days, Rp = ' + str(rplanet) + ", T0 = " + str(t0))
                    time_model, flux_model, flux_err_model = self.__make_model(time, flux0, flux_err, self.rstar,
                                                                               self.mstar, t0, period, rplanet)
                    file_name = os.path.join(inject_dir + '/P' + str(period) + '_R' + str(rplanet.value) + '_' + str(t0) +
                                             '.csv')
                    lc_df = pd.DataFrame(columns=['#time', 'flux', 'flux_err'])
                    lc_df['#time'] = time_model
                    lc_df['flux'] = flux_model
                    lc_df['flux_err'] = flux_err_model
                    lc_df.to_csv(file_name, index=False)
        return self

    def recovery(self, cores, sherlock_samples=0, known_transits=None, detrend_ws=0):
        assert detrend_ws is not None and isinstance(detrend_ws, (int, float))
        assert cores is not None
        assert known_transits is None or isinstance(known_transits, list)
        if self.object_info is None:
            self.retrieve_object_data()
        print('\n STELLAR PROPERTIES FOR THE SIGNAL SEARCH')
        print('================================================\n')
        print('limb-darkening estimates using quadratic LD (a,b)=', self.ab)
        print('mass =', format(self.mstar, '0.5f'))
        print('mass_min =', format(self.mstar_min, '0.5f'))
        print('mass_max =', format(self.mstar_max, '0.5f'))
        print('radius =', format(self.rstar, '0.5f'))
        print('radius_min =', format(self.rstar_min, '0.5f'))
        print('radius_max =', format(self.rstar_max, '0.5f'))
        reports_df = pd.DataFrame(columns=['period', 'radius', 'epoch', 'found', 'snr', 'sde', 'run'])
        inject_dir = self.build_inject_dir()
        for file in os.listdir(inject_dir):
            if file.endswith(".csv"):
                try:
                    period = float(re.search("P([0-9]+\\.[0-9]+)", file)[1])
                    r_planet = float(re.search("R([0-9]+\\.[0-9]+)", file)[1])
                    epoch = float(re.search("_([0-9]+\\.[0-9]+)\\.csv", file)[1])
                    df = pd.read_csv(inject_dir + file, float_precision='round_trip', sep=',',
                                     usecols=['#time', 'flux', 'flux_err'])
                    if len(df) == 0:
                        found = True
                        snr = 20
                        sde = 20
                        run = 1
                    else:
                        lc = lk.LightCurve(time=df['#time'], flux=df['flux'], flux_err=df['flux_err'])
                        clean = lc.remove_nans().remove_outliers(sigma_lower=float('inf'),
                                                                 sigma_upper=3)  # remove outliers over 3sigma
                        flux = clean.flux.value
                        time = clean.time.value
                        intransit = self.transit_masks(known_transits, time)
                        found, snr, sde, run = self.__tls_search(time, flux, self.radius, self.radiusmin,
                                                                 self.radiusmax, self.mass, self.massmin,
                                                                 self.massmax, self.ab, intransit, epoch, period, 0.5,
                                                                 time[len(time) - 1] - time[0], 5, cores, "default",
                                                                 detrend_ws, self.transits_min_count)
                    new_report = {"period": period, "radius": r_planet, "epoch": epoch, "found": found, "snr": snr,
                                  "sde": sde, "run": run}
                    reports_df = reports_df.append(new_report, ignore_index=True)
                    print("P=" + str(period) + ", R=" + str(r_planet) + ", T0=" + str(epoch) + ", FOUND WAS " + str(
                        found) +
                          " WITH SNR " + str(snr) + " AND SDE " + str(sde))
                    reports_df.to_csv(inject_dir + "a_tls_report.csv", index=False)
                except Exception as e:
                    traceback.print_exc()
                    print("File not valid: " + file)
        tls_report_df = pd.read_csv(inject_dir + "a_tls_report.csv", float_precision='round_trip', sep=',',
                                    usecols=['period', 'radius', 'epoch', 'found', 'snr', 'sde', 'run'])
        period_count = tls_report_df["period"].nunique()
        radius_count = tls_report_df["radius"].nunique()
        phases = len(tls_report_df) / (period_count * radius_count)
        period_count = period_count if period_count == 1 else period_count - 1
        radius_count = radius_count if radius_count == 1 else radius_count - 1
        step_period = (tls_report_df["period"].max() - tls_report_df["period"].min()) / period_count
        step_radius = (tls_report_df["radius"].max() - tls_report_df["radius"].min()) / radius_count
        self.plot_results(step_period, step_radius, phases)
        if sherlock_samples > 0:
            from sherlockpipe import sherlock
            from sherlockpipe.scoring.QuorumSnrBorderCorrectedSignalSelector import QuorumSnrBorderCorrectedSignalSelector

            class QuorumSnrBorderCorrectedStopWhenMatchSignalSelector(QuorumSnrBorderCorrectedSignalSelector):
                def __init__(self, strength=1, min_quorum=0, per=None, t0=None):
                    super().__init__()
                    self.strength = strength
                    self.min_quorum = min_quorum
                    self.per = per
                    self.t0 = t0

                def select(self, transit_results, snr_min, detrend_method, wl):
                    signal_selection = super(QuorumSnrBorderCorrectedStopWhenMatchSignalSelector, self) \
                        .select(transit_results, snr_min, detrend_method, wl)
                    if signal_selection.score == 0 or (
                            self.is_harmonic(signal_selection.transit_result.period, self.per) and
                            self.isRightEpoch(signal_selection.transit_result.t0, self.t0, self.per)):
                        signal_selection.score = 0
                    return signal_selection

                def is_harmonic(self, a, b, tolerance=0.05):
                    a = np.float(a)
                    b = np.float(b)
                    mod_ab = a % b
                    mod_ba = b % a
                    return (a > b and a < b * 3 + tolerance * 3 and (
                            abs(mod_ab % 1) <= tolerance or abs((b - mod_ab) % 1) <= tolerance)) or \
                           (b > a and a > b / 3 - tolerance / 3 and (
                                   abs(mod_ba % 1) <= tolerance or abs((a - mod_ba) % 1) <= tolerance))

                def isRightEpoch(self, t0, known_epoch, known_period):
                    right_epoch = False
                    for i in range(-5, 5):
                        right_epoch = right_epoch or (np.abs(t0 - known_epoch + i * known_period) < (
                                1. / 24.))
                    return right_epoch
            report = {}
            reports_df = pd.DataFrame(columns=['period', 'radius', 'epoch', 'found', 'snr', 'sde', 'run'])
            a = False
            samples_analysed = sherlock_samples
            for index, row in tls_report_df[::-1].iterrows():
                file = os.path.join(
                    'P' + str(row['period']) + '_R' + str(row['radius']) + '_' + str(row['epoch']) + '.csv')
                first_false = index > 0 and tls_report_df.iloc[index - 1]['found'] and \
                            not tls_report_df.iloc[index]['found']
                if first_false:
                    samples_analysed = 0
                elif tls_report_df.iloc[index]['found']:
                    samples_analysed = sherlock_samples
                if samples_analysed < sherlock_samples:
                    try:
                        samples_analysed = samples_analysed + 1
                        period = float(re.search("P([0-9]+\\.[0-9]+)", file)[1])
                        r_planet = float(re.search("R([0-9]+\\.[0-9]+)", file)[1])
                        epoch = float(re.search("_([0-9]+\\.[0-9]+)\\.csv", file)[1])
                        signal_selection_algorithm = QuorumSnrBorderCorrectedStopWhenMatchSignalSelector(1, 0, period,
                                                                                                         epoch)
                        df = pd.read_csv(inject_dir + file, float_precision='round_trip', sep=',',
                                         usecols=['#time', 'flux', 'flux_err'])
                        if len(df) == 0:
                            found = True
                            snr = 20
                            sde = 20
                            run = 1
                        else:
                            sherlock.Sherlock(False, object_infos=[MissionInputObjectInfo(self.id, inject_dir + file)]) \
                                .setup_detrend(True, True, 1.5, 4, 12, "biweight", 0.2, 1.0, 20, False,
                                               0.25, "cosine", None) \
                                .setup_transit_adjust_params(5, None, None, 10, None, None, 0.4, 14, 10,
                                                             20, 5, 5.5, 0.05, "mask", "quorum", 1, 0,
                                                             signal_selection_algorithm) \
                                .run()
                            df = pd.read_csv(self.id.replace(" ", "") + "_INP/candidates.csv", float_precision='round_trip', sep=',',
                                             usecols=['curve', 'period', 't0', 'run', 'snr', 'sde', 'rad_p',
                                                      'transits'])
                            snr = df["snr"].iloc[len(df) - 1]
                            run = df["run"].iloc[len(df) - 1]
                            sde = df["sde"].iloc[len(df) - 1]
                            per_run = 0
                            found_period = False
                            j = 0
                            for per in df["period"]:
                                if signal_selection_algorithm.is_harmonic(per, period / 2.):
                                    found_period = True
                                    t0 = df["t0"].iloc[j]
                                    break
                                j = j + 1
                            right_epoch = False
                            if found_period:
                                for i in range(-5, 5):
                                    right_epoch = right_epoch or (np.abs(t0 - epoch + i * period) < (
                                            1. / 24.))
                                    if right_epoch:
                                        snr = df["snr"].iloc[j]
                                        run = df["run"].iloc[j]
                                        sde = df["sde"].iloc[j]
                                        break
                            found = right_epoch
                        new_report = {"period": period, "radius": r_planet, "epoch": epoch, "found": found, "sde": sde,
                                      "snr": snr,
                                      "run": int(run)}
                        reports_df = reports_df.append(new_report, ignore_index=True)
                        reports_df.to_csv(inject_dir + "a_sherlock_report.csv", index=False)
                        print("P=" + str(period) + ", R=" + str(r_planet) + ", T0=" + str(epoch) + ", FOUND WAS " + str(
                            found) +
                              " WITH SNR " + str(snr) + "and SDE " + str(sde))
                    except Exception as e:
                        print(e)
                        print("File not valid: " + file)


    def __make_model(self, time, flux, flux_err, rstar, mstar, epoch, period, rplanet):
        # a = (7.495e-6 * period**2)**(1./3.)*u.au #in AU
        P1 = period * u.day
        a = np.cbrt((ac.G * mstar * P1 ** 2) / (4 * np.pi ** 2)).to(u.au)
        # print("radius_1 =", rstar.to(u.au) / a) #star radius convert from AU to in units of a
        # print("radius_2 =", rplanet.to(u.au) / a)
        texpo = 2. / 60. / 24.
        # print("T_expo = ", texpo,"dy")
        # tdur=t14(R_s=radius, M_s=mass,P=period,small_planet=False) #we define the typical duration of a small planet in this star
        # print("transit_duration= ", tdur*24*60,"min" )
        model = ellc.lc(
            t_obs=time,
            radius_1=rstar.to(u.au) / a,  # star radius convert from AU to in units of a
            radius_2=rplanet.to(u.au) / a,  # convert from Rearth (equatorial) into AU and then into units of a
            sbratio=0,
            incl=90,
            light_3=0,
            t_zero=epoch,
            period=period,
            a=None,
            q=1e-6,
            f_c=None, f_s=None,
            ldc_1=self.ab, ldc_2=None,
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
        flux_t = flux + model - 1.
        if model[0] > 0:
            result_flux = flux_t
            result_flux_err = flux_err
            result_time = time
        else:
            result_flux = []
            result_time = []
            result_flux_err = []
        return result_time, result_flux, result_flux_err

    def transit_masks(self, transit_masks, time):
        if transit_masks is None:
            transit_masks = []
        result = np.full(len(time), False)
        for mask in transit_masks:
            intransit = transit_mask(time, mask["P"], 2 * mask["D"], mask["T0"])
            result[intransit] = True
        return result

    def plot_results(self, step_period, step_radius, phases, plot_step_period=1, plot_step_radius=1):
        inject_dir = self.build_inject_dir()
        df = pd.read_csv(inject_dir + "a_tls_report.csv")
        min_period = df["period"].min()
        max_period = df["period"].max()
        min_rad = df["radius"].min()
        max_rad = df["radius"].max()
        if step_period <= 0:
            step_period = 0.1
        if step_radius <= 0:
            step_radius = 0.1
        period_grid = np.arange(min_period, max_period + 0.01, step_period)
        radius_grid = np.arange(min_rad, max_rad + 0.01, step_radius)
        result = np.zeros((len(period_grid), len(radius_grid)))
        for i in period_grid:
            ipos = int(round((i - min_period) * 1 / step_period))
            for j in radius_grid:
                jpos = int(round((j - min_rad) * 1 / step_radius))
                sel_df = df[self.__equal(df["period"], i)]
                sel_df = sel_df[self.__equal(sel_df["radius"], j)]
                found_count = len(sel_df[sel_df["found"]])
                result[ipos][jpos] = found_count / phases * 100
        result = np.transpose(result)
        fig, ax = plt.subplots()
        im = ax.imshow(result)
        ax.set_xticks(np.arange(len(period_grid), step=plot_step_period))
        ax.set_yticks(np.arange(len(radius_grid), step=plot_step_radius))
        ax.set_xticklabels(period_grid[0::plot_step_period])
        ax.set_yticklabels(radius_grid[0::plot_step_radius])
        ax.set_xlabel("Period")
        ax.set_ylabel("Radius")
        plt.setp(ax.get_xticklabels(), rotation=30, ha="right",
                 rotation_mode="anchor")
        cbar = ax.figure.colorbar(im, ax=ax, shrink=len(period_grid) / len(radius_grid))
        cbar.ax.set_ylabel("% of found transits", rotation=-90, va="bottom")
        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
                 rotation_mode="anchor")
        plt.gca().invert_yaxis()
        ax.set_title(self.id + " - TLS P/R recovery")
        fig.tight_layout()
        plt.savefig(inject_dir + "a_tls_report.png")
        plt.close()

    def __tls_search(self, time, flux, rstar, rstar_min, rstar_max, mass, mstar_min, mstar_max, ab, intransit, epoch,
                     period, min_period, max_period, min_snr, cores, transit_template, ws, transits_min_count):
        SNR = 1e12
        FOUND_SIGNAL = False
        time = time[~intransit]
        flux = flux[~intransit]
        time, flux = cleaned_array(time, flux)
        run = 0
        if ws > 0:
            flux = wotan.flatten(time, flux, window_length=ws, return_trend=False, method='biweight', break_tolerance=0.5)
        #::: search for the rest
        while (SNR >= min_snr) and (not FOUND_SIGNAL):
            model = transitleastsquares(time, flux)
            # R_starx = rstar / u.R_sun
            results = model.power(u=ab,
                                  R_star=rstar,  # rstar/u.R_sun,
                                  R_star_min=rstar_min,  # rstar_min/u.R_sun,
                                  R_star_max=rstar_max,  # rstar_max/u.R_sun,
                                  M_star=mass,  # mstar/u.M_sun,
                                  M_star_min=mstar_min,  # mstar_min/u.M_sun,
                                  M_star_max=mstar_max,  # mstar_max/u.M_sun,
                                  period_min=min_period,
                                  period_max=max_period,
                                  n_transits_min=transits_min_count,
                                  show_progress_bar=False,
                                  use_threads=cores,
                                  transit_template=transit_template
                                  )
            SNR = results.snr
            if results.snr >= min_snr:
                intransit_result = transit_mask(time, results.period, 2 * results.duration, results.T0)
                time = time[~intransit_result]
                flux = flux[~intransit_result]
                time, flux = cleaned_array(time, flux)
                right_period = self.__is_multiple_of(results.period, period / 2.)
                right_epoch = False
                for tt in results.transit_times:
                    for i in range(-5, 5):
                        right_epoch = right_epoch or (np.abs(tt - epoch + i * period) < (
                                1. / 24.))  # check if any epochs matches to within 1 hour
                #            right_depth   = (np.abs(np.sqrt(1.-results.depth)*rstar - rplanet)/rplanet < 0.05) #check if the depth matches
                if right_period and right_epoch:
                    FOUND_SIGNAL = True
                    break
            run = run + 1
        return FOUND_SIGNAL, results.snr, results.SDE, run

    def __equal(self, a, b, tolerance=0.01):
        return np.abs(a - b) < tolerance

    def __is_multiple_of(self, a, b, tolerance=0.05):
        a = np.float(a)
        b = np.float(b)
        mod_ab = a % b
        mod_ba = b % a
        return (a > b and a < b * 3 + tolerance * 3 and (
                    abs(mod_ab % 1) <= tolerance or abs((b - mod_ab) % 1) <= tolerance)) or \
               (b > a and a > b / 3 - tolerance / 3 and (
                           abs(mod_ba % 1) <= tolerance or abs((a - mod_ba) % 1) <= tolerance))