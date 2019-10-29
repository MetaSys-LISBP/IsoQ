import pandas as pd
import isocor as ic
import isoq
import isocor.ui.isocordb
from pathlib import Path
from os.path import expanduser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pkg_resources
import logging


class run():
    def __init__(self, folder, datafile, calibfile, tracer, resolution, mz_of_resolution, tracer_purity, correct_NA_tracer, resolution_formula_code, purity15N=None, exp_CID_IS=None, verbose = False):

        # create logger
        self.logger = logging.getLogger()
        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s', "%Y-%m-%d %H:%M:%S")

        # send logs to filehandler (to send logs to sys.stderr, use 'strm_hdlr = logging.StreamHandler()')
        strm_hdlr = logging.FileHandler(str(Path(folder, "results.log")), mode='w')
        strm_hdlr.setFormatter(formatter)
        self.logger.addHandler(strm_hdlr)

        # set log level (also catch isocor logs if verbose)
        if verbose:
            self.logger.setLevel(logging.DEBUG)
            self.logger.debug("Verbose mode.")
        else:
            self.logger.setLevel(logging.INFO)

        # log general information on the process
        self.logger.info('------------------------------------------------')
        self.logger.info('GENERAL INFORMATION')
        self.logger.info('------------------------------------------------')
        self.logger.info('IsoQ version: {}'.format(isoq.__version__))
        self.logger.info('15N-purity of 13C15N-internal standard: {}'.format(purity15N))
        self.purity15N = purity15N
        self.logger.info('exp. CID of 13C15N-internal standard: {}'.format(exp_CID_IS))
        self.exp_CID_IS = exp_CID_IS

        # maximal relative error allowed when constructing calibration curves (samples above this threeshold are removed)
        self.max_error = 0.2
        # minimum number of calibration points to keep
        self.min_cal_points = 5
        
        # create isocor environment
        self.display('create isocor environment...')
        baseenv = isocor.ui.isocordb.EnvComputing()
        home = Path(expanduser('~'), 'isocordb')
        baseenv.registerIsopotes(isotopesfile = Path(home, 'Isotopes.dat'))
        baseenv.registerDerivativesDB(derivativesfile = Path(home, 'Derivatives.dat'))
        baseenv.registerMetabolitesDB(metabolitesfile = Path(home, 'Metabolites.dat'))

        # load measurements
        datafile = Path(folder, datafile)
        if not datafile.is_file():
            raise ValueError(
                "Measurements file not found in:\n'{}'.".format(datafile))
        try:
            with open(str(datafile), 'r', encoding='utf-8') as fp:
                self.dfDatafile = pd.read_csv(fp, delimiter='\t', keep_default_na=False)
            self.dfDatafile = self.dfDatafile.replace("N/F", 0)
        except Exception as err:
            raise ValueError("An unknown error has occurred opening the measurements file ('{}').\n\nPlease check this file (details on the expected format can be found in the documentation), correct the issue and rerun IsoQ.\n\nTraceback for debugging:\n{}".format(datafile, err))
        self.logger.info('measurements file: {}'.format(datafile))

        # load concentrations
        calibfile = Path(folder, calibfile)
        if not calibfile.is_file():
            raise ValueError(
                "Calibration file not found in:\n'{}'.".format(calibfile))
        try:
            with open(str(calibfile), 'r', encoding='utf-8') as fp:
                self.dfCalibfile = pd.read_csv(fp, delimiter='\t', keep_default_na=False, na_values = "N/F")
            #self.dfCalibfile = self.dfCalibfile.replace("N/F", 0)
        except Exception as err:
            raise ValueError("An unknown error has occurred opening the calibration file ('{}').\n\nPlease check this file (details on the expected format can be found in the documentation), correct the issue and rerun IsoQ.\n\nTraceback for debugging:\n{}".format(calibfile, err))
        self.logger.info('calibration file: {}'.format(calibfile))
        
        # identify samples
        self.display('parse data files...')
        self.data = self.getSamplesList()
        self.samples, self.calib = self.filterSamplesList()
        self.samples.sort()

        # run calibration
        self.display('run calibration...')
        cal_data = self.parseCalib()
        self.pp = PdfPages(str(Path(folder, 'results.pdf')))
        for metabolite in cal_data.keys():
            cal_data[metabolite] = self.runCalib(cal_data, metabolite, 1, [])
        self.pp.close()

        # process measurements
        self.display('process samples...')
        df_iso, df_met = pd.DataFrame(), {'type':{}, 'comment':{}}
        for sample in self.samples:
            signals = self.getSignalsList(sample)
            metabolites = list(self.getMetabolitesList(signals))
            metabolites.sort()
            df_met[sample] = {}
            for metabolite in metabolites:
                stp = self.getSignals(metabolite, sample)
                isoclust = self.getIsoclust(stp, metabolite)
                isostandard = self.getIS(stp, metabolite)
                self.logger.info("*********")
                self.display("  sample, metabolite, mode: ({}, {}, {})".format(sample, metabolite, cal_data[metabolite]['mode']))
                self.logger.info("mass fractions (area): {}".format(isoclust))
                self.logger.info("13C15N-internal standard (area): {}".format(isostandard))
                try:
                    self.logger.debug("create IsoCor corrector")
                    corrector = ic.MetaboliteCorrectorFactory(
                                formula=baseenv.getMetaboliteFormula(metabolite), tracer=tracer, resolution=resolution, label=metabolite,
                                data_isotopes=baseenv.dictIsotopes, mz_of_resolution=mz_of_resolution,
                                derivative_formula='', tracer_purity=tracer_purity,
                                correct_NA_tracer=correct_NA_tracer, resolution_formula_code=resolution_formula_code,
                                charge=baseenv.getMetaboliteCharge(metabolite))
                    self.logger.debug("remove contribution of IS to CID")
                    if len(isostandard):
                        isoclustCorIS = self.removeIScontribution(isoclust, isostandard, corrector)
                    else:
                        isoclustCorIS = [isoclust, False]
                    self.logger.debug("correct for naturally occuring isotopes")
                    _, CID, _, enr = corrector.correct(isoclustCorIS[0])
                    tmp = isostandard
                    if len(tmp) == 1 and tmp:
                        CID_IS_ratio = sum(isoclust)/isostandard[0]
                        self.logger.info("total CID area / IS area: {}".format(CID_IS_ratio))
                        pool = self.getConc(cal_data[metabolite]['sim_fun'], CID_IS_ratio, cal_data[metabolite]['xlim'])
                        df_met['type'][metabolite] = "IS_ratio"
                        df_met['comment'][metabolite] = ""
                        self.logger.info("calc. concentration: {}".format(pool))
                    else:
                        CID_area = sum(isoclust)
                        self.logger.info("total CID area (no IS found): {}".format(CID_area))
                        if cal_data[metabolite]['mode'] == '12C area':
                            pool = self.getConc(cal_data[metabolite]['sim_fun'], CID_area, cal_data[metabolite]['xlim'])
                        else:
                            pool = CID_area
                        df_met['type'][metabolite] = "area"
                        df_met['comment'][metabolite] = ""
                        self.logger.info("pool (using calibration with total CID area, no IS found): {}".format(CID_area))
                except Exception as err:
                    self.logger.info("Error: {}".format(err))
                    df_met['comment'][metabolite] = err
                    pool = np.nan
                    isoclustCorIS = ([np.nan], None)

                # gather results
                df_met[sample][metabolite] = pool
                self.logger.debug(isoclust)
                self.logger.debug(CID)
                self.logger.debug(isoclustCorIS)
                for i, line in enumerate(zip(*(isoclust, isoclustCorIS[0], CID, [enr]*len(CID)))):
                    df_iso = pd.concat((df_iso, pd.DataFrame([line], index=pd.MultiIndex.from_tuples([[sample, metabolite, i]], names=[
                        'sample', 'metabolite', 'isotopologue']), columns=['area', 'area_corrected_IS', 'isotopologue_fraction', 'mean_enrichment'])))

        # save results
        self.display('save results...')
        df_iso.to_csv(str(Path(folder, 'results_CID.csv')), sep='\t')
        pd.DataFrame.from_dict(df_met).to_csv(str(Path(folder, 'results_conc.csv')), sep='\t')
        self.display('done')

    def getConc(self, fun, y, xlim):
        if not np.isfinite(y):
            return 'Error: exp. value ({}) must be a number.'.format(y)
        intconc = (fun - y).roots
        self.logger.info('polynomial roots: {}'.format(intconc))
        conc = [intconc[i] for i in range(len(intconc)) if xlim[0] <= intconc[i] <= xlim[1]]
        if len(conc) == 1:
            return conc[0]
        elif len(conc) == 0:
            return 'Error: out of calibration range ({}).'.format(xlim)
        else:
            return 'Error: invalid calibration curve, several values are returned ({}).'.format(conc)

    def display(self, m):
        """
        Send message to logs and console
        """
        self.logger.info(m)
        print(m)
        
    def getSamplesList(self):
        return [i[0] for i in self.dfDatafile[['Filename']].drop_duplicates().values]

    def filterSamplesList(self):
        s = [i for i in self.data if 'calib' not in i]
        c = [i for i in self.data if 'calib' in i]
        return s, c

    def getSignalsList(self, sample):
        return [i for i in self.dfDatafile[self.dfDatafile['Filename'] == sample]['Compound'].drop_duplicates().values]

    def getMetabolitesList(self, signals):
        return set([i.split(" ")[0] for i in signals])

    def getSignals(self, metabolite, sample):
        tmp = self.dfDatafile[self.dfDatafile['Filename'] == sample]
        return tmp[tmp['Compound'].str.match(metabolite)]

    def getIsoclust(self, dfSampleMet, metabolite):
        tmp = dfSampleMet[dfSampleMet['Compound'].str.match(metabolite + " M")]
        tmp = tmp.sort_values(by=['Compound'])
        return [float(i) for i in tmp['Area'].values]

    def getIS(self, dfSampleMet, metabolite):
        tmp = dfSampleMet[dfSampleMet['Compound'].str.match(metabolite + " 13C15N SI")]
        return [float(i) for i in tmp['Area'].values]

    def removeIScontribution(self, isoclust, isostandard, corrector):
        warning = False
        nN = corrector.formula['N']
        if not nN:
            raise ValueError('This metabolite does not contain N atoms.')
        nC = corrector.formula['C']
        if self.exp_CID_IS:
            try:
                corrected = isoclust - self.exp_CID_IS*isostandard
            except:
                raise ValueError('Experimental CID of IS does not comply with requirements (e.g. same length as mass fractions vector).')
        elif self.purity15N:
            contrib = [1.]
            for _ in range(nN):
                contrib = np.convolve(contrib, self.purity15N)
            contrib = np.concatenate([[0.] * nC, contrib*isostandard])
            corrected = isoclust - contrib[0:len(isoclust)]
        else:
            raise ValueError('Not enough information provided to correct CID for the contribution of IS.')
        if any(corrected < 0):
            self.logger.info('Warning: some fractions are negative after correction for IS contribution ({}), they are flatten to 0.'.format(corrected))
            corrected[corrected < 0] = 0
        return corrected, warning
        
    def runCalib(self, cal_data, metabolite, npass, iexcluded):
        res = {}
        datam = pd.DataFrame(cal_data[metabolite])
        if datam['CID_area'].isnull().values.all() or datam['concentration'].isnull().values.all() or datam['CID_area'].isnull().values.all() or datam['concentration'].isnull().values.all():
            res['mode'] = 'no calibration'
            res['x'] = np.nan
            res['xlim'] = np.nan
            res['y'] = np.nan
            res['r2'] = np.nan
            res['coeffs'] = np.nan
            res['sim_fun'] = lambda x: x
            res['relative_residuals'] = np.nan
        else:
            if datam['IS_area'].isnull().values.all():
                y = datam['CID_area']
                res['mode'] = '12C area'
            else:
                y = datam['CID_area']/datam['IS_area']
                res['mode'] = '12C/13C ratio'
            datap = pd.concat([datam['concentration'], y], axis=1, ignore_index=True).dropna()
            datap = datap.drop(index=iexcluded)
            res['x'] = np.array(datap[0].tolist())
            res['y'] = np.array(datap[1].tolist())
            res['xlim'] = [min(res['x']), max(res['x'])]
            # if you prefer to use the polyfit function (np doc suggests using Polynomial.fit for numerical stability):
            #fres = np.polyfit(res['x'][idx], res['y'][idx], 2, w=1/res['x'][idx])
            # Note: 'fres' object must replace 'coe' and must be passed directly to poly1d (coefficients must not be reversed!!)
            fres = np.polynomial.polynomial.Polynomial.fit(res['x'], res['y'], 2, w=1.0/res['x'])
            r2 = round(np.corrcoef(res['x'], res['y'])[0,1]**2, 3)
            res['r2'] = r2
            coe = list(fres.convert().coef)
            coe.reverse()
            res['coeffs'] = coe
            res['sim_fun'] = np.poly1d(coe)
            res['relative_residuals'] = (res['sim_fun'](res['x'])-res['y'])/res['y']
        self.logger.info('{} - pass {}'.format(metabolite, npass))
        self.logger.info('mode: {}'.format(res['mode']))
        self.logger.info('concentrations: {}'.format(res['x']))
        self.logger.info('{}: {}'.format(res['mode'], res['y']))
        self.logger.info('calibration range: {}'.format(res['xlim']))
        self.logger.info('polynomial coefficients: {}'.format(res['coeffs']))
        self.logger.info('R2: {}'.format(res['r2']))
        self.logger.info('best fit: {}'.format(res['sim_fun'](res['x'])))
        self.logger.info('relative errors: {}'.format(res['relative_residuals']))
        self.plotCalib(res, metabolite, npass, iexcluded)
        # run another calibration pass excluding the calibration sample which show the highest error, if > 20%
        if isinstance(res['relative_residuals'], np.ndarray):
            tmp = list(abs(res['relative_residuals']))
            excluded = tmp.index(max(tmp)) if max(tmp) > self.max_error else None
            if excluded is not None and len(res['relative_residuals']) > self.min_cal_points:
                iexcluded.append(list(datap.index)[excluded])
                res = self.runCalib(cal_data, metabolite, npass+1, iexcluded)
        return res

    def plotCalib(self, res, metabolite, npass, iexcluded):
        plt.figure(1)
        if 'r2' in res.keys():
            if np.isnan(res['r2']):
                stitle = '{} - no calibration data'.format(metabolite)
            else:
                stitle = '{} (pass {}, R² = {}'.format(metabolite, npass, res['r2'])
                if len(iexcluded):
                    stitle += ', excluding cal. sample #{}'.format(iexcluded[-1]+1)
                stitle += ')'
            plt.suptitle(stitle)
            if res['mode'] != 'no calibration':
                plt.subplot(211)
                xp = np.linspace(0, max(res['x']), 100)
                _ = plt.plot(res['x'], res['y'], '.', xp, res['sim_fun'](xp), '-')
                plt.ylim(-max(res['y'])*0.05, max(res['y'])*1.05)
                plt.xlim(-max(res['x'])*0.05, max(res['x'])*1.05)
                yleg = max(res['y'])/2
                xleg = -max(res['x'])*0.18
                plt.text(xleg, yleg, res['mode'], horizontalalignment='center', verticalalignment='center', rotation=90)
                plt.grid(True)
                
                plt.subplot(212)
                _ = plt.plot(res['x'], res['relative_residuals'], '.')
                plt.ylim(min(-self.max_error*1.25, min(res['relative_residuals'])*1.1), max(self.max_error*1.25, max(res['relative_residuals'])*1.1))
                plt.xlim(-max(res['x'])*0.05, max(res['x'])*1.05)
                yleg = min(-self.max_error*1.25, min(res['relative_residuals'])*1.1) - (max(self.max_error*1.25, max(res['relative_residuals'])*1.1) - min(-self.max_error*1.25, min(res['relative_residuals'])*1.1))*0.2
                xleg = max(res['x'])/2
                plt.text(xleg, yleg, 'concentration', horizontalalignment='center', verticalalignment='center')
                yleg = max(self.max_error*1.25, max(res['relative_residuals'])*1.1) - (max(self.max_error*1.25, max(res['relative_residuals'])*1.1) - min(-self.max_error*1.25, min(res['relative_residuals'])*1.1))/2
                xleg = -max(res['x'])*0.18
                plt.text(xleg, yleg, 'relative error', horizontalalignment='center', verticalalignment='center', rotation=90)
                plt.grid(True)
                plt.axhline(y=-self.max_error, color='r', linestyle='-')
                plt.axhline(y=self.max_error, color='r', linestyle='-')
            plt.savefig(self.pp, format='pdf')
            plt.close()
        else:
            self.logger.error("Cannot construct calibration curve for metabolite '{}'".format(metabolite))

    def parseCalib(self):
        cal_data = {}
        for sample in self.calib:
            signals = self.getSignalsList(sample)
            metabolites = self.getMetabolitesList(signals)
            for metabolite in metabolites:
                stp = self.getSignals(metabolite, sample)
                isoclust = self.getIsoclust(stp, metabolite)
                isostandard = self.getIS(stp, metabolite)
                cal_data[metabolite] = cal_data.get(metabolite, {'CID_area':[], 'IS_area':[], 'concentration':[]})
                cal_data[metabolite]['CID_area'].append(sum(isoclust))
                if len(isostandard):
                    cal_data[metabolite]['IS_area'].append(isostandard[0])
                else:
                    cal_data[metabolite]['IS_area'].append(np.nan)
                try:
                    scol = sample.split('_')
                    idcol = scol.index('calib')+1
                    calpt = scol[idcol]
                    calconc = self.dfCalibfile[(self.dfCalibfile['Compound'] == metabolite)][calpt]
                    cal_data[metabolite]['concentration'].append(calconc.values[0])
                except:
                    cal_data[metabolite]['concentration'].append(np.nan)
        return cal_data


