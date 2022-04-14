import numpy as np
import pandas as pd
import photospline as psp
import matplotlib.pyplot as plt
import glob
import os

# class that can be used for flux calculation
class SeyfertFlux():
    def __init__(self, path_to_splinetables, log_xray_lumin):
        '''
        Args
            path_to_splinetables (str): directory that contains the photospline .fits files
            log_xray_lumin (float): log10(intrinsic luminosity 2-10 keV),
                                    i.e. value from BASS column "logL2-10-intr"
            distance_in_mpc (float): distance to source in Mpc
                                     default = 1 Mpc
        '''

        self.log_xray_lumin = log_xray_lumin
        self.path_to_splinetables = path_to_splinetables

        # defines the "safe" range of the energy flux spline
        self._crit_log_energy_flux = -50
        self._crit_log_nu_energy_upper = 7.0
        self._crit_log_nu_energy_lower = 2.0

        loaded_log_lumin = self._load_spline()

        # generate relative flux scale based on requested xray logL
        # and loaded xray_logL
        self._relative_flux_scale = 10**(log_xray_lumin - loaded_log_lumin)
        print(f"will correct neutrino flux by luminosity ratio: {self._relative_flux_scale:.2f}")

    def _load_spline(self):
        fs = glob.glob(os.path.join(self.path_to_splinetables+'*.fits'))
        lls = [float(f.split("_")[-1][:4]) for f in fs]
        dls = [np.abs(ll-self.log_xray_lumin) for ll in lls]

        # find splinetable that matches best the request luminosity
        select = np.argmin(dls)
        f = fs[select]
        logL_select = lls[select]
        self._loaded_lumin = logL_select

        print('selected spline:', f, 'for log-xray-luminosity:', self.log_xray_lumin)
        self.spline = psp.SplineTable(f)
        return float(lls[select])

    def get_loaded_lumin(self):
        return self._loaded_lumin

    def get_splinetable_lumin_scale(self):
        return self._relative_flux_scale

    def get_splinetable(self):
        return self.spline

    def get_flux(self, enu, distance_in_mpc=1.):
        '''
        Computes neutrino flux expected for given neutrino energy.
        Args
            enu (float): neutrino energy in GeV
            distance_in_mpc (float): distance to source in Mpc
                                                 default = 1 Mpc
        Returns
            neutrino flux (float) in units of 1/GeV 1/cm^2 1/s
        '''

        log_enu = np.log10(enu)
        log_energy_flux = self.spline.evaluate_simple([log_enu])

        # convert energy flux to particle flux. account for source distance.
        neutrino_flux = 10**(log_energy_flux - 2.0*log_enu - 2.0*np.log10(distance_in_mpc))

        # need take care of very small fluxes (set to 0 beyond critical energy)
        # or below critical flux
        if isinstance(log_energy_flux, np.ndarray):
            out_of_bounds1 = log_energy_flux < self._crit_log_energy_flux
            out_of_bounds2 = np.logical_or(log_enu < self._crit_log_nu_energy_lower,
                                           log_enu > self._crit_log_nu_energy_upper)
            neutrino_flux[np.logical_or(out_of_bounds1, out_of_bounds2)] = 0

        else:
            if (log_energy_flux < self._crit_log_energy_flux) or \
                (log_enu < self._crit_log_nu_energy_lower) or \
                (log_enu > self._crit_log_nu_energy_upper):
                return 0

        return neutrino_flux * self._relative_flux_scale


# function that computes event expectations
def get_nevents_from_source(mc, src_dec, neutrino_flux, ltime=8.7, dec_bw=0.5):
    '''
    Computes expected number of events for a source with known flux from MC file.
    Args
        mc (np.ndarray): a MC file that holds MC events with fields following  NuSources convention.
        src_dec (float): declination of the source location (in deg)
        neutrino flux (python function): a function that neutrino energy (in GeV) as argument
                                         and returns corresponding neutrino flux in 1/GeV 1/cm^2 1/s
        etime (float): live time in units of year
                        default = 8.7yr corresponding to Northern Tracks v005p00
        dec_bw (float): specifies declination band (src_dec-dec_bw, src_dec+dec_bw) that is used to
                        to compute event expectation (in deg)
                        default = 0.5deg

    Returns
        expected number of events (float)

    '''

    # convert ltime from year to seconds
    ltime *= 365.25 * 24. * 3600.

    # convert args from deg to radians
    src_dec = np.radians(src_dec)
    bw = np.radians(dec_bw)

    # define sin_dec range
    sin_dec_edges = [np.sin(src_dec - bw), np.sin(src_dec + bw)]
    delta_sin_dec = sin_dec_edges[1] - sin_dec_edges[0]

    # mask matching dec band
    #mc_true_sin_dec = np.sin(mc['trueDec'])
    #mask = np.logical_and(mc_true_sin_dec < sin_dec_edges[1], mc_true_sin_dec > sin_dec_edges[0])
    mc_true_dec = mc['trueDec']
    mask = np.logical_and(mc_true_dec > src_dec - bw, mc_true_dec < src_dec + bw)
    solid_angle = delta_sin_dec * 2 * np.pi
    smc = mc[mask]

    return np.sum(smc['ow']*neutrino_flux(smc['trueE'])) * ltime/solid_angle









