"""
Created on Wed Jul 14 14:26:35 2021

@author: jroquette
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib.colors import LinearSegmentedColormap
from astropy.table import Table
from scipy.interpolate import interp1d
from astropy import units as u

 # define the color-table and a normalization used for FUV
 cmapFUV = LinearSegmentedColormap.from_list('mycmap', [
     (0.0, 'mediumpurple'),
     (0.2, 'rebeccapurple'),
     (0.38, 'forestgreen'),
     (0.58, 'gold'),
     (0.75, 'orange'),
     (0.8, 'tab:red'),
     (1.0, 'maroon')])
 normFUV = mpl.colors.Normalize(vmin=np.log10(1.7), vmax=5)

def plotCmapFUV(**kargs):
    ax = None
    fig = None
    fontsize = None
    if 'fig' in kargs.keys():
        fig = kargs['fig']
    if 'ax' not in kargs.keys():
        fig, ax = plt.subplots(figsize=(0.5, 3), dpi=100)
        plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.97)
    else:
        ax = kargs['ax']
    if 'fontsize' in kargs.keys():
        fontsize = kargs['fontsize']
    else:
        fontsize = 14
    cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmapFUV, norm=normFUV,
                                    orientation='vertical')
    cb1.ax.set_ylabel(r'log(F$_\mathrm{FUV}$) - (log(G$_0$)', rotation=270,
                      labelpad=20, fontsize=fontsize)

class DiskWithFUV:
    """
    Interpolates timescales for disk dissipation as a function of the local
    FUV flux, base on the results of the disk-dissipation model by Winter et
    al. (2020) MNRAS 491 90
    https://ui.adsabs.harvard.edu/abs/2020MNRAS.491..903W/abstract

    WKC20 grid provides:
        - tauvisc_Myr
        - Mstar_Msol
        - FUV_G0
        - Mdisc_Msol
        - Rscale_au
        - alphavisc
        - Tdisp_Myr

    -----
    Version: 14th July 2021



                - input -
    mass: in Msun
    FUV:  in Go
    tau_vis=1. is the recommended value by Andrew Winter
    curves: default=False
    ---
    - output-
    if curves=False, returns tau_d which is the timescale for
        disk-dissipation by external photoevaporation under FUV
    if curves=True, returns a curve of FUV,tau_d for the mass given

    """
    def __init__(self):
        """
        Load tabulations for a fixed viscous timescales as a grid

        ---
        Usage:

            disk_model=DiskWithFUV()

        """
        import os
        datadir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'data/disk_model/WKC20/')
        self.WKC20 = Table.read(datadir +
                                'FUVdisp_R40_tauvisc.csv', format='csv')
        #masses for  which disk-models are calculated
        self.M = np.array([0.1, 0.3, 0.5, 0.8, 1., 1.3, 1.6, 1.9])
        self.FUV = np.unique(np.array(self.WKC20['FUV_G0']))
    def get_tauD_vs_FUV(self, mass, tau_vis=1.):
        """
        Select within the grid the tabulations for a mass and and
        viscous timescale.

        Note that WKC20 models are estimated such as under no FUV influence,
        disks are dissipated in 10 Myr due to internal processes. Hence
        I am appending tau_D=10 Myr at FUV=0 to these arrays.

        ____
        input:
            mass = mass in Msun
            tau_vis [1.0,2.0,5.0] = viscous timescale in Myrs

        """
        self.mass = mass #needs to keep track to know where I am in the grid
        select_tau_vis = self.WKC20['# tauvisc_Myr'] == tau_vis
        if mass in self.M: #just select table for the M and interpolate values
            select_mass = self.WKC20['Mstar_Msol'] == mass
            self.tauD_Myr = 1e6*np.array(self.WKC20['Tdisp_Myr'][select_mass \
                                                    & select_tau_vis])
        elif (mass > self.M.max()) or (mass < self.M.min()):
            print('outside mass limit')
            self.tauD_Myr = np.full(len(self.FUV), np.nan)
        else: # interpolate a new table
            j = np.argsort(abs(self.M - mass))[0]
            if self.M[j] > mass:
                x = [self.M[j - 1], self.M[j]]
            else:
                x = [self.M[j], self.M[j + 1]]
            self.tauD_Myr = \
               np.array([1e6*np.interp(mass, x, y) for y in list(zip(\
                               self.WKC20['Tdisp_Myr'][select_tau_vis &
                               (self.WKC20['Mstar_Msol'] == x[0])],
                               self.WKC20['Tdisp_Myr'][select_tau_vis &
                               (self.WKC20['Mstar_Msol'] ==  x[1])]))])

    def get_tauD(self, mass, FUV, tau_vis=1.):
        self.get_tauD_vs_FUV(mass, tau_vis=tau_vis)
        tauD_ = np.insert(self.tauD_Myr, 0, 10e6, axis=0)
        FUV_ = np.insert(self.FUV, 0, 1, axis=0)
        if FUV in self.FUV:
            return self.tauD_Myr[self.FUV == FUV][0]
        elif FUV > 1e4:
            return self.tauD_Myr[self.FUV == 1e4][0]
        elif FUV > 1:
            #note that
            tauD = 10**interp1d(np.log10(FUV_), np.log10(tauD_),
                                kind='linear', bounds_error=False,
                                fill_value='extrapolate')(np.log10(FUV))
        elif FUV > 0:
            tauD = 10e6
        elif FUV < 0:
            print("I can't be negative!")
            tauD = np.nan
        else:
            tauD = 10e6
        return tauD


class Parravano:
    """
    Based on Parravano et al. (2003) ApJ 584 797 [1], this class provides
    the tools to estimate the amount of far-ultraviolet radiation outputed
    by massive stars. The scalling relations used here are defined in their
    Table 1 and are available for stars in the mass range 1.8-120Msun

    [1] https://ui.adsabs.harvard.edu/abs/2003ApJ...584..797P/abstract
    """
    def __init__(self, mass):
        """
        """
        self.M = float(mass)
        self.G0 = u.def_unit(r'G_0', 1.6e-3 << u.erg/(u.cm**2)/u.s )
    def LEUV(self):
        """
        Gives extreme-ultraviolet (h\nu > 13.6eV) luminosities as a function
        of stellar mass.

        ----

        L_EUV has units photons/s
        """
        if self.M >= 5.:
            if self.M < 7.:
                return 2.23*1e34*(self.M**11.5) << u.ph/u.s
            elif self.M < 12:
                return 3.69*1e36*(self.M**8.87) << u.ph/u.s
            elif self.M < 20:
                return 4.8*1e38*(self.M**7.85) << u.ph/u.s
            elif self.M < 30:
                return 3.12*1e41*(self.M**4.91) << u.ph/u.s
            elif self.M < 40:
                return 2.8*1e44*(self.M**2.91) << u.ph/u.s
            elif self.M < 60:
                return 3.49*1e45*(self.M**2.23) << u.ph/u.s
            elif self.M < 120:
                return 2.39*1e46*(self.M**1.76) << u.ph/u.s
        else:
            return np.nan

    def LFUV_LH2(self):
        """
        Parravano et al. (2003) gives scalling relations for L_FUV-L_H2,
        which is basically the FUV band excluding the narrow "H_2 band".

        ----

        L_FUV-L_H2 has units of Lsun.
        """
        if self.M >= 1.8:
            if self.M < 2.:
                return 2.77*1e-4*(self.M**11.8) << u.L_sun
            elif self.M < 2.5:
                return 1.88*1e-3*(self.M**9.03) << u.L_sun
            elif self.M < 3.:
                return 1.19*1e-2*(self.M**7.03) << u.L_sun
            elif self.M < 6.:
                return 1.47*1e-1*(self.M**4.76) << u.L_sun
            elif self.M < 9.:
                return 8.22*1e-1*(self.M**3.78) << u.L_sun
            elif self.M < 12.:
                return 2.29*(self.M**3.31) << u.L_sun
            elif self.M < 30.:
                return 2.7*1e1*(self.M**2.32) << u.L_sun
            elif self.M < 120.:
                return 3.99*1e2*(self.M**1.54) << u.L_sun
        else:
            return np.nan

    def LH2(self):
        """
        Gives the luminosity in the narrow "H_2 band"
        This covers about 10% of the FUV band in the wavelength range
        912-1100\AA.

        ----

        L_H2 has units of Lsun
        """
        if self.M >= 1.8:
            if self.M < 3:
                return 1.98*1e-14*(self.M**26.6) << u.L_sun
            elif self.M < 4.:
                return 2.86*1e-8*(self.M**13.7) << u.L_sun
            elif self.M < 6.:
                return 1.35*1e-4*(self.M**7.61) << u.L_sun
            elif self.M < 9.:
                return 1.1*1e-2*(self.M**5.13) << u.L_sun
            elif self.M < 12.:
                return 1.07*1e-1*(self.M**4.09) << u.L_sun
            elif self.M < 15.:
                return 5.47*1e-1*(self.M**3.43) << u.L_sun
            elif self.M < 30.:
                return 9.07*(self.M**2.39) << u.L_sun
            elif self.M < 120.:
                return 9.91*1e1*(self.M**1.69) << u.L_sun
        else:
            return np.nan

    def LFUV(self):
        """
        Gives the luminosity in the whole FUV band, which covers
        6eV<=h\nu<=13.4eV, or 919-2070\AA

        ----

        L_FUV is given in units erg/s (cgs)
        """
        return (self.LFUV_LH2() + self.LH2()).cgs

    def localFlux(self,L, d):
        """
        Given any luminosity, return the flux at a distance d, in parsec
        """
        if not isinstance(d, type(u.pc)):
            d = d << u.pc
        buffer = 1e-10*u.pc
        return L/4./np.pi/(d-buffer)**2

    def localFUV(self, d):
        """
        Return the Local FUV flux at a distance d (in parsec), in G0 units
        """
        return self.localFlux(self.LFUV(), d).to(self.G0)

    def localEUV(self, d):
        """
        Return the Local EUV flux at a distance d (in parsec), in
        photons/s/cm2
        """
        return self.localFlux(self.LEUV(),
                              d).decompose(bases=[u.cm, u.ph, u.s])
