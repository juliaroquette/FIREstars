"""
Created on Mon Jul 12 14:17:00 2021

@author: juliaroquette

Contains functions for loading stellar evolution models that will be fed
to the SpinEvolution.py code.

Version of 13th July 2021:

    Only the Baraffe et al (2015, BHAC15) model is currently available as
    BHAC15_MassTrack

"""
import numpy as np
from astropy.table import Table
from astropy import units as u

class BHAC15_MassTrack:
    """
    Loads Baraffe et al. (2015, BHAC15) Mass Tracks

    Tracks are for Solar Metalicity stars.
    source: http://perso.ens-lyon.fr/isabelle.baraffe/BHAC15dir/

    ------------------------------------------------------------------------
    Version date 13th July 2021

    Usage:

        To load the grid:

            grid = BHAC15_MassTrack()

        This will return the object:

            grid.BHAC15

        which contains the whole grid. Parameter available are detailed in [1]
        To get mass tracks for a specific mass (after loading the grid!):

        grid.getMassTrack(mass)

        this will load grid.Params for the basic Params:

            grid.Mass = mass in Msun
            grid.Teff = effective temperature in K
            grid.Age = age in years
            grid.R = Radius in Rsun
            grid.k2rad = radiative gyration radius
            grid.k2conv = convective gyration radius
            grid.k2_2 = gyration radius - see note [2]
            grid.I = Momentum of Inertia (Isun) see note [2]

        If full=True, it will load the remaining parameters available with
        the BHAC15 grid:

                grid.Lum = Luminosity in Lsun
                grid.logg = surface gravity (logg)
                grid.log_Li = Li abundance
                grid.Tcentral = central temperature in K
                grid.Rhocentral = central density in gr/cc
                grid.Mrad = radiative core mass in Msun
                grid.Rrad =  radiative core in Rsun

    ------------------------------------------------------------------------

    [1] Extracted from BHAC15 Notes:
    BHAC15 tracks and internal structure for brown dwarfs and low mass star
    (0.075 Msun to 1.4 Msun)

    M/Ms:  mass of the star in units of solar mass
    log t: age of the star (in yr)
    Teff: effective temperature (in K)
    L/Ls: log luminosity in units of solar luminosity
          (value used Ls=3.839d+33)
    g: log g  (surface gravity)
    R/Rs : radius of the star in units of solar radius
           (value used Rs=6.96d10)
    log(Li/Li0): log of the ratio of surface lithium abundance to initial
                 abundance
    log Tc : log of central temprature
    log ROc: log of central density (in gr/cc)
    Mrad: mass of radiative core (in solar mass)
    Rrad: radius of radiatif core (in solar radius)
    k2conv: convective gyration radius
    k2rad: radiative gyration radius

    ------------------

    [2] NOTE: Adopted convention to calculate the gyration radii:

    k2conv=[ I/(R**2 * Mstar) ]**(1/2)

    with I the moment of inertia defined here by :

    I = 2/3 integral(r**2 * dm)
    with the integral over the mass of the convective region

    Same for k2rad, corresponding to the radiative zone.

    Example: To calculate the moment of inertia with the definition
    above for a 1 Msun star at t = 4.6 Gyr, the table gives

    k2conv= 8.993E-02
    k2rad= 2.519E-01

    k2**2 = k2conv**2 + k2rad**2 ~ 0.071
    (which is the typical value found for the Sun)

    and the moment of inertia is:
     ===> I = k2**2 * Mstar * R**2

    END OF NOTE

    ------------------------------------------------------------------------

    '*_res files are based on the original grid, but were resampled to
    guarantee a regular timesampling across different mass tracks


    """

    def __init__(self, full=False):
        import os
        datadir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'data/stellar_model/BHAC15/')
        self.M_o = 1.99e33 << u.g
        self.R_o = 6.96e10 << u.cm
        self.I_o = 7e53 << u.g*(u.cm**2)
        self.masses = np.array([0.075, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5,
                                0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3,1.4])
        self.BHAC15 = {}
        for key in ['M_Ms', 'log_t_yr', 'Teff',  'R_Rs', 'I_Is', 'k2conv',
                    'k2rad']:
            self.BHAC15[key] = []
        if bool(full):
            for key in ['L_Ls', 'g', 'Log_Li_Li0', 'logTc', 'logROc', 'Mrad',
                         'Rrad']:
                self.BHAC15[key] = []
        for mass in self.masses:
            filename = '{:5.3f}'.format(mass)[0] + 'p'
            filename += '{:5.3f}'.format(mass)[2:] + 'Msu_masstrack_res.txt'
            data = Table.read(datadir + filename, format='ascii')
            for key in ['M_Ms', 'log_t_yr', 'Teff',  'R_Rs', 'k2conv',
                        'k2rad']:
                self.BHAC15[key].extend(data[key])
            if bool(full):
                for key in ['L_Ls', 'g', 'Log_Li_Li0', 'logTc', 'logROc',
                            'Mrad', 'Rrad']:
                    self.BHAC15[key].extend(data[key])
            self.BHAC15['I_Is'].extend((data['k2conv']**2 + data['k2rad']**2)*
                                       data['M_Ms']*data['R_Rs']**2)
        valid=np.isfinite(self.BHAC15['M_Ms'])
        for key in self.BHAC15.keys():
            self.BHAC15[key]=np.array(self.BHAC15[key])[valid]
        self.full = full

    def getMassTrack(self, mass):
        if mass in self.masses:
            select_mass = (self.BHAC15['M_Ms'] == mass)
            self.Teff = self.BHAC15['Teff'][select_mass]
            self.Mass = self.BHAC15['M_Ms'][select_mass]
            self.Age = 10**self.BHAC15['log_t_yr'][select_mass]
            self.R = self.BHAC15['R_Rs'][select_mass]
            self.I = self.BHAC15['I_Is'][select_mass]
            if bool(self.full):
                self.log_Li = self.BHAC15['Log_Li_Li0'][select_mass]
                self.Tcentral = 10**self.BHAC15['logTc'][select_mass]
                self.Rhocentral = 10**self.BHAC15['logROc'][select_mass]
                self.Mrad = self.BHAC15['Mrad'][select_mass]
                self.Rrad = self.BHAC15['Rrad'][select_mass]
                self.Lum = self.BHAC15['L_Ls'][select_mass]
                self.logg = self.BHAC15['g'][select_mass]
                self.k2conv = self.BHAC15['k2conv'][select_mass]
                self.k2rad = self.BHAC15['k2rad'][select_mass]
                self.k2_2 = self.BHAC15['k2rad']**2 + self.BHAC15['k2conv']**2
        elif (mass < 0.075) or (mass > 1.4):
             print('outside mass range')
        else:
            i = np.argsort(abs(self.masses - mass))[0]
            if self.masses[i] > mass:
                mass_b = self.masses[i]
                mass_a = self.masses[i - 1]
            else:
                mass_a = self.masses[i]
                mass_b = self.masses[i + 1]
            x = [mass_a, mass_b]
            select_a = (self.BHAC15['M_Ms'] == mass_a)
            select_b = (self.BHAC15['M_Ms'] == mass_b)
            def interpolate(key, select_a = select_a, select_b = select_b):
                return np.array([np.interp(mass, x, y) for y in
                                 list(zip(self.BHAC15[key][select_a],
                                          self.BHAC15[key][select_b]))])
            self.Mass = mass
            self.Teff = interpolate('Teff')
            self.Age = 10**interpolate('log_t_yr')
            self.R = interpolate('R_Rs')
            self.k2rad = interpolate('k2rad')
            self.k2conv = interpolate('k2conv')
            self.k2_2=self.k2conv**2+self.k2rad**2
            self.I=self.k2_2*self.Mass*(self.R**2)
            if bool(full):
                self.Lum = interpolate('L_Ls')
                self.logg = interpolate('g')
                self.log_Li = interpolate('Log_Li_Li0')
                self.Tcentral = 10**interpolate('logTc')
                self.Rhocentral = 10**interpolate('logROc')
                self.Mrad = interpolate('Mrad')
                self.Rrad = interpolate('Rrad')
