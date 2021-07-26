#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 13:25:10 2021

@author: jroquette
"""

import numpy as np
from astropy.table import Table
import pandas as pd
import os
datadir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'data/observations/')     
 
class hPer:
    """
    hPer database 
    
    Reference table is Moraux+2013
    http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A+A/560/A13    
    Reference for Spectral types is Currie et al. 2010
    http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/ApJS/186/191
    Reference for Mass transformations: Roquette et al. 2021 Appendix A
    ---
    Usage:
        hPer = spin.ObservedDatasets.hPer()
    ---        
    LAST UPDATE: 
        21 July 2021
    """
    def __init__(self, filename='hPer.fit'):
        print('Loading data,RA,Dec,Prot,Mass,Amp,SpT,EBV')
        print('Binary stars were removed following the flags in Moraux+2013')
        self.data = Table.read(datadir + filename, format='fits') 
        valid = np.isfinite(self.data['Per']) & \
                (self.data['Mass_B15_median'] >0)
        self.RA = self.data['RAJ2000'][valid]
        self.Dec = self.data['DEJ2000'][valid]
        self.Prot = self.data['Per'][valid]
        self.Mass = self.data['Mass_B15_median'][valid]           
        self.Amp = self.data['Amp'][valid]
        self.SpT = self.data['SpectralType'][valid]
        self.EBV = self.data['E_B-V_'][valid]

        
class NGC2264:
    """
    Reference table is a combination of:
        1. Venuti et al. 2017: 
            http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A+A/599/A23
        4. Lamm et al. 2005: 
            http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A+A/430/1005
    Reference for disks: Alana Souza (Private comunication)
    References for Mass transformations:  Roquette et al. 2021 Appendix A
    
    ----
    
    Usage:
        ngc2264 = spin.ObservedDatasets.NGC2264()
  
    ----
    
    Notes:
        
        ngc2264.Disk:   = 0 if unknown
                        > 0 has a disk
                        < 0 diskless
    -----
    
    LAST UPDATE: 
        21 July 2021
                
    """
    def __init__(self, filename='NGC2264.fit'):
        self.data = Table.read(datadir + filename, format='fits') 
        print('Loading data, RA, Dec, Prot, Mass, SpT, Av, Disk')
        valid = (self.data['Adopted_Period'] > 0.) & \
                (self.data['Mass_B15_Adopted'] > 0.)
        self.RA = self.data['RA'][valid]
        self.Dec = self.data['DE'][valid]
        print("Period adopted from Venuti+17, or Lamm+05 otherwise")
        self.Prot = self.data['Adopted_Period'][valid]
        self.Mass = self.data['Mass_B15_Adopted'][valid]
        self.AV = self.data['Av_V14'][valid]
        self.SpT = self.data['spt_V14']  [valid]  
        self.Disk = self.data['Disked'][valid]
        
class USco:
    """
    Reference table:
        Rebull et al. 2018  
        http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/AJ/155/196

    Reference for disks: Rebull et al. 2018
    References for Mass: Roquette et al. 2021 Appendix A
    
    ----
    
    Usage:
        USco = spin.ObservedDatasets.USco()
  
    ----
        
    Notes:
        ngc2264.Disk:   = 0 if unknown
                        > 0 has a disk
                        < 0 diskless
    -----
    
    LAST UPDATE: 
        21 July 2021
                
    """
    def __init__(self,filename = 'USco.fit'):
        self.data = Table.read(datadir+filename, format='fits') 
        valid = (self.data['Per1'] > 0.) & (self.data['Mass_B15'] > 0.)
        print('Loading data,RA,Dec,Prot,Mass,EBV,Disk')
        self.RA = self.data['RAJ2000'][valid]
        self.Dec = self.data['DEJ2000'][valid]
        self.Prot = self.data['Per1'][valid]
        self.Mass = self.data['Mass_B15'][valid]             
        self.EBV = self.data['E_B-V'][valid]
        self.Disk = self.data['Disked'][valid]
 
class Praesepe:
    """
    Praesepe database 
    
    Reference table is Rebull+17
    https://cdsarc.unistra.fr/viz-bin/cat/J/ApJ/839/92
    References for Mass: Roquette et al. 2021 Appendix A
    
    ----
    
    Usage:
        praesepe = spin.ObservedDatasets.Praesepe()
  
    -----
    
    LAST UPDATE: 
        21 July 2021
    """
    def __init__(self, filename='Praesepe.fit'):
        print('Loading data,RA,Dec,Prot,Mass')
        
        self.data = Table.read(datadir+filename, format='fits') 
        valid = (self.data['PPer_R17'] > 0.) & \
            (self.data['Mass_Ks_B15_solar_med'] > 0)
        self.RA = self.data['RAJ2000_R17'][valid]
        self.Dec = self.data['DEJ2000_R17'][valid]
        self.Prot = self.data['PPer_R17'][valid]
        self.Mass = self.data['Mass_Ks_B15_solar_med'][valid]        
                                         
class Pleiades:
    """
    Pleiades database 
    
    Reference table is Rebull+16
    https://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=J/AJ/152/113/table2
    References for Mass: Roquette et al. 2021 Appendix A
    
    ----
    
    Usage:
        praesepe = spin.ObservedDatasets.Praesepe()
  
    -----
    
    LAST UPDATE: 
        21 July 2021
    """
    def __init__(self, filename='Pleiades.fit'):
        print('Loading data,RA,Dec,Prot,Mass,Amp')
        self.data = Table.read(datadir+filename, format='fits') 
        valid = (self.data['Prot_R17'] > 0.) & \
                (self.data['Mass_Ks_B15_solar_med'] > 0.)
        self.RA = self.data['RAJ2000_R17'][valid]
        self.Dec = self.data['DEJ2000_R17'][valid]
        self.Prot = self.data['Prot_R17'][valid]
        self.Amp=self.data['Amp_R17'][valid]
        self.Mass = self.data['Mass_Ks_B15_solar_med'][valid]
                                
class NGC6811:
    """
    NGC6811 database 
    
    Reference table is Curtis+19 
    https://iopscience.iop.org/article/10.3847/1538-4357/ab2393/pdf
    But also available in Meibom+11

     References for Mass: Roquette et al. 2021 Appendix A
    
    ----
    
    Usage:
        praesepe = spin.ObservedDatasets.Praesepe()
  
    -----
    
    LAST UPDATE: 
        21 July 2021
    """
    def __init__(self,filename='NGC6811.fit'):
        print('Loading data,RA,Dec,Prot,Mass,Amp')
        self.data = Table.read(datadir+filename, format='fits') 
        valid = (self.data['Per_C19'] > 0.) & \
                (self.data['Mass_Gaia_B15_solar_med'] > 0.)
        self.RA = self.data['ra_C19'][valid]
        self.Dec = self.data['dec_C19'][valid]
        self.Prot = self.data['Per_C19'][valid]
        self.Mass = self.data['Mass_Gaia_B15_solar_med'][valid]
        self.Teff = self.data['Teff_C19'][valid]
        self.SpT = self.data['SpT_C19'][valid]
 
    
def PeriodMass_RollingPercentile(x, y, percentile, window_size=0.14, 
                                 secondary_window_size=0.05, min_data=12, 
                                 secondary_min_data=4):
    """
    Estimates the rolling percentile of a period mass distribution

    ----
    
    INPUT:
        x:                      mass
        y:                      periods
        percentile:            [0.,1.] desired percentile
        window_size:           [0.,1.] percentage of the total amount of 
                               datapoints  to be used as window size for 
                               estimating the rolling statistics
        secondary_window_size: secondary window size for calculating rolling 
                               statistics at higher masses, where typically 
                               there are less datapoints available
        min_data:              minimum datapoints to be used for calculating 
                               the rolling statistics close to the extremity
        secondary_min_data:    same as min_data, but for the portion of the 
                               distrbution at higher masses
    ----
    OUTPUT:
        returns a curve of estimated percentiles at each x
         x_: masses - it's effectivelly the same as x, but ordered for 
             ascending masses
         out: periods - percentile at each point
    """
    valid = np.isfinite(x) & np.isfinite(y)
    x = np.array(x[valid])
    y = np.array(y[valid])
    sort_ = np.argsort(x)
    x_ = x[sort_]
    out=np.full((len(x)),np.nan)
    series = pd.Series(y[sort_])
    wd = int(window_size*len(x_))
    wd_ = int(secondary_window_size*len(x_))
    out[:-wd] = series.rolling(window=wd, min_periods=min_data).quantile\
                (percentile, interpolation='lower').values[:-wd]
    out[-wd:] = series.rolling(window=wd_, min_periods=secondary_min_data)\
              .quantile(percentile, interpolation='higher').values[-wd:]
    return x_, out