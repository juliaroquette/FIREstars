"""
Created on Tue July 13 12:01 2021

@author: juliaroquette


Spin evolution model implemented in the development of the study by
Roquette et al. 2021.

The class SpinEvolutionCode performs the spin evolution modeling of stars
in the mass range 0.1-1.3 Msun. For a details on the theoretical background
in the model, see setailed notes in SpinEvolutionModel.ipynb


"""
import numpy as np
import astropy.units as u

from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d
from astropy.constants import G

from FIREstars.FIREstars.StarEvolution import BHAC15_MassTrack

class SpinEvolutionCode:
    """
    My spin evolution code

    """
    def __init__(self, t0):
        """
        M should be in solar masses
        to should be in years - the initial time step of the model
        Omegao should be in days
        """
        #define constants
        self.CHI = 10.
        self.P = 2.
        self.M_o = 1.99e33 << u.g
        self.OMEGA_o = u.def_unit(r'\Omega_\odot',  2.6*1e-6*u.Hz)
        self.R_o = 6.96e10 << u.cm
        self.I_o = 7e53 << u.g*(u.cm**2)
        self.TIME_o = 4.55e9 << u.yr
        self.TAU_CZ_o = 13.7798 << u.d # Value updated in the context of B15
                        # 12.9 << u.d
        self.baraffe = BHAC15_MassTrack() #loads the Baraffe+15 grid.
        self.t0 = t0 << u.yr
        # Makes sure that the initial time isn't before the initial time in
        #the BHAC15 grid
        if self.t0.value < 10**np.nanmin(self.baraffe.BHAC15['log_t_yr']):
            self.t0 = 10**np.nanmin(self.baraffe.BHAC15['log_t_yr']) << u.yr
            print('Minimum value must be inside the age range of the Baraffe \
                   Model, I reset it to {0}'.format(self.t0.value/1e6))

    def interpolateBaraffe(self,t):
        """
        given a time t in yrs, interpolate parameters
        of interest from Baraffe+2015 models
        and update local variables for R, I and T
        --
        input:
        t: time in yrs
        """
        #bouds_error=False will returns np.nan if out of bounds
        self.R, self.I, self.T = interp1d(self.baraffe.Age,
                                          (self.baraffe.Rarius, 
                                           self.baraffe.InertiaMomentum,
                                            self.baraffe.Teff),
                                              kind='linear', bounds_error=False
                                            )(t)

    def update_dIdt(self, t):
        """
        calculates the derivative dI/dt at a position t
        and updated the self.dIdt object
        """
        y = self.baraffe.InertiaMomentum*self.M_o*self.R_o**2
        x = (self.baraffe.Age*u.yr).to(u.s)
        spl = UnivariateSpline(x, y, k=3)
        self.dIdt = spl.derivative(n = 1)(t.to(u.s).value)*(u.cm**2)*u.g/u.s

    def update_To(self):
        """
        updates the Matt+15 torque scalling term `self.To`
        """
        self.To = 6.3e30*((self.R)**3.1)*((self.M)**0.5) << u.erg

    def tau_w(self, Omega):
        """
        Derives the wind torque based on the current values of the objects for
        M, R and T, while the rotation Omega is given as input
        """
        if Omega < self.saturation_limit():
            return - self.To*((self.tau_cz/self.TAU_CZ_o)**self.P)*\
                ((Omega.value)**(self.P + 1))
        else:
            return - self.To*(self.CHI**self.P)*(Omega.value)

    def tau_cz_CS11(self):
        """
        Derives the convective turnover timescale as in Cranmer & Saar 2011
        """
        self.tau_cz = (314.24*np.exp( -(self.T/1952.5) - (self.T/6250.)**18) +
                       0.002) << u.d

    def time_update(self,t):
        """
        updates parameters at any time t
        """
        self.interpolateBaraffe(t.value) #this will set the values of R,I&T
        self.tau_cz_CS11() # set tau_cz
        self.update_To() # set To
        self.update_dIdt(t) #get dIdt

    def saturation_limit(self):
        """
        Gives the saturation limit based on the current
        local variables.
        """
        return (self.CHI*self.TAU_CZ_o/self.tau_cz).value << self.OMEGA_o


    def Euler_(self, f, Omega, t0, dt, wind=True, structure=True,
               breakup=True):
        """
        Estimates one step of the Euler method

        ----
        if breakup = True: saturates Omega at the critical value
        """
        O = Omega + dt*self.f(Omega, wind = wind, structure = structure)
        if bool(breakup):
            Ocrit = self.BreakUp()
            if (O > Ocrit):
                O = Ocrit
        return O, t0 + dt

    def f_W(self, Omega):
        """
        Wind torque term of the rotation evolution equation
        Omega is the rotation rate in solar units
        """
        return (self.tau_w(Omega)/(self.I*self.M_o*self.R_o**2))

    def f_dIdt(self, Omega):
        """
        Term including of the rotational evolution including the
        variation in the Momentum of Ineria

        """
        return - Omega*self.dIdt/(self.I*self.M_o*self.R_o**2)


    def f(self, Omega, wind=True, structure=True):
        """
        Set up a function that will be feed to the Euler method
        """
        return bool(wind).real*self.f_W(Omega) + \
               bool(structure).real*self.f_dIdt(Omega)

    def get_dt(self, Omega, e, structure=True, wind=True):
        """
        estimates the best timestep by comparing how much the wind and the
        structure will contribute to changing omega, given an efficiency term
        e
        """
        if bool(wind) and bool(structure):
            dt_s = - e*self.I*self.M_o*(self.R_o**2)/self.dIdt
            dt_w = e*self.I*self.M_o*(self.R_o**2)*(Omega)/self.tau_w(Omega)
            return min(abs(dt_s), abs(dt_w)).to(u.yr)
        elif not bool(wind) and bool(structure):
            return abs(e*self.I*self.M_o*(self.R_o**2)/self.dIdt)
        elif bool(wind) and not bool(structure):
            return abs(e*self.I*self.M_o*(self.R_o**2)*\
                       (Omega)/self.tau_w(Omega))

    def BreakUp(self):
        """
        estimates the breakup rotation
        for a star of given mass and age.
        Omega=sqrt(GM/Re**3)
        Re=1.5Rp
        Ro used in Baraffe models: Rs=6.96d10cm
        Mo=1.99e33g
        ---
        Returns the Breakup limit in OmegaSun units
        """
        return np.sqrt((self.M*self.M_o*G.cgs)/(1.5*self.R*self.R_o)**3)\
            << self.OMEGA_o


    def dOmegadt(self, M, Omega0, t, tau_d=0, e=0.1, wind=True,
                 structure=True, snapshot=False, breakup=True):
        """
        __
        input
        _
        M: stellar mass in solar masses
        Omega0: initial rotation
        t: vector with the key timesteps
        e: [0,1] tolerance of minimal variation of Omega do define a
           intermediary timestep
        wind: [True] for activate wind term
        structure: [True] for activate structure term
        snapshot: [True] the model will return only data for the timesteps
                  listed in t
        breakup:  [True] saturates the rotation at the break-up speed
        __
        output
        Omega and t with no units
        _

        """
        #define a local variable for the Mass
        self.M = M
        if M == 0.1:
            self.M += 0.000001
        # locally load the grid for M
        self.baraffe.getMassTrack(self.M)
        # define all initial parameters relevant to the spin evolution
        self.time_update(self.t0)
        #
        try:
            len(t)
        except TypeError:
            t = np.array([t])
        else:
            t = np.array(t)
        Omega0 = Omega0 << self.OMEGA_o
        if Omega0 >= self.BreakUp():
            print('Initial rotation faster than Break-Up!')
        t = t << u.yr
        if bool(snapshot):
            t_out = np.full(len(t) + 1, np.nan)
            t_out[0] = self.t0.value
            t_out[1:] = 1.*t
        if bool(snapshot): t_mas = []
        for i, T in enumerate(t):
            if T.value > np.nanmax(self.baraffe.Age): #makes sure I am not
                                                #going beyond Baraffe's models
                if bool(snapshot): t_mas.append(i + 1) # save the index
                t[i] = (np.nanmax(self.baraffe.Age) - 1e6) << u.yr
                print('Maximum value must be inside the age range of the \
                      Baraffe Model')
        tau_d= tau_d << u.yr
        if (tau_d != 0) & (tau_d.value < 1e5):
            print(r'Is $\tau_D$ in the right units?')
        t_ = []
        Omega_ = []
        t_.append(self.t0.value)
        Omega_.append(Omega0.value)
        tk_o = 1.*self.t0
        # test if spin-evolution needs to be calculated at all:
        if t.max() > tau_d:
            n = 0
            # test if tau_D is before t0:
            if tau_d > self.t0:

                while t[n] < tau_d:
                    t_.append(t[n].value)
                    Omega_.append(Omega0.value)
                    tk_o = t[n]
                    n += 1
                #next register the moment disk was lost:
                if not bool(snapshot):
                    t_.append(tau_d.value)
                    Omega_.append(Omega0.value)
                tk_o = tau_d
            #initiate the parameters
            self.time_update(tk_o)
            dt=self.get_dt(Omega0, e, wind = wind, structure=structure)
            for T in t[n:]:
                while tk_o + dt < T:
                    O, tk = self.Euler_(self.f, Omega0, tk_o, dt, wind=wind,
                                        structure=structure, breakup=breakup)
                    if not bool(snapshot):
                        t_.append(tk.value)
                        Omega_.append(O.value)
                    Omega0 = 1.*O
                    tk_o = tk
                    self.time_update(tk_o)
                    dt = self.get_dt(Omega0, e, wind=wind,
                                     structure=structure)
                self.time_update(tk_o)
                dt = T - tk_o
                O,tk = self.Euler_(self.f, Omega0, tk_o, dt, wind=wind,
                                   structure=structure, breakup=breakup)
                t_.append(tk.value)
                Omega_.append(O.value)
                Omega0 = 1.*O
                tk_o = tk
                self.time_update(tk_o)
                dt = self.get_dt(Omega0, e, wind=wind, structure=structure)
            t_ = np.array(t_)
            Omega_ = np.array(Omega_)
            if bool(snapshot):
                t_ = t_out
                Omega_[t_mas] = np.nan
            return t_, Omega_
        else:
            return np.insert(t.value, 0, self.t0.value, axis=0), \
                             np.full(len(t) + 1, Omega0.value)

    def get_BreakUp(self, M, t):
        """
        input:
            M in Msun
            t in years

        Returns the Breakup limit as a function of time for a star with
        mass M in Msun
        """
        self.M = M
        if M == 0.1:
            self.M += 0.000001
        # locally load the grid for M
        self.baraffe.getMassTrack(self.M)
        try:
            len(t)
        except TypeError:
            t = np.array([t])
        else:
            t = np.array(t)
        t = t << u.yr
        O_crit = []
        for t_ in t:
            self.time_update(t_)
            O_crit.append(self.BreakUp().value)
        return np.array(O_crit)

    def get_SaturationLimit(self, M, t):
        """
        input:
            M in Msun
            t in years

        Returns the saturation limit as a function of time for a star with
        mass M in Msun
        """
        self.M = M
        if M == 0.1:
            self.M += 0.000001
        # locally load the grid for M
        self.baraffe.getMassTrack(self.M)
        t = np.array(t)
        t = t << u.yr
        O_sat = []
        for t_ in t:
            self.time_update(t_)
            O_sat.append(self.saturation_limit().value)
        return np.array(O_sat)

    def isogyrochrone(self, initial_period, time, fuv=False, tau_d=False,
                      dm=0.025, e=0.01, tau_vis=1.0, wind=True,
                      structure=True, breakup=True, initial_age=False,
                      get_breakup_limit=False):
        """
        Calculate "isogyrochrones" by running the spin evolution model for the
        range of masses 0.1-1.3Mo, considering coeval stars sharing the same
        initial condition

        ----
        input:

        initial_period:            Initial Period for the isogyrochrones
        time:                      Isogyrochrone ages in years. It can be
                                   either a single value or an array of values
        fuv [False]:               Set as False if no FUV level is considered.
                                   Otherwise, give a value in units of G0 for
                                   deriving disk-lifetimes from the Winter+20
                                   models
        tau_d [False]:             disk-locking duration in years if want to
                                   inform this parameter by hand, otherwise
                                   leave as False and isogyrochrones will be
                                   run without any disk-locking
        dm [0.025]:                Sets how many steps of mass will be
                                   calculated (linearily) between the values
                                   0.1 and 1.3Mo
        e [0.01]:                  tolerance parameter for the spin evolution
                                   code. Sets by how much Omega will have to
                                   change for a time-step to be calculated
                                   with the Euler method
        wind [True]:               True for a Spin Evolution including the
                                   wind-torque
        structure [True]:          True for a spin evolution including the
                                   stellar structure variation terms
        breakup [True]:            True to prevent stars from ever rotating
                                   faster than the break-up limit
        initial_age [False]:       If True, output will include an
                                   isogyrochrone for the initial age which was
                                   defined when calling the
                                   SpinEvolutionCode(t0)
        get_breakup_limit [False]: If True, returns the break-up limit as a
                                    function of mass and age

        ----
        Usage:

        mass, period = isogyrochrone(initial_period, time, fuv=1000.)

        ----
        output:

        mass:   array with masses in the isogyrochrone
        period: array with rotational periods in the isogyrochrone. Dimensions
                are [dimensions of time array, mass dimension]
        """
        try:
            len(time)
        except TypeError:
            time = np.array([time])
        else:
            time = np.array(time)
        if bool(fuv):
            from fire.FUVfunctions import DiskWithFUV
            disk = DiskWithFUV()
            get_tauD = lambda x: disk.get_tauD(x, fuv, tau_vis=tau_vis)
        elif bool(tau_d):
            get_tauD = lambda x: tau_d
        else: get_tauD = lambda x: 0.
        mass = np.linspace(0.1, 1.3, int((1.3 - 0.1)/dm))
        period_out = np.full((bool(initial_age).real + len(time), len(mass)),
                             np.nan)
        if not bool(get_breakup_limit):
            period_out[0,:] = bool(initial_age).real*initial_period
            for i,m in enumerate(mass):
                period_out[:, i] = omega2period(self.dOmegadt(m,
                                    period2omega(initial_period), time,
                                    tau_d=get_tauD(m), e=e, wind=wind,
                                    structure=structure, snapshot=True,
                                    breakup=breakup)[1]\
                                            [1 - 1*bool(initial_age).real:])
        else:
            for i, m in enumerate(mass):
                period_out[:, i] = \
                    omega2period(self.get_BreakUp(m, np.insert(time, 0,
                                                    self.t0.value, axis=0)))
        return mass, period_out


def period2omega(period):
    """
    Convert Period (days) to Omega (OmegaSun)
    """
    period = period << u.d
    OmegaSun = 2.6e-6 << u.Hz
    return (2.*np.pi/(period.to(u.s))/OmegaSun).value

def omega2period(Omega):
    """
    Convert Omega (OmegaSun) to Period (days)
    """
    OmegaSun = 2.6e-6 << u.Hz
    Omega = Omega << OmegaSun
    return 2.*np.pi*(u.s).to(u.d)/(Omega.value)


