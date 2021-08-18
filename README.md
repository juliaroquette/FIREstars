# Far-ultraviolet Irradiated Rotational Evolution model for low mass stars (FIREstars)

**Author: Julia Roquette**


This package, which is part of the research presented in Roquette et al. (2021), includes a compilation of tools for computing rotational evolution models of stars evolving under the irradiation of far-ultraviolet radiation. The package can be used to derive the rotational evolution of stars in the mass range 0.1-1.3<img src="https://render.githubusercontent.com/render/math?math=\rm{M}_\odot">.

-----------
## Installation

This package requires a `Python3` installation with [`astropy`](https://www.astropy.org/), [`matplotlib`](https://matplotlib.org/), [`numpy`](https://numpy.org/), [`pandas`](https://pandas.pydata.org/) and [`scipy`](https://www.scipy.org/).

To install, first clone this repository:

```bat
git clone https://github.com/juliaroquette/FIREstars.git
```

Next, use a terminal to access the repository directory and use `pip` to install:

```bat
pip install .
```

-----------

## Rotational Evolution model: `SpinEvolution.py`

The module [`SpinEvolution.py`](https://github.com/juliaroquette/FIREstars/blob/main/FIREstars/SpinEvolution.py) contains a class `SpinEvolutionCode` which comprises all tools to be used for modelling the rotational evolution of stars. Research notes on the implementation of this class can be found in the jupyter-notebook [jupyter/SpinEvolutionModel.ipynb](https://github.com/juliaroquette/FIREstars/blob/main/jupyter/SpinEvolutionModel.ipynb).
The model implemented here is based on three assumptions:

1. *Disk-locking:* During the early-PMS phase, stars with disks are locked to their disks and remain with constant rotation. The duration of the disk-locking phase is parametrised using FUV-irradiated disk-dissipation models by [Winter et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020MNRAS.491..903W/abstract) and a tool for reading and using disk-models is implemented as part of [`FUVfunctions.py`](https://github.com/juliaroquette/FIREstars/blob/main/FIREstars/FUVfunctions.py). See notes on [jupyter/FUV_TauD.ipynb](https://github.com/juliaroquette/FIREstars/blob/main/jupyter/FUV_TauD.ipynb).

2. *Internal Structure:* The internal structure of stars is described by the stellar evolution models of [Baraffe et al. (2015)](https://ui.adsabs.harvard.edu/abs/2015A%26A...577A..42B/abstract). Tools for reading these models are implemented in [`StarEvolution.py`](https://github.com/juliaroquette/FIREstars/blob/main/FIREstars/StarEvolution.py). See notes on REFERENCE.

3. *Magnetised Winds*: At all ages, stars are subject to a wind-torque. The wind-torque adopted is the one by [Matt et al. (2015)](https://ui.adsabs.harvard.edu/abs/2015ApJ...799L..23M/abstract).

### Usage:


```python
from FIREstars.SpinEvolution import SpinEvolutionCode
```

Initialize the class by providing an initial time for the model, `t0`, in years.

```python
spin = SpinEvolutionCode(t0)
```

Following initialization, the rotational evolution is estimated by `spin.dOmegadt`:

```python
time, omega = spin.dOmegadt(M, Omega0, t, tau_d=0, e=0.01, wind=True, structure=True, snapshot=False, breakup=True)
```

Where:
  - `M`: is the stellar mass in <img src="https://render.githubusercontent.com/render/math?math=\rm{M}_\odot">.
  - `Omega0` is the initial rotation rate, in <img src="https://render.githubusercontent.com/render/math?math=\Omega_\odot">.
  - `t` is one or multiple ages for model calculation:
    - If `snapshot=True`, `t` must be an array or list of ages at which snapshot values of the model will be saved.
    - Otherwise if `snapshot=False`, `t` must be the final age in the model.
  - `tau_d` is the disk-locking duration (in years) set as zero as default.
  - `e` is a tolerance factor which defines the size of the time-step in foward-step Euler-method. For example, in the default `e=0.01`, the model is estimating using time-steps large enough to increase <img src="https://render.githubusercontent.com/render/math?math=\Omega(t)"> in 1% in each step.
  -  `wind` turns on/off the wind-torque.
  -  `structure` turns on/off the structure term.
  -  `snapshot` sets the way in which the model's outputs are returned:
    - If `snapshot=False` (default) the model will return pairs of <img src="https://render.githubusercontent.com/render/math?math=\Omega(t)"> and <img src="https://render.githubusercontent.com/render/math?math=t"> at each time-step required for running the model with a tolerance factor `e`.
    - If `snapshot=True`, <img src="https://render.githubusercontent.com/render/math?math=\Omega(t)"> and <img src="https://render.githubusercontent.com/render/math?math=t"> are returned only for pre-defined timesteps provided inputed in the array `t`.
- `breakup` if set to `True`, the spin evolution model has a condition that prevents stars from rotating faster than the break-up limit.

**Example 1**: Spin evolution of a <img src="https://render.githubusercontent.com/render/math?math=1\rm{M}_\odot"> star from 0.5 Myr to 4.5 Gyr, with a disk-lifetime of 5 Myr and an initial rotational period of 8 d.

First, to transform the period of 8 d to <img src="https://render.githubusercontent.com/render/math?math=\Omega"> and vice-versa, [`SpinEvolution.py`](https://github.com/juliaroquette/FIREstars/blob/main/FIREstars/SpinEvolution.py)  includes the functions `period2omega` and `omega2period`:

```python
from astroFIREstars.SpinEvolution import  SpinEvolutionCode, period2omega
spin = SpinEvolutionCode(0.5e6)
time, omega = spin.dOmegadt(1.0, period2omega(8.), [4.5e9], tau_d=5e6)
```

**Example 2**: Rotational period of a <img src="https://render.githubusercontent.com/render/math?math=1\rm{M}_\odot"> star at the ages 10 Myr, 120 Myr and 4.5 Gyr, considering a disk-lifetime of 5 Myr and an initial rotational period of 8 d at the age 1 Myr.

```python
from FIREstars.SpinEvolution import  SpinEvolutionCode, period2omega
spin = SpinEvolutionCode(1e6)
time, omega = spin.dOmegadt(1.0, period2omega(8.), [10e6, 120e6, 4.5e9], tau_d=5e6, snapshot=true)
```

### Further usage

#### Break-up Limit:
`SpinEvolution` can also be used to estimate the break-up limit as a function of age.

**Example 3:** breakup limit for a <img src="https://render.githubusercontent.com/render/math?math=1\rm{M}_\odot"> at the ages 10 Myr and 1Gyr.


 ```python
 spin = SpinEvolutionCode(1e6)
 omega_crit = spin.get_BreakUp(1., [10e6, 1e9])
```

#### Saturation limit:

Similarly, the saturation limit used in the wind-torque can also be tracked.

**Example 4:** saturation limit for a <img src="https://render.githubusercontent.com/render/math?math=1\rm{M}_\odot"> at the ages 10 Myr and 1Gyr.

```python
 spin = SpinEvolutionCode(1e6)
 omega_sat = spin.get_SaturationLimit1., [10e6, 1e9])
```

## Rotational Evolution in Period-Mass space

The class `SpinEvolutionCode` also includes a tool for calculating "isogyrochrones" (see notes in [jupyter/PeriodMassDiagrams.ipynb](https://github.com/juliaroquette/FIREstars/blob/main/jupyter/PeriodMassDiagrams.ipynb)). Isogyrochrones are rotation tracks in period-mass space, in which the rotation of stars in the whole mass range 0.1-1.3 <img src="https://render.githubusercontent.com/render/math?math=1\rm{M}_\odot"> is evolved from a common set of initial conditions to a given age.

```python
spin = SpinEvolutionCode(t0)
mass, period = spin.isogyrochrone(initial_period, time, fuv=False, tau_d=False,
                                  dm=0.025, e=0.01, tau_vis=1.0, wind=True,
                                  structure=True, breakup=True, initial_age=False,
                                  get_breakup_limit=False)
```

Where:

- `initial_period` is the rotational period at the time `t0`.
- `time` provide a single or multiple ages (in years) at which the isogyrochrones will be estimated.
- `fuv` set as `False` for a model ignoring the environment. Otherwise, `fuv` is the local far-ultraviolet flux in Habing flux units (<img src="https://render.githubusercontent.com/render/math?math=\rm{G}_0">).
- `tau_vis` sets the viscous timescale in the disk-dissipation model (see notes in [jupyter/FUV_TauD.ipynb](https://github.com/juliaroquette/FIREstars/blob/main/jupyter/FUV_TauD.ipynb)). Possible values are 1, 2 or 5 Myr.
- `initial_age` if `True`, includes the isogyrochrone ate the age `t0`.
- `get_breakup_limit` if `True` returns an isogyrochrone with the break-up limit.
- `tau_d` is set as `False` if using the keyword `fuv`. Otherwise, `tau_d` sets a common disk-lockig duration for all the stars in the isogyrochrone.
- The remaining keywords are the same as in `dOmegadt`

**Example 5:** Creates an isogyrochrone at the age of 15 Myrs, for stars with initial rotation of 8 d at an age of 1 Myr, which are evolving under an FUV-flux of 1,000<img src="https://render.githubusercontent.com/render/math?math=\rm{G}_0">. Using steps of <img src="https://render.githubusercontent.com/render/math?math=0.05\rm{M}_\odot"> and considering a disk with viscosity timescale of 1 Myr.

```python
spin = SpinEvolutionCode(1e6)
mass, period = spin.isogyrochrone(8., 15e6, fuv=1000., dm=0.05, e=0.01, tau_vis=1.0)
```

------

## Observational datasets:

Roquette et al. (2021) compare the models implemented in this package with rolling-percentiles of observed datasets of rotation rates in open clusters. These datasets are processed in the module [ObservedDatasets.py](https://github.com/juliaroquette/FIREstars/blob/main/FIREstars/ObservedDatasets.py). Each of the datasets is processed in a class with the name of the observed cluster.

- `hPer`



https://github.com/juliaroquette/FIREstars/tree/main/data/observations

[jupyter/PeriodMassDiagrams.ipynb](https://github.com/juliaroquette/FIREstars/blob/main/jupyter/PeriodMassDiagrams.ipynb)

------

## `FUVfunctions.py` module

The module

------

## `StarEvolution.py` module


