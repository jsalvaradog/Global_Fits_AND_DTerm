import shelve, logging, sys, copy
logging.basicConfig(level=logging.ERROR)
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.stats import norm
import importlib_resources
import matplotlib
import matplotlib.pyplot as plt

import gepard as g
from gepard.fits import th_KM15, KM15
from gepard.fits import th_AFKM12
from gepard.fits import GLO15b
from gepard import cff, data, dvcs, eff, gpd

import os
import json

if not os.path.exists('figs'):
    os.makedirs('figs')

import numpy as np
import matplotlib.pyplot as plt
import gepard as g


class CFFPlotter:
    def __init__(self, theory):
        """
        Plot ImH and ReH as functions of xi for fixed t and Q2.

        Parameters
        ----------
        theory : Gepard theory object
            Fitted Gepard theory (e.g. th)
        """
        self.theory = theory

    @staticmethod
    def Dterm(t, D0=-1.6, M=0.7, p=3):
        """D-term parametrization"""
        return D0/(1-t/M**2)**p

    @staticmethod
    def pressure(r, Dterm, **args):
        """
        pressureN function
        Pressure distribution r^2*p(r) in GeV/fm for given Dterm(t) function
        """
        GeVfm = 0.197
        Mp = 0.938
        rGeV = r/GeVfm
        # variable change x=-t^2, so QUADPACK Fourier integral method can be used
        intg = quad(lambda x: 2*x*(-x**2)*Dterm(-x**2, **args), 0, np.inf, weight='sin', wvar=rGeV)
        return intg[0]/(24*np.pi**2*Mp*rGeV) * rGeV**2 / GeVfm

    def plot(self, ax=None, label_prefix=""):
        """
        Produce proton pressure with uncertainty bands.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Existing axes to plot on. If None, a new figure and axes are created.
        label_prefix : str, optional
            Prefix added to legend labels (e.g. theory name or fit identifier).

        Returns
        -------
        fig : matplotlib.figure.Figure
            The figure containing the plot.
        """
        
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 6))
        else:
            fig = ax.figure

        plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)

        # Sampling from Gaussian/normal distribution using scipy's stats routines
        ss = 200 #  sample size
        print("D(t=0)=", self.theory.parameters['C'], " +- ", self.theory.parameters_errors['C'])
        print("M=", self.theory.parameters['mC2']**0.5, " +- ", self.theory.parameters_errors['mC2']**0.5)
        d0, M, p = -1.0*2.607434678826921*(18./25.)*(5./9.), self.theory.parameters['mC2']**0.5, 3
        d0err = np.sqrt(0.14**2 + 0.33**2) * (18./25.)*(5./9.)
        M2err = np.sqrt(0.1**2 + 0.15**2)
        Merr = M2err/2/M

        d0MC = norm.rvs(d0, d0err, size=ss)
        MMC = norm.rvs(M, Merr, size=ss)

        rs = np.linspace(0.01,2,500)
        psmean = []
        psup = []
        psdown = []
        for r in rs:
            pressMC = []
            for k in range(ss):
                pressMC.append(self.pressure(r, self.Dterm, D0=d0MC[k], M=MMC[k], p=3))
            mean, std = 1000*np.mean(pressMC), 1000*np.std(pressMC)
            psmean.append(mean)
            psup.append(mean+std)
            psdown.append(mean-std)
            
        ax.plot(rs, psmean, 'k-')
        ax.fill_between(rs, psup, psdown,
                            edgecolor='black', linewidth=1, zorder=0, label=label_prefix)
        ax.set_ylabel(r'$r^2 p(r) \quad[\times 10^{-3}\; {\rm GeV}/{\rm fm}]$', fontsize=25)
        ax.set_xlabel(r'$r\quad[{\rm fm}]$', fontsize=25)
        ax.tick_params(axis='both', labelsize=20)
        
        ax.axhline(linestyle='--', color='g', lw=1)
        ax.set_xlim(0,2)
        #ax.set_ylim(-8, 17)
        ax.legend(loc="upper right", fontsize=15)
        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(5))
        ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1))    
        return fig
    
            
class KMA(KM15):  # noqa: D101, E501
    def subtraction(self, pt):
        #Dispersion relations subtraction constant.
        return self.parameters['C']/(1.-pt.t/self.parameters['mC2'])**3

class KMA_H(KMA):
  def subtraction(self, pt):
      #Dispersion relations subtraction constant.
      return self.parameters['C']/(1.-pt.t/self.parameters['mC2'])**3
  def ImE(self, pt: data.DataPoint):
      return 0
  def ImHt(self, pt: data.DataPoint):
      return 0
  def ImEt(self, pt: data.DataPoint):
      return 0
  def ReE(self, pt: data.DataPoint):
      return 0
  def ReHt(self, pt: data.DataPoint):
      return 0
  def ReEt(self, pt: data.DataPoint):
      return 0

#########################
# initial configuration
#########################

plt.rc('text', usetex=True)
params = {'text.latex.preamble' : [r'\usepackage{amssymb}', r'\usepackage{amsmath}']}
# Join the list of packages into a single string
params['text.latex.preamble'] = '\n'.join(params['text.latex.preamble'])
plt.rcParams.update(params)

#########################
# Create model and fit
#########################

par_KMA = {'tmv2': 15.94293227053628, 'rS': 1.0, 'alv': 0.43, 'tal': 0.43,
            'mpi2': 15.999998816676918, 'Nv': 1.35, 'rv': 0.918393047884448,
            'Nsea': 0.0, 'alS': 1.13, 'rpi': 2.6463144464701536,  'alpS': 0.15,
            'C': 2.7678681812890016, 'tNv': 0.6, 'bS': 2.0, 'tbv': 0.4000000003259146,
            'bv': 0.40000206775282354, 'mv2': 0.6228490959696493, 'alpv': 0.85,
            'talp': 0.85, 'mC2': 1.4497411858248308, 'trv': 0.881085721967267,
            'pows': 2.0, 'Ealpg': 0.15, 'md2': 1.0, 'mg2': 0.7, 'Ealps': 0.15,
            'powg': 2.0,
            'al0d': 0.5, 'delms2': 0.0,  'delmg2': 0.0, 'nd': 1.0, 'alpd': 1.0,
            'ns':  0.15203911208796006, 'al0s': 1.1575060246398083, 'alps': 0.15,
            'ms2': 0.4818827240886959, 'secs': 1.0707825621025808, 'ng': 0.5,
            'this': -0.36618269477432946, 'al0g': 1.247316701070471, 'alpg': 0.15,
            'mg2': 0.7, 'secg': -2.990809378821039, 'thig': 0.9052207712570559,
            'kaps': 0.0, 'kapg': 0.0}

th1 = th_KM15;
th2 = KMA()
th3 = KMA_H()

type = "Global_Fit" # "Global_Fit" or "DTerm"

model_name_1 = 'KM15'
model_name_2 = 'KMA'
model_name_3 = 'KMA'

th1.name = model_name_1
th2.name = model_name_2
th3.name = model_name_3

th1.parameters.update(par_KMA)
th2.parameters.update(par_KMA)
th3.parameters.update(par_KMA)

f1 = g.MinuitFitter(GLO15b, th1)
f2 = g.MinuitFitter(GLO15b, th2)
f3 = g.MinuitFitter(GLO15b, th3)

f1.fix_parameters('ALL')
f2.fix_parameters('ALL')
f3.fix_parameters('ALL')

f1.release_parameters('mv2', 'rv', 'bv', 'C', 'mC2', 'tmv2', 'trv', 'tbv', 'rpi', 'mpi2', 'ms2', 'secs', 'this', 'secg', 'thig')
f2.release_parameters('mv2', 'rv', 'bv', 'C', 'mC2', 'tmv2', 'trv', 'tbv', 'rpi', 'mpi2', 'ms2', 'secs', 'this', 'secg', 'thig')
f3.release_parameters('mv2', 'rv', 'bv', 'C', 'mC2', 'tmv2', 'trv', 'tbv', 'rpi', 'mpi2', 'ms2', 'secs', 'this', 'secg', 'thig')

# Load saved fit results
with open(f"{type}_{model_name_1}.json", "r") as fitres:
    fit_results_1 = json.load(fitres)
with open(f"{type}_{model_name_2}.json", "r") as fitres:
    fit_results_2 = json.load(fitres)
with open(f"{type}_{model_name_3}.json", "r") as fitres:
    fit_results_3 = json.load(fitres)

# Restore Minuit values
for par, val in fit_results_1["values"].items():
    f1.minuit.values[par] = val
for par, val in fit_results_2["values"].items():
    f2.minuit.values[par] = val
for par, val in fit_results_3["values"].items():
    f3.minuit.values[par] = val

# Restore Minuit errors
for par, err in fit_results_1["errors"].items():
    f1.minuit.errors[par] = err
for par, err in fit_results_2["errors"].items():
    f2.minuit.errors[par] = err
for par, err in fit_results_3["errors"].items():
    f3.minuit.errors[par] = err

# Synchronize values/errors back to the Gepard theory object
f1.theory.parameters.update(f1.minuit.values.to_dict())
f2.theory.parameters.update(f2.minuit.values.to_dict())
f3.theory.parameters.update(f3.minuit.values.to_dict())

f1.theory.parameters_errors = f1.minuit.errors.to_dict()
f2.theory.parameters_errors = f2.minuit.errors.to_dict()
f3.theory.parameters_errors = f3.minuit.errors.to_dict()

# restore theory parameters
th1.parameters.update(f1.minuit.values.to_dict())
th2.parameters.update(f2.minuit.values.to_dict())
th3.parameters.update(f3.minuit.values.to_dict())


# restore errors
th1.parameters_errors = f1.minuit.errors.to_dict()
th2.parameters_errors = f2.minuit.errors.to_dict()
th3.parameters_errors = f3.minuit.errors.to_dict()

#f.print_parameters()

plotter1 = CFFPlotter(th1)
plotter2 = CFFPlotter(th2)
plotter3 = CFFPlotter(th3)

fig, ax = plt.subplots(figsize=(8, 6))

plotter1.plot(ax=ax, label_prefix="KM15")
plotter2.plot(ax=ax, label_prefix="KMA")
plotter3.plot(ax=ax, label_prefix="KMA_H")

plt.show()

plt.show()
plt.savefig(f'figs/Comparison_DTerm_{type}.pdf')
