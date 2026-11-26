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
import math as m
import os
import json

if not os.path.exists('figs'):
    os.makedirs('figs')

import numpy as np
import matplotlib.pyplot as plt
import gepard as g

def Multipole_D(t, D0=-1.6, M=0.7, p=3):
        """D-term parametrization"""
        return D0/(1-t/M**2)**p


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

    def plot(self, ax=None, label_prefix="", style=''):
        """
        Produce proton D-term with uncertainty bands.

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

        plt.subplots_adjust(left=0.17, right=0.95, top=0.95, bottom=0.17)

        # Sampling from Gaussian/normal distribution using scipy's stats routines
        ss = 500 #  sample size
        
        d0, M, p = -1.0*self.theory.parameters['C']*18/25, self.theory.parameters['mC2']**0.5, 3.0
        d0err = self.theory.parameters_errors['C']*18/25
        M2err = self.theory.parameters_errors['mC2']
        perr = 0.5 #np.sqrt(0.25**2 + 0.5**2)
        Merr = M2err/(2*M)
        print("\n"+label_prefix)
        print("D(t=0)=", d0, " +- ", d0err)
        print("M^2=", M*M, " +- ", M2err)

        d0MC = norm.rvs(d0, d0err, size=ss)
        MMC = norm.rvs(M, Merr, size=ss)
        pMC = norm.rvs(p, perr, size=ss)
        
        ts = np.linspace(-0.45,0,100)
        psmean = []
        psup = []
        psdown = []
        for t in ts:
            dtermMC = []
            for k in range(ss):
                dtermMC.append(self.Dterm(t, D0=d0MC[k], M=MMC[k], p=pMC[k]))
            mean, std = np.mean(dtermMC), np.std(dtermMC)
            psmean.append(mean)
            psup.append(mean+std)
            psdown.append(mean-std)
            
        line = ax.plot(-ts, psmean)
        color = line[0].get_color()
        
        ax.fill_between(-ts, psup, psdown,edgecolor=color, linewidth=1, label=label_prefix, alpha=0.3, hatch=style)
        ax.set_ylabel(r'$D^{Q}(t)$', fontsize=30)
        ax.set_xlabel(r'$|t| [\text{GeV}^{2}/c^{4}]$', fontsize=30)
        ax.tick_params(axis='both', labelsize=25)
        
        ax.axhline(linestyle='--', color='gray', lw=1)
        ax.set_xlim(0,0.45)
        ax.set_ylim(-2.4, 1.0)
        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.1))    
        ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.02))    
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

th1 = KMA()
th2 = KMA()
th3 = KMA()

type = "Global_Fit" # "Global_Fit" or "DTerm"

model_name_1 = 'KMA'
model_name_2 = 'Nature'
model_name_3 = 'KMA_ME_FTn_correlated'

#model_name_1 = 'Nature'
#model_name_3 = 'Nature_Plus'

th1.name = model_name_1
th2.name = model_name_2
th3.name = model_name_3

th1.parameters.update(par_KMA)
th2.parameters.update(par_KMA)
th3.parameters.update(par_KMA)

f1 = g.MinuitFitter(GLO15b, th1)
f2 = g.MinuitFitter(GLO15b, th2)
f3 = g.MinuitFitter(GLO15b, th3)

# Load saved fit results
with open(f"Fits/{type}_{model_name_1}.json", "r") as fitres:
    fit_results_1 = json.load(fitres)
with open(f"Fits/{type}_{model_name_2}.json", "r") as fitres:
    fit_results_2 = json.load(fitres)
with open(f"Fits/{type}_{model_name_3}.json", "r") as fitres:
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
f1.fix_parameters('ALL')
f2.fix_parameters('ALL')
f3.fix_parameters('ALL')

f1.release_parameters('mv2', 'rv', 'bv', 'C', 'mC2', 'tmv2', 'trv', 'tbv', 'rpi', 'mpi2', 'ms2', 'secs', 'this', 'secg', 'thig')
f2.release_parameters('mv2', 'rv', 'bv', 'C', 'mC2', 'tmv2', 'trv', 'tbv', 'rpi', 'mpi2', 'ms2', 'secs', 'this', 'secg', 'thig')
f3.release_parameters('mv2', 'rv', 'bv', 'C', 'mC2', 'tmv2', 'trv', 'tbv', 'rpi', 'mpi2', 'ms2', 'secs', 'this', 'secg', 'thig')


plotter1 = CFFPlotter(th1)
plotter2 = CFFPlotter(th2)
plotter3 = CFFPlotter(th3)

fig, ax = plt.subplots(figsize=(8, 6))

plotter1.plot(ax=ax, label_prefix="KM15")
#plotter2.plot(ax=ax, label_prefix='CLAS data [10]', style='xx')
plotter3.plot(ax=ax, label_prefix=r'KM15-B', style='//')

# Plot Multipole_D for randomly initialized parameter values.
t_vals = np.linspace(-0.45, 0.0, 100)
t_CLAS = np.array([0.11,0.15,0.2,0.26,0.34])
D_CLAS = np.array([-1.57,-1.36,-1.28,-1.10,-0.92])
D_CLAS_err = np.array([0.1,0.08,0.08,0.07,0.07])

plt.plot(-t_vals, Multipole_D(t_vals, D0=-1.59*(4/5), M= 0.487, p=0.841), color='red', linestyle='--', label='Dispersive analysis [27]')
plt.plot(-t_vals, Multipole_D(t_vals, D0=-2.35*(4/5), M= 0.646, p=2), color='blue', linestyle=':', label='Chiral Quark Soliton Model [28]')
plt.plot(-t_vals, Multipole_D(t_vals, D0=-2.27*18/25, M=m.sqrt(1.02), p=2.76), color='black', linestyle='-.', label='CLAS data fit [4,26]')
plt.errorbar(t_CLAS, D_CLAS*4/5, yerr=D_CLAS_err*4/5, fmt='o', color='black', capsize=3, label='CLAS data [4,26]')

plt.xlabel(r'$|t|\ [\text{GeV}^{2}/c^{4}]$')
plt.ylabel(r'$D^{Q}(t)$')
plt.grid(alpha=0.3)
plt.legend(loc="upper left", fontsize=17, facecolor='white')

plt.show()
plt.savefig(f'figs/Comparison_DTerm_{type}.pdf')
