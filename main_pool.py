import shelve, logging, sys, copy
logging.basicConfig(level=logging.ERROR)
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.stats import norm
import importlib_resources
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, "./gepard")

import gepard as g
from gepard.fits import GLO15b
from gepard import cff, data, dvcs, eff, gpd

import os
import json
from multiprocessing import Pool
import atexit


if not os.path.exists('figs'):
    os.makedirs('figs')


class CFFPlotter:
    def __init__(self, theory, t=-0.2, Q2=10.0, xi_min=0.02, xi_max=0.25, npts=100):
        """
        Plot ImH and ReH as functions of xi for fixed t and Q2.

        Parameters
        ----------
        theory : Gepard theory object
            Fitted Gepard theory (e.g. th)
        t : float
            Momentum transfer
        Q2 : float
            Virtuality
        xi_min, xi_max : float
            xi range
        npts : int
            Number of points
        """
        self.theory = theory
        self.t = t
        self.Q2 = Q2
        self.xi = np.geomspace(xi_min, xi_max, npts)

    def predict_cff(self, observable):
        """
        Predict CFF values and uncertainties.
        """
        values = []
        errors = []

        for xi in self.xi:

            # Gepard DataPoint expects xB, not xi.
            # Convert xi -> xB: xi = xB/(2-xB)
            xB = 2 * xi / (1 + xi)

            pt = g.DataPoint(
                xB=xB,
                t=self.t,
                Q2=self.Q2
            )

            result = self.theory.predict(
                pt,
                observable=observable,
                uncertainty=True
            )

            values.append(pt.xi*result[0])
            errors.append(pt.xi*result[1])

        return np.array(values), np.array(errors)

    def plot(self):
        """
        Produce ImH and ReH plots with uncertainty bands.
        """

        fig, ax = plt.subplots(figsize=(8, 6))

        for cff, color in zip(["ImH", "ReH"], ["C0", "C1"]):

            y, dy = self.predict_cff(cff)

            ax.plot(
                self.xi,
                y,
                label=cff
            )

            ax.fill_between(
                self.xi,
                y - dy,
                y + dy,
                alpha=0.3
            )

        #ax.set_xscale("log")
        ax.set_xlabel(r"$\xi$")
        ax.set_ylabel(r"$\xi*$CFF")

        ax.legend()
        ax.grid(True)

        return fig
        
class KMA(g.fits.KM15):  # noqa: D101, E501
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
# datasets
#########################
# + g.dset[98]
pts_Volker = g.select(g.dset[7] + g.dset[100], criteria=['val != 0'])
print('\nVOLKER DATAPOINTS')
g.describe_data(pts_Volker)
#pts_Volker

print('\nKM15 DATAPOINTS')
g.describe_data(GLO15b)
GLO15b

# My measurements without proton
data={}
data_file_path = "./BSA-CLAS12-JSAG.dat"
with open(data_file_path, 'r', encoding='utf-8') as f:
    file_content = f.read()
dataset = g.DataSet(datafile=file_content)
for pt in dataset:
    pt.to_conventions()
data[dataset.id] = dataset
g.dset.update(data)
print('\nJSAG')
g.describe_data(g.dset[999])
print('\nJSAG Small |t|')
pts_Me = g.select(g.dset[999], criteria=['tm < 0.12', 'val != 0'])
g.describe_data(pts_Me)

#########################
# Create model and fit
#########################

#th = KMA_H()
th = KMA()
#th = KM15();

_worker_theory=th

model_name = 'KMA_ME_smallt_pool'
th.name = model_name

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
th.parameters.update(par_KMA)
#th._release_parameters('mv2', 'rv', 'bv', 'C', 'mC2', 'tmv2', 'trv', 'tbv', 'rpi', 'mpi2', 'ms2', 'secs', 'this', 'secg', 'thig')

# To save the fitted model
#f = g.MinuitFitter(pts_Volker + pts_Me, th)
f = g.MinuitFitter(GLO15b + pts_Me, th)
f.fix_parameters('ALL')
#f.release_parameters('rv','bS','bv', 'mC2','C')
f.release_parameters('mv2', 'rv', 'bv', 'C', 'mC2', 'tmv2', 'trv', 'tbv', 'rpi', 'mpi2', 'ms2', 'secs', 'this', 'secg', 'thig')
if(True):
	f.fit()

	fit_results = {
	    "values": f.minuit.values.to_dict(),
	    "errors": f.minuit.errors.to_dict(),
	    "fval": f.minuit.fval,
	    "valid": f.minuit.valid,
	}

	# Save to a JSON file
	with open(f"Global_Fit_{model_name}.json", "w") as fitres:
	    json.dump(fit_results, fitres, indent=4)
else:
	# Load saved fit results
	with open(f"Global_Fit_{model_name}.json", "r") as fitres:
	    fit_results = json.load(fitres)

	# Restore Minuit values
	for par, val in fit_results["values"].items():
	    f.minuit.values[par] = val

	# Restore Minuit errors
	for par, err in fit_results["errors"].items():
	    f.minuit.errors[par] = err

	# Synchronize values/errors back to the Gepard theory object
	f.theory.parameters.update(f.minuit.values.to_dict())
	f.theory.parameters_errors = f.minuit.errors.to_dict()

	# restore theory parameters
	th.parameters.update(f.minuit.values.to_dict())

	# restore errors
	th.parameters_errors = f.minuit.errors.to_dict()

f.print_parameters()

plotter = CFFPlotter(
    th,
    t=-0.2,
    Q2=10
)

fig = plotter.plot()
plt.show()
plt.savefig(f'figs/{model_name}.pdf')

#th.predict(pts[0], uncertainty=True)
#pt = g.DataPoint(xB=0.01, t=-0.2, Q2=10)
#th.predict(pt, observable='ImH', uncertainty=True)
