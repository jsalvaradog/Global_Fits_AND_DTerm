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

if not os.path.exists('Fits'):
    os.makedirs('Fits')
        
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


def compute_chi2_and_pull(th, data, dataset_id):
    chi2 = 0.0
    pull = 0.0
    npoints = 0
    for pt in data:
        prediction = th.predict(pt) if hasattr(th, "predict") else th(pt)
        error = getattr(pt, "err", getattr(pt, "error", None))
        if error is None or not np.isfinite(error) or error <= 0:
            continue
        chi2 += ((pt.val - prediction) / error) ** 2
        pull += (pt.val - prediction) / error
        npoints += 1

    ndf = npoints
    print(f"Dataset {dataset_id} chi2/N: {chi2:.6g}/{npoints} = {chi2 / ndf:.6g}")
    #print(f"Dataset {dataset_id} pull: {pull:.6g}/sqrt({npoints}) = {pull / np.sqrt(npoints):.6g}\n")



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
# My measurements without proton
data={}
data_file_path = "./BSA-CLAS12-JSAG_FTn.dat"
with open(data_file_path, 'r', encoding='utf-8') as f:
    file_content = f.read()
dataset = g.DataSet(datafile=file_content)
for pt in dataset:
    pt.to_conventions()
data[dataset.id] = dataset
g.dset.update(data)

# Maxime measurements
data={}
data_file_path = "./BSA-CLAS12-Defurne_FTn.dat"
with open(data_file_path, 'r', encoding='utf-8') as f:
    file_content = f.read()
dataset = g.DataSet(datafile=file_content)
for pt in dataset:
    pt.to_conventions()
data[dataset.id] = dataset
g.dset.update(data)
pts_Maxime = g.select(g.dset[998], criteria=['tm < 0.8', 'val != 0'])

#########################
# Create model and fit
#########################

th = KMA()
#th = KM15()
#model_name = 'KMA_HallA_Defurne_FTn'
#model_name = 'KMA_HallA_ME_FTn'
#model_name = 'KMA_Defurne_FTn_correlated'
#model_name = 'KMA_ME_FTn_correlated'
#model_name = 'KMA_ME_FTn_smallt_correlated'
#model_name = 'KMA_Maxime_FTn_correlated_tlt0p8'
model_name = 'KMA'

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
th._release_parameters('mv2', 'rv', 'bv', 'C', 'mC2', 'tmv2', 'trv', 'tbv', 'rpi', 'mpi2', 'ms2', 'secs', 'this', 'secg', 'thig')

# Load saved fit results
with open(f"Fits/Global_Fit_{model_name}.json", "r") as fitres:
    fit_results = json.load(fitres)

# Restore the fitted theory parameters.
th.parameters.update(fit_results["values"])
th.parameters_errors = fit_results.get("errors", {})

g.describe_data(GLO15b)
g.describe_data(g.dset[999])
g.describe_data(g.select(g.dset[999], criteria=['tm < 0.12', 'val != 0']))
g.describe_data(pts_Maxime)
g.describe_data(g.dset[143] + g.dset[144])

# Compute chi2/ndf and pull for the predictions on dataset 
print("\nChi2/ndf:")
data_ids = [67, 94, 81, 68, 66, 96, 53, 54, 95, 52, 39, 63, 45, 46, 47, 101, 117, 102, 116, 143, 144, 999]
for id in data_ids:
    compute_chi2_and_pull(th, g.dset[id], id)

compute_chi2_and_pull(th, pts_Maxime, 998)

