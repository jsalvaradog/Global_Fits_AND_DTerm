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

t_NN = np.linspace(0.01, 1.01, 100)

ReH_NN = np.array([
    -1.67, -1.60, -1.55, -1.50, -1.46, -1.42, -1.38, -1.35, -1.33, -1.31,
    -1.29, -1.27, -1.25, -1.24, -1.23, -1.22, -1.22, -1.21, -1.21, -1.20,
    -1.20, -1.20, -1.20, -1.19, -1.19, -1.19, -1.19, -1.19, -1.19, -1.19,
    -1.18, -1.17, -1.17, -1.16, -1.15, -1.14, -1.13, -1.13, -1.13, -1.14,
    -1.16, -1.18, -1.22, -1.27, -1.33, -1.41, -1.51, -1.63, -1.77, -1.94,
    -2.13, -2.34, -2.60, -2.88, -3.21, -3.59, -4.01, -4.48, -5.01, -5.60,
    -6.24, -6.94, -7.70, -8.51, -9.37, -10.28, -11.24, -12.23, -13.26, -14.33,
    -15.41, -16.52, -17.65, -18.77, -19.90, -21.01, -22.10, -23.16, -24.18, -25.16,
    -26.08, -26.94, -27.74, -28.47, -29.13, -29.73, -30.26, -30.72, -31.11, -31.44,
    -31.71, -31.91, -32.06, -32.16, -32.20, -32.20, -32.15, -32.08, -31.97, -31.84
])

ImH_NN = np.array([
    7.78, 7.40, 7.09, 6.83, 6.62, 6.42, 6.25, 6.08, 5.93, 5.78,
    5.64, 5.50, 5.37, 5.25, 5.13, 5.01, 4.91, 4.81, 4.72, 4.64,
    4.56, 4.49, 4.42, 4.36, 4.30, 4.24, 4.18, 4.13, 4.07, 4.02,
    3.96, 3.91, 3.86, 3.81, 3.75, 3.70, 3.65, 3.60, 3.55, 3.50,
    3.45, 3.40, 3.36, 3.31, 3.26, 3.21, 3.17, 3.12, 3.07, 3.02,
    2.97, 2.91, 2.85, 2.79, 2.72, 2.65, 2.57, 2.48, 2.38, 2.27,
    2.15, 2.02, 1.88, 1.74, 1.60, 1.46, 1.36, 1.29, 1.29, 1.37,
    1.55, 1.84, 2.25, 2.77, 3.39, 4.09, 4.84, 5.61, 6.37, 7.09,
    7.74, 8.31, 8.77, 9.14, 9.41, 9.59, 9.72, 9.81, 9.89, 9.98,
    10.09, 10.23, 10.42, 10.64, 10.88, 11.15, 11.44, 11.73, 12.02,
    12.30
])

ReH_NN_err = np.array([
    4.92, 3.98, 3.20, 2.57, 2.07, 1.68, 1.39, 1.19, 1.05, 0.96,
    0.90, 0.85, 0.81, 0.79, 0.77, 0.76, 0.77, 0.78, 0.80, 0.82,
    0.84, 0.86, 0.88, 0.89, 0.89, 0.90, 0.90, 0.91, 0.92, 0.95,
    0.98, 1.03, 1.09, 1.16, 1.24, 1.32, 1.40, 1.47, 1.55, 1.62,
    1.69, 1.76, 1.83, 1.91, 2.02, 2.15, 2.32, 2.54, 2.81, 3.14,
    3.52, 3.97, 4.49, 5.09, 5.80, 6.61, 7.55, 8.62, 9.84, 11.19,
    12.67, 14.28, 16.00, 17.80, 19.68, 21.60, 23.57, 25.56, 27.57, 29.59,
    31.62, 33.64, 35.67, 37.70, 39.72, 41.73, 43.74, 45.74, 47.72, 49.68,
    51.62, 53.51, 55.38, 57.20, 58.99, 60.74, 62.44, 64.11, 65.72, 67.27,
    68.74, 70.13, 71.44, 72.66, 73.80, 74.87, 75.88, 76.84, 77.78, 78.69
])

ImH_NN_err = np.array([
    14.12, 11.40, 9.10, 7.17, 5.59, 4.30, 3.29, 2.52, 1.96, 1.59,
    1.37, 1.24, 1.17, 1.13, 1.11, 1.11, 1.11, 1.13, 1.15, 1.18,
    1.22, 1.25, 1.28, 1.32, 1.36, 1.40, 1.44, 1.48, 1.53, 1.58,
    1.64, 1.69, 1.75, 1.81, 1.88, 1.94, 2.01, 2.07, 2.14, 2.22,
    2.30, 2.38, 2.47, 2.57, 2.67, 2.80, 2.94, 3.10, 3.29, 3.52,
    3.81, 4.15, 4.58, 5.10, 5.76, 6.58, 7.61, 8.88, 10.46, 12.40,
    14.76, 17.62, 21.03, 25.05, 29.69, 34.98, 40.88, 47.33, 54.26, 61.56,
    69.16, 76.96, 84.91, 92.93, 100.97, 109.00, 116.95, 124.77, 132.44, 139.91,
    147.18, 154.25, 161.13, 167.83, 174.37, 180.74, 186.95, 192.98, 198.81, 204.42,
    209.81, 214.95, 219.86, 224.55, 229.04, 233.35, 237.51, 241.54, 245.47, 249.31
])

ReHt_NN = np.array([
    2.90, 2.69, 2.51, 2.33, 2.18, 2.05, 1.92, 1.82, 1.72, 1.63,
    1.55, 1.47, 1.40, 1.32, 1.25, 1.17, 1.09, 1.02, 0.94, 0.86,
    0.78, 0.71, 0.64, 0.58, 0.52, 0.46, 0.41, 0.36, 0.32, 0.29,
    0.25, 0.22, 0.20, 0.18, 0.17, 0.17, 0.18, 0.19, 0.21, 0.24,
    0.28, 0.32, 0.37, 0.42, 0.48, 0.55, 0.62, 0.70, 0.79, 0.88,
    0.99, 1.10, 1.22, 1.35, 1.48, 1.62, 1.76, 1.90, 2.04, 2.18,
    2.31, 2.43, 2.53, 2.63, 2.71, 2.77, 2.83, 2.87, 2.92, 2.96,
    3.02, 3.09, 3.18, 3.29, 3.42, 3.57, 3.74, 3.93, 4.12, 4.33,
    4.53, 4.73, 4.92, 5.11, 5.28, 5.44, 5.59, 5.73, 5.86, 5.99,
    6.10, 6.21, 6.32, 6.42, 6.51, 6.59, 6.67, 6.75, 6.81, 6.87
])

ImHt_NN = np.array([
    2.00, 1.97, 1.95, 1.93, 1.92, 1.91, 1.91, 1.92, 1.93, 1.95,
    1.98, 2.00, 2.03, 2.06, 2.08, 2.11, 2.13, 2.15, 2.16, 2.16,
    2.17, 2.16, 2.16, 2.15, 2.13, 2.12, 2.10, 2.08, 2.06, 2.03,
    2.00, 1.98, 1.94, 1.91, 1.87, 1.83, 1.79, 1.74, 1.68, 1.62,
    1.56, 1.49, 1.42, 1.35, 1.27, 1.18, 1.10, 1.01, 0.92, 0.83,
    0.73, 0.62, 0.51, 0.40, 0.28, 0.15, 0.01, -0.14, -0.30, -0.47,
    -0.66, -0.86, -1.06, -1.26, -1.47, -1.66, -1.85, -2.03, -2.20, -2.36,
    -2.52, -2.68, -2.85, -3.04, -3.25, -3.48, -3.76, -4.07, -4.43, -4.84,
    -5.30, -5.82, -6.38, -6.99, -7.64, -8.32, -9.02, -9.73, -10.45, -11.16,
    -11.86, -12.54, -13.21, -13.86, -14.49, -15.09, -15.67, -16.23, -16.76, -17.27
])

ReHt_NN_err = np.array([
    6.94, 5.85, 4.92, 4.14, 3.49, 2.96, 2.54, 2.20, 1.94, 1.74,
    1.60, 1.50, 1.45, 1.43, 1.42, 1.43, 1.44, 1.45, 1.45, 1.45,
    1.45, 1.44, 1.44, 1.43, 1.42, 1.42, 1.42, 1.42, 1.42, 1.43,
    1.44, 1.46, 1.48, 1.51, 1.54, 1.57, 1.60, 1.63, 1.67, 1.71,
    1.78, 1.86, 1.98, 2.13, 2.32, 2.54, 2.81, 3.11, 3.46, 3.85,
    4.31, 4.82, 5.40, 6.03, 6.73, 7.48, 8.27, 9.12, 10.01, 10.96,
    11.98, 13.06, 14.22, 15.46, 16.77, 18.16, 19.60, 21.10, 22.62, 24.16,
    25.70, 27.24, 28.76, 30.26, 31.75, 33.22, 34.66, 36.08, 37.48, 38.84,
    40.17, 41.44, 42.67, 43.83, 44.92, 45.94, 46.89, 47.77, 48.59, 49.36,
    50.07, 50.76, 51.41, 52.05, 52.68, 53.30, 53.93, 54.56, 55.19, 55.82
])

ImHt_NN_err = np.array([
    6.68, 5.71, 4.91, 4.27, 3.78, 3.42, 3.14, 2.92, 2.74, 2.56,
    2.38, 2.19, 2.00, 1.82, 1.67, 1.56, 1.50, 1.49, 1.52, 1.57,
    1.63, 1.67, 1.69, 1.70, 1.69, 1.66, 1.63, 1.59, 1.55, 1.53,
    1.50, 1.49, 1.49, 1.49, 1.50, 1.51, 1.52, 1.53, 1.54, 1.56,
    1.59, 1.63, 1.70, 1.79, 1.90, 2.04, 2.22, 2.43, 2.67, 2.95,
    3.29, 3.68, 4.14, 4.71, 5.41, 6.27, 7.33, 8.63, 10.20, 12.04,
    14.18, 16.58, 19.22, 22.04, 24.99, 28.00, 31.02, 34.00, 36.92, 39.77,
    42.56, 45.28, 47.96, 50.61, 53.25, 55.86, 58.46, 61.02, 63.53, 65.97,
    68.32, 70.58, 72.74, 74.80, 76.77, 78.66, 80.48, 82.21, 83.88, 85.47,
    86.99, 88.45, 89.85, 91.19, 92.47, 93.70, 94.87, 96.00, 97.08, 98.12
])


def NN_CFF_Prediction(cffs):
    if "ReH" in cffs:
        plt.plot(t_NN, ReH_NN*xi_val, '--', color='black')
        plt.fill_between(
            t_NN,
            (ReH_NN - ReH_NN_err)*xi_val,
            (ReH_NN + ReH_NN_err)*xi_val,
            color='black',
            alpha=0.1, 
            hatch='.',
            label=r'NN $\mathfrak{Re}[\mathcal{H}]$ [24]'
        )
        plt.ylim(-0.5, 0.75)

    if "ImH" in cffs:
        plt.plot(t_NN, ImH_NN*xi_val, '--', color='red')
        plt.fill_between(
            t_NN,
            (ImH_NN - ImH_NN_err)*xi_val,
            (ImH_NN + ImH_NN_err)*xi_val,
            color='red',
            alpha=0.1, 
            hatch='.',
            label=r'NN $\mathfrak{Im}[\mathcal{H}]$ [24]'
        )
        plt.ylim(-0.5, 1.5)

    if "ReHt" in cffs:
        plt.plot(t_NN, ReHt_NN*xi_val, '--', color='black')
        plt.fill_between(
            t_NN,
            (ReHt_NN - ReHt_NN_err)*xi_val,
            (ReHt_NN + ReHt_NN_err)*xi_val,
            color='black',
            alpha=0.1, 
            hatch='.',
            label=r'NN $\mathfrak{Re}[\widetilde{\mathcal{H}}]$ [24]'
        )
        plt.ylim(-0.25, 1.0)

    if "ImHt" in cffs:
        plt.plot(t_NN, ImHt_NN*xi_val, '--', color='red')
        plt.fill_between(
            t_NN,
            (ImHt_NN - ImHt_NN_err)*xi_val,
            (ImHt_NN + ImHt_NN_err)*xi_val,
            color='red',
            alpha=0.1, 
            hatch='.',
            label=r'NN $\mathfrak{Im}[\widetilde{\mathcal{H}}]$ [24]'
        )
        plt.ylim(-0.25, 1.0)


class CFFPlotter:
    def __init__(self, theory, xi=0.1, Q2=2.0, mt_min=0.02, mt_max=1.0, npts=100):
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
        mt_min, mt_max : float
            |t| range
        npts : int
            Number of points
        """
        self.theory = theory
        self.xi = xi
        self.Q2 = Q2
        self.mt = np.geomspace(mt_min, mt_max, npts)

    def predict_cff(self, observable):
        """
        Predict CFF values and uncertainties.
        """
        values = []
        errors = []

        for mt in self.mt:
            # Gepard DataPoint expects xB, not xi.
            # Convert xi -> xB: xi = xB/(2-xB)
            xi=self.xi
            xB = 2 * xi / (1 + xi)

            pt = g.DataPoint(
                xB=xB,
                t=-1.0*mt,
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

    def plot(self, ax=None, label_prefix="", style='', linestyle='-', cffs=["ImH", "ReH"]):
        """
        Produce ImH and ReH plots with uncertainty bands.

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

        cff_labels = {
            "ImH": r"$\mathfrak{Im}[\mathcal{H}]$",
            "ReH": r"$\mathfrak{Re}[\mathcal{H}]$",
            "ImE": r"$\mathfrak{Im}[\mathcal{E}]$",
            "ReE": r"$\mathfrak{Re}[\mathcal{E}]$",
            "ImHt": r"$\mathfrak{Im}[\mathcal{\widetilde{H}}]$",
            "ReHt": r"$\mathfrak{Re}[\mathcal{\widetilde{H}}]$",
            "ImEt": r"$\mathfrak{Im}[\mathcal{\widetilde{E}}]$",
            "ReEt": r"$\mathfrak{Re}[\mathcal{\widetilde{E}}]$"
        }
        cff_title = {
            "ImH": r"$\mathcal{H}$",
            "ReH": r"$\mathcal{H}$",
            "ImE": r"$\mathcal{E}$",
            "ReE": r"$\mathcal{E}$",
            "ImHt": r"$\mathcal{\widetilde{H}}$",
            "ReHt": r"$\mathcal{\widetilde{H}}$",
            "ImEt": r"$\mathcal{\widetilde{E}}$",
            "ReEt": r"$\mathcal{\widetilde{E}}$"
        }
        
        for cff in cffs:
            y, dy = self.predict_cff(cff)
            print(f"Max uncertainty for {self.theory.name} {cff}: {np.max(dy):.4f}, Min uncertainty: {np.min(dy):.4f}")
            label = f"{label_prefix} {cff_labels[cff]}".strip()

            line = ax.plot(
                self.mt,
                y,
                linestyle=linestyle
            )
            color = line[0].get_color()
            edgecolor_value = 'black' if label_prefix == "KM15 data," else color

            ax.fill_between(
                self.mt,
                y - dy,
                y + dy,
                hatch=style,
                alpha=0.4,
                label=label,
                edgecolor=color
            )

            ax.set_xlabel(r"$|t| [\text{GeV}^2/c^4]$", fontsize=25)
            ax.set_ylabel(rf"$\xi${cff_title[cff]}$(\xi,t)$", fontsize=25)
        ax.tick_params(axis='both', labelsize=20)        
        ax.grid(True)

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
model_name_2 = 'KMA_ME_FTn_smallt_correlated'
model_name_3 = 'KMA_ME_FTn_correlated'

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
xi_val=0.1

plotter1 = CFFPlotter(th1, xi=xi_val)
plotter2 = CFFPlotter(th2, xi=xi_val)
plotter3 = CFFPlotter(th3, xi=xi_val)


# Normal comparison
for cff in ["ReH", "ImH", "ReHt", "ImHt"]:    
    cffs_choice=[cff]
    fig, ax = plt.subplots(figsize=(8, 6))

    plotter1.plot(ax=ax, label_prefix='KM15,', style='', linestyle='-.', cffs=cffs_choice)
    #plotter2.plot(ax=ax, label_prefix=r'KM15-A', style='', cffs=cffs_choice)
    plotter3.plot(ax=ax, label_prefix=r'KM15-B,', style='//', cffs=cffs_choice)
    NN_CFF_Prediction(cffs_choice)

    plt.legend(loc="upper left",fontsize=15)
    plt.xlim(0.02, 1.0)
    plt.show()
    plt.savefig(f'figs/CFF_vs_t/Comparison_{type}_vs_t_{cff}.pdf')


"""
# Small t comparison
for cff in ["ReH", "ImH", "ReHt", "ImHt"]:    
    cffs_choice=[cff]
    fig, ax = plt.subplots(figsize=(8, 6))

    plotter1.plot(ax=ax, label_prefix='KM15,', style='\\\\', linestyle=':', cffs=cffs_choice)
    plotter2.plot(ax=ax, label_prefix=r'KM15-A', style='', linestyle='-.', cffs=cffs_choice)
    plotter3.plot(ax=ax, label_prefix=r'KM15-B,', style='//', cffs=cffs_choice)
    #NN_CFF_Prediction(cffs_choice)

    plt.legend(loc="upper right",fontsize=15)
    plt.xlim(0.02, 1.0)
    plt.show()
    plt.savefig(f'figs/Comparison_smallt/Comparison_{type}_vs_t_smallt_comp_{cff}.pdf')
"""