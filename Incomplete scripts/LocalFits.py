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

class LocalFit_Model(eff.KellyEFF, gpd.PWNormGPD, cff.CFF, dvcs.BM10tw2):
    def __init__(self, **kwargs):
        # initial values of parameters and limits on their values
        self.add_parameters({'ReH': 1.0, 'ReE': 1.0, 'ReHt': 1.0, 'ReEt': 1.0,
                                'ImH': 1.0, 'ImE': 1.0, 'ImHt': 1.0, 'ImEt': 1.0,})
        super().__init__(**kwargs)

    def ImH(self, pt: data.DataPoint):
        return self.parameters['ImH']
    def ImE(self, pt: data.DataPoint):
        return self.parameters['ImE']
    def ImHt(self, pt: data.DataPoint):
        return self.parameters['ImHt']
    def ImEt(self, pt: data.DataPoint):
        return self.parameters['ImEt']
    def ReH(self, pt: data.DataPoint):
        return self.parameters['ReH']
    def ReE(self, pt: data.DataPoint):
        return self.parameters['ReE']
    def ReHt(self, pt: data.DataPoint):
        return self.parameters['ReHt']
    def ReEt(self, pt: data.DataPoint):
        return self.parameters['ReEt']

class LocalFit:
    def __init__(self, theory, t=-0.2, xi_min=0.02, xi_max=0.25, npts=100):
        """
        Plot ImH and ReH as functions of xi for fixed t and Q2.

        Parameters
        ----------
        theory : Gepard theory object
            Fitted Gepard theory (e.g. th)
        t : float
            Momentum transfer
       xi_min, xi_max : float
            xi range
        npts : int
            Number of points
        """
        self.theory = theory
        self.t = t
        self.xi = np.geomspace(xi_min, xi_max, npts)
        self.xi_min = xi_min
        self.xi_max = xi_max
        
        
    def tmin(self, xB, Q2):
        """BMK Eq. (31)"""

        eps2 = (4. * xB**2 * Mp2) / Q2
        return (-Q2 * ( 2. * (1.-xB)*(1. - np.sqrt(1.+eps2)) + eps2 ) / (
                4. * xB * (1.-xB) + eps2 ))
        
        
    def makefit(self, tbin, nxis):
        """ Perform simultaneous fit of both ImH and ReH."""
        binD = {}
        binS = {}
        binM = {}
        binDS = {}
        for xbin in range(nxis):
            binD[(tbin,xbin)] = BSD.loc[(tmlims[tbin] < BSD.tm) &
                                        (BSD.tm < tmlims[tbin+1]) &
                                (xilims[xbin] < BSD.xi) & (BSD.xi < xilims[xbin+1])]
            binS[(tbin,xbin)] = BSS.loc[(tmlims[tbin] < BSS.tm) &
                                        (BSS.tm < tmlims[tbin+1]) &
                                (xilims[xbin] < BSS.xi) & (BSS.xi < xilims[xbin+1])]
            binM[(tbin,xbin)] = BSM.loc[(tmlims[tbin] < BSM.tm) &
                                        (BSM.tm < tmlims[tbin+1]) &
                                (xilims[xbin] < BSM.xi) & (BSM.xi < xilims[xbin+1])]
            """
            print(BSD.loc[(tmlims[tbin] < BSD.tm) &
                                        (BSD.tm < tmlims[tbin+1]) &
                                (xilims[xbin] < BSD.xi) & (BSD.xi < xilims[xbin+1])].head())
            print(BSS.loc[(tmlims[tbin] < BSS.tm) &
                                        (BSS.tm < tmlims[tbin+1]) &
                                (xilims[xbin] < BSS.xi) & (BSS.xi < xilims[xbin+1])].head())
            print(BSM.loc[(tmlims[tbin] < BSM.tm) &
                                        (BSM.tm < tmlims[tbin+1]) &
                                (xilims[xbin] < BSM.xi) & (BSM.xi < xilims[xbin+1])].head())
            """

            binDS[(tbin,xbin)] = pd.concat([binD[(tbin,xbin)],
                                            binS[(tbin,xbin)]], ignore_index=True)
        totchi = 0
        ndof = 0
        CFFs = ['ReH', 'ImH']
        CFFvals = []
        CFFerrs = []
        xis = []
        pts = []
        for xbin in range(nxis):
            # Model creation
            th = LocalFit_Model()
            th.parameters.update({'ReH': 1.0, 'ImH': 1.0, 'ReE': 0.0, 'ImE': 0.0, 
                                'ReHt': 0.0, 'ImHt': 0.0, 'ReEt': 0.0, 'ImEt': 0.0})
            th._release_parameters('ReH', 'ImH')#, 'ReE', 'ImE', 'ReHt', 'ImHt', 'ReEt', 'ImEt')
            
            # Fit to data in given bin
            f = g.MinuitFitter(binDS[(tbin,xbin)].pt.values, th)
            f.fix_parameters('ALL')
            f.release_parameters('ReH', 'ImH')#, 'ReE', 'ImE', 'ReHt', 'ImHt', 'ReEt', 'ImEt')
            
            xis.append(binDS[(tbin,xbin)].xi.mean())
            pts.append(binDS[(tbin,xbin)].pt.values[0])
            
            #f.printMode = 2
            f.fit()
            f.print_parameters()
            totchi += f.minuit.fval
            ndof += len(f.fitpoints)-len(CFFs)

            vals= []
            errs = []
            for p in CFFs:
                    vals.append(th.m.parameters[p])
                    errs.append(np.sqrt(th.m.covariance[p,p]))
            CFFvals.append(vals)
            CFFerrs.append(errs)

        CFFvals = np.array(CFFvals)
        CFFerrs = np.array(CFFerrs)

        print("\n chisq/ndof = {:.1f}/{}".format(totchi,ndof))
        return xis, pts, CFFvals, CFFerrs

    def makefitstep(self, tbin, nxis):
        """ Perform step-wise fit, first ImH, then ReH ... not significant"""
        binD = {}
        binS = {}
        binM = {}
        binDS = {}
        for xbin in range(nxis):
            binD[(tbin,xbin)] = BSD.loc[(tmlims[tbin] < BSD.tm) &
                                        (BSD.tm < tmlims[tbin+1]) &
                                (xilims[xbin] < BSD.xi) & (BSD.xi < xilims[xbin+1])]
            binS[(tbin,xbin)] = BSS.loc[(tmlims[tbin] < BSS.tm) &
                                        (BSS.tm < tmlims[tbin+1]) &
                                (xilims[xbin] < BSS.xi) & (BSS.xi < xilims[xbin+1])]
            binM[(tbin,xbin)] = BSM.loc[(tmlims[tbin] < BSM.tm) &
                                        (BSM.tm < tmlims[tbin+1]) &
                                (xilims[xbin] < BSM.xi) & (BSM.xi < xilims[xbin+1])]
            """
            print(BSD.loc[(tmlims[tbin] < BSD.tm) &
                                        (BSD.tm < tmlims[tbin+1]) &
                                (xilims[xbin] < BSD.xi) & (BSD.xi < xilims[xbin+1])].head())
            print(BSS.loc[(tmlims[tbin] < BSS.tm) &
                                        (BSS.tm < tmlims[tbin+1]) &
                                (xilims[xbin] < BSS.xi) & (BSS.xi < xilims[xbin+1])].head())
            print(BSM.loc[(tmlims[tbin] < BSM.tm) &
                                        (BSM.tm < tmlims[tbin+1]) &
                                (xilims[xbin] < BSM.xi) & (BSM.xi < xilims[xbin+1])].head())
            """

            binDS[(tbin,xbin)] = pd.concat([binD[(tbin,xbin)],
                                            binS[(tbin,xbin)]], ignore_index=True)

        totchi = 0
        ndof = 0
        CFFs = ['ReH', 'ImH']
        CFFvals = []
        CFFerrs = []
        xis = []
        pts = []
        for xbin in range(nxis):
            # Model creation
            th = LocalFit_Model()
            th.parameters.update({'ReH': 1.0, 'ImH': 1.0, 'ReE': 0.0, 'ImE': 0.0, 
                                'ReHt': 0.0, 'ImHt': 0.0, 'ReEt': 0.0, 'ImEt': 0.0})

            xis.append(binM[(tbin,xbin)].xi.mean())
            pts.append(binM[(tbin,xbin)].pt.values[0])

            # Fit to BSD data in given bin
            f = g.MinuitFitter(binM[(tbin,xbin)].pt.values, th)
            f.fix_parameters('ALL')
            f.release_parameters('ImH')
            f.fit()
            totchi += f.minuit.fval
            ndof += len(f.fitpoints)-1
            valImH = th.m.parameters['ImH']
            errImH = np.sqrt(th.m.covariance['ImH','ImH'])
            f.print_parameters()
            
            # Now fit to BSS data
            f = g.MinuitFitter(binM[(tbin,xbin)].pt.values, th)
            f.fix_parameters('ALL')
            f.release_parameters('ReH')
            f.fit()
            totchi += f.minuit.fval
            ndof += len(f.fitpoints)-1
            valReH = th.m.parameters['ReH']
            errReH = np.sqrt(th.m.covariance['ReH','ReH'])
            f.print_parameters()
            
            vals= [valImH, valReH]
            errs = [errImH, errReH]
            CFFvals.append(vals)
            CFFerrs.append(errs)

        CFFvals = np.array(CFFvals)
        CFFerrs = np.array(CFFerrs)

        print("\n chisq/ndof = {:.1f}/{}".format(totchi,ndof))
        return xis, pts, CFFvals, CFFerrs


    def PlotLocalFit(self,i,j, ax=None):
        xis, pts, CFFvals, CFFerrs = self.makefit(i,j)
        """Plot the results of the local fits together with NNet and KM fits"""
        xirange = np.linspace(self.xi_min, self.xi_max, 20)
        CFFs = ['ImH', 'ReH']
        ex_pt = copy.deepcopy(pts[0])

        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 6))
        else:
            fig = ax.figure

        plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)

        ax.errorbar(xis, CFFvals[:,0], CFFerrs[:,0], linestyle='None', marker='s',
                    label='ReH LocalFit')
        ax.errorbar(xis, CFFvals[:,1], CFFerrs[:,1], linestyle='None', marker='s',
                    label='ImH LocalFit')

        ax.set_xlabel(r'$\xi$', fontsize=14)
        ax.set_xlim(0.03, 0.31)

        ax.set_xlabel(r"$\xi$", fontsize=25)
        ax.set_ylabel(r"CFF", fontsize=25)
        ax.set_ylim(-10., 40.)
        ax.set_xlim(self.xi_min, self.xi_max)
        ax.tick_params(axis='both', labelsize=20)
        ax.axhline(y=0, linewidth=0.8, linestyle=':', color='k')  # y=0 thin line
        ax.text(0.05, 1, r'$t = {}\,{{\rm GeV}}^2$'.format(self.t), fontsize=15,
                    bbox=dict(facecolor='gold', alpha=0.3))
        ax.legend(fontsize=16)
        return fig    



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

            values.append(result[0])
            errors.append(result[1])

        return np.array(values), np.array(errors)

    def plot(self, ax=None, label_prefix=""):
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
            "ReH": r"$\mathfrak{Re}[\mathcal{H}]$"
        }
        
        for cff in ["ImH", "ReH"]:
            y, dy = self.predict_cff(cff)

            label = f"{label_prefix} {cff_labels[cff]}".strip()

            ax.plot(
                self.xi,
                y
            )

            ax.fill_between(
                self.xi,
                y - dy,
                y + dy,
                alpha=0.3,
                label=label
            )

        ax.set_xlabel(r"$\xi$", fontsize=25)
        ax.set_ylabel(r"CFF", fontsize=25)
        ax.tick_params(axis='both', labelsize=20)

        ax.legend(fontsize=16)
        
        ax.grid(True)

        return fig

tmlims = [0, 0.13, 0.18, 0.22, 0.3, 0.4, 0.5]
xilims = [0, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.26]

BSA = g.dset[7].df()   # ALU(90 deg) by CLAS 2007
BSS = g.dset[102].df()  # CLAS 2015
BSD = g.dset[101].df()  # CLAS 2015

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
BSM = g.dset[999].df()  # CLAS12 small |t|

Volker = (g.dset[7] + g.dset[98] + g.dset[100]).df()

"""
xis, pts, CFFvals, CFFerrs = makefit(0,3)
fig = makefig()
plt.savefig('figs/LocalFits_tbin0.pdf', bbox_inches='tight')

xis, pts, CFFvals, CFFerrs = makefit(1,3)
fig = makefig()

xis, pts, CFFvals, CFFerrs = makefit(2,6)
fig = makefig()

xis, pts, CFFvals, CFFerrs = makefit(3,7)
fig = makefig()

xis, pts, CFFvals, CFFerrs = makefit(4,8)
fig = makefig()
"""


th1 = th_KM15
type = "Global_Fit" # "Global_Fit" or "DTerm"
model_name_1 = 'KM15'
f1 = g.MinuitFitter(GLO15b, th1)
f1.fix_parameters('ALL')
f1.release_parameters('mv2', 'rv', 'bv', 'C', 'mC2', 'tmv2', 'trv', 'tbv', 'rpi', 'mpi2', 'ms2', 'secs', 'this', 'secg', 'thig')
# Load saved fit results
with open(f"{type}_{model_name_1}.json", "r") as fitres:
    fit_results_1 = json.load(fitres)
# Restore Minuit values
for par, val in fit_results_1["values"].items():
    f1.minuit.values[par] = val
# Restore Minuit errors
for par, err in fit_results_1["errors"].items():
    f1.minuit.errors[par] = err
# Synchronize values/errors back to the Gepard theory object
f1.theory.parameters.update(f1.minuit.values.to_dict())
f1.theory.parameters_errors = f1.minuit.errors.to_dict()
# restore theory parameters
th1.parameters.update(f1.minuit.values.to_dict())
# restore errors
th1.parameters_errors = f1.minuit.errors.to_dict()

plot_KM15 = CFFPlotter(th1, t=-0.11)
LocalFit = LocalFit(th1, t=-0.11)

fig, ax = plt.subplots(figsize=(8, 6))
plot_KM15.plot(ax=ax, label_prefix="KM15")
LocalFit_Plot = LocalFit.PlotLocalFit(0,3, ax=ax)
plt.show()
plt.savefig(f'figs/Comparison_LocalFits.pdf')