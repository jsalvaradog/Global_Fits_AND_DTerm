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
    def __init__(self, theory, theory_DR, theory_MB, t=-0.2, xi_min=0.02, xi_max=0.25, npts=100):
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
        self.theory_DR = theory_DR
        self.theory_MB = theory_MB
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
        binDS = {}
        for xbin in range(nxis):
            binDS[(tbin,xbin)] = BSM.loc[(tmlims[tbin] < BSM.tm) &
                                        (BSM.tm < tmlims[tbin+1]) &
                                (xilims[xbin] < BSM.xi) & (BSM.xi < xilims[xbin+1]) &
                                (1.0 < BSM.Q2) & (BSM.Q2 < 1.6)] # [1, 1.6, 2.4, 3.25]
                                #(1.6 < BSM.Q2) & (BSM.Q2 < 2.4)] # [1, 1.6, 2.4, 3.25]
                                #(2.4 < BSM.Q2) & (BSM.Q2 < 3.24)] # [1, 1.6, 2.4, 3.25]
            #print(binDS[(tbin,xbin)][['tm', 'xi', 'xB', 'FTn']].head())
        CFFs = ['ReH', 'ImH']
        CFFvals = []
        CFFerrs = []
        xis = []
        pts = []
        for xbin in range(nxis):
            # Gepard DataPoint expects xB, not xi.
            pt_fix = g.DataPoint(
                xB=binDS[(tbin,xbin)].xB.mean(),
                t=-1.0*binDS[(tbin,xbin)].tm.mean(),
                Q2=binDS[(tbin,xbin)].Q2.mean()
            )
            print("Q2= {:.3f}, tm={:.3f}, xi={:.3f}, xB={:.3f}".format(binDS[(tbin,xbin)].Q2.mean(), binDS[(tbin,xbin)].tm.mean(), binDS[(tbin,xbin)].xi.mean(), binDS[(tbin,xbin)].xB.mean()))
            th_int_ImH = self.theory_DR.predict(
                pt_fix,
                observable="ReH",
                uncertainty=True
            )

            th_ReH_MB = self.theory_MB.predict(
                pt_fix,
                observable="ReH",
                uncertainty=True
            )

            th_ImH = self.theory.predict(
                pt_fix,
                observable="ImH",
                uncertainty=True
            )

            th_ImE = self.theory.predict(
                pt_fix,
                observable="ImE",
                uncertainty=True
            )

            th_ImHt = self.theory.predict(
                pt_fix,
                observable="ImHt",
                uncertainty=True
            )

            th_ImEt = self.theory.predict(
                pt_fix,
                observable="ImEt",
                uncertainty=True
            )

            # Model creation
            th = LocalFit_Model()
            th.parameters.update({'ReH': th_ReH_MB[0], 'ImH': th_ImH[0], 'ReE': 0.0, 'ImE': th_ImE[0], 
                                'ReHt': 0.0, 'ImHt': th_ImHt[0], 'ReEt': 0.0, 'ImEt': th_ImEt[0]})
            th._release_parameters(*CFFs)#, 'ReE', 'ImE', 'ReHt', 'ImHt', 'ReEt', 'ImEt')
            
            # Fit to data in given bin
            f = g.MinuitFitter(binDS[(tbin,xbin)].pt.values, th)
            f.fix_parameters('ALL')
            f.release_parameters(*CFFs)#, 'ReE', 'ImE', 'ReHt', 'ImHt', 'ReEt', 'ImEt')
            
            xis.append(binDS[(tbin,xbin)].xi.mean())
            pts.append(binDS[(tbin,xbin)].pt.values[0])
            
            #f.printMode = 2
            f.fit()
            print("Valid minimum? ", f.minuit.valid, ", chi2=", f.minuit.fval)
            f.print_parameters()
            
            vals= []
            errs = []
            for p in CFFs:
                vals.append(th.m.parameters[p])
                errs.append(np.sqrt(th.m.covariance[p,p]))
                if 'ImH' not in CFFs:
                    vals.append(th_ImH[0])
                    errs.append(th_ImH[1])
            CFFvals.append(vals)
            CFFerrs.append(errs)

            print("ImH = {:.3f} +- {:.3f}, ImE = {:.3f} +- {:.3f}, ImHt = {:.3f} +- {:.3f}, ImEt = {:.3f} +- {:.3f}".format(th_ImH[0], th_ImH[1], th_ImE[0], th_ImE[1], th_ImHt[0], th_ImHt[1], th_ImEt[0], th_ImEt[1]))
            print("int ImH = {:.3f}".format(th_int_ImH[0]))
            print("ReH (MB) = {:.3f}".format(th_ReH_MB[0]))
            print(" DTerm = {:.3f} +- {:.3f}".format(th.m.parameters['ReH']-th_ReH_MB[0]-th_int_ImH[0], np.sqrt(th.m.covariance['ReH','ReH'])))
            #print("DTerm = {:.3f} +- {:.3f}".format(th.m.parameters['ReH']-th_int_ImH[0], np.sqrt(th.m.covariance['ReH','ReH'])))
            tm_vals.append(binDS[(tbin,xbin)].tm.mean())
            xi_vals.append(binDS[(tbin,xbin)].xi.mean())
            Dterm_vals.append(th.m.parameters['ReH']-th_ReH_MB[0]-th_int_ImH[0])
            Dterm_errs.append(np.sqrt(th.m.covariance['ReH','ReH']))
            

        CFFvals = np.array(CFFvals)
        CFFerrs = np.array(CFFerrs)

        print("\n")
        return xis, pts, CFFvals, CFFerrs


    def PlotLocalFit(self,i,j, ax=None):
        xis, pts, CFFvals, CFFerrs = self.makefit(i,j)
        """Plot the results of the local fits together with NNet and KM fits"""
        xirange = np.linspace(self.xi_min, self.xi_max, 20)
        ex_pt = copy.deepcopy(pts[0])

        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 6))
        else:
            fig = ax.figure

        plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)

        ax.errorbar(xis, CFFvals[:,0], CFFerrs[:,0], linestyle='None', marker='s')
                    #,label='ReH LocalFit')
        ax.errorbar(xis, CFFvals[:,1], CFFerrs[:,1], linestyle='None', marker='s')
                    #,label='ImH LocalFit')

        ax.set_xlabel(r'$\xi$', fontsize=14)
        ax.set_xlim(0.03, 0.31)

        ax.set_xlabel(r"$\xi$", fontsize=25)
        ax.set_ylabel(r"CFF", fontsize=25)
        ax.set_ylim(-10., 40.)
        ax.set_xlim(self.xi_min, self.xi_max)
        ax.tick_params(axis='both', labelsize=20)
        ax.axhline(y=0, linewidth=0.8, linestyle=':', color='k')  # y=0 thin line
        #ax.text(0.05, 1, r'$t = {}\,{{\rm GeV}}^2$'.format(self.t), fontsize=15,
        #            bbox=dict(facecolor='gold', alpha=0.3))
        ax.legend(fontsize=16)
        return fig    



class CFFPlotter_vs_xi:
    def __init__(self, theory, t=-0.2, Q2=1.35, xi_min=0.02, xi_max=0.25, npts=100):
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
    
class CFFPlotter_vs_t:
    def __init__(self, theory, Q2=1.35, tm_min=0.001, tm_max=0.4, npts=100):
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
        tm_min, tm_max : float
            t range
        npts : int
            Number of points
        """
        self.theory = theory
        self.Q2 = Q2
        self.tm = np.geomspace(tm_min, tm_max, npts)

    @staticmethod
    def Dterm(t, D0=-2.726203077588369*(18./25.), M=np.sqrt(1.6372939157165587), p=3):
        """D-term parametrization"""
        return D0/(1.-t/M)**p

    def plot(self, ax=None):
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
        
        d0=-2.726203077588369*(18./25.)
        M=np.sqrt(1.6372939157165587)
        p=3
        
        d0err = 0.4775*18/25
        M2err = 0.573
        perr = 0.5 #np.sqrt(0.25**2 + 0.5**2)
        Merr = M2err/(2*M)

        d0MC = norm.rvs(d0, d0err, size=ss)
        MMC = norm.rvs(M, Merr, size=ss)
        pMC = norm.rvs(p, perr, size=ss)
        
        psmean = []
        psup = []
        psdown = []
        for t in self.tm:
            dtermMC = []
            for k in range(ss):
                dtermMC.append(self.Dterm(-1.0*t, D0=d0MC[k], M=MMC[k], p=pMC[k]))
            mean, std = np.mean(dtermMC), np.std(dtermMC)
            psmean.append(mean)
            psup.append(mean+std)
            psdown.append(mean-std)
            
        line = ax.plot(self.tm, psmean)
        color = line[0].get_color()
        
        ax.fill_between(self.tm, psup, psdown,edgecolor=color, linewidth=1, label="Global Fit", alpha=0.3)
        ax.set_ylabel(r'$D^{Q}(t)$', fontsize=25)
        ax.set_xlabel(r'$|t| [\text{GeV}^{2}/c^{4}]$', fontsize=25)
        ax.tick_params(axis='both', labelsize=20)
        
        ax.axhline(linestyle='--', color='gray', lw=1)
        #ax.set_xlim(0,0.45)
        #ax.set_ylim(-2.4, 1.0)
        #ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.5))
        #ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.1))    
        #ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.02))    
        return fig
        
class KMA(KM15):  # noqa: D101, E501
    def subtraction(self, pt):
        #Dispersion relations subtraction constant.
        return self.parameters['C']/(1.-pt.t/self.parameters['mC2'])**3
        #return 0
    #def ImH(self, pt, xi=0):
    #    """Imaginary part of CFF H."""
    #    return g.DispersionFixedPoleCFF.ImH(self, pt, xi)

class KM_DR(KM15):
    def subtraction(self, pt):
        #Dispersion relations subtraction constant.
        # set to zero so ReH is directly the integral
        return 0
    #def ImH(self, pt, xi=0):
    #    """Imaginary part of CFF H."""
    #    return g.DispersionFixedPoleCFF.ImH(self, pt, xi)
    def ReH(self, pt):
        """Real part of CFF H."""
        # removed MB component for the same reason
        return g.DispersionFixedPoleCFF.ReH(self, pt, imfun=g.DispersionFixedPoleCFF.ImH)

class KM_MB(KM15):
    def subtraction(self, pt):
        #Dispersion relations subtraction constant.
        # set to zero so ReH is directly the integral
        return 0
    def ImH(self, pt, xi=0):
        """Imaginary part of CFF H."""
        return g.MellinBarnesCFF.ImH(self, pt)
    def ReH(self, pt):
        """Real part of CFF H."""
        # removed MB component for the same reason
        return g.MellinBarnesCFF.ReH(self, pt)


tmlims = [0, 0.04, 0.07, 0.12, 0.17, 0.25, 0.40, 1.0]
xilims = [0, 0.06, 0.09, 0.16, 0.29, 0.43]
tm_vals = []
xi_vals = []
Dterm_vals = []
Dterm_errs = []


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
BSM = g.select(g.dset[999], criteria=['val != 0']).df()
g.describe_data(g.select(g.dset[999], criteria=['val != 0']))

th1 = KMA()
th2 = KM_DR()
th3 = KM_MB()

type = "Global_Fit" # "Global_Fit" or "DTerm"
model_name_1 = 'KMA_ME_FTn_correlated'
model_name_2 = 'KMA_ME_FTn_correlated'
model_name_3 = 'KMA_ME_FTn_correlated'

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
# Restore Minuit values
for par, val in fit_results_1["values"].items():
    f1.minuit.values[par] = val
    f2.minuit.values[par] = val
    f3.minuit.values[par] = val
# Restore Minuit errors
for par, err in fit_results_1["errors"].items():
    f1.minuit.errors[par] = err
    f2.minuit.errors[par] = err
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

plot_KM15 = CFFPlotter_vs_xi(th1, t=-0.11)
LocalFit = LocalFit(th1, th2, th3, t=-0.11)

fig, ax = plt.subplots(figsize=(8, 6))
plot_KM15.plot(ax=ax, label_prefix="KM15")
LocalFit_Plot = LocalFit.PlotLocalFit(0,3, ax=ax)
LocalFit_Plot = LocalFit.PlotLocalFit(1,3, ax=ax)
LocalFit_Plot = LocalFit.PlotLocalFit(2,4, ax=ax)
LocalFit_Plot = LocalFit.PlotLocalFit(3,4, ax=ax)
LocalFit_Plot = LocalFit.PlotLocalFit(4,4, ax=ax)
LocalFit_Plot = LocalFit.PlotLocalFit(5,4, ax=ax)
#LocalFit_Plot = LocalFit.PlotLocalFit(6,5, ax=ax)

    
plt.show()
plt.savefig(f'figs/Comparison_LocalFits.pdf')

print("\n")
for i in range(len(tm_vals)):
    print("tm = {:.3f}, xi = {:.3f}, DTerm = {:.3f} +- {:.3f}".format(tm_vals[i], xi_vals[i], Dterm_vals[i], Dterm_errs[i]))

fig, ax = plt.subplots(figsize=(8, 6))
plot_KM15_2 = CFFPlotter_vs_t(th1)
plot_KM15_2.plot(ax=ax)

plt.errorbar(tm_vals, Dterm_vals, yerr=Dterm_errs, fmt='o', label='Local Fits')

# Plot Multipole_D for randomly initialized parameter values.
t_vals = np.linspace(-0.45, 0.0, 100)
t_CLAS = np.array([0.11,0.15,0.2,0.26,0.34])
D_CLAS = np.array([-1.57,-1.36,-1.28,-1.10,-0.92])
D_CLAS_err = np.array([0.1,0.08,0.08,0.07,0.07])
plt.errorbar(t_CLAS, D_CLAS*4/5, yerr=D_CLAS_err*4/5, fmt='o', color='black', capsize=3, label='CLAS data')

plt.legend(fontsize=16)
plt.show()
plt.savefig(f'figs/Comparison_LocalFits_vs_t.pdf')
