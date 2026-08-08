import shelve, logging, sys, copy
logging.basicConfig(level=logging.ERROR)
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.stats import norm
from scipy.optimize import curve_fit
import importlib_resources
import matplotlib
import matplotlib.pyplot as plt
import math as m

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
        
        
    def makefit(self):
        """ Perform simultaneous fit of both ImH and ReH."""
        binD = {}
        
        binD = BSD.loc[(0 < BSD.xi) & (BSD.xi < 1) &
                        (0 < BSD.tm) & (BSD.tm < 10)]                        

        CFFvals = []
        CFFerrs = []

        xs_vals = []
        xs_vals_err = []
        xs_fit = []
        xs_fit_errs = []
        phi_vals = []

        phi_prev=0
        bin_count=1
        iprev=0

        dterm = pd.DataFrame(columns=['xB', 'Q2', 't', 'xi', 'Dterm']) #create dictionary or df with xB, Q2, t, xi, Dterm, Dterm_err, ReH, ReH_err, ImH, ImH_err, DR_H, DR_H_err
        # Model creation
        th = LocalFit_Model()
        th.parameters.update({'ReH': 1.0, 'ImH': 0.0, 'ReE': 0.0, 'ImE': 0.0, 
                            'ReHt': 0.0, 'ImHt': 0.0, 'ReEt': 0.0, 'ImEt': 0.0})
                
        CFFs = ['ReH', 'ImH']
        for iter, row in binD.iterrows():
            phi_val=(m.pi - row.phi)
            pt = g.data.DataPoint(Q2=row.Q2, xB=row.xB, t=-1.0*row.tm, phi=row.phi, xi=row.xi
                                    , observable=row.observable, process="ep2epgamma", exptype="fixed target"
                                    , in1energy = 5.75, in1charge=-1,in1polarization=0,in2polarization='U',frame='Trento')
            pt.prepare()
            if(phi_val < phi_prev):
                # select data
                XS_Data = binD.iloc[iprev:iter-1]
                print(XS_Data[['xB', 'Q2', 'tm', 'xi', 'phi', 'val', 'err']].head())
                bin_count+=1
                iprev=iter       
                
                # Fit to data in given bin
                f = g.MinuitFitter(XS_Data.pt.values, th)
                f.fix_parameters('ALL')

                th._release_parameters('ReH', 'ImH')
                f.release_parameters('ReH', 'ImH')

                # More accurate (and slower) minimization
                f.minuit.strategy = 2

                # Allow more function evaluations
                f.minuit.migrad(ncall=100000)

                # Stage 1: fit dominant CFF only
                f.minuit.fixed["ReH"] = False
                f.minuit.fixed["ImH"] = False
                f.minuit.fixed["ReE"] = True
                f.minuit.fixed["ImE"] = True
                f.minuit.fixed["ReHt"] = True
                f.minuit.fixed["ImHt"] = True
                f.minuit.fixed["ReEt"] = True
                f.minuit.fixed["ImEt"] = True

                # Set initial values (optional but often very important)
                f.minuit.values["ReH"] = 1.0
                f.minuit.values["ImH"] = 1.0
                f.minuit.values["ReE"] = 0.0
                f.minuit.values["ImE"] = 0.0
                f.minuit.values["ReHt"] = 0.0
                f.minuit.values["ImHt"] = 0.0
                f.minuit.values["ReEt"] = 0.0
                f.minuit.values["ImEt"] = 0.0

                # Set initial step sizes
                f.minuit.errors["ReH"] = 0.9
                f.minuit.errors["ImH"] = 0.
                f.minuit.errors["ReE"] = 0.
                f.minuit.errors["ImE"] = 0.
                f.minuit.errors["ReHt"] = 0.
                f.minuit.errors["ImHt"] = 0.
                f.minuit.errors["ReEt"] = 0.
                f.minuit.errors["ImEt"] = 0.

                # Run minimization
                f.minuit.migrad()
                f.minuit.hesse()
                #f.minuit.minos()

                # Sync covariance back to GEPARD
                f.covsync()
                print(f.minuit.valid)     # True if Minuit considers the minimum valid
                chi2 = f.minuit.fval
                ndf = len(f.fitpoints)-len(CFFs)
                if ndf and ndf > 0:
                    print(f"chi2/ndf = {chi2:.6f}/{ndf} = {chi2/ndf:.6f}")
                else:
                    print(f"chi2 = {chi2:.6f}, ndf = {ndf}")
                
                # Set ReH to 10
                #th.m.parameters["ReH"] = 2.0
                

                # Compute XS predictions for the local fit result 
                for s, iphi in enumerate(phi_vals):
                    pt1 = g.data.DataPoint(Q2=row.Q2, xB=row.xB, t=-1.0*row.tm, phi=m.pi -iphi*m.pi/180
                                            , observable=row.observable, process="ep2epgamma", exptype="fixed target"
                                            , in1energy = 5.75, in1charge=-1,in1polarization=0,in2polarization='U',frame='Trento')
                    pt1.prepare()
                    xs_pred = th.predict(pt1, observable=pt1.observable, uncertainty=True)
                    xs_fit.append(xs_pred[0])
                    xs_fit_errs.append(xs_pred[1])
                    print(f"{iphi:.4f}  {xs_pred[0]:.4f}  {xs_pred[1]:.4f}  {xs_vals[s]:.4f}  {xs_vals_err[s]:.4f}")

                # Store Local Fit results
                vals= []
                errs = []
                for p in CFFs:
                        vals.append(th.m.parameters[p])
                        errs.append(np.sqrt(th.m.covariance[p,p]))
                CFFvals.append(vals)
                CFFerrs.append(errs)
                print('(', vals[0], vals[1], '), (', errs[0], errs[1], ')')

                # Plot fitted datapoints and fit for XS observable
                fig_xs, ax_xs = plt.subplots(figsize=(8, 6))
                plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)

                xs_fit = np.array(xs_fit)
                xs_fit_errs = np.array(xs_fit_errs)
                xs_vals = np.array(xs_vals)
                xs_vals_err = np.array(xs_vals_err)
                phi_vals = np.array(phi_vals)
                
                ax_xs.fill_between(phi_vals, xs_fit - xs_fit_errs, xs_fit + xs_fit_errs,
                                    alpha=0.25, label='XS Local Fit')
                ax_xs.plot(phi_vals, xs_fit, linestyle='-')
                ax_xs.errorbar(phi_vals, xs_vals, yerr=xs_vals_err, fmt='o', linestyle='None', label='XS Measurements', markersize=8)
                ax_xs.set_xlabel(r"$\phi$", fontsize=25)
                ax_xs.set_ylabel(r"XS", fontsize=25)
                ax_xs.tick_params(axis='both', labelsize=20)
                ax_xs.legend(fontsize=16)
                ax_xs.grid(True, alpha=0.3)

                # Add xB and Q2 values to plot
                textstr = f'$x_B$ = {pt.xB:.3f}\n$Q^2$ = {pt.Q2:.3f} GeV$^{2}/c^{4}$\n$t$ = {pt.t:.3f} GeV$^{2}/c^{4}$'
                ax_xs.text(0.98, 0.97, textstr, transform=ax_xs.transAxes, fontsize=14,
                          verticalalignment='top', horizontalalignment='right',
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
                
                #plt.show()

                # Compute the Dterm locally
                DR = self.theory.predict(pt, observable='ReH', uncertainty=True)
                DR_H = DR[0]
                DR_H_err = DR[1]
                DTerm = DR_H - CFFvals[-1][0]
                DTerm_err = np.sqrt(DR_H_err**2 + CFFerrs[-1][0]**2)

                # Summarize the results in the DataFrame
                dterm.loc[len(dterm)] = {
                    'xB': pt.xB,
                    'Q2': pt.Q2,
                    't': pt.t,
                    'xi': pt.xi,
                    'Dterm': DR_H - CFFvals[-1][0],
                    'Dterm_err': np.sqrt(DR_H_err**2 + CFFerrs[-1][0]**2),
                    'ReH': CFFvals[-1][0],
                    'ReH_err': CFFerrs[-1][0],
                    'ImH': CFFvals[-1][1],
                    'ImH_err': CFFerrs[-1][1],
                    'DR_H': DR_H,
                    'DR_H_err': DR_H_err                    
                }
                
                xs_fit = []
                xs_fit_errs = []
                xs_vals = []
                xs_vals_err = []
                phi_vals = []
            else:
                xs_vals.append(row.val)
                xs_vals_err.append(row.err)
                phi_vals.append(round(phi_val*180/m.pi,4))
            phi_prev=phi_val

        return dterm


    def PlotLocalFit(self, ax=None):
        xilims = [0, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.26]
        Dterm_Data = self.makefit()

        # Tripole function
        def tripole(t, D0, m2):
            return D0 / (1.0 - t / m2) ** 3

        # Create figure with 8 subplots (4x2)
        fig, axes = plt.subplots(4, 2, figsize=(14, 16))
        axes = axes.flatten()
        plt.subplots_adjust(left=0.1, right=0.95, top=0.97, bottom=0.08, hspace=0.3, wspace=0.3)
        
        # Loop over each xi bin
        for i in range(len(xilims) - 1):
            xi_min = xilims[i]
            xi_max = xilims[i + 1]
            
            # Filter data for this xi bin
            mask = (Dterm_Data['xi'] >= xi_min) & (Dterm_Data['xi'] < xi_max)
            data_bin = Dterm_Data[mask]
            
            if len(data_bin) == 0:
                continue
            
            tvals = data_bin['t'].values
            DTerm = data_bin['Dterm'].values
            DTerm_err = data_bin['Dterm_err'].values
            
            ax = axes[i]
            
            # Fit tripole model
            if len(tvals) > 2:
                p0 = (DTerm[0], 1.0)
                bounds = ([-10, 0], [0, 10])
                try:
                    popt, pcov = curve_fit(tripole, tvals, DTerm, p0=p0, bounds=bounds, maxfev=10000)
                    perr = np.sqrt(np.diag(pcov))
                except Exception:
                    popt = p0
                    perr = (np.nan, np.nan)
            else:
                popt = (np.nan, np.nan)
                perr = (np.nan, np.nan)

            print(rf'Tripole fit: $D_0={popt[0]:.3f}\pm {perr[0]:.3f}$, $m^2={popt[1]:.3f}\pm {perr[1]:.3f}$')
            
            # Plot data points
            ax.errorbar(np.abs(tvals), DTerm, yerr=DTerm_err, fmt='s', linestyle='None', 
                       label='D-term points', markersize=6, capsize=4)
            
            # Plot fit curve
            if not np.isnan(popt[0]):
                tfit = np.linspace(np.min(tvals), np.max(tvals), 200)
                ax.plot(np.abs(tfit), tripole(tfit, *popt), '-', linewidth=2,
                       label=rf'Tripole: $D_0={popt[0]:.3f}\pm {perr[0]:.3f}$, $m^2={popt[1]:.3f}\pm {perr[1]:.3f}$')
            
            ax.set_xlabel(r"$|t|$ [GeV$^{2}/c^{4}$]", fontsize=12)
            ax.set_ylabel(r"D-term", fontsize=12)
            ax.set_title(rf'$\xi \in [{xi_min:.2f}, {xi_max:.2f})$', fontsize=13)
            ax.tick_params(axis='both', labelsize=10)
            ax.axhline(y=0, linewidth=0.8, linestyle=':', color='k')
            ax.legend(fontsize=10, loc='best')
            ax.grid(True, alpha=0.3)    
        return fig    

class KM_DR(KM15):
    def ReH(self, pt):
        """Real part of CFF H."""
        return g.DispersionFixedPoleCFF.ReH(self, pt, imfun=g.DispersionFixedPoleCFF.ImH)

tmlims = [0, 0.13, 0.18, 0.22, 0.3, 0.4, 0.5]
xilims = [0, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.26]

BSA = g.dset[7].df()   # ALU(90 deg) by CLAS 2007
#BSS = g.dset[102].df()  # CLAS 2015
#BSD = g.dset[101].df()  # CLAS 2015
BSD = g.select(g.dset[98], criteria=['val != 0', 'phi < 130*3.14159265359/180', 'phi > -130*3.14159265359/180']).df()  # BSS 1504.02009 CLAS data base E145M1
BSS = g.select(g.dset[100], criteria=['val != 0']).df()  # BSS 1504.02009 CLAS data base E145M1, restricted kinematics
# &(BSD.phi > 50*m.pi/180.) & (BSD.phi < 300*m.pi/180.)
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


th1 = KM_DR()
par_KM15 = {'tmv2': 15.94293227053628, 'rS': 1.0, 'alv': 0.43, 'tal': 0.43,
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
th1.parameters.update(par_KM15)

type = "Global_Fit" # "Global_Fit" or "DTerm"
model_name_1 = 'KM15'
#model_name_1 = 'KM15_ME_smallt'

f1 = g.MinuitFitter(GLO15b, th1)
f1.fix_parameters('ALL')
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
th1.parameters.update({'C': 0.0})

LocalFit_CLAS = LocalFit(th1, t=-0.11)
fig, ax = plt.subplots(figsize=(8, 6))
LocalFit_Plot = LocalFit_CLAS.PlotLocalFit(ax=ax)
plt.show()

#plt.savefig(f'figs/Comparison_DTerm_Local_Fit.pdf')
