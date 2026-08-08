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
        
        
    def makefit(self, xibin, ntms):
        """ Perform simultaneous fit of both ImH and ReH."""
        binD = {}
        
        for tbin in range(ntms):
            binD[(xibin,tbin)] = BSD.loc[(xilims[xibin] < BSD.xi) &
                        (BSD.xi < xilims[xibin+1]) &
                        (tmlims[tbin] < BSD.tm) & (BSD.tm < tmlims[tbin+1])]            
            print(binD[(xibin,tbin)].head())
            

        totchi = 0
        ndof = 0
        CFFs = ['ReH', 'ImH']
        CFFvals = []
        CFFerrs = []
        xis = []
        pts = []
        
        pts_fit = []
        xs_vals = []
        xs_vals_err = []
        for tbin in range(ntms):
            nrows = binD[(xibin, tbin)].shape[0]
            print("nrows =", nrows)
            # Model creation
            th = LocalFit_Model()
            th.parameters.update({'ReH': 1.0, 'ImH': 1.0, 'ReE': 0.0, 'ImE': 0.0, 
                                'ReHt': 0.0, 'ImHt': 0.0, 'ReEt': 0.0, 'ImEt': 0.0})
            
            # Fit to data in given bin
            f = g.MinuitFitter(binD[(xibin,tbin)].pt.values, th)
            f.fix_parameters('ALL')

            th._release_parameters('ReH', 'ImH')
            f.release_parameters('ReH', 'ImH')
            
            xis.append(binD[(xibin,tbin)].xi.mean())
            pts.append(binD[(xibin,tbin)].pt.values[0])
                        
            f.fit()
            totchi += f.minuit.fval
            ndof += len(f.fitpoints)-len(CFFs)

            
            phi_prev=0
            bin_count=1
            # Store datapoints in pt objects            
            for idx, row in binD[(xibin,tbin)].iterrows():
                phi_val=(m.pi - row.phi)
                if(phi_val < phi_prev):
                    bin_count+=1
                pt = g.data.DataPoint(Q2=row.Q2, xB=row.xB, t=-1.0*row.tm, phi=row.phi, xi=row.xi
                                      , observable=row.observable, process="ep2epgamma", exptype="fixed target"
                                      , in1energy = 5.75, in1polarization = 1, in1charge=1)
                pt.prepare()
                pts_fit.append(pt)
                xs_vals.append(row.val)
                xs_vals_err.append(row.err)
                phi_prev=phi_val
        
            vals= []
            errs = []
            for k in range(bin_count):
                for p in CFFs:
                        vals.append(th.m.parameters[p])
                        errs.append(np.sqrt(th.m.covariance[p,p]))
                CFFvals.append(vals)
                CFFerrs.append(errs)
            print(binD[(xibin,tbin)].pt.values[0], ' (', vals[0], vals[1], '), (', errs[0], errs[1], ')')


        
        # Plot fitted datapoints and fit for XS observable
        xs_fit = []
        xs_fit_errs = []
        xs_meas = []
        xs_meas_err = []
        phi_vals = []
        phi_prev=0
        ifig=0
        for l, pt in enumerate(pts_fit):
            phi_val=(m.pi - pt.phi)
                
            if phi_val < phi_prev or l==len(pts_fit)-1:
                fig_xs, ax_xs = plt.subplots(figsize=(8, 6))
                plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)

                xs_fit = np.array(xs_fit)
                xs_fit_errs = np.array(xs_fit_errs)
                xs_meas = np.array(xs_meas)
                xs_meas_err = np.array(xs_meas_err)
                phi_vals = np.array(phi_vals)
                
                ax_xs.fill_between(phi_vals, xs_fit - xs_fit_errs, xs_fit + xs_fit_errs,
                                    alpha=0.25, label='XS Local Fit')
                ax_xs.plot(phi_vals, xs_fit, linestyle='-')
                ax_xs.errorbar(phi_vals, xs_meas, yerr=xs_meas_err, fmt='o', linestyle='None', label='XS Measurements', markersize=8)
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
                
                plt.show()

                xs_fit = []
                xs_fit_errs = []
                xs_meas = []
                xs_meas_err = []
                phi_vals = []
                
                ifig+=1
            else:
                #print("lala", CFFvals[ifig][0], CFFerrs[ifig][0])
                #print({param: CFFvals[ifig][k-1] for k, param in enumerate(CFFs)})
                print(ifig, l, '/',len(pts_fit)-1)
                th.parameters.update({param: CFFvals[ifig][k] for k, param in enumerate(CFFs)})
                th.parameters_errors = {param: CFFerrs[ifig][k] for k, param in enumerate(CFFs)}
                xs_pred = th.predict(pt, phi=pt.phi, observable=pt.observable, uncertainty=True)
                xs_fit.append(xs_pred[0])
                xs_fit_errs.append(xs_pred[1])
                xs_meas.append(xs_vals[l])
                xs_meas_err.append(xs_vals_err[l])
                print(f"{phi_val*180/np.pi:.4f}  {xs_pred[0]:.4f}  {xs_pred[1]:.4f}  {xs_vals[l]:.4f}  {xs_vals_err[l]:.4f}")
                phi_vals.append(round(phi_val*180/np.pi,4))

            phi_prev=phi_val

        return xis, pts, pts_fit, CFFvals, CFFerrs


    def PlotLocalFit(self,i,j, ax=None):
        xis, pts, pts_fit, CFFvals, CFFerrs = self.makefit(i,j)
        """Plot the results of the local fits together with NNet and KM fits"""
        xirange = np.linspace(self.xi_min, self.xi_max, 20)
        CFFs = ['ImH', 'ReH']
        ex_pt = copy.deepcopy(pts[0])

        DTerm = []
        DTerm_err = []
        tvals=[]
        l=0
        for pt in pts:
            DR = self.theory.predict(pt, observable='ReH', uncertainty=True)
            DR_H = DR[0]
            DR_H_err = DR[1]
            print(f"t={pt.t:.4f}, xi={pt.xi:.4f}; DR ReH={DR_H:.4f} +- {DR_H_err:.4f}; Local Fit H={CFFvals[l][0]:.4f} +- {CFFerrs[l][0]:.4f}, DTerm={DR_H-CFFvals[l][0]:.4f}")
            DTerm.append(DR_H-CFFvals[l][0])
            DTerm_err.append(np.sqrt(CFFerrs[l][0]**2 + DR_H_err**2))
            tvals.append(pt.t)
            l+=1

        DTerm = np.asarray(DTerm)
        tvals = np.asarray(tvals)

        def tripole(t, D0, m2):
            return D0 / (1.0 - t / m2) ** 3

        p0 = (DTerm[0], 1.0)
        bounds = ([-10, 0], [0, 10])
        try:
            popt, pcov = curve_fit(tripole, tvals, DTerm, p0=p0, bounds=bounds, maxfev=10000)
            perr = np.sqrt(np.diag(pcov))
        except Exception:
            popt = p0
            perr = (np.nan, np.nan)



        fig, ax = plt.subplots(figsize=(8, 6))
        plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)        
    
        ax.errorbar(np.abs(tvals), DTerm, yerr=DTerm_err, fmt='s', linestyle='None', label='D-term points')
        tfit = np.linspace(np.min(tvals), np.max(tvals), 200)
        ax.plot(np.abs(tfit), tripole(tfit, *popt), '-',
            label=rf'Tripole fit: $D_0={popt[0]:.3f}\pm {perr[0]:.3f}$, $m^2={popt[1]:.3f}\pm {perr[1]:.3f}$')
        print(rf'Tripole fit: $D_0={popt[0]:.3f}\pm {perr[0]:.3f}$, $m^2={popt[1]:.3f}\pm {perr[1]:.3f}$')
        ax.set_xlabel(r"$ |t| [\mathrm{GeV}^{2}/c^{4}]$", fontsize=25)
        ax.set_ylabel(r"D-term", fontsize=25)
        #ax.set_xlim(0.0, 0.55)
        #ax.set_ylim(-5., 1)
        ax.tick_params(axis='both', labelsize=20)
        ax.axhline(y=0, linewidth=0.8, linestyle=':', color='k')  # y=0 thin line
        ax.legend(fontsize=16)
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
BSD = g.select(g.dset[98], criteria=['val != 0', 'phi < 2.26892', 'phi > -2.26892']).df()  # BSS 1504.02009 CLAS data base E145M1
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
#LocalFit_Plot = LocalFit_CLAS.PlotLocalFit(0,6, ax=ax)
LocalFit_Plot = LocalFit_CLAS.PlotLocalFit(1,6, ax=ax)
#LocalFit_Plot = LocalFit_CLAS.PlotLocalFit(2,6, ax=ax)
#LocalFit_Plot = LocalFit_CLAS.PlotLocalFit(3,6, ax=ax)
#LocalFit_Plot = LocalFit_CLAS.PlotLocalFit(4,6, ax=ax)
##LocalFit_Plot = LocalFit_CLAS.PlotLocalFit(5,6, ax=ax)
plt.show()

#LocalFit_Plot = LocalFit_CLAS.PlotLocalFit(3,4, ax=ax)
#LocalFit_Plot = LocalFit_CLAS.PlotLocalFit(4,6, ax=ax)
#LocalFit_Plot = LocalFit_CLAS.PlotLocalFit(5,1, ax=ax)


#plt.savefig(f'figs/Comparison_DTerm_Local_Fit.pdf')
#plt.savefig(f'figs/Comparison_Pressure_Local_Fit.pdf')