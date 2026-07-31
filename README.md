# VGG DDVCS

Python routine to do a global fits of DVCS experimental data to the KM parametrization using Gepard.

**main.py perform the fits.** 
**plot.py creates a plot comparing the effect of the DVCS data at small |t|**
- Configuration must be chosen on the code itself
- Fit result stores in a json file (as it takes a few days to be done)
- This work is based on KK Matters Arising work (see notebook on github)


**Used Datasets are**
- GLO15b : Dataset used to fit the KM15 model
- pts_Volker = g.dset[7] + g.dset[98] + g.dset[100] : Dataset used for the proton-pressure extraction of [Volker](https://misportal.jlab.org/sti/publications/14962/attachments/1228/JLAB-PHY-17-2536.pdf)
	- Notice the repply from [Kresimir](https://wwwcompass.cern.ch/compass/gpd/meetings/2022/Pressure_Nature2019_Kumericki_s41586-019-1211-6.pdf) where he cannot reproduce the curve and wonders if there is a factor missing. [See Jupyter Notebook](https://github.com/openhep/dterm18)
- ME = g.dset[999] : Defined in this scope, are my DVCS BSA measurements without proton detection

**Theory objects**
- th_KM15 : KM15 model
- KMA : KM15 model but subtraction constant has p=3
- KMA_H : KM15 model but p=3 in the subtraction constant. Only GPD H


**Fitted configurations are:**
- GLO15b / KM15 -> Global_Fit_KM15.json
- GLO15b / KMA -> Global_Fit_KMA.json
- GLO15b / KMA_H -> Global_Fit_KMA_H.json

- GLO15b + ME / KM15 -> Global_Fit_ME_KM15.json
- GLO15b + ME / KMA -> Global_Fit_ME_KMA.json
- GLO15b + ME / KMA_H  -> Global_Fit_ME_KMA_H.json

- pts_Volker / KM15 -> DTerm_KM15.json
- pts_Volker / KMA -> DTerm_KMA.json
- pts_Volker / KMA_H -> DTerm_KMA_H.json

- pts_Volker + ME / KM15 -> DTerm_ME_KM15.json
- pts_Volker + ME / KMA -> DTerm_ME_KMA.json
- pts_Volker + ME / KMA_H -> DTerm_ME_KMA_H.json


