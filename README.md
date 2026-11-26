# VGG DDVCS

Python routine to do a global fits of DVCS experimental data to the KM parametrization using Gepard.

**main.py perform the fits.** 

**main_pool.py perform the fits but minuit minimization is parallelized with 4 cores. Requires the modified version of gepard included in this repository** 

**CFF_vs_t.py creates a plot comparing the effect of the DVCS data at small |t| on CFF H**

**CFF_vs_xi.py creates a plot comparing the effect of the DVCS data at small |t| on CFF H**

**DTerm.py creates a plot comparing different D-term predictions**

**Pressure.py creates a plot comparing the effect of the DVCS data at small |t| on the Pressure distribution**

**Data2FTn.py Transforms my BSA measurements as a function of phi, to the first two harmonics**

**Data2FTn_2023.py Transforms 2023 BSA measurements as a function of phi, to the first two harmonics**

**CFF_vs_t_KM15_xcheck.py creates a plot comparing KM15 and the re-fit of this work**

**CFF_vs_t_KM15_xcheck_2.py creates a plot comparing KM15, KM15-A and KM15-B**

**Get_Pull.py computes the pull of each dataset**

**LocalFits.py (unverified) Perform CFF local fits and Dterm extraction**

- Configuration must be chosen on the code itself
- Fit result stores in a json file (as it takes a few days to be done)
- This work is based on KK Matters Arising work (see notebook on github)


**Used Datasets are**
- GLO15b : Dataset used to fit the KM15 model
- pts_Volker = g.dset[7] + g.dset[98] + g.dset[100] : Dataset used for the proton-pressure extraction of [Volker](https://misportal.jlab.org/sti/publications/14962/attachments/1228/JLAB-PHY-17-2536.pdf)
	- Notice the reply from [Kresimir](https://wwwcompass.cern.ch/compass/gpd/meetings/2022/Pressure_Nature2019_Kumericki_s41586-019-1211-6.pdf) where he cannot reproduce the curve and wonders if there is a factor missing. [See Jupyter Notebook](https://github.com/openhep/dterm18)
- ME = g.dset[999] : Defined in this scope, are my DVCS BSA measurements without proton detection
- Defurne = g.dset[998] : Defined in this scope, are the first CLAS12 BSA measurements (2023) (Christiaens)
- HallA =  g.dset[143] + g.dset[144] : Defined in this scope, are the 2022 HallA measurements (Georges)

- FTn means that I used the dataset after doing a FT analysis (extraction of first two harmonics)

**Theory objects**
- th_KM15 : KM15 model
- KMA : KM15 model but subtraction constant has p=3


**Fitted configurations are:**
- GLO15b / KM15 -> Global_Fit_KM15.json
- GLO15b + ME / KM15 -> Global_Fit_KM15_ME_FTn.json
- GLO15b + ME@|t|<0.12 / KM15 -> Global_Fit_KM15_ME_FTn_smallt.json

- GLO15b / KMA -> Global_Fit_KMA.json
- GLO15b + ME / KMA -> Global_Fit_KMA_ME_FTn.json
- GLO15b + ME@|t|<0.12 GeV² / KMA -> Global_Fit_KMA_ME_FTn_smallt.json

- GLO15b + Defurne / KMA -> Global_Fit_KMA_Defurne_FTn.json
- GLO15b + Defurne@|t|<0.8 GeV² / KMA -> Global_Fit_KMA_Defurne_FTn_tlt0p8.json

- GLO15b + HallA + Defurne / KMA -> Global_Fit_KMA_HallA_Defurne_FTn.json
- GLO15b + HallA + ME / KMA -> Global_Fit_KMA_HallA_ME_FTn.json

**Usefull links**

- (KM15 fit) https://arxiv.org/pdf/1512.09014
- (KM15 parametrization, page 59) https://arxiv.org/pdf/0904.0458
- (pressure formula Eq. 35b) https://arxiv.org/pdf/1805.06596
- (Grips QCD workshop) https://www.jlab.org/Hall-B/shifts/csc/slides/2018/burkert-Primosten2018-2018-09-27.pdf
https://github.com/openhep/dterm18/blob/master/dterm.ipynb

