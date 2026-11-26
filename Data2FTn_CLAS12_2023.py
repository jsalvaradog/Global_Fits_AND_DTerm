import numpy as np
import pandas as pd

correlations = {
    1: 0.81408,
    2: 0.480194,
    3: 0.864481,
    4: 0.941704,
    5: 0.720709,
    6: 0.876433,
    7: 0.904457,
    8: 0.824154,
    9: -0.904354,
    10: 0.000679844,
    11: -0.0477064,
    13: 0.717769,
    14: 0.428317,
    15: -0.8934,
    16: 0.919677,
    17: 0.38172,
    18: 0.51825,
    19: -0.949111,
    20: -0.717797,
    21: 0.884613,
    22: -0.000176444,
    23: 0.874923,
    24: 0.764778,
    26: 0.360701,
    27: 0.915343,
    28: 0.870463,
    29: 0.425838,
    30: 0.960717,
    31: 0.767984,
    32: 0.907246,
    33: -0.794821,
    34: 0.891706,
    35: -0.0718554,
    36: -0.77373,
    37: 0.51143,
    38: 0.682444,
    39: -0.851594,
    40: 0.771712,
    41: 0.592385,
    42: 0,
    43: 0.937812,
    44: 0.647621,
    45: 0.958253,
    46: 0.727142,
    47: 0.858664,
    48: 0.68195,
    49: 0.266034,
    50: 0.830827,
    51: 0.290027,
    52: 0.737991,
    53: 0.892358,
    54: -0.56761,
    55: -0.828446,
    56: 0.885185,
    57: 0.275814,
    59: 0.239413,
    61: 0.470089,
    63: 0.932679,
    64: 0.683938,
}

# ------------------------------------------------------------
# Fourier coefficients
# ------------------------------------------------------------
def C(a, b, n):
    q = b / (1.0 + np.sqrt(1.0 - b**2))
    return (2.0 * a / b) * (-1)**(n-1) * q**n
# ------------------------------------------------------------
# Monte Carlo propagation
# ------------------------------------------------------------
def mc_uncertainty(a, b, da, db, corr, nsamples=100000):
    """
    Propagate uncertainties on (a, b) to (C1, C2) using Monte Carlo.
    a and b are sampled from a correlated 2D Gaussian.

    Parameters
    ----------
    a, b : float
        Central values.
    da, db : float
        Standard deviations of a and b.
    corr : float
        Correlation coefficient between a and b, in [-1, 1].
    nsamples : int
        Number of Monte Carlo samples.

    Returns
    -------
    sigma_C1, sigma_C2 : float
        Standard deviations of C1 and C2.
    """
    rng = np.random.default_rng()

    cov = np.array([
        [da**2, corr * da * db],
        [corr * da * db, db**2]
    ])

    samples = rng.multivariate_normal(
        mean=[a, b],
        cov=cov,
        size=nsamples
    )

    a_mc = samples[:, 0]
    b_mc = samples[:, 1]

    # Keep only physical values |a| < 1 and |b| < 1
    mask = (np.abs(a_mc) < 1.0) & (np.abs(b_mc) < 1.0)
    a_mc = a_mc[mask]
    b_mc = b_mc[mask]

    c1 = C(a_mc, b_mc, 1)
    c2 = C(a_mc, b_mc, 2)

    return np.std(c1, ddof=1), np.std(c2, ddof=1)

# ------------------------------------------------------------
# Read input file
# ------------------------------------------------------------
# Columns:
# Bin Q2 xB t a stat sys b stat sys
df = pd.read_csv(
    "NP_Data_Defurne.dat",
    sep=r"\s+",
    comment="#",
    header=0
)
df["t"] = -df["t"]

# Rename duplicate column names
df.columns = [
    "Bin", "Q2", "xB", "t",
    "a", "a_stat", "a_sys",
    "b", "b_stat", "b_sys"
]

kin_df_1 = df[["Q2", "xB", "t"]].copy()
kin_df_2 = df[["Q2", "xB", "t"]].copy()
# ------------------------------------------------------------
# Compute coefficients and uncertainties
# ------------------------------------------------------------
C1_vals = []
C2_vals = []
C1_stat = []
C2_stat = []
C1_sys = []
C2_sys = []

for _, row in df.iterrows():
    a = row["a"]
    b = row["b"]
    
    if(b==0 or b==-1): #b==-1 are non convergent fits on a2. There are only two points
        C1_vals.append(a)
        C2_vals.append(0.0)
        # Statistical propagation
        C1_stat.append(row["a_stat"])
        C2_stat.append(0.0)

        # Systematic propagation
        C1_sys.append(row["a_sys"])
        C2_sys.append(0.0)
    else:
        C1_vals.append(C(a, b, 1))
        C2_vals.append(C(a, b, 2))

        # Statistical propagation
        s1, s2 = mc_uncertainty(a, b, row["a_stat"], row["b_stat"], correlations[row["Bin"]])
        C1_stat.append(s1)
        C2_stat.append(s2)

        # Systematic propagation
        s1, s2 = mc_uncertainty(a, b, row["a_sys"], row["b_sys"], correlations[row["Bin"]])
        C1_sys.append(s1)
        C2_sys.append(s2)

# ------------------------------------------------------------
# Append results
# ------------------------------------------------------------

kin_df_1["FTn"] = -1
kin_df_1["Cn"] = C1_vals
kin_df_1["Cn_stat"] = C1_stat
kin_df_1["Cn_sys"] = C1_sys


kin_df_2["FTn"] = -2
kin_df_2["Cn"] = C2_vals
kin_df_2["Cn_stat"] = C2_stat
kin_df_2["Cn_sys"] = C2_sys

df = pd.concat([kin_df_1, kin_df_2], ignore_index=True, axis=0)
# remove points with FTn=-2 and Cn=0, which are non-convergent fits. Where fitted to only a1
df = df[~((df["FTn"] == -2) & (df["Cn"] == 0))]

# ------------------------------------------------------------
# Save output
# ------------------------------------------------------------
header_text = """id = 998
editor = JSAG

### Experiment

collaboration = CLAS
process = ep2epgamma
exptype = fixed target
year = 2023
reference = Phys. Rev. Lett. 130, 211902
inspiretex = 10.1103/PhysRevLett.130.211902


### Scattering Process

frame = Trento

in1particle = em
in1energy = 10.6
in1energyunit = GeV
in1polarizationvector = L
in1polarization = +1

in2particle = p

out1particle = em
out2particle = p
out3particle = gamma


### Observable 
# A_{LU}


y1name = ALU
y1namelong = Beam Spin Asymmetry
y1unit = 1
y1value = column5
y1errorstatistic = column6 
y1errorsystematic = column7 


### x-axes

x1name = Q2
x1unit = GeV^2
x1value = column1

x2name = xB
x2unit = 1
x2value = column2

x3name = tm
x3unit = GeV^2
x3value = column3

x4name = FTn
x4unit = 1
x4value = column4

### Data
table=
#
#  CLAS 2026 ALU FTn
#Q2  xB  tm  FTn  Cn  Cn_stat  Cn_sys
#########################################
"""

with open("BSA-CLAS12-Defurne_FTn.dat", "w") as f:
    f.write(header_text)

df.to_csv(
    "BSA-CLAS12-Defurne_FTn.dat",
    sep=" ",
    index=False,
    float_format="%.4f",
    mode="a",
    header=True
)

print("Done. Results written to BSA-CLAS12-Defurne_FTn.dat")
