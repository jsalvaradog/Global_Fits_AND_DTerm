import numpy as np
import pandas as pd

correlations = {
    1: 0.833868,
    2: 0.918328,
    3: 0.815894,
    4: 0.806549,
    5: 0.590322,
    6: 0.481608,
    7: 0.118247,
    8: 0.857132,
    9: 0.710517,
    10: 0.441741,
    11: 0.29826,
    12: 0.144065,
    13: 0.392206,
    14: 0.159177,
    15: 0.842064,
    16: 0.750704,
    17: 0.751817,
    18: 0.685861,
    19: 0.533421,
    20: 0.484316,
    21: -0.0352973,
    22: -0.221437,
    23: 0.86927,
    24: 0.496285,
    25: 0.666892,
    26: -0.0392934,
    27: -0.00130214,
    28: 0.91969,
    29: 0.939444,
    30: -0.00749357,
    31: 0.853128,
    32: -0.00232756,
    33: 0.47191,
    34: 0.726359,
    35: 0.706929,
    36: -0.000142848,
    37: -0.000071898,
    38: -0.0000255706,
    39: -0.0000479402,
    40: -0.000175708,
    41: 0.816769,
    42: 0.358942,
    43: 0.19334,
    44: 0.258163,
    45: 0.164185,
    46: 0.0860699,
    47: 0.000105765,
    48: 0.905195,
    49: 0.718564,
    50: 0.649431,
    51: 0.529332,
    52: -0.390294,
    53: 0.90734,
    54: 0.939251,
    55: 0.78251,
    56: -0.000141515,
    57: -0.000507702,
    58: 0.864737,
    59: 0.826528,
    60: 0.396668,
    61: -0.000241236,
    62: -0.00000532797,
    63: -0.00001391,
    64: -0.0000100183,
    65: -0.00000415747,
    66: 0.0000596179,
    67: 0.714229,
    68: 0.71061,
    69: 0.348146,
    70: 0.45204,
    71: 0.0584522,
    74: 0.493089,
    75: 0.42045,
    76: -0.0000371531,
    77: -0.0000262345,
    78: 0.137744,
    79: 0.00050437,
    81: 0.0952768,
    82: 0.0765559,
    83: 0.236029,
    84: 0.293729,
    85: 0.00000360676,
    86: 0.716337,
    87: 0.430398,
    88: -0.00010565,
    89: 0.0354334,
    90: -0.0000278554,
    91: 0.664696,
    93: 0.268571,
    94: 0.411531,
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
import numpy as np

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
    "NP_Data_JSAG.dat",
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
# ------------------------------------------------------------
# Save output
# ------------------------------------------------------------
header_text = """id = 999
editor = JSAG

### Experiment

collaboration = CLAS
process = ep2epgamma
exptype = fixed target
year = 2025
reference = arXiv:9999.99999
inspiretex = Alvarado:2026saf


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

with open("BSA-CLAS12-JSAG_FTn.dat", "w") as f:
    f.write(header_text)

df.to_csv(
    "BSA-CLAS12-JSAG_FTn.dat",
    sep=" ",
    index=False,
    float_format="%.4f",
    mode="a",
    header=True
)

print("Done. Results written to BSA-CLAS12-JSAG_FTn.dat")