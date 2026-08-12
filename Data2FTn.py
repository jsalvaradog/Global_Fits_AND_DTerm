import numpy as np
import pandas as pd

# ------------------------------------------------------------
# Fourier coefficients
# ------------------------------------------------------------
def C(a, b, n):
    q = b / (1.0 + np.sqrt(1.0 - b**2))
    return (2.0 * a / b) * (-1)**(n-1) * q**n
# ------------------------------------------------------------
# Monte Carlo propagation
# ------------------------------------------------------------
def mc_uncertainty(a, b, da, db, nsamples=100000):
    """
    Propagate uncertainties on (a,b) to (C1,C2) using Monte Carlo.
    Assumes a and b are independent Gaussian variables.
    Returns sigma(C1), sigma(C2).
    """
    rng = np.random.default_rng()

    a_mc = rng.normal(a, da, nsamples)
    b_mc = rng.normal(b, db, nsamples)

    # Keep only physical values |b| < 1
    mask = np.abs(b_mc) < 1.0
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
    s1, s2 = mc_uncertainty(a, b, row["a_stat"], row["b_stat"])
    C1_stat.append(s1)
    C2_stat.append(s2)

    # Systematic propagation
    s1, s2 = mc_uncertainty(a, b, row["a_sys"], row["b_sys"])
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