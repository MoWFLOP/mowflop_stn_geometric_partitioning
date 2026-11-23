"""Convert the entropy curve written by diagnose_entropy into the H.npy read by entropy_curve.py.

Usage:
    python paper/entropy_csv_to_npy.py entropy_curve_ns509_p100_i50.csv H.npy
"""
import sys

import numpy as np
import pandas as pd

np.save(sys.argv[2], pd.read_csv(sys.argv[1])["entropy"].to_numpy())
