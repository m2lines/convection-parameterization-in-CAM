"""A smoke test for the ANN model.

This test checks that the model can be loaded from a weights file in both pt format and
netcdf format and that they produce the expected output when given an input of all ones.
This ensures that it is equivalent to the Fortran NN model.
"""

import os
from pathlib import Path
import torch
import numpy as np
from models import ANN, load_from_netcdf_params


os.chdir(Path(__file__).parent)

expected = np.loadtxt("nn_ones.txt").astype(np.float32)

model1 = ANN().load("nn_state.pt")  # load from pt file
model2 = load_from_netcdf_params(
    "qobsTTFFFFFTF30FFTFTF30TTFTFTFFF80FFTFTTF2699FFFF_X01_no_qp_no_adv_"
    "surf_F_Tin_qin_disteq_O_Trad_rest_Tadv_qadv_qout_qsed_RESCALED_7epochs"
    "_no_drop_REAL_NN_layers5in61out148_BN_F_te70.nc"
)  # load from netcdf file

x = torch.ones(61)

actual1 = model1.forward(x).detach().numpy()
actual2 = model2.forward(x).detach().numpy()

assert np.all(actual1 == actual2)
assert np.allclose(expected, actual1, atol=1e-6, rtol=2e-6)
