"""A smoke test for the ANN model.

This test checks that the model can be loaded and that it produces
the expected output when given an input of all ones.
It ensures that it is equivalent to the Fortran NN model.
"""

import torch
import numpy as np
from models import ANN, load_from_netcdf_params


expected = np.loadtxt("nn_ones.txt")

models = [
    ANN().load("nn_state.pt"),
    load_from_netcdf_params(
        "qobsTTFFFFFTF30FFTFTF30TTFTFTFFF80FFTFTTF2699FFFF_X01_no_qp_no_adv_"
        "surf_F_Tin_qin_disteq_O_Trad_rest_Tadv_qadv_qout_qsed_RESCALED_7epochs"
        "_no_drop_REAL_NN_layers5in61out148_BN_F_te70.nc"
    )
]

x = torch.ones(61)

for model in models:
    actual = model.forward(x).detach().numpy()
    assert np.allclose(expected, actual, atol=1e-6, rtol=1e-5)
