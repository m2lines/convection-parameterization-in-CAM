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
# nn_ones.txt is the output of the Fortran NN model given an input of all ones.

model1 = ANN().load("nn_state.pt")  # load from the pytorch weights
model2 = load_from_netcdf_params(
    "qobsTTFFFFFTF30FFTFTF30TTFTFTFFF80FFTFTTF2699FFFF_X01_no_qp_no_adv_"
    "surf_F_Tin_qin_disteq_O_Trad_rest_Tadv_qadv_qout_qsed_RESCALED_7epochs"
    "_no_drop_REAL_NN_layers5in61out148_BN_F_te70.nc"
)  # load from the NetCDF weights of the pretrained Fortran NN model
# file created at https://github.com/yaniyuval/Neural_nework_parameterization/blob/f81f5f695297888f0bd1e0e61524590b4566bf03/NN_training/src/ml_train_nn.py#L417

x = torch.ones(61)

actual1 = model1.forward(x).detach().numpy()
actual2 = model2.forward(x).detach().numpy()

assert np.all(actual1 == actual2)
assert np.allclose(expected, actual1, atol=3e-8, rtol=2e-6)
# Values of atol and rtol are chosen to be the lowest that still pass the test.

print("Smoke tests passed")
