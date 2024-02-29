"""A smoke test for the ANN model.

This test checks that the model can be loaded and that it produces
the expected output when given an input of all ones.
It ensures that it is equivalent to the Fortran NN model.
"""

import torch
import numpy as np
from models import ANN

expected = np.loadtxt("nn_ones.txt")

model = ANN().load("nn_state.pt")

x = torch.ones(61)
actual = super(ANN, model).forward(x).detach().numpy()

assert np.allclose(expected, actual)
