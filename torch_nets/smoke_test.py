import torch
import numpy as np
from models import ANN

expected = np.loadtxt("nn_ones.txt")

model = ANN().load("nn_state.pt")

x = torch.ones(61)
actual = super(ANN, model).forward(x).detach().numpy()

assert np.allclose(expected, actual)
