"""Neural network architectures."""
import netCDF4 as nc

from torch.nn import Module, Linear, Dropout
from torch.nn import functional as F


class ANN(Module): #pylint: disable=too-many-instance-attributes
    """This seems to be the model used in the paper.

    Paper: https://doi.org/10.1029/2020GL091363


    Parameters
    ----------
    n_in : int
        Number of input features.
    n_out : int
        Number of output features.
    neurons : int
        The number of neurons in the bhidden layers.
    dropout : float
        The dropout probability to apply in the hidden layers.

    """

    def __init__(self, n_in, n_out, neurons=128, dropout=0.0):
        """Build ``ANN``."""
        super().__init__()
        self.linear1 = Linear(n_in, neurons)
        self.linear2 = Linear(neurons, neurons)
        self.linear3 = Linear(neurons, neurons)
        self.linear4 = Linear(neurons, neurons)
        self.linear5 = Linear(neurons, n_out)

        self.lin_drop1 = Dropout(dropout)
        self.lin_drop2 = Dropout(dropout)
        self.lin_drop3 = Dropout(dropout)
        self.lin_drop4 = Dropout(dropout)

    def _load_params_from_netcdf(self, nc_file: str):
        """Load weights and biases from ``nc_file``.

        Parameters
        ----------
        nc_file : str
            Path to the file containing the weights.

        Returns
        -------
        nc.Dataset
            A netcdf dataset containing the weights.

        """
        return nc.Dataset(nc_file)  # pylint: disable=no-member

    def forward(self, batch):
        """Pass ``batch`` through the model.

        Parameters
        ----------
        batch : Tensor
            A mini-batch of inputs.

        Returns
        -------
        Tensor
            The result of passing ``batch`` through the model.

        """
        batch = F.relu(self.linear1(batch))
        batch = self.lin_drop1(batch)

        batch = F.relu(self.linear2(batch))
        batch = self.lin_drop2(batch)

        batch = F.relu(self.linear3(batch))
        batch = self.lin_drop3(batch)

        batch = F.relu(self.linear4(batch))
        batch = self.lin_drop4(batch)

        batch = self.linear5(batch)

        return batch
