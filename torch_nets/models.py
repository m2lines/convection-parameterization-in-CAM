"""Neural network architectures."""
import netCDF4 as nc

from torch import as_tensor, no_grad  # pylint: disable=no-name-in-module

from torch.nn import Module, Linear, Dropout
from torch.nn import functional as F


class ANN(Module):  # pylint: disable=too-many-instance-attributes
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

    def __init__(
        self,
        n_in: int = 61,
        n_out: int = 148,
        neurons=128,
        dropout=0.0,
    ):
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


@no_grad()
def endow_with_netcdf_params(model: Module, nc_file: str):
    """Endow the model with weights and biases in the netcdf file.

    Parameters
    ----------
    model : Module
        The model to add the parameters to.
    nc_file : str
        The netcdf file containing the parameters.

    """
    data_set = nc.Dataset(nc_file)  # pylint: disable=no-member

    for name, layer in model.named_children():
        if not isinstance(layer, Linear):
            continue

        layer_num = int("".join(filter(lambda x: x.isdigit(), name)))

        weight = as_tensor(
            data_set.variables[f"w{layer_num}"][:],
            dtype=layer.weight.dtype,
            device=layer.weight.device,
        )

        bias = as_tensor(
            data_set.variables[f"b{layer_num}"][:],
            dtype=layer.bias.dtype,
            device=layer.bias.device,
        )

        layer.weight[:] = weight[:]
        layer.bias[:] = bias[:]


if __name__ == "__main__":
    model = ANN(61, 148)

    nc_file = "qobsTTFFFFFTF30FFTFTF30TTFTFTFFF80FFTFTTF2699FFFF_X01_no_qp_no_adv_surf_F_Tin_qin_disteq_O_Trad_rest_Tadv_qadv_qout_qsed_RESCALED_7epochs_no_drop_REAL_NN_layers5in61out148_BN_F_te70.nc"

    endow_with_netcdf_params(model, nc_file)
