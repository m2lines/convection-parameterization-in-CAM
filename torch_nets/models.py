"""Neural network architectures."""
import netCDF4 as nc  # type: ignore

import torch
from torch import nn
from torch.nn import functional as F


class ANN(nn.Module):  # pylint: disable=too-many-instance-attributes
    """Model used in the paper.

    Paper: https://doi.org/10.1029/2020GL091363


    Parameters
    ----------
    n_in : int
        Number of input features.
    n_out : int
        Number of output features.
    neurons : int
        The number of neurons in the hidden layers.
    dropout : float
        The dropout probability to apply in the hidden layers.

    Notes
    -----
    If you are doing inference, always remember to put the model in eval model,
    by using ``model.eval()``, so the dropout layers are turned off.

    """

    def __init__(
        self,
        n_in: int = 61,
        n_out: int = 148,
        n_layers: int = 5,
        neurons=128,
        dropout=0.0,
        device="cpu",
        features_mean=None,
        features_std=None,
        outputs_mean=None,
        outputs_std=None,
        output_groups=None
    ):
        """Build ``ANN``."""
        super().__init__()
        
        n_units = [n_in] + [neurons] * (n_layers - 1) + [n_out]
        
        self.layers = [nn.Linear(n_units[i], n_units[i + 1])
                       for i in range(n_layers)]
        self.dropout = nn.Dropout(dropout)
        self.features_mean = features_mean
        self.features_std = features_std
        self.outputs_mean = outputs_mean
        self.outputs_std = outputs_std
        
        if output_groups is not None:
            assert len(output_groups) == len(outputs_mean)
            self.outputs_mean = torch.cat([
                torch.tensor([x] * g, dtype=torch.float32)
                for x, g in zip(outputs_mean, output_groups)
            ])
            self.outputs_std = torch.cat([
                torch.tensor([x] * g, dtype=torch.float32)
                for x, g in zip(outputs_std, output_groups)
            ])

        self.to(torch.device("cpu"))


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
        
        if self.features_mean is not None:
            batch = (batch - self.features_mean) / self.features_std
        
        for layer in self.layers[:-1]:
            batch = self.dropout(F.relu(layer(batch)))

        batch = self.layers[-1](batch)
        
        if self.outputs_mean is not None:
            batch = batch * self.outputs_std + self.outputs_mean

        return batch


    def initialize(self, weights_file: str, use_pkl: bool = True):
        """Initialise model with with saved weights.

        Parameters
        ----------
        weights_file : str
            File name for saved weights
        use_pkl : bool
            Whether to read weights from a .pkl file (if `True`)
            or netCDF file (if `False`). Defaults is `True`.
        """
        import pickle
        if use_pkl:
            with open(weights_file, "rb") as f:
                weights = pickle.load(f)
            self.load_state_dict(weights)
        else:
            endow_with_netcdf_params(self, weights_file)
        return self


@torch.no_grad()
def endow_with_netcdf_params(model: Module, nc_file: str):
    """Endow the model with weights and biases in the netcdf file.

    Parameters
    ----------
    model : Module
        The model to add the parameters to.
    nc_file : str
        The netcdf file containing the parameters.

    Notes
    -----
    This function is protected with ``no_grad`` because PyTorch doesn't like
    it when we edit a layer's weights/biases with gradients enabled.

    """
    data_set = nc.Dataset(nc_file)  # pylint: disable=no-member

    for i, layer in enumerate(model.layers):
        layer.weight[:] = data_set[f"w{i+1}"]
        layer.bias[:] = data_set[f"b{i+1}"]
