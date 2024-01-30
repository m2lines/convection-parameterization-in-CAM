"""Neural network architectures."""
from typing import Any

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
    n_layers : int
        Number of layers.
    neurons : int
        The number of neurons in the hidden layers.
    dropout : float
        The dropout probability to apply in the hidden layers.
    device : str
        The device to put the model on.
    features_mean : ndarray
        The mean of the input features.
    features_std : ndarray
        The standard deviation of the input features.
    outputs_mean : ndarray
        The mean of the output features.
    outputs_std : ndarray
        The standard deviation of the output features.
    output_groups : ndarray
        The number of output features in each group of the ouput.

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
        neurons: int = 128,
        dropout: float = 0.0,
        device: str = "cpu",
        features_mean: Any = None,
        features_std: Any = None,
        outputs_mean: Any = None,
        outputs_std: Any = None,
        output_groups: Any = None,
    ):
        """Build ``ANN``."""
        super().__init__()

        n_units = [n_in] + [neurons] * (n_layers - 1) + [n_out]

        self.layers = [nn.Linear(n_units[i], n_units[i + 1]) for i in range(n_layers)]
        self.dropout = dropout
        self.features_mean = features_mean
        self.features_std = features_std
        self.outputs_mean = outputs_mean
        self.outputs_std = outputs_std

        if features_mean is not None:
            assert features_std is not None
            assert features_mean.shape == features_std.shape
            self.features_mean = torch.tensor(features_mean, dtype=torch.float32)
            self.features_std = torch.tensor(features_std, dtype=torch.float32)

        if outputs_mean is not None:
            assert outputs_std is not None
            assert outputs_mean.shape == outputs_std.shape
            if output_groups is None:
                self.outputs_mean = torch.tensor(outputs_mean, dtype=torch.float32)
                self.outputs_std = torch.tensor(outputs_std, dtype=torch.float32)
            else:
                assert len(output_groups) == len(outputs_mean)
                self.outputs_mean = torch.cat(
                    [
                        torch.tensor([x] * g, dtype=torch.float32)
                        for x, g in zip(outputs_mean, output_groups)
                    ]
                )
                self.outputs_std = torch.cat(
                    [
                        torch.tensor([x] * g, dtype=torch.float32)
                        for x, g in zip(outputs_std, output_groups)
                    ]
                )

        self.to(torch.device(device))

    def forward(self, batch: torch.Tensor):
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
            batch = F.relu(layer(batch))
            batch = F.dropout(batch, p=self.dropout, training=self.training)

        batch = self.layers[-1](batch)

        if self.outputs_mean is not None:
            batch = batch * self.outputs_std + self.outputs_mean

        return batch


@torch.no_grad()
def endow_with_netcdf_params(model: nn.Module, nc_file: str):
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
