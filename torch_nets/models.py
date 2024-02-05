"""Neural network architectures."""

from typing import Optional

import netCDF4 as nc  # type: ignore
import torch
from torch import nn, Tensor


class ANN(nn.Sequential):
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

    def __init__(  # pylint: disable=too-many-arguments,too-many-locals
        self,
        n_in: int = 61,
        n_out: int = 148,
        n_layers: int = 5,
        neurons: int = 128,
        dropout: float = 0.0,
        device: str = "cpu",
        features_mean: Optional[Tensor] = None,
        features_std: Optional[Tensor] = None,
        outputs_mean: Optional[Tensor] = None,
        outputs_std: Optional[Tensor] = None,
        output_groups: Optional[Tensor] = None,
    ):
        """Initialize the ANN model."""
        dims = [n_in] + [neurons] * (n_layers - 1) + [n_out]
        layers = []

        for i in range(n_layers):
            layers.append(nn.Linear(dims[i], dims[i + 1]))
            if i < n_layers - 1:
                layers.append(nn.ReLU())
                layers.append(nn.Dropout(dropout))

        super().__init__(*layers)

        if features_mean is not None:
            assert features_std is not None
            assert len(features_mean) == len(features_std)
            fmean = torch.tensor(features_mean)
            fstd = torch.tensor(features_std)

        if outputs_mean is not None:
            assert outputs_std is not None
            assert len(outputs_mean) == len(outputs_std)
            if output_groups is None:
                omean = torch.tensor(outputs_mean)
                ostd = torch.tensor(outputs_std)
            else:
                assert len(output_groups) == len(outputs_mean)
                omean = torch.tensor(
                    [x for x, g in zip(outputs_mean, output_groups) for _ in range(g)]
                )
                ostd = torch.tensor(
                    [x for x, g in zip(outputs_std, output_groups) for _ in range(g)]
                )

        self.register_buffer("features_mean", fmean)
        self.register_buffer("features_std", fstd)
        self.register_buffer("outputs_mean", omean)
        self.register_buffer("outputs_std", ostd)

        self.to(torch.device(device))

    def forward(self, input: Tensor):  # pylint: disable=redefined-builtin
        """Pass the input through the model.

        Parameters
        ----------
        input : Tensor
            A mini-batch of inputs.

        Returns
        -------
        Tensor
            The model output.

        """
        if self.features_mean is not None:
            input = (input - self.features_mean) / self.features_std

        # pass the input through the layers
        output = super().forward(input)

        if self.outputs_mean is not None:
            output = output * self.outputs_std + self.outputs_mean

        return output

    def load(self, path: str):
        """Load the model from a checkpoint.

        Parameters
        ----------
        path : str
            The path to the checkpoint.

        """
        state = torch.load(path)
        for key in ["features_mean", "features_std", "outputs_mean", "outputs_std"]:
            if key in state and getattr(self, key) is None:
                setattr(self, key, state[key])
        self.load_state_dict(state)

    def save(self, path: str):
        """Save the model to a checkpoint.

        Parameters
        ----------
        path : str
            The path to save the checkpoint to.

        """
        torch.save(self.state_dict(), path)


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

    for i, layer in enumerate(l for l in model.modules() if isinstance(l, nn.Linear)):
        layer.weight.data = torch.tensor(data_set[f"w{i+1}"][:])
        layer.bias.data = torch.tensor(data_set[f"b{i+1}"][:])

    model.outputs_mean = torch.tensor(data_set["oscale_mean"][:])
    model.outputs_std = torch.tensor(data_set["oscale_stnd"][:])
    model.features_mean = torch.tensor(data_set["fscale_mean"][:])
    model.features_std = torch.tensor(data_set["fscale_stnd"][:])
