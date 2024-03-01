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
        output_groups: Optional[list] = None,
    ):
        """Initialize the ANN model."""
        dims = [n_in] + [neurons] * (n_layers - 1) + [n_out]
        layers = []

        for i in range(n_layers):
            layers.append(nn.Linear(dims[i], dims[i + 1]))
            if i < n_layers - 1:
                layers.append(nn.ReLU())  # type: ignore
                layers.append(nn.Dropout(dropout))  # type: ignore

        super().__init__(*layers)

        fmean = fstd = omean = ostd = None

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
        
        Override the forward method of nn.Sequential to add normalization
        to the input and denormalization to the output.

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

        # pass the input through the layers using nn.Sequential.forward
        output = super().forward(input)

        if self.outputs_mean is not None:
            output = output * self.outputs_std + self.outputs_mean

        return output

    def load(self, path: str) -> "ANN":
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
        return self

    def save(self, path: str):
        """Save the model to a checkpoint.

        Parameters
        ----------
        path : str
            The path to save the checkpoint to.

        """
        torch.save(self.state_dict(), path)


def load_from_netcdf_params(nc_file: str, dtype: str = "float32") -> ANN:
    """Load the model with weights and biases from the netcdf file.

    Parameters
    ----------
    nc_file : str
        The netcdf file containing the parameters.
    dtype : str
        The data type to cast the parameters to.

    """
    data_set = nc.Dataset(nc_file)  # pylint: disable=no-member

    model = ANN(
        features_mean=data_set["fscale_mean"][:].astype(dtype),
        features_std=data_set["fscale_stnd"][:].astype(dtype),
        outputs_mean=data_set["oscale_mean"][:].astype(dtype),
        outputs_std=data_set["oscale_stnd"][:].astype(dtype),
        output_groups=[30, 29, 29, 30, 30],
    )

    for i, layer in enumerate(l for l in model.modules() if isinstance(l, nn.Linear)):
        layer.weight.data = torch.tensor(data_set[f"w{i+1}"][:].astype(dtype))
        layer.bias.data = torch.tensor(data_set[f"b{i+1}"][:].astype(dtype))

    return model


if __name__ == "__main__":
    # Load the model from the netcdf file and save it to a checkpoint.
    net = load_from_netcdf_params(
        "qobsTTFFFFFTF30FFTFTF30TTFTFTFFF80FFTFTTF2699FFFF_X01_no_qp_no_adv_"
        "surf_F_Tin_qin_disteq_O_Trad_rest_Tadv_qadv_qout_qsed_RESCALED_7epochs"
        "_no_drop_REAL_NN_layers5in61out148_BN_F_te70.nc"
    )
    net.save("nn_state.pt")
    print("Model saved to nn_state.pt")
