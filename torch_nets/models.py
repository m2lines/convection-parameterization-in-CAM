"""Neural network architectures."""

from torch.nn import Module, Linear, Dropout
from torch.nn import functional as F


class ANN(Module):
    """This seems to be the model used in the paper.

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

    """

    def __init__(self, n_in, n_out, neurons=128, dropout=0.0):
        """Build ``ANN``."""
        super().__init__()
        self.linear1 = Linear(n_in, neurons)
        self.linear2 = Linear(neurons, neurons)
        self.linear3 = Linear(neurons, neurons)
        self.linear4 = Linear(neurons, neurons)
        self.linear5 = Linear(neurons, n_out)

        self.lin_drop = Dropout(
            dropout
        )  # regularization method to prevent overfitting.
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
        batch = self.lin_drop(batch)

        batch = F.relu(self.linear2(batch))
        batch = self.lin_drop2(batch)

        batch = F.relu(self.linear3(batch))
        batch = self.lin_drop3(batch)

        batch = F.relu(self.linear4(batch))
        batch = self.lin_drop4(batch)

        batch = self.linear5(batch)

        return batch
