"""Neural network architectures."""

from torch.nn import Module, Linear, Dropout, BatchNorm1d
from torch.nn import functional as F


class Net_ANN_no_BN(Module):
    def __init__(self, n_in, n_out, neurons=128, dropoff=0.0):
        super(Net_ANN_no_BN, self).__init__()
        self.linear1 = Linear(n_in, neurons)
        self.linear2 = Linear(neurons, n_out)
        self.lin_drop = Dropout(
            dropoff
        )  # regularization method to prevent overfitt1ing.

    def forward(self, x):
        x = F.relu(self.linear1(x))
        x = self.lin_drop(x)
        x = self.linear2(x)
        return x


class Net_ANN_5(Module):
    def __init__(self, n_in, n_out, dropoff=0.0):
        super(Net_ANN_5, self).__init__()
        self.linear1 = Linear(n_in, 256)
        self.linear2 = Linear(256, 256)
        self.linear3 = Linear(256, 256)
        self.linear4 = Linear(256, 256)
        self.linear5 = Linear(256, n_out)

        self.dense1_bn = BatchNorm1d(256)
        self.dense2_bn = BatchNorm1d(256)
        self.dense3_bn = BatchNorm1d(256)
        self.dense4_bn = BatchNorm1d(256)

        self.lin_drop = Dropout(
            dropoff
        )  # regularization method to prevent overfitting.
        self.lin_drop2 = Dropout(dropoff)
        self.lin_drop3 = Dropout(dropoff)
        self.lin_drop4 = Dropout(dropoff)

    def forward(self, x):
        x = self.dense1_bn(F.relu(self.linear1(x)))
        x = self.lin_drop(x)
        x = self.dense2_bn(F.relu(self.linear2(x)))
        x = self.lin_drop2(x)
        x = self.dense3_bn(F.relu(self.linear3(x)))
        x = self.lin_drop3(x)
        x = self.dense4_bn(F.relu(self.linear4(x)))
        x = self.lin_drop4(x)
        x = self.linear5(x)
        return x


class Net_ANN_3_no_BN(Module):
    def __init__(self, n_in, n_out, neurons=128, dropoff=0.0):
        super(Net_ANN_3_no_BN, self).__init__()
        self.linear1 = Linear(n_in, neurons)
        self.linear2 = Linear(neurons, neurons)
        self.linear3 = Linear(neurons, n_out)

        self.lin_drop = Dropout(
            dropoff
        )  # regularization method to prevent overfitting.
        self.lin_drop2 = Dropout(dropoff)

    def forward(self, x):
        x = F.relu(self.linear1(x))
        x = self.lin_drop(x)
        x = F.relu(self.linear2(x))
        x = self.lin_drop2(x)
        x = self.linear3(x)
        return x


class Net_ANN_4_no_BN(Module):
    def __init__(self, n_in, n_out, neurons=128, dropoff=0.0):
        super(Net_ANN_4_no_BN, self).__init__()
        self.linear1 = Linear(n_in, neurons)
        self.linear2 = Linear(neurons, neurons)
        self.linear3 = Linear(neurons, neurons)
        self.linear4 = Linear(neurons, n_out)

        self.lin_drop = Dropout(
            dropoff
        )  # regularization method to prevent overfitting.
        self.lin_drop2 = Dropout(dropoff)
        self.lin_drop3 = Dropout(dropoff)

    def forward(self, x):
        x = F.relu(self.linear1(x))
        x = self.lin_drop(x)
        x = F.relu(self.linear2(x))
        x = self.lin_drop2(x)
        x = F.relu(self.linear3(x))
        x = self.lin_drop3(x)
        x = self.linear4(x)
        return x


class Net_ANN_5_no_BN(Module):
    """This seems to be the model used in the paper.

    Paper: https://doi.org/10.1029/2020GL091363

    """

    def __init__(self, n_in, n_out, neurons=128, dropoff=0.0):
        super(Net_ANN_5_no_BN, self).__init__()
        self.linear1 = Linear(n_in, neurons)
        self.linear2 = Linear(neurons, neurons)
        self.linear3 = Linear(neurons, neurons)
        self.linear4 = Linear(neurons, neurons)
        self.linear5 = Linear(neurons, n_out)

        self.lin_drop = Dropout(
            dropoff
        )  # regularization method to prevent overfitting.
        self.lin_drop2 = Dropout(dropoff)
        self.lin_drop3 = Dropout(dropoff)
        self.lin_drop4 = Dropout(dropoff)

    def forward(self, x):
        x = F.relu(self.linear1(x))
        x = self.lin_drop(x)
        x = F.relu(self.linear2(x))
        x = self.lin_drop2(x)
        x = F.relu(self.linear3(x))
        x = self.lin_drop3(x)
        x = F.relu(self.linear4(x))
        x = self.lin_drop4(x)
        x = self.linear5(x)
        return x


class Net_ANN_6_no_BN(Module):
    def __init__(self, n_in, n_out, neurons=128, dropoff=0.0):
        super(Net_ANN_6_no_BN, self).__init__()
        self.linear1 = Linear(n_in, neurons)
        self.linear2 = Linear(neurons, neurons)
        self.linear3 = Linear(neurons, neurons)
        self.linear4 = Linear(neurons, neurons)
        self.linear5 = Linear(neurons, neurons)
        self.linear6 = Linear(neurons, n_out)

        self.lin_drop = Dropout(
            dropoff
        )  # regularization method to prevent overfitting.
        self.lin_drop2 = Dropout(dropoff)
        self.lin_drop3 = Dropout(dropoff)
        self.lin_drop4 = Dropout(dropoff)
        self.lin_drop5 = Dropout(dropoff)

    def forward(self, x):
        x = F.relu(self.linear1(x))
        x = self.lin_drop(x)
        x = F.relu(self.linear2(x))
        x = self.lin_drop2(x)
        x = F.relu(self.linear3(x))
        x = self.lin_drop3(x)
        x = F.relu(self.linear4(x))
        x = self.lin_drop4(x)
        x = F.relu(self.linear5(x))
        x = self.lin_drop5(x)
        x = self.linear6(x)
        return x
