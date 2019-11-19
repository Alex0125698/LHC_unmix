import proxmin
import proxmin.operators
import torch
import torch.nn as nn
from torch.optim.sgd import SGD
from torch.optim.optimizer import required
import numpy as np
from functools import partial
# import matplotlib.pyplot as plt

# two parameters A,S ; model is matrix product
class NMF(nn.Module):
    def __init__(self, B, N, K):
        super(NMF, self).__init__()
        self.A = nn.Parameter(torch.rand(B, K, requires_grad=True))
        self.S = nn.Parameter(torch.rand(K, N, requires_grad=True))

    def forward(self):
        return torch.matmul(self.A, self.S)

class PGM(SGD):
    def __init__(self, params, proxs, lr=required, momentum=0, dampening=0,
                 nesterov=False):
        kwargs = dict(lr=lr, momentum=momentum, dampening=dampening, weight_decay=0, nesterov=nesterov)
        super().__init__(params, **kwargs)

        if len(proxs) != len(self.param_groups):
            raise ValueError("Invalid length of argument proxs: {} instead of {}".format(len(proxs), len(self.param_groups)))

        for group, prox in zip(self.param_groups, list(proxs)):
            group.setdefault('prox', prox)

    def step(self, closure=None):
        # this performs a gradient step
        # optionally with momentum or nesterov acceleration
        super().step(closure=closure)

        for group in self.param_groups:
            prox = group['prox']

            # here we apply the proximal operator to each parameter in a group
            for p in group['params']:
                p.data = prox(p.data)

# Proximal Constraints
#prox_plus = torch.threshold(0,0)
# def prox_unity(X):
#    return X / X.sum(dim=0).unsqueeze(dim=0)

def prox_plus(X):
    """Projection onto non-negative numbers
    """
    below = X < 0
    X[below] = 0
    return X

def prox_unity0(X):
    axis = 0
    X[:] = X / np.sum(X, axis=axis, keepdims=True)
    return X

def prox_unity1(X):
    axis = 1
    X[:] = X / np.sum(X, axis=axis, keepdims=True)
    return X

# some data cube Y: B x N and we want to factor it into K components
n_epoch = 12
B = 100
N = 100
K = 12
Y = torch.rand(B,N)

nmf = NMF(B, N, K)
Y_ = nmf()

# define loss function
loss_fn = nn.MSELoss(reduction='sum')

pA = prox_plus
pS = prox_unity1

# set proximal constraints and learning-rates for each matrix
param_list = [{'params': nmf.A, 'lr': 0.01},
              {'params': nmf.S, 'lr': 10}]
prox_list = [pA,pS]
optimizer = PGM(param_list, prox_list, momentum=0.5)

for epoch in range(n_epoch):
    Y_ = nmf()
    loss = loss_fn(Y_, Y)
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()

print("done")