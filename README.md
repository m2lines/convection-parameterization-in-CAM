# Convection Parameterization in CAM

## Description
This repository contains code as part of an effort to deploy machine learning (ML) models of geophysical parameterisations into the [Community Earth System Model (CESM)](https://www.cesm.ucar.edu/).
This work is part of the [M<sup>2</sup>LInES](https://m2lines.github.io/) project aiming to improve performance of climate models using ML models for subgrid parameterizations.

A Neural Net providing a subgrid parameterization of atmospheric convection in a [single column model](https://www.arm.gov/publications/proceedings/conf04/extended_abs/randall_da.pdf) has been developed and successfully deployed as part of an atmospheric simulation.
The work is described in a [GRL paper](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020GL091363) with [accompanying code available](https://github.com/yaniyuval/Neural_nework_parameterization/tree/v.1.0.3). The repository contains the neural net and its implementation into a simple system for atmospheric modelling, [SAM](http://rossby.msrc.sunysb.edu/~marat/SAM.html).

The aims of this repository are to:
1. develop a standalone fortran module based on this neural net that can be used elsewhere,
2. deploy the module in another atmospheric model, and
3. evaluate its performance.

We may also perform an investigation into interfacing the pytorch implementation of the Neural Net using the [pytorch-fortran bridging code](https://github.com/Cambridge-ICCS/fortran-pytorch-lib) developed at the [Institute of Computing for Climate Science](https://cambridge-iccs.github.io/).

The model will first be deployed into the [Single Column Atmospheric Model (SCAM)](https://www.cesm.ucar.edu/models/simple/scam) - a single column version of the CESM.
We plan to evaluate performance using SCAM in the gateIII configuration for tropical convection in a similar manner described by the [SCAM6 pulication in JAMES](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2018MS001578).
This will compare model performance to data from an intense observation period (IOP) described in an [AMS publication](https://journals.ametsoc.org/view/journals/atsc/36/1/1520-0469_1979_036_0053_saposs_2_0_co_2.xml).

Long term developments of this project will seek to re-deploy more complex ML parameterizations into mode complex atmospheric models such as the [Community Atmospheric Model (CAM)](https://www.cesm.ucar.edu/models/cam) part of the CESM.


## Repository structure

```
├── NN_module
│   └── ...
└── torch_nets
    └── ...
```



## Contents

### `NN_module/`
This folder contains the fortran neural net extracted from the [code referenced above](https://github.com/yaniyuval/Neural_nework_parameterization/tree/v.1.0.3), along with any dependencies, that may be compiled as a standalone fortran module.

Currently there is code that can be built on CSD3 using the included shell script.

This now needs cleaning up, testing, and a proper makefile creating (see open issues #9 and #10).

### ``torch_nets/``
The directory contains the PyTorch versions of the neural networks we are interested in.


## Running tests
Instructions for running the tests to be added here.



## Contributing

This repository is currently private as it is new and work in progress.
Open tickets can be viewed at ['Issues'](https://github.com/m2lines/convection-parameterization-in-CAM/issues).

To contribute find a relevant issue or open a new one and assign yourself to work on it.
Then create a branch in which to add your contribution and open a pull request.
Once ready assign a reviewer and request a code review.
Merging should _only_ be performed once a code review has taken place.
