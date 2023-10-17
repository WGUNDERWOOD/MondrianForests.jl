# MondrianForests.jl <img src="replication/logo/logo.svg" alt="Mondrian forests logo" align="right" width=150 />

Mondrian random forests in Julia

[![CI](https://github.com/WGUNDERWOOD/MondrianForests.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/WGUNDERWOOD/MondrianForests.jl/actions/workflows/CI.yml)
[![license: GPLv3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![codecov](https://codecov.io/gh/WGUNDERWOOD/MondrianForests.jl/branch/main/graph/badge.svg?token=DbGSOsocw6)](https://codecov.io/gh/WGUNDERWOOD/MondrianForests.jl/)
[![docs](https://img.shields.io/readthedocs/MondrianForests.jl?label=Julia%20docs)](https://wgunderwood.github.io/MondrianForests.jl/stable/)

## Introduction

This repository provides implementations of Mondrian random forests in Julia.
This code is based on methods detailed in
[Cattaneo, Klusowski and Underwood, 2023, arXiv:2310:09702](https://arxiv.org/abs/2310.09702).
This package provides:

- Fitting Mondrian random forests
- Fitting debiased Mondrian random forests
- Selecting the lifetime parameter with polynomial estimation
- Selecting the lifetime parameter with generalized cross-validation

### Branches

The main branch contains stable versions.
Other branches may be unstable,
and are for development purposes only.

### License

This repository and its included Julia package are licensed under
[GPLv3](http://gplv3.fsf.org/).

## Julia package

The Julia package is named **MondrianForests.jl**

### Installation

From the Julia General registry:

```
using Pkg
Pkg.add("MondrianForests")
```

### Usage

TODO give some examples here

### Dependencies

- Distributions
- Random
- Suppressor
- Test

### Documentation
Documentation for the **MondrianForests** package is available on
[the web](https://wgunderwood.github.io/MondrianForests.jl/stable/).
