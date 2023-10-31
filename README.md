# MondrianForests.jl <img src="docs/src/assets/logo.svg" alt="Mondrian forests logo" align="right" width=220 />

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

```julia
using Pkg
Pkg.add("MondrianForests")
```

### Usage

```julia
using MondrianForests

# sample a two-dimensional Mondrian tree
d = 2
lambda = 2.0
tree = MondrianTree(d, lambda)
println()
show(tree)
println()

# generate some data
# covariates X are two-dimensional
# response Y is one-dimensional
# true regression function is zero
n_data = 100
data = MondrianForests.generate_uniform_data_uniform_errors(d, n_data)
X = data["X"]
Y = data["Y"]
println("covariates: ")
display(X[1:5])
println("\nresponses: ")
display(Y[1:5])

# select a lifetime parameter
# with generalized cross-validation
n_trees = 50
n_subsample = 30
debias_order = 0
lambdas = collect(range(0.5, 10.0, step=0.5))
lambda = select_lifetime_gcv(lambdas, n_trees, X, Y, debias_order, n_subsample)
println("\nlambda chosen by GCV: ", lambda)

# fit and evaluate a Mondrian random forest
x_evals = [(0.5, 0.5), (0.2, 0.8)]
estimate_var = true
forest = MondrianForest(lambda, n_trees, x_evals, X, Y, estimate_var)
println("\nestimated regression function:")
display(forest.mu_hat)
println("\nestimated estimator variance:")
display(forest.Sigma_hat)
println("\nestimated confidence band:")
display(forest.confidence_band)

# fit and evaluate a debiased Mondrian random forest
debiased_forest = DebiasedMondrianForest(lambda, n_trees, x_evals, debias_order, X, Y, estimate_var)
println("\ndebiased estimated regression function:")
display(debiased_forest.mu_hat)
println("\ndebiased estimated estimator variance:")
display(debiased_forest.Sigma_hat)
println("\ndebiased estimated confidence band:")
display(debiased_forest.confidence_band)
```

### Dependencies

- Distributions
- Random
- Suppressor
- Test

### Documentation
Documentation for the **MondrianForests** package is available on
[the web](https://wgunderwood.github.io/MondrianForests.jl/stable/).
