"""
    select_lifetime_gcv(lambdas::Vector{Float64}, n_trees::Int,
                        X_data::Vector{NTuple{d,Float64}}, Y_data::Vector{Float64},
                        debias_order::Int, n_subsample::Int) where {d}

Select the lifetime parameter for a (debiased) Mondrian random forest
using generalized cross-validation.

# Examples

```julia
lambdas = collect(range(0.5, 10.0, step=0.5))
n_trees = 20
debias_order = 0
n_subsample = 40
data = generate_uniform_data_uniform_errors(2, 50)
X_data = data["X"]
Y_data= data["Y"]
lambda = select_lifetime_gcv(lambdas, n_trees, X_data, Y_data, debias_order, n_subsample)
```
"""
function select_lifetime_gcv(lambdas::Vector{Float64}, n_trees::Int,
                             X_data::Vector{NTuple{d,Float64}}, Y_data::Vector{Float64},
                             debias_order::Int, n_subsample::Int) where {d}
    n_lambdas = length(lambdas)
    gcvs = [NaN for _ in 1:n_lambdas]
    for l in 1:n_lambdas
        lambda = lambdas[l]
        gcvs[l] = get_gcv(lambda, n_trees, X_data, Y_data, debias_order, n_subsample)
    end

    best_l = argmin(gcvs)
    best_lambda = lambdas[best_l]
    return best_lambda
end

"""
    get_gcv(lambda::Float64, n_trees::Int,
            X_data::Vector{NTuple{d,Float64}}, Y_data::Vector{Float64},
            debias_order::Int, n_subsample::Int) where {d}

Get the generalized cross-validation score of a lifetime parameter lambda.

# Examples

```julia
lambda = 3.0
n_trees = 20
debias_order = 0
n_subsample = 40
data = generate_uniform_data_uniform_errors(2, 50)
X_data = data["X"]
Y_data= data["Y"]
lambda = get_gcv(lambdas, n_trees, X_data, Y_data, debias_order, n_subsample)
```
"""
function get_gcv(lambda::Float64, n_trees::Int,
                 X_data::Vector{NTuple{d,Float64}}, Y_data::Vector{Float64},
                 debias_order::Int, n_subsample::Int) where {d}
    n_data = length(X_data)
    a_bar_d = get_a_bar_d(debias_order, d)
    if n_data <= a_bar_d * lambda^d
        return Inf
    else
        @assert 1 <= n_subsample <= n_data
        ids = sample(1:n_data, n_subsample, replace=false)
        X_evals = X_data[ids]
        Y_evals = Y_data[ids]
        forest = DebiasedMondrianForest(lambda, n_trees, X_evals, debias_order, X_data, Y_data)
        gcv = sum((Y_evals - forest.mu_hat) .^ 2) / n_data
        gcv /= (1 - a_bar_d * lambda^d / n_data)
        return gcv
    end
end

function get_a_bar_d(debias_order::Int, d::Int)
    debias_scaling = MondrianForests.get_debias_scaling(debias_order)
    J = debias_order
    return sum(debias_scaling .^ d) / (J + 1)
end
