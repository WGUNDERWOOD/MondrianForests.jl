"""
A Mondrian random forest is determined by:
- `lambda`: the non-negative lifetime parameter
- `n_trees`: the number of Mondrian trees in the forest
- `n_data`: the number of data points
- `n_evals`: the number of evaluation points
- `x_evals`: the evaluation points
- `significance_level`: the significance level for confidence intervals
- `X_data`: the covariate data
- `Y_data`: the response data
- `trees`: a list of the trees used in the forest
- `mu_hat`: the estimated regression function
- `sigma2_hat`: the estimated residual variance function
- `Sigma_hat`: the estimated forest variance function
- `confidence_band`: a confidence band for the regression function
"""
mutable struct MondrianForest{d}
    # parameters
    const lambda::Float64
    const n_trees::Int
    const n_data::Int
    const n_evals::Int
    const x_evals::Vector{NTuple{d,Float64}}
    const significance_level::Float64
    # data
    const X_data::Vector{NTuple{d,Float64}}
    const Y_data::Vector{Float64}
    # estimates
    trees::Vector{MondrianTree{d}}
    mu_hat::Vector{Float64}
    sigma2_hat::Vector{Float64}
    Sigma_hat::Vector{Float64}
    confidence_band::Vector{Tuple{Float64,Float64}}
end

"""
    MondrianForest(lambda::Float64, n_trees::Int, x_evals::Vector{NTuple{d,Float64}},
                   X_data::Vector{NTuple{d,Float64}}, Y_data::Vector{Float64},
                   estimate_var::Bool=false, significance_level::Float64=0.05) where {d}

Fit a Mondrian random forest to data.
If `estimate_var` is `false`, do not estimate the variance or construct confidence bands.
This can speed up computation significantly.

# Examples

```julia
lambda = 3.0
n_trees = 20
x_evals = [(0.5, 0.5), (0.2, 0.8)]
data = generate_uniform_data_uniform_errors(2, 50)
X_data = data["X"]
Y_data= data["Y"]
estimate_var = true
significance_level = 0.05
forest = MondrianForest(lambda, n_trees, x_evals, X_data, Y_data,
                        estimate_var, significance_level)
```
"""
function MondrianForest(lambda::Float64, n_trees::Int, x_evals::Vector{NTuple{d,Float64}},
                        X_data::Vector{NTuple{d,Float64}}, Y_data::Vector{Float64},
                        estimate_var::Bool=false, significance_level::Float64=0.05) where {d}
    n_data = length(X_data)
    n_evals = length(x_evals)
    forest = MondrianForest(lambda, n_trees, n_data, n_evals, x_evals, significance_level,
                            X_data, Y_data, MondrianTree{d}[], Float64[], Float64[],
                            Float64[], Tuple{Float64,Float64}[])
    forest.trees = [MondrianTree(d, lambda) for b in 1:n_trees]
    Ns = Matrix{Int}(undef, (n_trees, n_evals))
    @inbounds Threads.@threads for s in 1:n_evals
        @inbounds for b in 1:n_trees
            Ns[b, s] = sum(are_in_same_leaf(X, x_evals[s], forest.trees[b]) for X in forest.X_data)
        end
    end # COV_EXCL_LINE
    estimate_mu_hat(forest, Ns)
    if estimate_var
        estimate_sigma2_hat(forest, Ns)
        estimate_Sigma_hat(forest, Ns)
        construct_confidence_band(forest)
    end
    return forest
end

function estimate_mu_hat(forest::MondrianForest{d}, Ns::Matrix{Int}) where {d}
    mu_hat = [0.0 for _ in 1:(forest.n_evals)]
    Y_bar = sum(forest.Y_data) / forest.n_data

    @inbounds Threads.@threads for s in 1:(forest.n_evals)
        x_eval = forest.x_evals[s]
        @inbounds for b in 1:(forest.n_trees)
            if Ns[b, s] > 0
                I = sum(are_in_same_leaf(forest.X_data[i], x_eval, forest.trees[b])
                        .*
                        forest.Y_data[i] for i in 1:(forest.n_data))
                mu_hat[s] += I / Ns[b, s]
            else
                mu_hat[s] += Y_bar
            end
        end
    end # COV_EXCL_LINE

    forest.mu_hat = mu_hat / forest.n_trees
    return nothing
end

function estimate_sigma2_hat(forest::MondrianForest{d}, Ns::Matrix{Int}) where {d}
    sigma2_hat = [0.0 for _ in 1:(forest.n_evals)]

    @inbounds Threads.@threads for s in 1:(forest.n_evals)
        @inbounds for b in 1:(forest.n_trees)
            if Ns[b, s] > 0
                I = sum(are_in_same_leaf(forest.X_data[i], forest.x_evals[s], forest.trees[b])
                        .*
                        (forest.Y_data[i] - forest.mu_hat[s])^2 for i in 1:(forest.n_data))
                sigma2_hat[s] += I / Ns[b, s]
            end
        end
    end # COV_EXCL_LINE

    sigma2_hat ./= forest.n_trees
    forest.sigma2_hat = sigma2_hat
    return nothing
end

function estimate_Sigma_hat(forest::MondrianForest{d}, Ns::Matrix{Int}) where {d}
    Sigma_hat = [0.0 for _ in 1:(forest.n_evals)]

    @inbounds Threads.@threads for s in 1:(forest.n_evals)
        x_eval = forest.x_evals[s]
        @inbounds for i in 1:(forest.n_data)
            A = 0.0
            @inbounds for b in 1:(forest.n_trees)
                if are_in_same_leaf(forest.X_data[i], x_eval, forest.trees[b])
                    A += 1 / Ns[b, s]
                end
            end
            Sigma_hat[s] += (A / forest.n_trees)^2
        end
    end # COV_EXCL_LINE

    Sigma_hat .*= forest.sigma2_hat .* forest.n_data / forest.lambda^d
    forest.Sigma_hat = Sigma_hat
    return nothing
end

function construct_confidence_band(forest::MondrianForest{d}) where {d}
    q = quantile(Normal(0, 1), 1 - forest.significance_level / 2)
    width = q .* sqrt.(forest.Sigma_hat) .* sqrt(forest.lambda^d / forest.n_data)
    forest.confidence_band = [(forest.mu_hat[s] - width[s], forest.mu_hat[s] + width[s])
                              for s in 1:(forest.n_evals)]
    return nothing
end

"""
    Base.show(forest::MondrianForest{d}) where {d}

Show a Mondrian random forest.
"""
function Base.show(forest::MondrianForest{d}) where {d}
    println("lambda: ", forest.lambda)
    println("n_data: ", forest.n_data)
    println("n_trees: ", forest.n_trees)
    println("n_evals: ", length(forest.x_evals))
    println("x_evals: ", forest.x_evals)
    println("mu_hat: ", forest.mu_hat)
    println("sigma2_hat: ", forest.sigma2_hat)
    println("Sigma_hat: ", forest.Sigma_hat)
    return nothing
end
