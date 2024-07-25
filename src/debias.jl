"""
A debiased Mondrian random forest is determined by:
- `lambda`: the non-negative lifetime parameter
- `n_trees`: the number of Mondrian trees in the forest
- `n_data`: the number of data points
- `n_evals`: the number of evaluation points
- `x_evals`: the evaluation points
- `debias_order`: the order of debiasing to apply
- `significance_level`: the significance level for confidence intervals
- `X_data`: the covariate data
- `Y_data`: the response data
- `debias_scaling`: the lifetime scaling values
- `debias_coeffs`: the debiasing coefficients
- `trees`: a list of the trees used in the forest
- `mu_hat`: the estimated regression function
- `sigma2_hat`: the estimated residual variance function
- `Sigma_hat`: the estimated forest variance function
- `confidence_band`: a confidence band for the regression function
"""
mutable struct DebiasedMondrianForest{d}
    # parameters
    const lambda::Float64
    const n_trees::Int
    const n_data::Int
    const n_evals::Int
    const x_evals::Vector{NTuple{d,Float64}}
    const debias_order::Int
    const significance_level::Float64
    # data
    const X_data::Vector{NTuple{d,Float64}}
    const Y_data::Vector{Float64}
    # estimates
    debias_scaling::Vector{Float64}
    debias_coeffs::Vector{Float64}
    trees::Matrix{MondrianTree{d}}
    mu_hat::Vector{Float64}
    sigma2_hat::Vector{Float64}
    Sigma_hat::Vector{Float64}
    confidence_band::Vector{Tuple{Float64,Float64}}
end

"""
    DebiasedMondrianForest(lambda::Float64, n_trees::Int, x_evals::Vector{NTuple{d,Float64}},
                           debias_order::Int,
                           X_data::Vector{NTuple{d,Float64}}, Y_data::Vector{Float64},
                           estimate_var::Bool=false,
                           significance_level::Float64=0.05) where {d}

Fit a debiased Mondrian random forest to data.
If `estimate_var` is `false`, do not estimate the variance or construct confidence bands.
This can speed up computation significantly.

# Examples

```julia
lambda = 3.0
n_trees = 20
x_evals = [(0.5, 0.5), (0.2, 0.8)]
debias_order = 1
data = generate_uniform_data_uniform_errors(2, 50)
X_data = data["X"]
Y_data= data["Y"]
estimate_var = true
significance_level = 0.05
forest = DebiasedMondrianForest(lambda, n_trees, x_evals, debias_order, X_data, Y_data,
                                estimate_var, significance_level)
```
"""
function DebiasedMondrianForest(lambda::Float64, n_trees::Int, x_evals::Vector{NTuple{d,Float64}},
                                debias_order::Int,
                                X_data::Vector{NTuple{d,Float64}}, Y_data::Vector{Float64},
                                estimate_var::Bool=false,
                                significance_level::Float64=0.05) where {d}
    n_data = length(X_data)
    n_evals = length(x_evals)

    forest = DebiasedMondrianForest(lambda, n_trees, n_data, n_evals, x_evals,
                                    debias_order, significance_level, X_data, Y_data,
                                    Float64[], Float64[],
                                    #! format: off
                                    MondrianTree{d}[;;],
                                    #! format: on
                                    Float64[], Float64[], Float64[],
                                    Tuple{Float64,Float64}[])

    forest.debias_scaling = get_debias_scaling(debias_order)
    forest.debias_coeffs = get_debias_coeffs(debias_order)
    forest.trees = Matrix{MondrianTree{d}}(undef, (n_trees, debias_order + 1))
    for j in 0:debias_order
        new_lambda = lambda * forest.debias_scaling[j + 1]
        for b in 1:n_trees
            forest.trees[b, j + 1] = MondrianTree(d, new_lambda)
        end
    end

    Ns = Array{Int,3}(undef, (n_trees, debias_order + 1, n_evals))
    X_data = forest.X_data
    @inbounds Threads.@threads for s in 1:n_evals
        x_eval = x_evals[s]
        @inbounds for j in 0:debias_order
            @inbounds for b in 1:n_trees
                tree = forest.trees[b, j + 1]
                Ns[b, j + 1, s] = sum(are_in_same_leaf(X, x_eval, tree) for X in X_data)
            end
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

function get_debias_scaling(debias_order::Int)
    J = debias_order
    debias_scaling = [1.5^r for r in 0:J]
    return debias_scaling
end

function get_debias_coeffs(debias_order::Int)
    J = debias_order
    debias_scaling = get_debias_scaling(debias_order)
    A = zeros(J + 1, J + 1)
    for r in 1:(J + 1)
        for s in 1:(J + 1)
            A[r, s] = debias_scaling[r]^(2 - 2 * s)
        end
    end
    e0 = [[1]; [0 for _ in 1:J]]
    debias_coeffs = A \ e0
    return debias_coeffs
end

function estimate_mu_hat(forest::DebiasedMondrianForest{d}, Ns::Array{Int,3}) where {d}
    mu_hat = [0.0 for _ in 1:(forest.n_evals)]
    Y_bar = sum(forest.Y_data) / forest.n_data

    @inbounds Threads.@threads for s in 1:(forest.n_evals)
        x_eval = forest.x_evals[s]
        @inbounds for j in 0:(forest.debias_order)
            coeff = forest.debias_coeffs[j + 1]
            @inbounds for b in 1:(forest.n_trees)
                if Ns[b, j + 1, s] > 0
                    tree = forest.trees[b, j + 1]
                    I = sum(are_in_same_leaf(forest.X_data[i], x_eval, tree)
                            .*
                            forest.Y_data[i] for i in 1:(forest.n_data))
                    mu_hat[s] += coeff * I / Ns[b, j + 1, s]
                else
                    mu_hat[s] += coeff * Y_bar
                end
            end
        end
    end # COV_EXCL_LINE

    forest.mu_hat = mu_hat / forest.n_trees
    return nothing
end

function estimate_sigma2_hat(forest::DebiasedMondrianForest{d}, Ns::Array{Int,3}) where {d}
    n_data = forest.n_data
    sigma2_hat = [0.0 for _ in 1:(forest.n_evals)]
    j = 0
    @assert forest.debias_scaling[j + 1] == 1

    @inbounds Threads.@threads for s in 1:(forest.n_evals)
        x_eval = forest.x_evals[s]
        mu_hat = forest.mu_hat[s]
        @inbounds for b in 1:(forest.n_trees)
            if Ns[b, j + 1, s] > 0
                tree = forest.trees[b, j + 1]
                I = sum(are_in_same_leaf(forest.X_data[i], x_eval, tree)
                        .*
                        (forest.Y_data[i] - mu_hat)^2 for i in 1:n_data)
                sigma2_hat[s] += I / Ns[b, j + 1, s]
            end
        end
    end # COV_EXCL_LINE

    sigma2_hat ./= forest.n_trees
    forest.sigma2_hat = sigma2_hat
    return nothing
end

function estimate_Sigma_hat(forest::DebiasedMondrianForest{d}, Ns::Array{Int,3}) where {d}
    Sigma_hat = [0.0 for _ in 1:(forest.n_evals)]

    @inbounds Threads.@threads for s in 1:(forest.n_evals)
        x_eval = forest.x_evals[s]
        @inbounds for i in 1:(forest.n_data)
            X = forest.X_data[i]
            A = 0.0
            @inbounds for j in 0:(forest.debias_order)
                coeff = forest.debias_coeffs[j + 1]
                @inbounds for b in 1:(forest.n_trees)
                    tree = forest.trees[b, j + 1]
                    if are_in_same_leaf(X, x_eval, tree)
                        A += coeff / Ns[b, j + 1, s]
                    end
                end
            end
            Sigma_hat[s] += (A / forest.n_trees)^2
        end
    end # COV_EXCL_LINE

    n_data = forest.n_data
    lambda = forest.lambda
    Sigma_hat .*= forest.sigma2_hat .* n_data / lambda^d
    forest.Sigma_hat = Sigma_hat
    return nothing
end

function construct_confidence_band(forest::DebiasedMondrianForest{d}) where {d}
    n_data = forest.n_data
    n_evals = forest.n_evals
    lambda = forest.lambda
    mu_hat = forest.mu_hat
    q = quantile(Normal(0, 1), 1 - forest.significance_level / 2)
    width = q .* sqrt.(forest.Sigma_hat) .* sqrt(lambda^d / n_data)
    confidence_band = [(mu_hat[s] - width[s], mu_hat[s] + width[s]) for s in 1:n_evals]
    forest.confidence_band = confidence_band
    return nothing
end

"""
    Base.show(forest::DebiasedMondrianForest{d}) where {d}

Show a debiased Mondrian random forest.
"""
function Base.show(forest::DebiasedMondrianForest{d}) where {d}
    println("lambda: ", forest.lambda)
    println("n_data: ", forest.n_data)
    println("n_trees: ", forest.n_trees)
    println("n_evals: ", length(forest.x_evals))
    println("x_evals: ", forest.x_evals)
    println("debias_order: ", forest.debias_order)
    println("debias_scaling: ", forest.debias_scaling)
    println("debias_coeffs: ", forest.debias_coeffs)
    println("mu_hat: ", forest.mu_hat)
    println("sigma2_hat: ", forest.sigma2_hat)
    return println("Sigma_hat: ", forest.Sigma_hat)
end
