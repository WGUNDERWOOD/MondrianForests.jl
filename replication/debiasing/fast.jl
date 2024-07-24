using Distributions

mutable struct FastDebiasedMondrianForest
    # parameters
    const lambda::Float64
    const n_trees::Int
    const n_data::Int
    const x_eval::Vector{Float64}
    const debias_order::Int
    const significance_level::Float64
    # data
    const X_data::Vector{Vector{Float64}}
    const Y_data::Vector{Float64}
    # estimates
    debias_scaling::Vector{Float64}
    debias_coeffs::Vector{Float64}
    #trees::Matrix{MondrianTree{d}}
    mu_hat::Float64
    sigma2_hat::Float64
    Sigma_hat::Float64
    confidence_band::Vector{Tuple{Float64,Float64}}
end

function FastDebiasedMondrianForest(lambda::Float64, n_trees::Int,
        x_eval::Vector{Float64},
        debias_order::Int,
        X_data::Vector{Vector{Float64}}, Y_data::Vector{Float64},
        significance_level::Float64=0.05)

    n_data = length(X_data)

    forest = FastDebiasedMondrianForest(lambda, n_trees, n_data, x_eval,
                                        debias_order, significance_level, X_data, Y_data, Float64[], Float64[], NaN, NaN, NaN,
                                        Tuple{Float64,Float64}[])

    forest.debias_scaling = get_debias_scaling(debias_order)
    forest.debias_coeffs = get_debias_coeffs(debias_order)

    estimate_mu_hat(forest)
    estimate_sigma2_hat(forest)
    estimate_Sigma_hat(forest)
    #construct_confidence_band(forest)
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

function estimate_mu_hat(forest::FastDebiasedMondrianForest)
    mu_hat = 0.0
    Y_bar = sum(forest.Y_data) / forest.n_data
    x_eval = forest.x_eval
    lambda = forest.lambda

    @inbounds for j in 0:(forest.debias_order)
        scaling = forest.debias_scaling[j + 1]
        E_dist = Exponential(1 / (lambda * scaling))
        coeff = forest.debias_coeffs[j + 1]
        @inbounds for b in 1:(forest.n_trees)
            E_lower = rand(E_dist, d)
            E_upper = rand(E_dist, d)
            cell_lower = max.(0, x_eval .- E_lower)
            cell_upper = min.(1, x_eval .+ E_lower)
            denom = sum(all(cell_lower .<= forest.X_data[i] .<= cell_upper)
                        for i in 1:(forest.n_data))
            if denom > 0
                numer = sum(all(cell_lower .<= forest.X_data[i] .<= cell_upper)
                            * forest.Y_data[i] for i in 1:(forest.n_data))
                mu_hat += coeff * numer / denom
            else
                mu_hat += coeff * Y_bar
            end
        end
    end

    forest.mu_hat = mu_hat / forest.n_trees
    return nothing
end


function estimate_sigma2_hat(forest::FastDebiasedMondrianForest)
    sigma2_hat = 0.0
    x_eval = forest.x_eval
    lambda = forest.lambda
    E_dist = Exponential(1 / lambda)

    @inbounds for b in 1:(forest.n_trees)
        E_lower = rand(E_dist, d)
        E_upper = rand(E_dist, d)
        cell_lower = max.(0, x_eval .- E_lower)
        cell_upper = min.(1, x_eval .+ E_lower)
        denom = sum(all(cell_lower .<= forest.X_data[i] .<= cell_upper)
                    for i in 1:(forest.n_data))
        if denom > 0
            numer = sum(all(cell_lower .<= forest.X_data[i] .<= cell_upper)
                        * (forest.Y_data[i] - forest.mu_hat)^2
                        for i in 1:(forest.n_data))
            sigma2_hat += numer / denom
        end
    end

    forest.sigma2_hat = sigma2_hat / forest.n_trees
    return nothing
end

function estimate_Sigma_hat(forest::FastDebiasedMondrianForest)
    Sigma_hat = 0.0
    x_eval = forest.x_eval
    lambda = forest.lambda

    @inbounds for i in 1:(forest.n_data)
        X = forest.X_data[i]
        A = 0.0
        @inbounds for j in 0:(forest.debias_order)
            scaling = forest.debias_scaling[j + 1]
            E_dist = Exponential(1 / (lambda * scaling))
            coeff = forest.debias_coeffs[j + 1]
            @inbounds for b in 1:(forest.n_trees)
                E_lower = rand(E_dist, d)
                E_upper = rand(E_dist, d)
                cell_lower = max.(0, x_eval .- E_lower)
                cell_upper = min.(1, x_eval .+ E_lower)
                if all(cell_lower .<= X .<= cell_upper)
                    A += coeff / Ns[b, j + 1, s]
                end
            end
        end
        Sigma_hat[s] += (A / forest.n_trees)^2
    end

    @inbounds for b in 1:(forest.n_trees)
        E_lower = rand(E_dist, d)
        E_upper = rand(E_dist, d)
        cell_lower = max.(0, x_eval .- E_lower)
        cell_upper = min.(1, x_eval .+ E_lower)
        denom = sum(all(cell_lower .<= forest.X_data[i] .<= cell_upper)
                    for i in 1:(forest.n_data))
        if denom > 0
            numer = sum(all(cell_lower .<= forest.X_data[i] .<= cell_upper)
                        * (forest.Y_data[i] - forest.mu_hat)^2
                        for i in 1:(forest.n_data))
            sigma2_hat += numer / denom
        end
    end

    forest.sigma2_hat = sigma2_hat / forest.n_trees
    return nothing




end

#=

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
=#


lambda = 10.0
n_trees = 1000
d = 1
n = 1000
x_eval = [0.5 for _ in 1:d]
debias_order = 0
significance_level = 0.05

X_dist = Uniform(0, 1)
mu = (x -> sum(sin.(pi .* x)))
sigma = 0.1
eps_dist = Normal(0, sigma)
X = [[rand(X_dist) for j in 1:d] for i in 1:n]
Y = [mu(X[i]) + rand(eps_dist) for i in 1:n]

forest = FastDebiasedMondrianForest(lambda, n_trees,
                           x_eval,
                           debias_order,
                           X, Y,
                           significance_level)
println(forest.mu_hat)
println(forest.sigma2_hat)
