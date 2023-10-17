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

# TODO switch index order in Ns to be b, j, s

function DebiasedMondrianForest(lambda::Float64, n_trees::Int, x_evals::Vector{NTuple{d,Float64}},
        debias_order::Int, significance_level::Float64,
        X_data::Vector{NTuple{d,Float64}}, Y_data::Vector{Float64},
        estimate_var::Bool=false) where {d}
    n_data = length(X_data)
    n_evals = length(x_evals)

    debiased_forest = DebiasedMondrianForest(lambda, n_trees, n_data, n_evals, x_evals,
                                             debias_order, significance_level, X_data, Y_data,
                                             Float64[], Float64[], MondrianTree{d}[;;],
                                             Float64[], Float64[], Float64[],
                                             Tuple{Float64,Float64}[])

    get_debias_params(debiased_forest)
    debiased_forest.trees = Matrix{MondrianTree{d}}(undef, (n_trees, debias_order + 1))
    for j in 0:debias_order
        new_lambda = lambda * debiased_forest.debias_scaling[j + 1]
        for b in 1:n_trees
            debiased_forest.trees[b, j+1] = MondrianTree(d, new_lambda)
        end
    end

    Ns = Array{Int, 3}(undef, (n_trees, debias_order + 1, n_evals))
    X_data = debiased_forest.X_data
    @inbounds Threads.@threads for s in 1:n_evals
        x_eval = x_evals[s]
        @inbounds for j in 0:debias_order
            @inbounds for b in 1:n_trees
                tree = debiased_forest.trees[b, j+1]
                Ns[b,j+1,s] = sum(are_in_same_cell(X, x_eval, tree) for X in X_data)
            end
        end
    end

    estimate_mu_hat(debiased_forest, Ns)
    if estimate_var
        estimate_sigma2_hat(debiased_forest, Ns)
        estimate_Sigma_hat(debiased_forest, Ns)
        construct_confidence_band(debiased_forest)
    end
    return debiased_forest
end

function get_debias_params(debiased_forest)
    J = debiased_forest.debias_order
    debiased_forest.debias_scaling = [1.5^r for r in 0:J]
    A = zeros(J + 1, J + 1)
    for r in 1:(J+1)
        for s in 1:(J+1)
            A[r, s] = debiased_forest.debias_scaling[r]^(2 - 2 * s)
        end
    end
    e0 = [[1]; [0 for _ in 1:J]]
    debiased_forest.debias_coeffs = A \ e0
    return nothing
end

function estimate_mu_hat(debiased_forest::DebiasedMondrianForest{d}, Ns::Array{Int, 3}) where {d}
    mu_hat = [0.0 for _ in 1:(debiased_forest.n_evals)]
    Y_bar = sum(debiased_forest.Y_data) / debiased_forest.n_data

    @inbounds Threads.@threads for s in 1:(debiased_forest.n_evals)
        x_eval = debiased_forest.x_evals[s]
        @inbounds for j in 0:debiased_forest.debias_order
            coeff = debiased_forest.debias_coeffs[j+1]
            @inbounds for b in 1:(debiased_forest.n_trees)
                if Ns[b,j+1,s] > 0
                    tree = debiased_forest.trees[b,j+1]
                    I = sum(are_in_same_cell(debiased_forest.X_data[i], x_eval, tree)
                            .* debiased_forest.Y_data[i] for i in 1:(debiased_forest.n_data))
                    mu_hat[s] += coeff * I / Ns[b,j+1,s]
                else
                    mu_hat[s] += coeff * Y_bar
                end
            end
        end
    end

    debiased_forest.mu_hat = mu_hat / debiased_forest.n_trees
    return nothing
end

function estimate_sigma2_hat(debiased_forest::DebiasedMondrianForest{d}, Ns::Array{Int, 3}) where {d}
    n_data = debiased_forest.n_data
    sigma2_hat = [0.0 for _ in 1:(debiased_forest.n_evals)]
    j = 0
    @assert debiased_forest.debias_scaling[j+1] == 1

    @inbounds Threads.@threads for s in 1:(debiased_forest.n_evals)
        x_eval = debiased_forest.x_evals[s]
        mu_hat = debiased_forest.mu_hat[s]
        @inbounds for b in 1:(debiased_forest.n_trees)
            if Ns[b,j+1,s] > 0
                tree = debiased_forest.trees[b,j+1]
                I = sum(are_in_same_cell(debiased_forest.X_data[i], x_eval, tree)
                        .* (debiased_forest.Y_data[i] - mu_hat)^2 for i in 1:n_data)
                sigma2_hat[s] += I / Ns[b,j+1,s]
            end
        end
    end

    sigma2_hat ./= debiased_forest.n_trees
    debiased_forest.sigma2_hat = sigma2_hat
    return nothing
end

function estimate_Sigma_hat(debiased_forest::DebiasedMondrianForest{d}, Ns::Array{Int, 3}) where {d}
    Sigma_hat = [0.0 for _ in 1:(debiased_forest.n_evals)]

    @inbounds Threads.@threads for s in 1:(debiased_forest.n_evals)
        x_eval = debiased_forest.x_evals[s]
        @inbounds for i in 1:(debiased_forest.n_data)
            X = debiased_forest.X_data[i]
            A = 0.0
            @inbounds for j in 0:(debiased_forest.debias_order)
                coeff = debiased_forest.debias_coeffs[j+1]
                @inbounds for b in 1:(debiased_forest.n_trees)
                    tree = debiased_forest.trees[b,j+1]
                    if are_in_same_cell(X, x_eval, tree)
                        A += coeff / Ns[b,j+1,s]
                    end
                end
            end
            Sigma_hat[s] += (A / debiased_forest.n_trees)^2
        end
    end

    n_data = debiased_forest.n_data
    lambda = debiased_forest.lambda
    Sigma_hat .*= debiased_forest.sigma2_hat .* n_data / lambda^d
    debiased_forest.Sigma_hat = Sigma_hat
    return nothing
end


function construct_confidence_band(debiased_forest::DebiasedMondrianForest{d}) where {d}
    n_data = debiased_forest.n_data
    n_evals = debiased_forest.n_evals
    lambda = debiased_forest.lambda
    mu_hat = debiased_forest.mu_hat
    q = quantile(Normal(0, 1), 1 - debiased_forest.significance_level / 2)
    width = q .* sqrt.(debiased_forest.Sigma_hat) .* sqrt(lambda^d / n_data)
    confidence_band = [(mu_hat[s] - width[s], mu_hat[s] + width[s]) for s in 1:n_evals]
    debiased_forest.confidence_band = confidence_band
    return nothing
end

function Base.show(debiased_forest::DebiasedMondrianForest{d}) where {d}
    println("lambda: ", debiased_forest.lambda)
    println("n_data: ", debiased_forest.n_data)
    println("n_trees: ", debiased_forest.n_trees)
    println("n_evals: ", length(debiased_forest.x_evals))
    println("x_evals: ", debiased_forest.x_evals)
    println("debias_order: ", debiased_forest.debias_order)
    println("debias_scaling: ", debiased_forest.debias_scaling)
    println("debias_coeffs: ", debiased_forest.debias_coeffs)
    println("mu_hat: ", debiased_forest.mu_hat)
    println("sigma2_hat: ", debiased_forest.sigma2_hat)
    return println("Sigma_hat: ", debiased_forest.Sigma_hat)
end
