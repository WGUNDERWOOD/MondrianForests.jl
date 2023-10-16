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

function MondrianForest(lambda::Float64, n_trees::Int, x_evals::Vector{NTuple{d,Float64}},
                        significance_level::Float64, X_data::Vector{NTuple{d,Float64}},
                        Y_data::Vector{Float64}, estimate_var::Bool=false) where {d}
    n_data = length(X_data)
    n_evals = length(x_evals)
    forest = MondrianForest(lambda, n_trees, n_data, n_evals, x_evals, significance_level,
                            X_data, Y_data, MondrianTree{d}[], Float64[], Float64[],
                            Float64[], Tuple{Float64,Float64}[])
    forest.trees = [MondrianTree(d, lambda) for b in 1:n_trees]
    Ns = Matrix{Int}(undef, (n_trees, n_evals))
    @inbounds Threads.@threads for s in 1:n_evals
        @inbounds for b in 1:n_trees
            Ns[b,s] = sum(are_in_same_cell(X, x_evals[s], forest.trees[b]) for X in forest.X_data)
        end
    end
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
            if Ns[b,s] > 0
                I = sum(are_in_same_cell(forest.X_data[i], x_eval, forest.trees[b])
                        .* forest.Y_data[i] for i in 1:(forest.n_data))
                mu_hat[s] += I / Ns[b,s]
            else
                mu_hat[s] += Y_bar
            end
        end
    end

    forest.mu_hat = mu_hat / forest.n_trees
    return nothing
end

function estimate_sigma2_hat(forest::MondrianForest{d}, Ns::Matrix{Int}) where {d}
    sigma2_hat = [0.0 for _ in 1:(forest.n_evals)]

    @inbounds Threads.@threads for s in 1:(forest.n_evals)
        @inbounds for b in 1:(forest.n_trees)
            if Ns[b,s] > 0
                I = sum(are_in_same_cell(forest.X_data[i], forest.x_evals[s], forest.trees[b])
                        .* (forest.Y_data[i] - forest.mu_hat[s])^2 for i in 1:(forest.n_data))
                sigma2_hat[s] += I / Ns[b,s]
            end
        end
    end

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
                if are_in_same_cell(forest.X_data[i], x_eval, forest.trees[b])
                    A += 1 / Ns[b,s]
                end
            end
            Sigma_hat[s] += (A / forest.n_trees)^2
        end
    end

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
