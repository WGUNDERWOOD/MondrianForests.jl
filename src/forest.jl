mutable struct MondrianForest{d}
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
    trees::Vector{Vector{MondrianTree{d}}}
    mu_hat::Vector{Float64}
    sigma2_hat::Vector{Float64}
    Sigma_hat::Vector{Float64}
    confidence_band::Vector{Tuple{Float64,Float64}}
    gcv_dof::Union{Nothing, Float64}
end

function MondrianForest(lambda::Float64, n_trees::Int, x_evals::Vector{NTuple{d,Float64}},
                        debias_order::Int, significance_level::Float64,
                        X_data::Vector{NTuple{d,Float64}}, Y_data::Vector{Float64},
                        estimate_var::Bool=false, get_gcv::Bool=false) where {d}
    n_data = length(X_data)
    n_evals = length(x_evals)
    forest = MondrianForest(lambda, n_trees, n_data, n_evals, x_evals, debias_order,
                            significance_level, X_data, Y_data, Float64[], Float64[],
                            Vector{MondrianTree{d}}[], Float64[], Float64[], Float64[],
                            Tuple{Float64,Float64}[], nothing)
    get_debias_params(forest)
    forest.trees = [[MondrianTree(d, lambda / forest.debias_scaling[r + 1])
                     for b in 1:n_trees] for r in 0:debias_order]
    Ns = [Vector{Int}[] for x in x_evals]
    Threads.@threads for s in 1:n_evals
        Ns[s] = [[sum(are_in_same_cell(X, x_evals[s], forest.trees[r + 1][b])
                      for X in forest.X_data)
                  for b in 1:(forest.n_trees)]
                 for r in 0:(forest.debias_order)]
    end
    estimate_mu_hat(forest, Ns)
    if estimate_var
        estimate_sigma2_hat(forest, Ns)
        estimate_Sigma_hat(forest, Ns)
        construct_confidence_band(forest)
    end
    if get_gcv
        get_gcv_dof(forest)
    end
    return forest
end

function get_debias_params(forest)
    J = forest.debias_order
    forest.debias_scaling = [0.95^r for r in 0:J]
    A = zeros(J + 1, J + 1)
    for r in 1:(J + 1)
        for s in 1:(J + 1)
            A[r, s] = forest.debias_scaling[r]^(2 * s - 2)
        end
    end
    e0 = [[1]; [0 for _ in 1:J]]
    return forest.debias_coeffs = A \ e0
end

function estimate_mu_hat(forest::MondrianForest{d}, Ns::Vector{Vector{Vector{Int}}}) where {d}
    mu_hat = [0.0 for _ in 1:(forest.n_evals)]
    Y_bar = sum(forest.Y_data) / forest.n_data

    Threads.@threads for s in 1:(forest.n_evals)
        for r in 0:(forest.debias_order)
            for b in 1:(forest.n_trees)
                if Ns[s][r + 1][b] > 0
                    I = sum(are_in_same_cell(forest.X_data[i], forest.x_evals[s],
                                             forest.trees[r + 1][b])
                            .*
                            forest.Y_data[i] for i in 1:(forest.n_data))
                    mu_hat[s] += forest.debias_coeffs[r + 1] * I / Ns[s][r + 1][b]
                else
                    mu_hat[s] += Y_bar
                end
            end
        end
    end

    return forest.mu_hat = mu_hat / forest.n_trees
end

function estimate_sigma2_hat(forest::MondrianForest{d}, Ns::Vector{Vector{Vector{Int}}}) where {d}
    mu_hat = [0.0 for _ in 1:(forest.n_evals)]

    for s in 1:(forest.n_evals)
        for b in 1:(forest.n_trees)
            if Ns[s][1][b] > 0
                I = sum(are_in_same_cell(forest.X_data[i], forest.x_evals[s],
                                         forest.trees[1][b])
                        .*
                        forest.Y_data[i] for i in 1:(forest.n_data))
                mu_hat[s] += I / Ns[s][1][b]
            end
        end
    end

    mu_hat ./= forest.n_trees

    sigma2_hat = [0.0 for _ in 1:(forest.n_evals)]

    for s in 1:(forest.n_evals)
        for b in 1:(forest.n_trees)
            if Ns[s][1][b] > 0
                I = sum(are_in_same_cell(forest.X_data[i], forest.x_evals[s],
                                         forest.trees[1][b])
                        .*
                        forest.Y_data[i]^2 for i in 1:(forest.n_data))
                sigma2_hat[s] += I / Ns[s][1][b]
            end
        end
    end

    sigma2_hat ./= forest.n_trees
    sigma2_hat .-= mu_hat .^ 2
    return forest.sigma2_hat = sigma2_hat
end

function estimate_Sigma_hat(forest::MondrianForest{d}, Ns::Vector{Vector{Vector{Int}}}) where {d}
    Sigma_hat = [0.0 for _ in 1:(forest.n_evals)]

    for s in 1:(forest.n_evals)
        for i in 1:(forest.n_data)
            A = 0.0
            for r in 0:(forest.debias_order)
                for b in 1:(forest.n_trees)
                    if are_in_same_cell(forest.X_data[i], forest.x_evals[s], forest.trees[r + 1][b])
                        A += forest.debias_coeffs[r + 1] / Ns[s][r + 1][b]
                    end
                end
            end
            Sigma_hat[s] += (A / forest.n_trees)^2
        end
    end

    Sigma_hat .*= forest.sigma2_hat .* forest.n_data / forest.lambda^d
    return forest.Sigma_hat = Sigma_hat
end

function get_gcv_dof(forest::MondrianForest{d}) where {d}
    gcv_dof = 0.0
    for i in 1:forest.n_data
        for b in 1:forest.n_trees
            Xi = forest.X_data[i]
            Tb = forest.trees[1][b]
            S = sum(are_in_same_cell(Xi, forest.X_data[j], Tb) for j in 1:forest.n_data)
            gcv_dof += 1 / (S * forest.n_trees)
        end
    end
    forest.gcv_dof = gcv_dof
end

function construct_confidence_band(forest::MondrianForest{d}) where {d}
    q = quantile(Normal(0, 1), 1 - forest.significance_level / 2)
    width = q .* sqrt.(forest.Sigma_hat) .* sqrt(forest.lambda^d / forest.n_data)
    return forest.confidence_band = [(forest.mu_hat[s] - width[s], forest.mu_hat[s] + width[s])
                                     for s in 1:(forest.n_evals)]
end

function Base.show(forest::MondrianForest{d}) where {d}
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
