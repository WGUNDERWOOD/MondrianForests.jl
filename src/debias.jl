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

    Ns = Array{Int, 3}(undef, (n_trees, n_evals, debias_order + 1))
    @inbounds for j in 0:debias_order
        @inbounds Threads.@threads for s in 1:n_evals
            @inbounds for b in 1:n_trees
                Ns[b,s,j+1] = sum(are_in_same_cell(X, x_evals[s], debiased_forest.trees[b, j+1])
                                  for X in debiased_forest.X_data)
            end
        end
    end

    # TODO define these
    estimate_mu_hat(debiased_forest, Ns)
    #if estimate_var
        #estimate_sigma2_hat(debiased_forest, Ns)
        #estimate_Sigma_hat(debiased_forest, Ns)
        #construct_confidence_band(debiased_forest)
    #end
    #if get_gcv
        #get_gcv_dof(forest)
    #end
    return debiased_forest
end

function get_debias_params(debiased_forest)
    J = debiased_forest.debias_order
    debiased_forest.debias_scaling = [2.0^r for r in 0:J]
    A = zeros(J + 1, J + 1)
    for r in 1:(J + 1)
        for s in 1:(J + 1)
            A[r, s] = debiased_forest.debias_scaling[r]^(2 * s - 2)
        end
    end
    e0 = [[1]; [0 for _ in 1:J]]
    debiased_forest.debias_coeffs = A \ e0
    return nothing
end

function estimate_mu_hat(debiased_forest::DebiasedMondrianForest{d}, Ns::Array{Int, 3}) where {d}
    mu_hat = [0.0 for _ in 1:(debiased_forest.n_evals)]
    Y_bar = sum(debiased_forest.Y_data) / debiased_forest.n_data

    @inbounds for j in 0:debiased_forest.debias_order
        coeff = debiased_forest.debias_coeffs[j+1]
        @inbounds Threads.@threads for s in 1:(debiased_forest.n_evals)
            x_eval = debiased_forest.x_evals[s]
            @inbounds for b in 1:(debiased_forest.n_trees)
                if Ns[b,s,j+1] > 0
                    tree = debiased_forest.trees[b,j+1]
                    I = sum(are_in_same_cell(debiased_forest.X_data[i], x_eval, tree)
                            .* debiased_forest.Y_data[i] for i in 1:(debiased_forest.n_data))
                    mu_hat[s] += coeff * I / Ns[b,s,j+1]
                else
                    mu_hat[s] += coeff * Y_bar
                end
            end
        end
    end

    debiased_forest.mu_hat = mu_hat / debiased_forest.n_trees
    return nothing
end

#=
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
=#
