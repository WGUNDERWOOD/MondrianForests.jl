mutable struct MondrianForest{d}
    # parameters
    const lambda::Float64
    const n_trees::Int
    const n_data::Int
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
    mu_hat::Float64
    sigma2_hat::Float64
    Sigma_hat::Float64
    confidence_interval::Tuple{Float64, Float64}
end

function MondrianForest(lambda::Float64, n_trees::Int, x_evals::Vector{NTuple{d,Float64}},
                        debias_order::Int, significance_level::Float64,
                        X_data::Vector{NTuple{d,Float64}}, Y_data::Vector{Float64}) where {d}

    n_data = length(X_data)
    forest = MondrianForest(lambda, n_trees, n_data, x_evals, debias_order, significance_level,
                            X_data, Y_data, Float64[], Float64[], Vector{MondrianTree{d}}[],
                            NaN, NaN, NaN, ntuple(x -> NaN, 2))
    get_debias_params(forest)
    forest.trees = [[MondrianTree(d, lambda / forest.debias_scaling[r+1])
                     for b in 1:n_trees] for r in 0:debias_order]
    #[[show(t) for t in ts] for ts in forest.trees]
    # TODO think about how to do this
    Ns = [[sum(is_in(forest.X_data[i], forest.cells[r+1][b]) for i in 1:forest.n_data)
          for b in 1:forest.n_trees] for r in 0:forest.debias_order]
    #estimate_mu_hat(forest, Ns)
    #estimate_sigma2_hat(forest, Ns)
    #estimate_Sigma_hat(forest, Ns)
    #construct_confidence_interval(forest)
    return forest
end

function get_debias_params(forest)
    J = forest.debias_order
    forest.debias_scaling = [0.95^r for r in 0:J]
    A = zeros(J+1, J+1)
    for r in 1:J+1
        for s in 1:J+1
            A[r, s] = forest.debias_scaling[r]^(2 * s - 2)
        end
    end
    e0 = [[1]; [0 for _ in 1:J]]
    forest.debias_coeffs = A \ e0
end

function estimate_mu_hat(forest::MondrianForest{d}, Ns::Vector{Vector{Int}}) where {d}

    mu_hat = 0.0

    for r in 0:forest.debias_order
        for b in 1:forest.n_trees
            if Ns[r+1][b] > 0
                I = sum(is_in(forest.X_data[i], forest.cells[r+1][b]) .* forest.Y_data[i]
                        for i in 1:forest.n_data)
                mu_hat += forest.debias_coeffs[r+1] * I / Ns[r+1][b]
            end
        end
    end

    forest.mu_hat = mu_hat / forest.n_trees
end

function estimate_sigma2_hat(forest::MondrianForest{d}, Ns::Vector{Vector{Int}}) where {d}

    sigma2_hat = 0.0

    for r in 0:forest.debias_order
        for b in 1:forest.n_trees
            if Ns[r+1][b] > 0
                I = sum(is_in(forest.X_data[i], forest.cells[r+1][b]) .* forest.Y_data[i]^2
                        for i in 1:forest.n_data)
                sigma2_hat += forest.debias_coeffs[r+1] * I / Ns[r+1][b]
            end
        end
    end

    sigma2_hat /= forest.n_trees
    sigma2_hat -= forest.mu_hat^2
    forest.sigma2_hat = sigma2_hat
end

function estimate_Sigma_hat(forest::MondrianForest{d}, Ns::Vector{Vector{Int}}) where {d}

    Sigma_hat = 0.0

    for i in 1:forest.n_data
        A = 0.0
        for r in 0:forest.debias_order
            for b in 1:forest.n_trees
                if is_in(forest.X_data[i], forest.cells[r+1][b])
                    A += forest.debias_coeffs[r+1] / Ns[r+1][b]
                end
            end
        end
        Sigma_hat += (A / forest.n_trees)^2
    end

    Sigma_hat *= forest.sigma2_hat * forest.n_data / forest.lambda^d
    forest.Sigma_hat = Sigma_hat
end

function construct_confidence_interval(forest::MondrianForest{d}) where {d}
    q = quantile(Normal(0, 1), 1 - forest.significance_level / 2)
    width = q * sqrt(forest.Sigma_hat) * sqrt(forest.lambda^d / forest.n_data)
    forest.confidence_interval = (forest.mu_hat - width, forest.mu_hat + width)
end

function Base.show(forest::MondrianForest{d}) where {d}

    println("lambda: ", forest.lambda)
    println("n_data: ", forest.n_data)
    println("n_trees: ", forest.n_trees)
    println("x_eval: ", forest.x_eval)
    println("debias_order: ", forest.debias_order)
    println("debias_scaling: ", forest.debias_scaling)
    println("debias_coeffs: ", forest.debias_coeffs)
    println("mu_hat: ", forest.mu_hat)
    println("sigma2_hat: ", forest.sigma2_hat)
    println("Sigma_hat: ", forest.Sigma_hat)
end
