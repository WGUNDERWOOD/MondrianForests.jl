mutable struct MondrianForest{d}
    # parameters
    const lambda::Float64
    const n_trees::Int
    const n_data::Int
    const x_eval::NTuple{d,Float64}
    const debias_order::Int
    const significance_level::Float64
    # data
    const X_data::Vector{NTuple{d,Float64}}
    const Y_data::Vector{Float64}
    # estimates
    debias_scaling::Vector{Float64}
    debias_coeffs::Vector{Float64}
    cells::Vector{Vector{MondrianCell{d}}}
    mu_hat::Float64
    sigma2_hat::Float64
    Sigma_hat::Float64
end

function MondrianForest(lambda::Float64, n_trees::Int, x_eval::NTuple{d,Float64}, debias_order::Int,
                        significance_level::Float64,
                        X_data::Vector{NTuple{d,Float64}}, Y_data::Vector{Float64}) where {d}

    n_data = length(X_data)
    cells = [[sample_mondrian_cell(x_eval, lambda) for b in 1:n_trees] for r in 0:debias_order]
    forest = MondrianForest(lambda, n_trees, n_data, x_eval, debias_order, X_data, Y_data,
                            Float64[], Float64[], cells, NaN, NaN, NaN)
    get_debias_params(forest)
    Ns = [[sum(is_in(forest.X_data[i], forest.cells[r+1][b]) for i in 1:forest.n_data)
          for b in 1:forest.n_trees] for r in 0:forest.debias_order]
    estimate_mu_hat(forest, Ns)
    estimate_sigma2_hat(forest, Ns)
    estimate_Sigma_hat(forest, Ns)
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

    mu_hat = 0.0

    for b in 1:forest.n_trees
        if Ns[1][b] > 0
            I = sum(is_in(forest.X_data[i], forest.cells[1][b]) .* forest.Y_data[i]
                    for i in 1:forest.n_data)
            mu_hat += I / Ns[1][b]
        end
    end

    mu_hat /= forest.n_trees

    sigma2_hat = 0.0

    for b in 1:forest.n_trees
        if Ns[1][b] > 0
            I = sum(is_in(forest.X_data[i], forest.cells[1][b]) .* forest.Y_data[i]^2
                    for i in 1:forest.n_data)
            sigma2_hat += I / Ns[1][b]
        end
    end

    sigma2_hat /= forest.n_trees
    sigma2_hat -= mu_hat^2
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
