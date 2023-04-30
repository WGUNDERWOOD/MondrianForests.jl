struct MondrianForest{d}
    # parameters
    lambda::Float64
    n_trees::Int
    n_data::Int
    n_data_reduced::Int
    x_eval::NTuple{d,Float64}
    # trees
    cells::Vector{MondrianCell{d}}
    largest_cell::MondrianCell{d}
    # data
    X_data::Vector{NTuple{d,Float64}}
    Y_data::Vector{Float64}
    mu_hat::Float64
    sigma2_hat::Float64
    Sigma_hat::Float64
    X_data_reduced::Vector{NTuple{d,Float64}}
    Y_data_reduced::Vector{Float64}
    membership::Vector{Vector{Bool}}
end

function sample_mondrian_cell(x::NTuple{d,Float64}, lambda::Float64) where {d}
    E_lower = [rand(Exponential(1)) for _ in 1:d]
    E_upper = [rand(Exponential(1)) for _ in 1:d]
    lower = max.(x .- E_lower ./ lambda, 0)
    upper = min.(x .+ E_upper ./ lambda, 1)
    lower = ntuple(i -> lower[i], d)
    upper = ntuple(i -> upper[i], d)
    return MondrianCell(lower, upper)
end

function get_largest_cell(cells::Vector{MondrianCell{d}}) where {d}
    lower = [minimum(cell.lower[j] for cell in cells) for j in 1:d]
    upper = [maximum(cell.upper[j] for cell in cells) for j in 1:d]
    lower = ntuple(i -> lower[i], d)
    upper = ntuple(i -> upper[i], d)
    return MondrianCell(lower, upper)
end

function MondrianForest(lambda::Float64, n_trees::Int, x_eval::NTuple{d,Float64},
                        X_data::Vector{NTuple{d,Float64}},
                        Y_data::Vector{Float64}) where {d}
    n_data = length(X_data)
    cells = [sample_mondrian_cell(x_eval, lambda) for _ in 1:n_trees]
    largest_cell = get_largest_cell(cells)
    indices_reduced = [i for i in 1:n_data if is_in(X_data[i], largest_cell)]
    X_data_reduced = [X_data[i] for i in indices_reduced]
    Y_data_reduced = [Y_data[i] for i in indices_reduced]
    n_data_reduced = length(X_data_reduced)
    membership = [[is_in(X_data[i], cells[b]) for i in 1:n_data_reduced] for b in 1:n_trees]
    mu_hat = estimate_mu_hat(membership, Y_data_reduced)
    sigma2_hat = estimate_sigma2_hat(membership, Y_data_reduced, mu_hat)
    Sigma_hat = estimate_Sigma_hat(membership, Y_data_reduced, sigma2_hat, d, lambda)

    return MondrianForest(
                          # parameters
                          lambda,
                          n_trees,
                          n_data,
                          n_data_reduced,
                          x_eval,
                          # trees
                          cells,
                          largest_cell,
                          # data
                          X_data,
                          Y_data,
                          mu_hat,
                          sigma2_hat,
                          Sigma_hat,
                          X_data_reduced,
                          Y_data_reduced,
                          membership)
end

function estimate_mu_hat(membership::Vector{Vector{Bool}}, Y_data_reduced::Vector{Float64})
    mu_hat = 0.0
    n_trees = length(membership)
    Ns = [sum(m for m in membership[b]) for b in 1:n_trees]

    for b in 1:n_trees
        N = Ns[b]
        if N > 0
            I = sum(membership[b] .* Y_data_reduced)
            mu_hat += I / N
        end
    end

    mu_hat /= n_trees
    return mu_hat
end

function estimate_sigma2_hat(membership::Vector{Vector{Bool}}, Y_data_reduced::Vector{Float64},
                             mu_hat::Float64)
    sigma2_hat = 0.0
    n_trees = length(membership)
    Ns = [sum(m for m in membership[b]) for b in 1:n_trees]

    for b in 1:n_trees
        N = Ns[b]
        if N > 0
            I = sum(membership[b] .* (Y_data_reduced .^ 2))
            sigma2_hat += I / N
        end
    end

    sigma2_hat /= n_trees
    sigma2_hat -= mu_hat^2
    return sigma2_hat
end

function estimate_Sigma_hat(membership::Vector{Vector{Bool}}, Y_data_reduced::Vector{Float64},
                            sigma2_hat::Float64, d::Int, lambda::Float64)

    # TODO new faster method for computing this
    A = 0.0
    n_trees = length(membership)
    n_data = length(Y_data_reduced)
    Ns = [sum(m for m in membership[b]) for b in 1:n_trees]

    for b1 in 1:n_trees
        for b2 in 1:n_trees
            if b1 < b2
                N1 = Ns[b1]
                N2 = Ns[b2]
                if (N1 > 0) && (N2 > 0)
                    m1 = membership[b1]
                    m2 = membership[b2]
                    I = sum(m1[i] * m2[i] for i in 1:n_data)
                    A += I / (N1 * N2)
                end
            end
        end
    end

    Sigma_hat = A * (n_data / lambda^d)
    Sigma_hat *= 2 / (n_trees * (n_trees - 1))
    Sigma_hat *= sigma2_hat
    return Sigma_hat
end
