function select_lifetime_polynomial(X_data::Vector{NTuple{d,Float64}}, Y_data::Vector{Float64},
                                    debias_order::Int) where {d}
    n = length(X_data)
    derivative_estimates = get_derivative_estimates_polynomial(X_data, Y_data, debias_order)
    sigma2_hat = get_variance_estimate_polynomial(X_data, Y_data, debias_order)

    omega_bar = get_omega_bar(debias_order)
    numerator = (4 * debias_order + 4) * omega_bar^2 / (debias_order + 2)^2
    numerator *= sum(sum(derivative_estimates[i])^2 for i in 1:n)

    denominator = d * sigma2_hat * get_V_omega(debias_order, d)
    lambda_hat = (numerator / denominator)^(1 / (4 * debias_order + 4 + d))

    return lambda_hat
end

function make_design_matrix_polynomial(X_data::Vector{NTuple{d,Float64}},
                                       debias_order::Int) where {d}
    n = length(X_data)
    J = debias_order
    design_matrix = ones((n, 1))

    for j in 1:d
        for s in 1:(2 * J + 4)
            column = [x[j]^s for x in X_data]
            design_matrix = hcat(design_matrix, column)
        end
    end

    @assert size(design_matrix) == (n, (2 * J + 4) * d + 1)
    return design_matrix
end

function get_derivative_estimates_polynomial(X_data::Vector{NTuple{d,Float64}},
                                             Y_data::Vector{Float64},
                                             debias_order::Int) where {d}
    n = length(X_data)
    J = debias_order
    derivative_vectors = [Vector{Float64}[] for i in 1:n]

    for i in 1:n
        for j in 1:d
            derivative_vector = zeros(1 + (j - 1) * (2 * J + 4) + (2 * J + 1))
            Xij = X_data[i][j]
            append!(derivative_vector, [1, Xij, Xij^2 / 2])
            append!(derivative_vector, zeros((d - j) * (2 * J + 4)))
            push!(derivative_vectors[i], derivative_vector)
        end
    end

    design_matrix = make_design_matrix_polynomial(X_data, debias_order)
    regression_vector = (design_matrix' * design_matrix) \ (design_matrix' * Y_data)
    derivative_estimates = [Float64[] for i in 1:n]

    for i in 1:n
        for j in 1:d
            derivative_estimate = derivative_vectors[i][j]' * regression_vector
            derivative_estimate *= factorial(2 * debias_order + 2)
            push!(derivative_estimates[i], derivative_estimate)
        end
    end

    return derivative_estimates
end

function get_variance_estimate_polynomial(X_data::Vector{NTuple{d,Float64}},
                                          Y_data::Vector{Float64},
                                          debias_order::Int) where {d}
    n = length(X_data)
    design_matrix = make_design_matrix_polynomial(X_data, debias_order)
    regression_vector = (design_matrix' * design_matrix) \ (design_matrix' * Y_data)
    sigma2_hat = Y_data' * Y_data - Y_data' * design_matrix * regression_vector
    sigma2_hat /= n - (2 * debias_order + 4) * d - 1
    return sigma2_hat
end

function get_V_omega(debias_order::Int, d::Int)
    J = debias_order
    a = [0.95^r for r in 0:J]
    A = zeros(J + 1, J + 1)
    for r in 1:(J + 1)
        for s in 1:(J + 1)
            A[r, s] = a[r]^(2 * s - 2)
        end
    end
    e0 = [[1]; [0 for _ in 1:J]]
    omega = A \ e0

    V = zeros(J + 1, J + 1)
    for r in 1:(J + 1)
        for s in 1:(J + 1)
            V[r, s] = 2 / (3 * a[r]) * (1 - a[s] / a[r] * log(a[r] / a[s] + 1))
        end
    end

    V = V + V'

    return sum(sum(V[r, s]^d * omega[r] * omega[s] for r in 1:(J + 1)) for s in 1:(J + 1))
end

function get_omega_bar(debias_order::Int)
    debias_scaling = MondrianForests.get_debias_scaling(debias_order)
    debias_coeffs = MondrianForests.get_debias_coeffs(debias_order)
    J = debias_order
    return sum(debias_coeffs .* debias_scaling .^ (-2 * J - 2))
end
