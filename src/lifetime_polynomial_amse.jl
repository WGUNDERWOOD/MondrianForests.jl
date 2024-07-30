"""
    select_lifetime_polynomial(X_data::Vector{NTuple{d,Float64}}, Y_data::Vector{Float64},
                               debias_order::Int=0) where {d}

Select the lifetime parameter for a (debiased) Mondrian random forest
using polynomial estimation.

# Examples

```julia
data = generate_uniform_data_uniform_errors(2, 50)
X_data = data["X"]
Y_data= data["Y"]
debias_order = 0
lambda = select_lifetime_polynomial(X_data, Y_data, debias_order)
```
"""
function select_lifetime_polynomial_amse(X_data::Vector{NTuple{d,Float64}}, Y_data::Vector{Float64},
        x_eval::NTuple{d,Float64},
        debias_order::Int=0) where {d}
    n = length(X_data)
    derivative_estimates = get_derivative_estimates_polynomial_amse(X_data, Y_data, x_eval, debias_order)
    sigma2_hat = get_variance_estimate_polynomial(X_data, Y_data, debias_order)

    omega_bar = get_omega_bar(debias_order)
    numerator = (4 * debias_order + 4) * omega_bar^2 / (debias_order + 2)^2
    numerator *= n * sum(derivative_estimates)^2

    denominator = d * sigma2_hat * get_V_omega(debias_order, d)
    lambda_hat = (numerator / denominator)^(1 / (4 * debias_order + 4 + d))

    return lambda_hat
end

function get_derivative_estimates_polynomial_amse(X_data::Vector{NTuple{d,Float64}},
                                            Y_data::Vector{Float64},
                                            x_eval::NTuple{d,Float64},
                                            debias_order::Int) where {d}
    n = length(X_data)
    J = debias_order
    derivative_vectors = Vector{Float64}[]

    for j in 1:d
        derivative_vector = zeros(1 + (j - 1) * (2 * J + 4) + (2 * J + 1))
        append!(derivative_vector, [1, x_eval[j], x_eval[j]^2 / 2])
        append!(derivative_vector, zeros((d - j) * (2 * J + 4)))
        push!(derivative_vectors, derivative_vector)
    end

    design_matrix = make_design_matrix_polynomial(X_data, debias_order)
    regression_vector = (design_matrix' * design_matrix) \ (design_matrix' * Y_data)
    derivative_estimates = Float64[]

    for j in 1:d
        derivative_estimate = derivative_vectors[j]' * regression_vector
        derivative_estimate *= factorial(2 * debias_order + 2)
        push!(derivative_estimates, derivative_estimate)
    end

    return derivative_estimates
end
