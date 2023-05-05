function select_lifetime_global_polynomial(X_data::Vector{NTuple{d,Float64}},
                                           debias_order::Int) where {d}

    design_matrix = make_design_matrix_global_polynomial(X_data, debias_order)
    derivative_vectors = get_derivative_vectors_global_polynomial(X_data, debias_order)
    display(derivative_vectors)

end

function make_design_matrix_global_polynomial(X_data::Vector{NTuple{d,Float64}},
                                              debias_order::Int) where {d}

    n = length(X_data)
    design_matrix = ones((n, 1))

    for j in 1:d
        for s in 1:(2 * debias_order + 4)
            column = [x[j]^s for x in X_data]
            design_matrix = hcat(design_matrix, column)
        end
    end

    return design_matrix
end


function get_derivative_vectors_global_polynomial(X_data::Vector{NTuple{d,Float64}},
                                                  debias_order::Int) where {d}

    n = length(X_data)
    derivative_vectors = [[Float64[] for j in 1:d] for i in 1:n]

    for i in 1:n
        for j in 1:d
            derivative_vector = zeros(1 + (j-1) * (2 * debias_order + 4) + (2 * debias_order + 1))
            Xij = X_data[i][j]
            append!(derivative_vector, [1, Xij, Xij^2 / 2])
            append!(derivative_vector, zeros((d - j) * (2 * debias_order + 4)))
            push!(derivative_vectors[i], derivative_vector)
        end
    end

    return derivative_vectors
end
