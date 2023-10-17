# TODO significance_level can be nothing if not estimate_var

function select_lifetime_gcv(lambdas::Vector{Float64},
        n_trees::Int,
        n_subsample::Int,
        X_data::Vector{NTuple{d,Float64}}, Y_data::Vector{Float64},
        debias_order::Int) where {d}

    n_lambdas = length(lambdas)
    gcvs = [NaN for _ in 1:n_lambdas]
    for l in 1:n_lambdas
        lambda = lambdas[l]
        gcvs[l] = get_gcv(lambda, n_trees, n_subsample, debias_order, X_data, Y_data)
    end

    best_l = argmin(gcvs)
    best_lambda = lambdas[best_l]
    return best_lambda
end


function get_gcv(lambda::Float64, n_trees::Int, n_subsample::Int, debias_order::Int,
        X_data::Vector{NTuple{d,Float64}}, Y_data::Vector{Float64}) where {d}

    n_data = length(X_data)
    a_bar_d = get_a_bar_d(debias_order, d)
    if n_data <= a_bar_d * lambda^d
        return Inf
    else
        @assert 1 <= n_subsample <= n_data
        significance_level = 0.0
        estimate_var = false
        ids = sample(1:n_data, n_subsample, replace=false)
        X_evals = X_data[ids]
        Y_evals = Y_data[ids]
        forest = DebiasedMondrianForest(lambda, n_trees, X_evals, debias_order,
                                        significance_level, X_data, Y_data, estimate_var)
        gcv = sum((Y_evals - forest.mu_hat).^2) / n_data
        gcv /= (1 - a_bar_d * lambda^d / n_data)
        return gcv
    end
end


function get_a_bar_d(debias_order::Int, d::Int)
    debias_scaling = MondrianForests.get_debias_scaling(debias_order)
    J = debias_order
    return sum(debias_scaling .^ d) / (J + 1)
end
