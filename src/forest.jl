mutable struct MondrianForest
    const lambda::Float64
    const n_trees::Int
    x::Vector{Float64}
    mu_hat::Float64
    sigma2_hat::Float64
    Sigma_hat::Float64
end

function MondrianForest(d::Int, lambda::Float64, n_trees::Int)
    return MondrianForest(lambda, n_trees, fill(NaN, d), NaN, NaN, NaN)
end

function sample_mondrian_cell(x::Vector{Float64}, lambda::Float64)
    d = length(x)
    E_lower = [rand(Exponential(1)) for _ in 1:d]
    E_upper = [rand(Exponential(1)) for _ in 1:d]
    lower = max.(x .- E_lower ./ lambda, 0)
    upper = min.(x .+ E_upper ./ lambda, 1)
    return MondrianCell(lower, upper)
end

function fit(forest::MondrianForest, x::Vector{Float64}, X::Vector{Vector{Float64}}, Y::Vector{Float64})
    forest.x = x
    estimate_mu_hat(forest, x, X, Y)
    estimate_sigma2_hat(forest, x, X, Y)
    estimate_Sigma_hat(forest, x, X, Y)
end

function estimate_mu_hat(forest::MondrianForest, x::Vector{Float64},
                         X::Vector{Vector{Float64}}, Y::Vector{Float64})

    n_data = length(X)
    forest.mu_hat = 0.0

    for b in 1:forest.n_trees
        N = 0.0
        I = 0.0
        cell = sample_mondrian_cell(forest.x, forest.lambda)
        for i in 1:n_data
            if is_in(X[i], cell)
                N += 1
                I += Y[i]
            end
        end

        if N > 0
            forest.mu_hat += I / N
        end
    end

    forest.mu_hat /= forest.n_trees
end

function estimate_sigma2_hat(forest::MondrianForest, x::Vector{Float64},
                             X::Vector{Vector{Float64}}, Y::Vector{Float64})

    n_data = length(X)
    forest.sigma2_hat = 0.0

    for b in 1:forest.n_trees
        N = 0.0
        I = 0.0
        cell = sample_mondrian_cell(forest.x, forest.lambda)
        for i in 1:n_data
            if is_in(X[i], cell)
                N += 1
                I += Y[i]^2
            end
        end

        if N > 0
            forest.sigma2_hat += I / N
        end
    end

    forest.sigma2_hat /= forest.n_trees
    forest.sigma2_hat -= forest.mu_hat^2
end

function estimate_Sigma_hat(forest::MondrianForest, x::Vector{Float64},
                            X::Vector{Vector{Float64}}, Y::Vector{Float64})

    n_data = length(X)
    forest.Sigma_hat = 0.0
    A = 0.0

    # TODO reuse the cells
    for b1 in 1:forest.n_trees
        println(b1)
        for b2 in 1:forest.n_trees
            if b1 != b2
                N1 = 0.0
                N2 = 0.0
                I = 0.0
                cell1 = sample_mondrian_cell(forest.x, forest.lambda)
                cell2 = sample_mondrian_cell(forest.x, forest.lambda)
                for i in 1:n_data
                    if is_in(X[i], cell1)
                        N1 += 1
                    end
                    if is_in(X[i], cell2)
                        N2 += 1
                    end
                    if is_in(X[i], cell1) && is_in(X[i], cell2)
                        I += 1
                    end
                end

                if (N1 > 0) && (N2 > 0)
                    A += I / (N1 * N2)
                end
            end
        end
    end

    d = length(forest.x)
    forest.Sigma_hat = A * (n_data / forest.lambda^d)
    forest.Sigma_hat /= forest.n_trees * (forest.n_trees-1)
    forest.Sigma_hat *= forest.sigma2_hat
end
