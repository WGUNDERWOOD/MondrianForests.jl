struct MondrianForest
    X::Vector{Vector{Float64}}
    Y::Vector{Float64}
    trees::Vector{MondrianTree}
end

function MondrianForest(d::Int, lambda::Float64, B::Int,
                        X::Vector{Vector{Float64}}, Y::Vector{Float64})

    trees = [MondrianTree(d, lambda) for _ in 1:B]
    return MondrianForest(X, Y, trees)
end
