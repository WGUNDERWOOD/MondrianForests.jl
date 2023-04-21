struct MondrianForest
    X::Vector{Vector{Float64}}
    Y::Vector{Float64}
    trees::Vector{MondrianTree}
    ids::Vector{Vector{String}}
end

function MondrianForest(d::Int, lambda::Float64, B::Int,
                        X::Vector{Vector{Float64}}, Y::Vector{Float64})

    trees = [MondrianTree(d, lambda) for _ in 1:B]
    ids = [[get_id(x, tree) for x in X] for tree in trees]
    return MondrianForest(X, Y, trees, ids)
end
