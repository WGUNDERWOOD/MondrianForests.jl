using Random
using Distributions
using StaticArrays

struct MondrianCell{d}
    lower::SVector{d,Float64}
    upper::SVector{d,Float64}
end

function MondrianCell(d::Int)
    lower = zeros(SVector{d,Float64})
    upper = ones(SVector{d,Float64})
    return MondrianCell(lower, upper)
end

struct MondrianTree{d}
    lambda::Float64
    creation_time::Float64
    cell::MondrianCell{d}
    split_axis::Union{Int,Nothing}
    split_location::Union{Float64,Nothing}
    tree_left::Union{MondrianTree,Nothing}
    tree_right::Union{MondrianTree,Nothing}
end

function MondrianTree(lambda::Float64, creation_time::Float64, cell::MondrianCell)
    d = length(cell.lower)
    size_cell = sum(cell.upper .- cell.lower)
    E = randexp() / size_cell
    if creation_time + E <= lambda
        split_axis = rand(DiscreteNonParametric(1:d, (cell.upper .- cell.lower) ./ size_cell))
        split_location = rand(Uniform(cell.lower[split_axis], cell.upper[split_axis]))
        left_upper = MArray(cell.upper)
        left_upper[split_axis] = split_location
        cell_left = MondrianCell(cell.lower, SVector(left_upper))
        right_lower = MArray(cell.lower)
        right_lower[split_axis] = split_location
        cell_right = MondrianCell(SVector(right_lower), cell.upper)
        tree_left = MondrianTree(lambda, creation_time + E, cell_left)
        tree_right = MondrianTree(lambda, creation_time + E, cell_right)
        return MondrianTree(lambda, creation_time, cell, split_axis, split_location, tree_left, tree_right)
    else
        return MondrianTree(lambda, creation_time, cell, nothing, nothing, nothing, nothing)
    end
end

function MondrianTree(d::Int, lambda::Float64)
    @assert lambda >= 0
    @assert d >= 1
    return MondrianTree(lambda, 0.0, MondrianCell(d))
end
