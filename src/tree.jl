using Random
using Distributions

struct MondrianCell{d}
    lower::NTuple{d, Float64}
    upper::NTuple{d, Float64}
end

function MondrianCell(d::Int)
    lower = zeros(d)
    upper = ones(d)
    return MondrianCell(lower, upper)
end

struct MondrianTree
    lambda::Float64
    id::String
    creation_time::Float64
    cell::MondrianCell
    split_axis::Union{Int,Nothing}
    split_location::Union{Float64,Nothing}
    tree_left::Union{MondrianTree,Nothing}
    tree_right::Union{MondrianTree,Nothing}
end

function MondrianTree(lambda::Float64, id::String, creation_time::Float64,
                      cell::MondrianCell)
    d = length(cell.lower)
    size_cell = sum(cell.upper .- cell.lower)
    E = randexp() / size_cell
    if creation_time + E <= lambda
        split_probabilities = (cell.upper .- cell.lower) ./ size_cell
        split_axis = rand(DiscreteNonParametric(1:d, split_probabilities))
        split_location = rand(Uniform(cell.lower[split_axis], cell.upper[split_axis]))
        left_upper = copy(cell.upper)
        left_upper[split_axis] = split_location
        cell_left = MondrianCell(cell.lower, left_upper)
        right_lower = copy(cell.lower)
        right_lower[split_axis] = split_location
        cell_right = MondrianCell(right_lower, cell.upper)
        tree_left = MondrianTree(lambda, id * "0", creation_time + E, cell_left)
        tree_right = MondrianTree(lambda, id * "1", creation_time + E, cell_right)
        return MondrianTree(lambda, id, creation_time, cell, split_axis, split_location,
                            tree_left, tree_right)
    else
        return MondrianTree(lambda, id, creation_time, cell, nothing, nothing,
                            nothing, nothing)
    end
end

function MondrianTree(d::Int, lambda::Float64)
    @assert lambda >= 0
    @assert d >= 1
    return MondrianTree(lambda, "", 0.0, MondrianCell(d))
end

function is_in(x::NTuple{d, Float64}, cell::MondrianCell) where {d}
    for i in 1:d
        if (cell.lower[i] >= x[i]) || (x[i] > cell.upper[i])
            return false
        end
    end
    return true
end
