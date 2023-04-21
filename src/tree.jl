using Random
using Distributions

struct MondrianCell
    lower::Vector{Float64}
    upper::Vector{Float64}
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

function is_leaf(tree::MondrianTree)
    return isnothing(tree.split_axis)
end

function get_id(x::Vector{Float64}, tree::MondrianTree)
    if is_leaf(tree)
        return tree.id
    else
        if x[tree.split_axis] <= tree.split_location
            return get_id(x, tree.tree_left)
        else
            return get_id(x, tree.tree_right)
        end
    end
end

function MondrianTree(d::Int, lambda::Float64)
    @assert lambda >= 0
    @assert d >= 1
    return MondrianTree(lambda, "", 0.0, MondrianCell(d))
end