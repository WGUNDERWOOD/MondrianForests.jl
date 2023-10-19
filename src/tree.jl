using Random
using Distributions

# TODO rename some functions to all or leaf
# TODO rewrite functions based on subtrees
# TODO use bool to check splits, not isnothing
# TODO rewrite replication files to use new functions
# TODO docs

"""
A Mondrian tree is determined by:
- `lambda`: the non-negative lifetime parameter
- `creation_time`: the time when the root cell was created during sampling
- `cell`: the root cell of the tree
- `is_split`: whether the root cell is split
- `split_axis`: the direction in which the root cell was split, if any
- `split_location`: the location on `split_axis` at which the root cell was split, if any
- `tree_left`: the left child tree of the root cell, if any
- `tree_right`: the right child tree of the root cell, if any
"""
struct MondrianTree{d}
    lambda::Float64
    creation_time::Float64
    cell::MondrianCell{d}
    is_split::Bool
    split_axis::Union{Int,Nothing}
    split_location::Union{Float64,Nothing}
    tree_left::Union{MondrianTree{d},Nothing}
    tree_right::Union{MondrianTree{d},Nothing}
end

"""Sample a Mondrian tree with a given lifetime, root cell and creation time."""
function MondrianTree(lambda::Float64, creation_time::Float64,
                      cell::MondrianCell{d}) where {d}
    size_cell = sum(cell.upper .- cell.lower)
    E = randexp() / size_cell
    if creation_time + E <= lambda
        split_probabilities = collect(cell.upper .- cell.lower) ./ size_cell
        split_axis = rand(DiscreteNonParametric(1:d, split_probabilities))
        split_location = rand(Uniform(cell.lower[split_axis], cell.upper[split_axis]))
        left_upper = collect(cell.upper)
        left_upper[split_axis] = split_location
        left_upper = ntuple(i -> left_upper[i], d)
        cell_left = MondrianCell(cell.id * "L", cell.lower, left_upper)
        right_lower = collect(cell.lower)
        right_lower[split_axis] = split_location
        right_lower = ntuple(i -> right_lower[i], d)
        cell_right = MondrianCell(cell.id * "R", right_lower, cell.upper)
        tree_left = MondrianTree(lambda, creation_time + E, cell_left)
        tree_right = MondrianTree(lambda, creation_time + E, cell_right)
        tree = MondrianTree(lambda, creation_time, cell, true, split_axis,
                            split_location, tree_left, tree_right)
    else
        tree = MondrianTree(lambda, creation_time, cell, false, nothing,
                            nothing, nothing, nothing)
    end
    return tree
end

"""Sample a Mondrian tree `M([0,1]^d, lambda)`."""
function MondrianTree(d::Int, lambda::Float64)
    if lambda < 0
        throw(DomainError(lambda, "lambda must be non-negative"))
    else
        return MondrianTree(lambda, 0.0, MondrianCell(d))
    end
end

"""
Check if two points are in the same leaf cell of a Mondrian tree.
"""
function are_in_same_cell(x1::Vector{Float64}, x2::Vector{Float64}, tree::MondrianTree)
    if tree.is_split
        if is_in(x1, tree.left.cell) && is_in(x2, tree.left.cell)
            return are_in_same_cell(x1, x2, tree.left)
        elseif is_in(x1, tree.right.cell) && is_in(x2, tree.right.cell)
            return are_in_same_cell(x1, x2, tree.right)
        else
            return false
        end
    else
        return is_in(x1, tree.cell) && is_in(x2, tree.cell)
    end
end

"""Get a list of the subtrees contained in a Mondrian tree."""
function get_subtrees(tree::MondrianTree{d}) where {d}
    if tree.is_split
        subtrees_left = get_subtrees(tree.tree_left)
        subtrees_right = get_subtrees(tree.tree_right)
        return [tree; subtrees_left; subtrees_right]
    else
        return [tree]
    end
end

"""Get a list of the leaves in a Mondrian tree."""
function get_leaves(tree::MondrianTree{d}) where {d}
    return [t for t in get_subtrees(tree) if !t.is_split]
end

"""Get the leaf of a Mondrian tree containing a point `x`."""
function get_leaf_containing(x::NTuple{d,Float64}, tree::MondrianTree{d}) where {d}
    return [t for t in get_leaves(tree) if is_in(x, t.cell)][]
end

"""Count the leaves of a Mondrian tree."""
function count_leaves(tree::MondrianTree{d}) where {d}
    return length(get_leaves(tree))
end

"""Restrict a Mondrian tree to a stopping time."""
function restrict(tree::MondrianTree{d}, time::Float64) where {d}
    if tree.is_split && tree.tree_left.creation_time <= time
        tree_left = restrict(tree.tree_left, time)
        tree_right = restrict(tree.tree_right, time)
        return MondrianTree{d}(tree.lambda, tree.creation_time, tree.cell, true,
                               tree.split_axis, tree.split_location, tree_left, tree_right)
    else
        return MondrianTree{d}(tree.lambda, tree.creation_time, tree.cell, false,
                               nothing, nothing, nothing, nothing)
    end
end

"""Show a Mondrian tree."""
function Base.show(tree::MondrianTree{d}) where {d}
    depth = length(tree.cell.id)
    if depth >= 1
        print(repeat("-", length(tree.cell.id)))
        print(" ")
    end
    if depth >= 1
        printstyled("$(tree.cell.id) ", bold=true, color=:light_magenta)
    else
        printstyled("MondrianTree ", bold=true, color=:yellow)
        print("in ")
        printstyled("dimension $d ", color=:cyan)
        print("with ")
        printstyled("lambda = $(round(tree.lambda, digits=4)) \n", color=:cyan)
        printstyled("Root ", bold=true, color=:light_magenta)
    end
    if tree.is_split
        print("split on axis $(tree.split_axis) ")
        print("at location $(round(tree.split_location, digits=4)) ")
        print("at time $(round(tree.tree_left.creation_time, digits=4)) ")
    else
        print("leaf at ")
        lower = round.(tree.cell.lower, digits=4)
        upper = round.(tree.cell.upper, digits=4)
        printstyled("$lower -- $upper ", color=:green)
    end
    print("\n")
    if tree.is_split
        show(tree.tree_left)
        show(tree.tree_right)
    end
    return nothing
end
