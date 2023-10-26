using Random
using Distributions

# TODO rewrite functions based on subtrees
# TODO rewrite replication files to use new functions
# TODO docs

"""
A Mondrian tree is determined by:
- `id`: a string to identify the tree
- `lambda`: the non-negative lifetime parameter
- `lower`: the lower coordinate of the cell
- `upper`: the upper coordinate of the cell
- `creation_time`: the time when the cell was created during sampling
- `is_split`: whether the cell is split
- `split_axis`: the direction in which the cell is split, if any
- `split_location`: the location on `split_axis` at which the cell is split, if any
- `tree_left`: the left child tree of the cell, if any
- `tree_right`: the right child tree of the cell, if any
"""
struct MondrianTree{d}
    id::String
    lambda::Float64
    lower::NTuple{d,Float64}
    upper::NTuple{d,Float64}
    creation_time::Float64
    is_split::Bool
    split_axis::Union{Int,Nothing}
    split_location::Union{Float64,Nothing}
    tree_left::Union{MondrianTree{d},Nothing}
    tree_right::Union{MondrianTree{d},Nothing}
end

"""
Sample a Mondrian tree with a given lifetime, lower and upper cell coordinates, and creation time.
"""
function MondrianTree(id::String, lambda::Float64, lower::NTuple{d,Float64},
                      upper::NTuple{d,Float64}, creation_time::Float64) where {d}
    size_cell = sum(upper .- lower)
    E = randexp() / size_cell
    if creation_time + E <= lambda
        split_probabilities = collect(upper .- lower) ./ size_cell
        split_axis = rand(DiscreteNonParametric(1:d, split_probabilities))
        split_location = rand(Uniform(lower[split_axis], upper[split_axis]))
        left_upper = collect(upper)
        left_upper[split_axis] = split_location
        left_upper = ntuple(i -> left_upper[i], d)
        right_lower = collect(lower)
        right_lower[split_axis] = split_location
        right_lower = ntuple(i -> right_lower[i], d)
        tree_left = MondrianTree(id * "L", lambda, lower, left_upper, creation_time + E)
        tree_right = MondrianTree(id * "R", lambda, right_lower, upper, creation_time + E)
        tree = MondrianTree(id, lambda, lower, upper, creation_time, true, split_axis,
                            split_location, tree_left, tree_right)
    else
        tree = MondrianTree(id, lambda, lower, upper, creation_time, false,
                            nothing, nothing, nothing, nothing)
    end
    return tree
end

"""Sample a Mondrian tree `M([0,1]^d, lambda)`."""
function MondrianTree(d::Int, lambda::Float64)
    if lambda < 0
        throw(DomainError(lambda, "lambda must be non-negative"))
    else
        lower = ntuple(i -> 0.0, d)
        upper = ntuple(i -> 1.0, d)
        return MondrianTree("", lambda, lower, upper, 0.0)
    end
end

"""
    is_in(x::NTuple{d,Float64}, tree::MondrianTree{d}) where {d}

Check if a point `x` is contained in a Mondrian tree.

# Examples
```jldoctest
x = ntuple(i -> 0.2, 2)
tree = MondrianTree(2, 3.0)
is_in(x, tree)

# output
true
```
"""
function is_in(x::NTuple{d,Float64}, tree::MondrianTree{d}) where {d}
    return all(tree.lower .<= x .<= tree.upper)
end

"""
    get_center(tree::MondrianTree{d}) where {d}

Get the center point of a Mondrian tree.

# Examples
```jldoctest
tree = MondrianTree(2, 3.0)
get_center(tree)

# output
(0.5, 0.5)
```
"""
function get_center(tree::MondrianTree{d}) where {d}
    return (tree.lower .+ tree.upper) ./ 2
end

"""
    get_volume(tree::MondrianTree{d}) where {d}

Get the d-dimensional volume of a Mondrian tree.

# Examples
```jldoctest
tree = MondrianTree(2, 3.0)
get_volume(tree)

# output
1.0
```
"""
function get_volume(tree::MondrianTree{d}) where {d}
    return prod(tree.upper .- tree.lower)
end

function apply_split(tree::MondrianTree{d}, split_lower::NTuple{d,Float64},
                     split_upper::NTuple{d,Float64}, split_time::Float64,
                     split_axis::Int, split_location::Float64) where {d}
    if tree.is_split
        tree_left = apply_split(tree.tree_left, split_lower, split_upper,
                                split_time, split_axis, split_location)
        tree_right = apply_split(tree.tree_right, split_lower, split_upper,
                                split_time, split_axis, split_location)
        return MondrianTree(tree.id, tree.lambda, tree.lower, tree.upper, tree.creation_time,
                            true, tree.split_axis, tree.split_location, tree_left, tree_right)
    else
        if all(tree.lower .< split_upper) && all(split_lower .< tree.upper)
            left_upper = min.(split_upper, tree.upper)
            right_lower = max.(split_lower, tree.lower)
            tree_left = MondrianTree(tree.id * "L", tree.lambda, tree.lower, left_upper,
                                     split_time, false, nothing, nothing, nothing, nothing)
            tree_right = MondrianTree(tree.id * "R", tree.lambda, right_lower, tree.upper,
                                      split_time, false, nothing, nothing, nothing, nothing)
            return MondrianTree(tree.id, tree.lambda, tree.lower, tree.upper, tree.creation_time,
                                true, split_axis, split_location, tree_left, tree_right)
        else
            return tree
        end
    end
end

function get_common_refinement(tree1::MondrianTree{d}, tree2::MondrianTree{d}) where {d}
    @assert tree1.id == tree2.id
    @assert tree1.lambda == tree2.lambda
    @assert tree1.lower == tree2.lower
    @assert tree1.upper == tree2.upper
    @assert tree1.creation_time == tree2.creation_time
    subtrees1 = get_subtrees(tree1)
    subtrees2 = get_subtrees(tree2)
    subtrees = [subtrees1; subtrees2]
    subtrees = [t for t in subtrees if t.is_split]
    subtrees = sort(subtrees, by=(x -> x.tree_left.creation_time))
    tree = MondrianTree(tree1.id, tree1.lambda, tree1.lower, tree1.upper, tree1.creation_time,
                        false, nothing, nothing, nothing, nothing)

    for subtree in subtrees
        split_location = subtree.split_location
        split_axis = subtree.split_axis
        lower = subtree.lower
        upper = subtree.upper
        split_lower = ntuple(j -> (j == split_axis ? split_location : lower[j]), d)
        split_upper = ntuple(j -> (j == split_axis ? split_location : upper[j]), d)
        split_time = subtree.tree_left.creation_time
        tree = apply_split(tree, split_lower, split_upper, split_time,
                           split_axis, split_location)
    end

    return tree
end

function get_common_refinement(trees::Vector{MondrianTree{d}}) where {d}
    @assert !isempty(trees)
    if length(trees) == 1
        return trees[]
    else
        refinement = get_common_refinement(trees[1], trees[2])
        println(length(get_leaves(trees[1])))
        println(length(get_leaves(refinement)))
        println()
        return get_common_refinement([refinement; trees[3:end]])
    end
end

"""
Check if two points are in the same leaf of a Mondrian tree.
"""
function are_in_same_leaf(x1::Vector{Float64}, x2::Vector{Float64}, tree::MondrianTree)
    if tree.is_split
        if is_in(x1, tree.tree_left) && is_in(x2, tree.tree_right)
            return are_in_same_leaf(x1, x2, tree.left)
        elseif is_in(x1, tree.tree_right) && is_in(x2, tree.tree_right)
            return are_in_same_leaf(x1, x2, tree.right)
        else
            return false
        end
    else
        return is_in(x1, tree) && is_in(x2, tree)
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
    return [t for t in get_leaves(tree) if is_in(x, t)][]
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
        return MondrianTree{d}(tree.lambda, tree.lower, tree.upper, tree.creation_time, true,
                               tree.split_axis, tree.split_location, tree_left, tree_right)
    else
        return MondrianTree{d}(tree.lambda, tree.lower, tree.upper, tree.creation_time, false,
                               nothing, nothing, nothing, nothing)
    end
end

"""Show a Mondrian tree."""
function Base.show(tree::MondrianTree{d}) where {d}
    depth = length(tree.id)
    if depth >= 1
        print(repeat("-", length(tree.id)))
        print(" ")
    end
    if depth >= 1
        printstyled("$(tree.id) ", bold=true, color=:light_magenta)
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
        lower = round.(tree.lower, digits=4)
        upper = round.(tree.upper, digits=4)
        printstyled("$lower -- $upper ", color=:green)
    end
    print("\n")
    if tree.is_split
        show(tree.tree_left)
        show(tree.tree_right)
    end
    return nothing
end
