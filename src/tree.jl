using Random
using Distributions

# TODO rename some functions to all or leaf
# TODO rewrite functions based on subtrees
# TODO use bool to check splits, not isnothing
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

#=
function get_intersection(tree1::MondrianTree{d}, tree2::MondrianTree{d}) where {d}
    @assert tree1.lambda == tree2.lambda
    @assert !tree1.is_split && !tree2.is_split
    lower = max.(tree1.lower, tree2.lower)
    upper = min.(tree1.upper, tree2.upper)
    if all(lower .< upper)
        return MondrianTree("", tree1.lambda, lower, upper, 0.0, false,
                            nothing, nothing, nothing, nothing)
    else
        return nothing
    end
end
=#

function get_common_refinement(tree1::MondrianTree{d}, tree2::MondrianTree{d}) where {d}
    @assert tree1.id == tree2.id
    @assert tree1.lambda == tree2.lambda
    if !tree1.is_split && !tree2.is_split
        lower = max.(tree1.lower, tree2.lower)
        upper = min.(tree1.upper, tree2.upper)
        if all(lower .< upper)
            creation_time = max(tree1.creation_time, tree2.creation_time)
            return MondrianTree(tree1.id, tree1.lambda, lower, upper, creation_time,
                                false, nothing, nothing, nothing, nothing)
        else
            return nothing
        end
    elseif !tree1.is_split
        first_split_time = tree2.tree_left.creation_time
        refinement_left = get_common_refinement(tree1, tree2.tree_left)
        refinement_right = get_common_refinement(tree1, tree2.tree_right)
        return MondrianTree(tree1.id, tree1.lambda, tree1.lower, tree1.upper,
                            tree1.creation_time, false, nothing, nothing, nothing, nothing)
    elseif !tree2.is_split
        first_split_time = tree1.tree_left.creation_time
    else
        first_split_time = min(tree1.tree_left.creation_time, tree2.tree_left.creation_time)
    end

    println(first_split_time)
        #if first_split_time == tree1.tree_left.creation_time
            #refinement_left =


    #leaves1 = get_leaves(tree1)
    #leaves2 = get_leaves(tree2)
    #common_refinement = MondrianTree{d}[]
    #for c1 in leaves1
        #for c2 in leaves2
            #c = get_intersection(c1, c2)
            #if isa(c, MondrianTree{d})
                #push!(common_refinement, c)
            #end
        #end
    #end
    #return unique(common_refinement)
end

#=
function get_common_refinement(trees::Vector{MondrianTree{d}}) where {d}
    if length(trees) == 1
        return trees[1]
    else
        return get_common_refinement(trees[1], get_common_refinement(trees[2:end]))
    end
end
=#


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
