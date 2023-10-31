using Random
using Distributions

# TODO rewrite replication files to use new functions

"""
A Mondrian tree is determined by:
- `id`: a string to identify the tree
- `lambda`: the non-negative lifetime parameter
- `lower`: the lower coordinate of the root cell
- `upper`: the upper coordinate of the root cell
- `creation_time`: the time when the root cell was created during sampling
- `is_split`: whether the root cell is split
- `split_axis`: the direction in which the root cell is split, if any
- `split_location`: the location on `split_axis` at which the root cell is split, if any
- `tree_left`: the left child tree of the root cell, if any
- `tree_right`: the right child tree of the root cell, if any
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
    MondrianTree(id::String, lambda::Float64, lower::NTuple{d,Float64},
                 upper::NTuple{d,Float64}, creation_time::Float64) where {d}

Sample a Mondrian tree with a given id, lifetime,
lower and upper cell coordinates, and creation time.
To be used in internal construction methods.

# Examples
```julia
id = ""
lambda = 3.0
lower = ntuple(i -> 0.2, d)
upper = ntuple(i -> 0.7, d)
creation_time = 0.0
tree = MondrianTree(id, lambda, lower, upper, creation_time)
```
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

"""
    MondrianTree(d::Int, lambda::Float64)

Sample a Mondrian tree \$\\mathcal{M}([0,1]^d, \\lambda)\$.

# Examples
```julia
d = 2
lambda = 3.0
tree = MondrianTree(d, lambda);
```
"""
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

Check if a point `x` is contained in the root cell of a Mondrian tree.

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

Get the center point of the root cell of a Mondrian tree.

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

Get the d-dimensional volume of the root cell of a Mondrian tree.

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

function apply_split(tree::MondrianTree{d}, split_tree::MondrianTree{d}) where {d}
    if tree.is_split
        tree_left = apply_split(tree.tree_left, split_tree)
        tree_right = apply_split(tree.tree_right, split_tree)
        return MondrianTree(tree.id, tree.lambda, tree.lower, tree.upper, tree.creation_time,
                            true, split_tree.split_axis, split_tree.split_location,
                            tree_left, tree_right)
    else
        if all(tree.lower .< split_tree.tree_left.upper) &&
           all(split_tree.tree_right.lower .< tree.upper)
            left_upper = min.(split_tree.tree_left.upper, tree.upper)
            right_lower = max.(split_tree.tree_right.lower, tree.lower)
            tree_left = MondrianTree(tree.id * "L", tree.lambda, tree.lower, left_upper,
                                     split_tree.tree_left.creation_time, false,
                                     nothing, nothing, nothing, nothing)
            tree_right = MondrianTree(tree.id * "R", tree.lambda, right_lower, tree.upper,
                                      split_tree.tree_left.creation_time, false,
                                      nothing, nothing, nothing, nothing)
            return MondrianTree(tree.id, tree.lambda, tree.lower, tree.upper, tree.creation_time,
                                true, split_tree.split_axis, split_tree.split_location,
                                tree_left, tree_right)
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
        tree = apply_split(tree, subtree)
    end

    return tree
end

"""
    get_common_refinement(trees::Vector{MondrianTree{d}}) where {d}

Get the common refinement of several Mondrian trees,
whose leaf cells are the intersections of all leaf cells in `trees`.
Preserves the split times and
returns a new equivalent `MondrianTree`.

# Examples
```julia
trees = [MondrianTree(2, 3.0) for _ in 1:3]
refined_tree = get_common_refinement(trees)
```
"""
function get_common_refinement(trees::Vector{MondrianTree{d}}) where {d}
    @assert !isempty(trees)
    if length(trees) == 1
        return trees[]
    else
        refinement = get_common_refinement(trees[1], trees[2])
        return get_common_refinement([trees[3:end]; refinement])
    end
end

"""
    are_in_same_leaf(x1::NTuple{d,Float64}, x2::NTuple{d,Float64}, tree::MondrianTree)

Check if two points are in the same leaf cell of a Mondrian tree.

# Examples
```jldoctest
d = 2
tree = MondrianTree(d, 0.0)
x1 = ntuple(i -> 0.2, d)
x2 = ntuple(i -> 0.7, d)
are_in_same_leaf(x1, x2, tree)

# output
true
```
"""
function are_in_same_leaf(x1::NTuple{d,Float64}, x2::NTuple{d,Float64},
                          tree::MondrianTree) where {d}
    if tree.is_split
        if is_in(x1, tree.tree_left) && is_in(x2, tree.tree_left)
            return are_in_same_leaf(x1, x2, tree.tree_left)
        elseif is_in(x1, tree.tree_right) && is_in(x2, tree.tree_right)
            return are_in_same_leaf(x1, x2, tree.tree_right)
        else
            return false
        end
    else
        return is_in(x1, tree) && is_in(x2, tree)
    end
end

"""
    get_subtrees(tree::MondrianTree{d}) where {d}

Get a list of the subtrees contained in a Mondrian tree.

# Examples
```julia
tree = MondrianTree(2, 3.0)
get_subtrees(tree)
```
"""
function get_subtrees(tree::MondrianTree{d}) where {d}
    if tree.is_split
        subtrees_left = get_subtrees(tree.tree_left)
        subtrees_right = get_subtrees(tree.tree_right)
        return [tree; subtrees_left; subtrees_right]
    else
        return [tree]
    end
end

"""
    get_leaves(tree::MondrianTree{d}) where {d}

Get a list of the leaves in a Mondrian tree.

# Examples
```julia
tree = MondrianTree(2, 3.0)
get_leaves(tree)
```
"""
function get_leaves(tree::MondrianTree{d}) where {d}
    return [t for t in get_subtrees(tree) if !t.is_split]
end

"""
    get_leaf_containing(x::NTuple{d,Float64}, tree::MondrianTree{d}) where {d}

Get the leaf of a Mondrian tree containing a point `x`.

# Examples
```julia
d = 2
x = ntuple(i -> 0.2, d)
tree = MondrianTree(d, 3.0)
get_leaf_containing(x, tree)
```
"""
function get_leaf_containing(x::NTuple{d,Float64}, tree::MondrianTree{d}) where {d}
    return [t for t in get_leaves(tree) if is_in(x, t)][]
end

"""
    count_leaves(tree::MondrianTree{d}) where {d}

Count the leaves of a Mondrian tree.

# Examples
```julia
tree = MondrianTree(2, 3.0)
count_leaves(tree)
```
"""
function count_leaves(tree::MondrianTree{d}) where {d}
    return length(get_leaves(tree))
end

"""
    restrict(tree::MondrianTree{d}, time::Float64) where {d}

Restrict a Mondrian tree to a stopping time.

# Examples
```julia
tree = MondrianTree(2, 3.0)
restrict(tree, 2.0)
```
"""
function restrict(tree::MondrianTree{d}, time::Float64) where {d}
    if tree.is_split && tree.tree_left.creation_time <= time
        tree_left = restrict(tree.tree_left, time)
        tree_right = restrict(tree.tree_right, time)
        return MondrianTree{d}(tree.id, tree.lambda, tree.lower, tree.upper, tree.creation_time,
                               true, tree.split_axis, tree.split_location, tree_left, tree_right)
    else
        return MondrianTree{d}(tree.id, tree.lambda, tree.lower, tree.upper, tree.creation_time,
                               false, nothing, nothing, nothing, nothing)
    end
end

"""
    Base.show(tree::MondrianTree{d}) where {d}

Show the recursive structure of a Mondrian tree.
"""
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
