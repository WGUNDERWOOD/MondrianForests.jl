using Random
using Distributions

struct MondrianTree{d}
    lambda::Float64
    id::String
    creation_time::Float64
    cell::MondrianCell{d}
    split_axis::Union{Int,Nothing}
    split_location::Union{Float64,Nothing}
    tree_left::Union{MondrianTree{d},Nothing}
    tree_right::Union{MondrianTree{d},Nothing}
end

function MondrianTree(lambda::Float64, id::String, creation_time::Float64,
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
        cell_left = MondrianCell(cell.lower, left_upper)
        right_lower = collect(cell.lower)
        right_lower[split_axis] = split_location
        right_lower = ntuple(i -> right_lower[i], d)
        cell_right = MondrianCell(right_lower, cell.upper)
        tree_left = MondrianTree(lambda, id * "L", creation_time + E, cell_left)
        tree_right = MondrianTree(lambda, id * "R", creation_time + E, cell_right)
        tree = MondrianTree(lambda, id, creation_time, cell, split_axis,
                            split_location, tree_left, tree_right)
    else
        tree = MondrianTree(lambda, id, creation_time, cell, nothing,
                            nothing, nothing, nothing)
    end
    return tree
end

function MondrianTree(d::Int, lambda::Float64)
    if lambda < 0
        throw(DomainError(lambda, "lambda must be non-negative"))
    else
        return MondrianTree(lambda, "", 0.0, MondrianCell(d))
    end
end

function get_cell_id(x::NTuple{d,Float64}, tree::MondrianTree{d}) where {d}
    if isnothing(tree.split_axis)
        return tree.id
    else
        if x[tree.split_axis] <= tree.split_location
            return get_cell_id(x, tree.tree_left)
        else
            return get_cell_id(x, tree.tree_right)
        end
    end
end

function are_in_same_cell(x1::NTuple{d,Float64}, x2::NTuple{d,Float64},
                          tree::MondrianTree{d}) where {d}
    if isnothing(tree.split_axis)
        return true
    else
        axis = tree.split_axis
        location = tree.split_location
        x1_left = (x1[axis] <= location)
        x2_left = (x2[axis] <= location)

        if x1_left != x2_left
            return false
        elseif x1_left
            return are_in_same_cell(x1, x2, tree.tree_left)
        else
            return are_in_same_cell(x1, x2, tree.tree_right)
        end
    end
end

function count_cells(tree::MondrianTree{d}) where {d}
    if isnothing(tree.split_axis)
        return 1
    else
        return count_cells(tree.tree_left) + count_cells(tree.tree_right)
    end
end

function get_splits(tree::MondrianTree{d}) where {d}
    if !isnothing(tree.split_axis)
        lower = tree.tree_right.cell.lower
        upper = tree.tree_left.cell.upper
        return [(lower, upper); get_splits(tree.tree_left); get_splits(tree.tree_right)]
    else
        return Tuple{NTuple{d,Float64},NTuple{d,Float64}}[]
    end
end

function get_cells(tree::MondrianTree{d}) where {d}
    if !isnothing(tree.split_axis)
        return [get_cells(tree.tree_left); get_cells(tree.tree_right)]
    else
        return [tree.cell]
    end
end

function get_ids(tree::MondrianTree{d}) where {d}
    if !isnothing(tree.split_axis)
        return [get_ids(tree.tree_left); get_ids(tree.tree_right)]
    else
        return [tree.id]
    end
end

function get_split_times(tree::MondrianTree{d}) where {d}
    if !isnothing(tree.split_axis)
        return [tree.creation_time; get_split_times(tree.tree_left); get_split_times(tree.tree_right)]
    else
        return [tree.creation_time]
    end
end

function restrict(tree::MondrianTree{d}, time::Float64) where {d}
    if !isnothing(tree.split_axis)
        if tree.tree_left.creation_time > time
            new_tree = MondrianTree{d}(tree.lambda, tree.id, tree.creation_time, tree.cell,
                                       nothing, nothing, nothing, nothing)
            return new_tree
        else
            new_tree = MondrianTree{d}(tree.lambda, tree.id, tree.creation_time, tree.cell,
                                       tree.split_axis, tree.split_location,
                                       restrict(tree.tree_left, time),
                                       restrict(tree.tree_right, time))
            return new_tree
        end
    else
        return tree
    end
end

function Base.show(tree::MondrianTree{d}) where {d}
    depth = length(tree.id)
    has_split = !isnothing(tree.split_axis)
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
    if has_split
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
    if has_split
        show(tree.tree_left)
        show(tree.tree_right)
    end
    return nothing
end
