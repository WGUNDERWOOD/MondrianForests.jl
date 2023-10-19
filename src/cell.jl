"""
    MondrianCell(id::String, lower::NTuple{d,Float64}, upper::NTuple{d,Float64}) where {d}

A Mondrian cell is determined by the coordinates of its lower and upper corner points.
The dimension `d` may be any positive integer.

# Examples
```jldoctest
lower = ntuple(i -> 0.2, 2)
upper = ntuple(i -> 0.7, 2)
MondrianCell("", lower, upper)

# output
MondrianCell{2}("", (0.2, 0.2), (0.7, 0.7))
```
"""
struct MondrianCell{d}
    id::String
    lower::NTuple{d,Float64}
    upper::NTuple{d,Float64}

    function MondrianCell(id::String, lower::NTuple{d,Float64}, upper::NTuple{d,Float64}) where {d}
        if !all(0 .<= lower .<= 1)
            throw(DomainError(lower, "coordinates of lower must be in [0, 1]"))
        elseif !all(0 .<= upper .<= 1)
            throw(DomainError(upper, "coordinates of upper must be in [0, 1]"))
        elseif !all(lower .<= upper)
            throw(ArgumentError("lower must be pointwise no greater than upper"))
        else
            return new{d}(id, lower, upper)
        end
    end
end

"""
    MondrianCell(d::Int)

Construct the d-dimensional unit hypercube Mondrian cell `[0,1]^d`.
Uses an empty string for the `id` field.

```jldoctest
MondrianCell(2)

# output
MondrianCell{2}("", (0.0, 0.0), (1.0, 1.0))
```
"""
function MondrianCell(d::Int)
    if d < 0
        throw(DomainError(d, "dimension must be non-negative"))
    else
        lower = ntuple(i -> 0.0, d)
        upper = ntuple(i -> 1.0, d)
        return MondrianCell("", lower, upper)
    end
end

"""
    is_in(x::NTuple{d,Float64}, cell::MondrianCell{d}) where {d}

Check if a point `x` is contained in a Mondrian cell.

# Examples
```jldoctest
x = ntuple(i -> 0.2, 2)
cell = MondrianCell(2)
is_in(x, cell)

# output
true
```
"""
function is_in(x::NTuple{d,Float64}, cell::MondrianCell{d}) where {d}
    return all(cell.lower .<= x .<= cell.upper)
end

"""
    get_center(cell::MondrianCell{d}) where {d}

Get the center point of a Mondrian cell.

# Examples
```jldoctest
cell = MondrianCell(2)
get_center(cell)

# output
(0.5, 0.5)
```
"""
function get_center(cell::MondrianCell{d}) where {d}
    return (cell.lower .+ cell.upper) ./ 2
end

"""
    get_volume(cell::MondrianCell{d}) where {d}

Get the d-dimensional volume of a Mondrian cell.

# Examples
```jldoctest
cell = MondrianCell(2)
get_volume(cell)

# output
1.0
```
"""
function get_volume(cell::MondrianCell{d}) where {d}
    return prod(cell.upper .- cell.lower)
end

function get_intersection(cell1::MondrianCell{d}, cell2::MondrianCell{d}) where {d}
    lower = max.(cell1.lower, cell2.lower)
    upper = min.(cell1.upper, cell2.upper)
    if all(lower .< upper)
        return MondrianCell(lower, upper)
    else
        return nothing
    end
end

function get_common_refinement(cells1::Vector{MondrianCell{d}},
                               cells2::Vector{MondrianCell{d}}) where {d}
    cells = MondrianCell{d}[]
    for c1 in cells1
        for c2 in cells2
            c = get_intersection(c1, c2)
            if isa(c, MondrianCell{d})
                push!(cells, c)
            end
        end
    end
    return unique(cells)
end

function get_common_refinement(cells::Vector{Vector{MondrianCell{d}}}) where {d}
    if length(cells) == 1
        return cells[1]
    else
        return get_common_refinement(cells[1], get_common_refinement(cells[2:end]))
    end
end

"""
    show(cell::MondrianCell{d}) where {d}

Show a Mondrian cell.
"""
function Base.show(cell::MondrianCell{d}) where {d}
    lower = round.(cell.lower, digits=4)
    upper = round.(cell.upper, digits=4)
    printstyled("$lower -- $upper ", color=:green)
    print("\n")
    return nothing
end
