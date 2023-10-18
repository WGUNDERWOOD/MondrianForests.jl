# TODO document all functions in the package

"""Struct to represent a Mondrian cell."""
struct MondrianCell{d}
    lower::NTuple{d,Float64}
    upper::NTuple{d,Float64}

    function MondrianCell(lower::NTuple{d,Float64}, upper::NTuple{d,Float64}) where {d}
        if !all(0 .<= lower .<= 1)
            throw(DomainError(lower, "coordinates of lower must be in [0, 1]"))
        elseif !all(0 .<= upper .<= 1)
            throw(DomainError(upper, "coordinates of upper must be in [0, 1]"))
        elseif !all(lower .<= upper)
            throw(ArgumentError("lower must be pointwise no greater than upper"))
        else
            return new{d}(lower, upper)
        end
    end
end

"""Construct the cell `[0,1]^d`."""
function MondrianCell(d::Int)
    if d < 0
        throw(DomainError(d, "dimension must be non-negative"))
    else
        lower = ntuple(i -> 0.0, d)
        upper = ntuple(i -> 1.0, d)
        return MondrianCell(lower, upper)
    end
end

"""Construct the cell `[0,1]^d`."""
function is_in(x::NTuple{d,Float64}, cell::MondrianCell) where {d}
    return all(cell.lower .<= x .<= cell.upper)
end

"""Get the center point of a cell."""
function get_center(cell::MondrianCell)
    return (cell.lower .+ cell.upper) ./ 2
end

"""Get the volume of a cell."""
function get_volume(cell::MondrianCell)
    return prod(cell.upper .- cell.lower)
end

"""Get the intersection of two cells."""
function get_intersection(cell1::MondrianCell{d}, cell2::MondrianCell{d}) where {d}
    lower = max.(cell1.lower, cell2.lower)
    upper = min.(cell1.upper, cell2.upper)
    if all(lower .< upper)
        return MondrianCell(lower, upper)
    else
        return nothing
    end
end

"""Get the common refinement of two sets of cells."""
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

"""Get the common refinement of many sets of cells."""
function get_common_refinement(cells::Vector{Vector{MondrianCell{d}}}) where {d}
    if length(cells) == 1
        return cells[1]
    else
        return get_common_refinement(cells[1], get_common_refinement(cells[2:end]))
    end
end

"""Show a cell."""
function Base.show(cell::MondrianCell{d}) where {d}
    lower = round.(cell.lower, digits=4)
    upper = round.(cell.upper, digits=4)
    printstyled("$lower -- $upper ", color=:green)
    print("\n")
    return nothing
end
