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

function MondrianCell(d::Int)
    if d < 0
        throw(DomainError(d, "dimension must be non-negative"))
    else
        lower = ntuple(i -> 0.0, d)
        upper = ntuple(i -> 1.0, d)
        return MondrianCell(lower, upper)
    end
end

function is_in(x::NTuple{d,Float64}, cell::MondrianCell) where {d}
    return all(cell.lower .<= x .<= cell.upper)
end

function get_center(cell::MondrianCell)
    return (cell.lower .+ cell.upper) ./ 2
end

function get_intersection(cell1::MondrianCell, cell2::MondrianCell)
    lower = max.(cell1.lower, cell2.lower)
    upper = min.(cell1.upper, cell2.upper)
    if all(lower .< upper)
        return MondrianCell(lower, upper)
    else
        return nothing
    end
end

function Base.show(cell::MondrianCell{d}) where {d}
    lower = round.(cell.lower, digits=4)
    upper = round.(cell.upper, digits=4)
    printstyled("$lower -- $upper ", color=:green)
    print("\n")
    return nothing
end
