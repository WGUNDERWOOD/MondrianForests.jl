struct MondrianCell{d}
    lower::NTuple{d,Float64}
    upper::NTuple{d,Float64}
end

function MondrianCell(d::Int)
    lower = ntuple(i -> 0.0, d)
    upper = ntuple(i -> 1.0, d)
    return MondrianCell(lower, upper)
end

function is_in(x::NTuple{d,Float64}, cell::MondrianCell) where {d}
    for i in 1:d
        if (cell.lower[i] >= x[i]) || (x[i] > cell.upper[i])
            return false
        end
    end
    return true
end

function sample_mondrian_cell(x::NTuple{d,Float64}, lambda::Float64) where {d}
    E_lower = [rand(Exponential(1)) for _ in 1:d]
    E_upper = [rand(Exponential(1)) for _ in 1:d]
    lower = max.(x .- E_lower ./ lambda, 0)
    upper = min.(x .+ E_upper ./ lambda, 1)
    lower = ntuple(i -> lower[i], d)
    upper = ntuple(i -> upper[i], d)
    return MondrianCell(lower, upper)
end

function Base.show(cell::MondrianCell{d}) where {d}
    lower = round.(cell.lower, digits=4)
    upper = round.(cell.upper, digits=4)
    printstyled("$lower -- $upper ", color=:green)
    print("\n")
    return nothing
end
