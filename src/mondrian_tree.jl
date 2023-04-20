using Random

struct MondrianCell{d}
    left::NTuple{d, Float64}
    right::NTuple{d, Float64}
end

function MondrianCell(d::Int)
    left = ntuple(i -> 0.0, d)
    right = ntuple(i -> 1.0, d)
    return MondrianCell(left, right)
end

struct MondrianTree{d}
    lambda::Float64
    cell::MondrianCell
    left::Union{MondrianTree, Nothing}
    right::Union{MondrianTree, Nothing}
    time::Float64
end

function MondrianTree(cell::MondrianCell, time::Float64, lambda::Float64)
    size_cell = sum(cell.right .- cell.left)
    E = randexp() / size_cell

    if time + E <= lambda
        j = rand(DiscreteNonParametric(1:d, (c.b .- c.a) / dim_c))
        sj = rand(Uniform(c.a[j], c.b[j]))
        b0 = copy(c.b)
        b0[j] = sj
        c0 = Cell(c.a, b0)
        a1 = copy(c.a)
        a1[j] = sj
        c1 = Cell(a1, c.b)
        m0 = sample_mondrian(c0, time + E, lambda)
        m1 = sample_mondrian(c1, time + E, lambda)
        return [m0; m1]
    else
        return [c]
    end
end

