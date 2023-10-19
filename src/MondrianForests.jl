module MondrianForests

# TODO better docs for exported functions

export MondrianCell
export is_in
export get_center
export get_volume
export show

# TODO tree
export MondrianTree
export get_cell_id
export are_in_same_cell
export count_cells
export restrict

# TODO forest
export MondrianForest
export fit

# TODO debias
export DebiasedMondrianForest

# TODO lifetime_gcv
export select_lifetime_gcv

# TODO lifetime_polynomial
export select_lifetime_polynomial

# TODO include source files
include("cell.jl")
include("tree.jl")
include("data.jl")
include("forest.jl")
include("debias.jl")
include("lifetime_polynomial.jl")
include("lifetime_gcv.jl")

end
