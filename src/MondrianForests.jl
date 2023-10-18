module MondrianForests

# TODO better docs for exported functions

# cell
export MondrianCell
export is_in
export get_center
export get_volume
export show

# tree
export MondrianTree
export get_cell_id
export are_in_same_cell
export count_cells
export restrict

# forest
export MondrianForest
export fit

# debias
export DebiasedMondrianForest

# lifetime_gcv
export select_lifetime_gcv

# lifetime_polynomial
export select_lifetime_polynomial

# include source files
include("cell.jl")
include("tree.jl")
include("data.jl")
include("forest.jl")
include("debias.jl")
include("lifetime_polynomial.jl")
include("lifetime_gcv.jl")

end
