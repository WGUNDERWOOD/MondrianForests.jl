module MondrianForests

# TODO better docs for exported functions

# TODO tree
export MondrianTree
export is_in
export get_center
export get_volume
export show
export get_subtrees
export get_leaves
export get_common_refinement
#export get_cell_id
#export are_in_same_cell
#export count_cells
#export restrict

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
include("tree.jl")
include("data.jl")
include("forest.jl")
include("debias.jl")
include("lifetime_polynomial.jl")
include("lifetime_gcv.jl")

end
