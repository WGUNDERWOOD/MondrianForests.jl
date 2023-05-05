module MondrianForests

include("cell.jl")
include("data.jl")
include("forest.jl")
include("tree.jl")
include("lifetime.jl")

export generate_data
export generate_uniform_data
export MondrianCell
export MondrianTree
export MondrianForest
export fit
export show
export select_lifetime_global_polynomial

end
