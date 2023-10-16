module MondrianForests

export MondrianCell
export MondrianTree
export generate_data
export generate_uniform_data
export MondrianForest
export DebiasedMondrianForest
export fit
export show
export select_lifetime_global_polynomial

include("cell.jl")
include("tree.jl")
include("data.jl")
include("forest.jl")
include("debias.jl")
include("lifetime.jl")

end
