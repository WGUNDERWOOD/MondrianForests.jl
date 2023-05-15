module MondrianForests

include("cell.jl")
export MondrianCell

include("tree.jl")
export MondrianTree

include("data.jl")
export generate_data
export generate_uniform_data

include("forest.jl")
export MondrianForest
export fit
export show

include("lifetime.jl")
export select_lifetime_global_polynomial

end
