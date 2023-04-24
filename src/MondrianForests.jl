module MondrianForests

include("cell.jl")
include("data.jl")
include("forest.jl")
include("tree.jl")

export generate_data
export generate_uniform_data
export MondrianCell
export MondrianTree
export MondrianForest
export fit

end
