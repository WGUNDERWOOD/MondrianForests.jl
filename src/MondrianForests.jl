module MondrianForests

include("data.jl")
include("tree.jl")
include("forest.jl")

export generate_data
export MondrianCell
export MondrianTree
export MondrianForest
export fit

end
