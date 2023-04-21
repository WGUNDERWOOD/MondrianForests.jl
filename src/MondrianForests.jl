module MondrianForests

include("cell.jl")
include("data.jl")
include("tree.jl")
include("forest.jl")
include("plot.jl")

export generate_data
export MondrianCell
export MondrianTree
export MondrianForest
export fit
export show

end
