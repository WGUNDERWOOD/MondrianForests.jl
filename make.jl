using JuliaFormatter
using MondrianForests
using Pkg
using Coverage

# format
format(MondrianForests)

# test
Pkg.test(coverage=true)

# coverage
coverage = process_folder()
covered_lines, total_lines = get_summary(coverage)
println("Coverage: $covered_lines / $total_lines")

# docs
include("docs/make.jl")

# replication
include("replication/construction_diagrams/construction_diagrams.jl")
include("replication/logo/logo.jl")
include("replication/partition_plots/partition_plots.jl")
include("replication/piet_diagram/piet_diagram.jl")
include("replication/readme_examples/readme_examples.jl")
include("replication/theorem_diagrams/theorem_diagrams.jl")
include("replication/weather/weather.jl")
include("replication/weather/weather_cv.jl")
