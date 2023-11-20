using MondrianForests
using Test
using Random
using Distributions
using Suppressor
using Aqua
using JuliaFormatter

Aqua.test_ambiguities(MondrianForests)
Aqua.test_unbound_args(MondrianForests)
Aqua.test_undefined_exports(MondrianForests)
Aqua.test_project_extras(MondrianForests)
Aqua.test_stale_deps(MondrianForests, ignore=[:Aqua, :Suppressor, :JuliaFormatter])
Aqua.test_deps_compat(MondrianForests)

@testset verbose = true "JuliaFormatter" begin
    @test JuliaFormatter.format(MondrianForests, overwrite=false)
end

include("test_tree.jl")
include("test_forest.jl")
include("test_debias.jl")
include("test_lifetime_polynomial.jl")
include("test_lifetime_gcv.jl")
include("test_data.jl")
