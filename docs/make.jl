#push!(LOAD_PATH,"../src/")

using Documenter
using MondrianForests

makedocs(sitename="MondrianForests.jl",
         modules=MondrianForests,
         pages=["Home" => "index.md",
                "Documentation" => "documentation.md"])

#deploydocs(repo = "github.com/WGUNDERWOOD/MondrianForests.jl.git")
