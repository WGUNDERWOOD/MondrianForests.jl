push!(LOAD_PATH, "../src/")

using Documenter
using MondrianForests

DocMeta.setdocmeta!(MondrianForests, :DocTestSetup,
                    :(using MondrianForests); recursive=true)

makedocs(sitename="MondrianForests.jl",
         modules=MondrianForests,
         pages=["Home" => "index.md",
                "Documentation" => "documentation.md"])

deploydocs(repo="github.com/WGUNDERWOOD/MondrianForests.jl.git")
