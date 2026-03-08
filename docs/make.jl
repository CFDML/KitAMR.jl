using Documenter
push!(LOAD_PATH, "..")
using KitAMR
type_page = ["Configuration"=>"type_configure.md"]
makedocs(;
    sitename = "KitAMR.jl",
    modules = [KitAMR],
    pages = [
        "Home"=>"index.md",
        "Types"=>type_page
        ],
    checkdocs = :none
)

deploydocs(
    repo = "github.com/CFDML/KitAMR.jl.git",
)