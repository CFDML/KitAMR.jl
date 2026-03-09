using Documenter
push!(LOAD_PATH, "..")
using KitAMR
type_page = [
    "Configuration"=>"type_configure.md",
    "Condition"=>"type_condition.md",
    "Gas"=>"type_gas.md",
    "Output"=>"type_output.md"
]
makedocs(;
    sitename = "KitAMR.jl",
    modules = [KitAMR],
    pages = [
        "Home"=>"index.md",
        "Types"=>type_page
        ],
    checkdocs = :none
)

deploydocs(;devbranch = "ge",
    repo = "github.com/CFDML/KitAMR.jl.git",
)