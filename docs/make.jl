using Documenter
push!(LOAD_PATH, "..")
using KitAMR
type_page = [
    "Configuration"=>"type_configure.md",
    "Data"=>"type_data.md",
    "Condition"=>"type_condition.md",
    "Gas"=>"type_gas.md",
    "Physical space"=>"type_physical_space.md",
    "Velocity space"=>"type_velocity_space.md",
    "Output"=>"type_output.md",
]
methods_page = [
    "Running a simulation"=>"methods_solve.md",
    "Initialization"=>"methods_initialize.md",
    "Reconstruction"=>"methods_reconstruct.md",
    "Flux"=>"methods_flux.md",
    "Iteration"=>"methods_iterate.md",
    "AMR"=>"methods_AMR.md",
    "Parallel"=>"methods_parallel.md",
    "Theory"=>"methods_theory.md",
    "Finalization"=>"methods_finalize.md",
    "IO"=>"methods_io.md",
]
makedocs(;
    sitename = "KitAMR.jl",
    modules = [KitAMR],
    pages = [
        "Home"=>"index.md",
        "Types"=>type_page,
        "Functions and Methods"=>methods_page,
        "User-defined functions"=>"user_defined.md",
        "Tutorial"=>"tutorial.md",
        "Cluster deployment"=>"cluster_deployment.md",
        "Index"=>"index_internal.md",
        "Limitations"=>"limitations.md",
    ],
    # format = Documenter.HTML(; collapselevel=1),
    checkdocs = :none,
    format = Documenter.HTML(),
)

deploydocs(; devbranch = "main", repo = "github.com/CFDML/KitAMR.jl.git")
