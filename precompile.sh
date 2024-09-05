echo Triggering precompilation
julia --project=. -e 'using Pkg; Pkg.precompile(); println("Precompilation complete");'