# Short representative workload for building an effi_study sysimage.
#
# Keep this small: PackageCompiler needs the method/type paths, not a physically
# meaningful run.  The workload covers the shared vacuum-jet setup plus both
# efficiency-study both_on mode:
#   PS-AMR + VS-AMR + partition scheduler path.
# The same concrete solver/configuration types cover both_off, so avoid running
# the dense-reference mode during sysimage construction.

using KitAMR
using MPI

setdefault_env(k, v) = haskey(ENV, k) ? nothing : (ENV[k] = v)

setdefault_env("VJE_PS_BASE_X", "2")
setdefault_env("VJE_PS_BASE_Y", "4")
setdefault_env("VJE_PS_LMAX", "1")
setdefault_env("VJE_VS_BASE", "8")
setdefault_env("VJE_VS_LMAX", "3")
setdefault_env("VJE_MAX_SIM_TIME", "0.01")
setdefault_env("VJE_MAX_STEPS", "1")
setdefault_env("VJE_ANIM_DT", "0")
setdefault_env("VJE_PS_INTERVAL", "1")
setdefault_env("VJE_VS_INTERVAL", "1")
setdefault_env("VJE_PARTITION_INTERVAL", "1")
setdefault_env("VJE_OUTDIR", mktempdir())

include("common.jl")

function run_precompile_mode(mode)
    spec = R2DVJEffi.mode_spec(mode)
    config = R2DVJEffi.build_config(spec)
    p4est = nothing
    ka = nothing
    try
        p4est, ka = initialize(config;
            prerefine_steps = spec.prerefine_steps,
            prerefine_reinit_ic = true)
        solve!(p4est, ka;
            ps_interval = R2DVJEffi.PS_INTERVAL,
            vs_interval = R2DVJEffi.VS_INTERVAL,
            partition_interval = R2DVJEffi.PARTITION_INTERVAL,
            vs_balance = true,
            break_on_convergence = false,
            listen_for_save = false,
            animation = false,
            status_check = false,
            progress = false,
            max_steps = 1)
    finally
        if p4est !== nothing && ka !== nothing
            KitAMR.finalize!(p4est, ka)
        end
    end
    return nothing
end

R2DVJEffi.with_mpi() do
    run_precompile_mode("both_on")
    MPI.Comm_rank(MPI.COMM_WORLD) == 0 &&
        println("effi_study precompile workload complete.")
end
