# This file was modified from P4est.jl
using Pkg
Pkg.add("MPIPreferences")
Pkg.add("Preferences")
Pkg.add("UUIDs")

@static if VERSION >= v"1.8"
    Pkg.compat("MPIPreferences", "0.1")
    Pkg.compat("Preferences", "1")
    Pkg.compat("UUIDs", "1")
end

const P4ESTTYPES_TEST = get(ENV, "P4ESTTYPES_TEST", "P4ESTTYPES_JLL_MPI_DEFAULT")
const P4ESTTYPES_TEST_LIBP4EST = get(ENV, "P4ESTTYPES_TEST_LIBP4EST", "")
const P4ESTTYPES_TEST_LIBSC = get(ENV, "P4ESTTYPES_TEST_LIBSC", "")

@static if P4ESTTYPES_TEST == "P4ESTTYPES_CUSTOM_MPI_CUSTOM"
    import MPIPreferences
    MPIPreferences.use_system_binary()
end

@static if P4ESTTYPES_TEST == "P4ESTTYPES_CUSTOM_MPI_CUSTOM"
    import UUIDs, Preferences
    Preferences.set_preferences!(
        UUIDs.UUID("7d669430-f675-4ae7-b43e-fab78ec5a902"), # UUID of P4est.jl
        "libp4est" => P4ESTTYPES_TEST_LIBP4EST,
        force = true,
    )
    Preferences.set_preferences!(
        UUIDs.UUID("7d669430-f675-4ae7-b43e-fab78ec5a902"), # UUID of P4est.jl
        "libsc" => P4ESTTYPES_TEST_LIBSC,
        force = true,
    )
end

@info "P4estTypes.jl tests configured" P4ESTTYPES_TEST P4ESTTYPES_TEST_LIBP4EST P4ESTTYPES_TEST_LIBSC
