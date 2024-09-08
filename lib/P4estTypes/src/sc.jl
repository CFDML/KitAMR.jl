module SC

using MPI
using ...P4est

sc_package_id() = unsafe_load(cglobal((:sc_package_id, P4est.LibP4est.libsc), Cint))
initialized() = sc_package_id() >= 0
finalize() = sc_finalize_noabort()

module LP
using ....P4est

"""
    P4estTypes.SC.LP.LogPriority

An abstract type for p4est and libsc log priorities.  The follow
priorities are available:

  - [`P4estTypes.SC.LP.Default`](@ref Default): the libsc default.
  - [`P4estTypes.SC.LP.Always`](@ref Always): log absolutely everything.
  - [`P4estTypes.SC.LP.Trace`](@ref Trace): prefix file and line number.
  - [`P4estTypes.SC.LP.Debug`](@ref Debug): any information on the internal state.
  - [`P4estTypes.SC.LP.Verbose`](@ref Verbose): information on conditions, decisions.
  - [`P4estTypes.SC.LP.Info`](@ref Info): most relevant things a function is doing.
  - [`P4estTypes.SC.LP.Statistics`](@ref Statistics): important for consistency/performance.
  - [`P4estTypes.SC.LP.Essential`](@ref Essential): a few lines at most for a major api function.
  - [`P4estTypes.SC.LP.Production`](@ref Production): log a few lines max per program.
  - [`P4estTypes.SC.LP.Error`](@ref Error): log errors only.
  - [`P4estTypes.SC.LP.Silent`](@ref Silent): never log anything.
"""
abstract type LogPriority end


"""
    P4estTypes.SC.LP.Always

Log priority indicating to log absolutely everything.
"""
struct Always <: LogPriority end
"""
    P4estTypes.SC.LP.Debug

Log priority indicating to log any information on the internal state.
"""
struct Debug <: LogPriority end
"""
    P4estTypes.SC.LP.Default

The libsc default log priority.
"""
struct Default <: LogPriority end
"""
    P4estTypes.SC.LP.Error

Log priority indicating to log errors only.
"""
struct Error <: LogPriority end
"""
    P4estTypes.SC.LP.Essential

Log priority indicating to log a few lines at most for a major api function.
"""
struct Essential <: LogPriority end
"""
    P4estTypes.SC.LP.Info

Log priority indicating to log most relevant things a function is doing.
"""
struct Info <: LogPriority end
"""
    P4estTypes.SC.LP.Production

Log priority indicating to log a few lines max per program.
"""
struct Production <: LogPriority end
"""
    P4estTypes.SC.LP.Silent

Log priority indicating to never log anything.
"""
struct Silent <: LogPriority end
"""
    P4estTypes.SC.LP.Statistics

Log priority indicating to log important for consistency/performance.
"""
struct Statistics <: LogPriority end
"""
    P4estTypes.SC.LP.Trace

Log priority indicating to log trace information with prefix file and
line number.
"""
struct Trace <: LogPriority end
"""
    P4estTypes.SC.LP.Verbose

Log priority indicating to log information on conditions, decisions.
"""
struct Verbose <: LogPriority end

Base.convert(::Type{Cint}, ::Always) = Cint(SC_LP_ALWAYS)
Base.convert(::Type{Cint}, ::Debug) = Cint(SC_LP_DEBUG)
Base.convert(::Type{Cint}, ::Default) = Cint(SC_LP_DEFAULT)
Base.convert(::Type{Cint}, ::Error) = Cint(SC_LP_ERROR)
Base.convert(::Type{Cint}, ::Essential) = Cint(SC_LP_ESSENTIAL)
Base.convert(::Type{Cint}, ::Info) = Cint(SC_LP_INFO)
Base.convert(::Type{Cint}, ::Production) = Cint(SC_LP_PRODUCTION)
Base.convert(::Type{Cint}, ::Silent) = Cint(SC_LP_SILENT)
Base.convert(::Type{Cint}, ::Statistics) = Cint(SC_LP_STATISTICS)
Base.convert(::Type{Cint}, ::Trace) = Cint(SC_LP_TRACE)
Base.convert(::Type{Cint}, ::Verbose) = Cint(SC_LP_VERBOSE)

# struct SCArray{T} <: AbstractArray{T,1}
#     pointer::Ptr{sc_array}
# end
# Base.IndexStyle(::SCArray) = IndexLinear()
# function Base.getindex(a::SCArray{T}, i::Int) where {T}
#     return GC.@preserve a unsafe_load(T(unsafe_load(a.pointer).array), i)
# end
# Base.size(a::SCArray) = GC.@preserve a (unsafe_load(a.pointer).elem_count,)

end # module

"""
    P4estTypes.SC.setverbosity(logpriority::P4estTypes.SC.LP.LogPriority)

Sets the verbosity of libsc. See [`P4estTypes.SC.LP.LogPriority`](@ref P4estTypes.SC.LP.LogPriority)
for a list of valid priorities.
"""
setverbosity(lp::SC.LP.LogPriority) = sc_package_set_verbosity(sc_package_id(), lp)

end # module
