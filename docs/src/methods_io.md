# IO
## Input
```@docs
read_config
```
```@docs
listen_for_save!
```
```@docs
check_for_save!
```
The periodic status report (and the `save` hook) is driven once per step by
```@docs
check!
```
!!! todo
    A bug exists here, which causes the `save` input will not be aware of, until one stdout is invoked.

## Output
```@docs
save_result
```
Animation frames are written incrementally during the run by
```@docs
check_for_animsave!
```

## Restart / checkpoint
Unlike [`save_result`](@ref) (which writes macroscopic fields for post-processing), a checkpoint
stores the full phase-space state so a run can be resumed — possibly on a different number of MPI
ranks. Write one with
```@docs
save_for_restart
```
and resume from it with
```@docs
restart
```
The configuration (including the initial-condition and user-defined functions) is saved with the
checkpoint; the same user-defined-function script must be loaded in the session that calls
`restart`, or a live `config` passed via the `config` keyword.
