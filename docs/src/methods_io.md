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