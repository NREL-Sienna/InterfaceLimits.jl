# interface-limits

A sandbox for developing models/methods to determine inter-regional transfer limits

## installation

```julia
julia> ]
(tmp) pkg> add https://github.nrel.gov/SIIP/InterfaceLimits.jl
```

## demo

```julia
using InterfaceLimits
include(joinpath(dirname(dirname(pathof(InterfaceLimits))),"5-bus-interface-demo.jl"))
```
