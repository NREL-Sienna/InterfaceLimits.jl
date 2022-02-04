# interface-limits

A sandbox for developing models/methods to determine inter-regional transfer limits

## installation

```julia
using Pkg
Pkg.add("https://github.nrel.gov/SIIP/interface-limits")
```

## demo

```julia
using InterfaceLimits
include(joinpath(dirname(dirname(pathof(InterfaceLimits))),"5-bus-interface-demo.jl"))
```
