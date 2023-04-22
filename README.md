# InterfaceLimits.jl

A basic Julia module that constructs a JuMP optimization model to find the maximum active
power transfer limit between interconnected areas in a power system. The module depends on
[PowerSystems.jl](https://github.com/nrel-siip/powersystems.jl) and accepts data in the form
of a `System`. Interface limits are determined for every pair of `Area`s that are directly
connected with an `ACBranch`.
## installation

```julia
julia> ]
(tmp) pkg> add https://github.nrel.gov/SIIP/InterfaceLimits.jl
```

## usage

```julia
using PowerSystems
using InterfaceLimits
sys = System("matpower.m")
limits = find_interface_limits(sys)
```

## demo

See the demo script in `examples/rts/rts-interface-demo.jl`

```julia
julia> using InterfaceLimits
julia> include(joinpath(dirname(dirname(pathof(InterfaceLimits))), "examples", "rts", "rts-interface-demo.jl"))
┌ Info: Solving interface limit problem with
└   solver = HiGHS.Optimizer
Running HiGHS 1.5.1 [date: 1970-01-01, git hash: 93f1876e4]
Copyright (c) 2023 HiGHS under MIT licence terms
Presolving model
708 rows, 1086 cols, 45312 nonzeros
708 rows, 1074 cols, 43896 nonzeros
Presolve : Reductions: rows 708(-924); columns 1074(-36); elements 43896(-2362)
Solving the presolved LP
Using EKK dual simplex solver - serial
  Iteration        Objective     Infeasibilities num(sum)
          0    -4.3502751027e+03 Pr: 708(173423); Du: 0(3.31781e-09) 0s
        706    -4.3500000000e+03 Pr: 0(0) 0s
Solving the original LP from the solution after postsolve
Model   status      : Optimal
Simplex   iterations: 706
Objective value     :  4.3500000000e+03
HiGHS run time      :          0.03
┌ Info: Interface limits calculated
│   df =
│    6×3 DataFrame
│     Row │ interface  transfer_limit  sum_capacity
│         │ String     Float64         Float64?
│    ─────┼─────────────────────────────────────────
│       1 │ 1_3                 500.0         500.0
│      ⋮  │     ⋮            ⋮              ⋮
│       6 │ 3_2                 500.0         500.0
└                                     4 rows omitted
[ Info: calculating n-1 interface limits
[ Info: Building interface limit optimization model
┌ Info: Solving interface limit problem with
└   solver = HiGHS.Optimizer
Running HiGHS 1.5.1 [date: 1970-01-01, git hash: 93f1876e4]
Copyright (c) 2023 HiGHS under MIT licence terms
Presolving model
83544 rows, 83922 cols, 293820 nonzeros
83544 rows, 83910 cols, 292404 nonzeros
Presolve : Reductions: rows 83544(-90888); columns 83910(-3600); elements 292404(-95170)
Solving the presolved LP
Using EKK dual simplex solver - serial
  Iteration        Objective     Infeasibilities num(sum)
          0    -4.3501483614e+03 Pr: 83544(21371.6); Du: 0(3.5465e-09) 0s
      20469    -4.3501161239e+03 Pr: 67221(594593); Du: 0(1.70483e-06) 5s
      53855    -3.6281666180e+03 Pr: 33777(15917.3); Du: 0(4.40279e-07) 10s
      84487    -3.4112624969e+03 Pr: 0(0); Du: 0(1.23329e-10) 12s
Solving the original LP from the solution after postsolve
Model   status      : Optimal
Simplex   iterations: 84487
Objective value     :  3.4112624969e+03
HiGHS run time      :         12.25
┌ Info: Interface limits calculated
│   df =
│    6×3 DataFrame
│     Row │ interface  transfer_limit  sum_capacity
│         │ String     Float64         Float64?
│    ─────┼─────────────────────────────────────────
│       1 │ 1_3               491.993         500.0
│      ⋮  │     ⋮            ⋮              ⋮
│       6 │ 3_2               489.511         500.0
└                                     4 rows omitted
6×3 DataFrame
 Row │ interface  transfer_limit  sum_capacity
     │ String     Float64         Float64?
─────┼─────────────────────────────────────────
   1 │ 1_3               491.993         500.0
   2 │ 1_2               725.548        1175.0
   3 │ 2_3               489.369         500.0
   4 │ 3_1               489.293         500.0
   5 │ 2_1               725.548        1175.0
   6 │ 3_2               489.511         500.0
```
