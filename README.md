# InterfaceLimits.jl

A basic Julia module that constructs a JuMP optimization model to find the maximum active
power transfer limit between interconnected areas in a power system. The module depends on
[PowerSystems.jl](https://github.com/nrel-siip/powersystems.jl) and accepts data in the form
of a `System`. Interface limits are determined for every pair of `Area`s that are directly
connected with an `ACBranch`.

## installation

```julia
julia> ]
pkg> add https://github.com/NREL-Sienna/InterfaceLimits.jl
```

## usage

```julia
julia> using PowerSystems
julia> using InterfaceLimits
julia> using HiGHS
julia> solver = optimizer_with_attributes(HiGHS.Optimizer);

julia> sys = System("matpower.m") # load a PowerSystems.jl System

julia> limits = find_interface_limits(sys, solver)
┌ Info:  1/3 Building interface limit optimization model
└   interface_key = "1" => "3"
Running HiGHS 1.5.1 [date: 1970-01-01, git hash: 93f1876e4]
Copyright (c) 2023 HiGHS under MIT licence terms
┌ Info: Solving interface limit problem with
└   solver = MathOptInterface.OptimizerWithAttributes(HiGHS.Optimizer, Pair{MathOptInterface.AbstractOptimizerAttribute, Any}[])
Presolving model
236 rows, 362 cols, 15104 nonzeros
236 rows, 358 cols, 14632 nonzeros
Presolve : Reductions: rows 236(-308); columns 358(-12); elements 14632(-786)
Solving the presolved LP
Using EKK dual simplex solver - serial
  Iteration        Objective     Infeasibilities num(sum)
          0    -1.0000323866e+03 Pr: 236(54595.1); Du: 0(1.03133e-09) 0s
        236    -1.0000000000e+03 Pr: 0(0) 0s
Solving the original LP from the solution after postsolve
Model   status      : Optimal
Simplex   iterations: 236
Objective value     :  1.0000000000e+03
HiGHS run time      :          0.01
df = leftjoin(df, interface_cap, on = :interface) = 2×3 DataFrame
 Row │ interface  transfer_limit  sum_capacity
     │ String     Float64         Float64?
─────┼─────────────────────────────────────────
   1 │ 1_3                 500.0         500.0
   2 │ 3_1                 500.0         500.0
┌ Info:  2/3 Building interface limit optimization model
└   interface_key = "1" => "2"
Running HiGHS 1.5.1 [date: 1970-01-01, git hash: 93f1876e4]
Copyright (c) 2023 HiGHS under MIT licence terms
┌ Info: Solving interface limit problem with
└   solver = MathOptInterface.OptimizerWithAttributes(HiGHS.Optimizer, Pair{MathOptInterface.AbstractOptimizerAttribute, Any}[])
Presolving model
236 rows, 362 cols, 15104 nonzeros
236 rows, 358 cols, 14632 nonzeros
Presolve : Reductions: rows 236(-308); columns 358(-12); elements 14632(-790)
Solving the presolved LP
Using EKK dual simplex solver - serial
  Iteration        Objective     Infeasibilities num(sum)
          0    -2.3501410037e+03 Pr: 236(54595.1); Du: 0(1.03133e-09) 0s
        243    -2.3500000000e+03 Pr: 0(0) 0s
Solving the original LP from the solution after postsolve
Model   status      : Optimal
Simplex   iterations: 243
Objective value     :  2.3500000000e+03
HiGHS run time      :          0.01
df = leftjoin(df, interface_cap, on = :interface) = 2×3 DataFrame
 Row │ interface  transfer_limit  sum_capacity
     │ String     Float64         Float64?
─────┼─────────────────────────────────────────
   1 │ 1_2                1175.0        1175.0
   2 │ 2_1                1175.0        1175.0
┌ Info:  3/3 Building interface limit optimization model
└   interface_key = "2" => "3"
Running HiGHS 1.5.1 [date: 1970-01-01, git hash: 93f1876e4]
Copyright (c) 2023 HiGHS under MIT licence terms
┌ Info: Solving interface limit problem with
└   solver = MathOptInterface.OptimizerWithAttributes(HiGHS.Optimizer, Pair{MathOptInterface.AbstractOptimizerAttribute, Any}[])
Presolving model
236 rows, 362 cols, 15104 nonzeros
236 rows, 358 cols, 14632 nonzeros
Presolve : Reductions: rows 236(-308); columns 358(-12); elements 14632(-786)
Solving the presolved LP
Using EKK dual simplex solver - serial
  Iteration        Objective     Infeasibilities num(sum)
          0    -1.0000329296e+03 Pr: 236(54595.1); Du: 0(1.03133e-09) 0s
        251    -1.0000000000e+03 Pr: 0(0) 0s
Solving the original LP from the solution after postsolve
Model   status      : Optimal
Simplex   iterations: 251
Objective value     :  1.0000000000e+03
HiGHS run time      :          0.01
df = leftjoin(df, interface_cap, on = :interface) = 2×3 DataFrame
 Row │ interface  transfer_limit  sum_capacity
     │ String     Float64         Float64?
─────┼─────────────────────────────────────────
   1 │ 2_3                 500.0         500.0
   2 │ 3_2                 500.0         500.0
┌ Info: Interface limits calculated
│   df =
│    6×3 DataFrame
│     Row │ interface  transfer_limit  sum_capacity
│         │ String     Float64         Float64?
│    ─────┼─────────────────────────────────────────
│       1 │ 1_3                 500.0         500.0
│      ⋮  │     ⋮            ⋮              ⋮
└                                     5 rows omitted
6×3 DataFrame
 Row │ interface  transfer_limit  sum_capacity
     │ String     Float64         Float64?
─────┼─────────────────────────────────────────
   1 │ 1_3                 500.0         500.0
   2 │ 3_1                 500.0         500.0
   3 │ 1_2                1175.0        1175.0
   4 │ 2_1                1175.0        1175.0
   5 │ 2_3                 500.0         500.0
   6 │ 3_2                 500.0         500.0
```

## Examples

See the example scripts in the `examples` folder.
