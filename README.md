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

```julia
julia> using InterfaceLimits
julia> include(joinpath(dirname(dirname(pathof(InterfaceLimits))), "examples", "5_bus", "5-bus-interface-demo.jl"))

[ Info: Unit System changed to UnitSystem.NATURAL_UNITS = 2
This is Ipopt version 3.14.4, running with linear solver MUMPS 5.4.1.

Number of nonzeros in equality constraint Jacobian...:       97
Number of nonzeros in inequality constraint Jacobian.:       36
Number of nonzeros in Lagrangian Hessian.............:        0

Total number of variables............................:       36
                     variables with only lower bounds:        0
                variables with lower and upper bounds:        0
                     variables with only upper bounds:        0
Total number of equality constraints.................:       21
Total number of inequality constraints...............:       36
        inequality constraints with only lower bounds:       18
   inequality constraints with lower and upper bounds:        0
        inequality constraints with only upper bounds:       18

iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
   0  0.0000000e+00 0.00e+00 2.69e-01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
   1  4.5747846e+02 2.84e-14 4.55e-01  -1.0 2.25e+02  -4.0 6.94e-01 1.00e+00f  1
   2  1.1661066e+03 1.14e-13 1.61e-01  -1.0 5.50e+02  -4.5 6.09e-01 6.64e-01f  1
   3  1.3063865e+03 1.42e-14 7.65e-02  -1.0 2.91e+02  -5.0 7.00e-01 5.69e-01f  1
   4  1.3573258e+03 5.68e-14 2.11e-02  -1.0 3.28e+02  -5.4 8.58e-01 8.28e-01f  1
   5  1.3582843e+03 5.68e-14 3.81e-06  -1.0 3.08e+00  -5.9 1.00e+00 1.00e+00f  1
   6  1.3586702e+03 1.14e-13 1.64e-07  -2.5 3.99e-01  -6.4 1.00e+00 1.00e+00f  1
   7  1.3586835e+03 1.14e-13 1.25e-08  -3.8 9.11e-02  -6.9 1.00e+00 1.00e+00f  1
   8  1.3586842e+03 5.68e-14 1.66e-10  -5.7 3.64e-03  -7.3 1.00e+00 1.00e+00f  1
   9  1.3586842e+03 1.14e-13 2.47e-13  -8.6 1.62e-05  -7.8 1.00e+00 1.00e+00f  1

Number of Iterations....: 9

                                   (scaled)                 (unscaled)
Objective...............:  -1.3586842241006284e+03    1.3586842241006284e+03
Dual infeasibility......:   2.4731437577565306e-13    2.4731437577565306e-13
Constraint violation....:   1.1368683772161603e-13    1.1368683772161603e-13
Variable bound violation:   0.0000000000000000e+00    0.0000000000000000e+00
Complementarity.........:   2.5059870666583991e-09    2.5059870666583991e-09
Overall NLP error.......:   2.5059870666583991e-09    2.5059870666583991e-09


Number of objective function evaluations             = 10
Number of objective gradient evaluations             = 10
Number of equality constraint evaluations            = 10
Number of inequality constraint evaluations          = 10
Number of equality constraint Jacobian evaluations   = 1
Number of inequality constraint Jacobian evaluations = 1
Number of Lagrangian Hessian evaluations             = 1
Total seconds in IPOPT                               = 0.009

EXIT: Optimal Solution Found.
┌ Info: Interface Limits Calculated
│   df =
│    3×3 DataFrame
│     Row │ interface  transfer_limit  sum_capacity
│         │ String     Float64         Float64?
│    ─────┼─────────────────────────────────────────
│       1 │ two_one           718.684         800.0
│       2 │ two_three         240.0           240.0
└       3 │ one_three         400.0           400.0
```
