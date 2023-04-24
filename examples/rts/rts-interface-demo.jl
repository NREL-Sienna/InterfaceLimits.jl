using InterfaceLimits
using PowerSystems
using HiGHS
# use a simple rts system
sys = System(
    joinpath(dirname(dirname(pathof(InterfaceLimits))), "examples", "rts", "RTS-GMLC.RAW"),
);

# set units to MW
set_units_base_system!(sys, "natural_units")

# calculate interfafce limits
@info "calculating n-0 interface limits"
@time interface_lims = find_interface_limits(sys, HiGHS.Optimizer);

@info "calculating n-0 interface limits with single problem"
@time interface_lims = find_monolithic_interface_limits(sys, HiGHS.Optimizer);

@info "calculating n-1 interface limits"
interface_lims = find_interface_limits(sys, HiGHS.Optimizer, security = true)

@info "calculating n-1 interface limits with single problem"
interface_lims = find_monolithic_interface_limits(sys, HiGHS.Optimizer, security = true)
