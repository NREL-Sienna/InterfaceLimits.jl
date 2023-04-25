using InterfaceLimits
using PowerSystems
using PowerNetworkMatrices
using HiGHS
solver = HiGHS.Optimizer
# use a simple rts system
sys = System(
    joinpath(dirname(dirname(pathof(InterfaceLimits))), "examples", "rts", "RTS-GMLC.RAW"),
);

# set units to MW
set_units_base_system!(sys, "natural_units")

# calculate interfafce limits
@info "calculating n-0 interface limits"
@time interface_lims = find_interface_limits(sys, solver);

@info "calculating n-0 interface limits with single problem"
@time interface_lims = find_monolithic_interface_limits(sys, solver);

lodf = LODF(sys)
@info "calculating n-1 interface limits"
interface_lims = find_interface_limits(sys, solver, lodf = lodf, security = true)

@info "calculating n-1 interface limits with single problem"
interface_lims = find_monolithic_interface_limits(sys, solver, lodf = lodf, security = true)


@info "calculating n-0 interface limits for just 1 interface"
interfaces = find_interfaces(sys)
interface_key = first(collect(keys(interfaces)))
interface = interfaces[interface_key]
@time interface_lims = find_interface_limits(sys, solver, interface_key, interface, interfaces);
@time interface_lims = find_interface_limits(sys, solver, interface_key, interface, interfaces, lodf = lodf, security = true);
