using InterfaceLimits
using PowerSystems
using PowerNetworkMatrices
using HiGHS
solver = optimizer_with_attributes(HiGHS.Optimizer)

# use a simple rts system
sys = System(
    joinpath(dirname(dirname(pathof(InterfaceLimits))), "examples", "rts", "RTS-GMLC.RAW"),
);

# set units to MW
set_units_base_system!(sys, "natural_units")

# calculate interfafce limits
@info "calculating n-0 interface limits"
@time interface_lims = find_interface_limits(sys, solver);

@info "calculating n-0 interface limits with ptdf rounding"
@time interface_lims =
    find_interface_limits(sys, solver, ptdf = VirtualPTDF(sys, tol = 10e-3))

lodf = VirtualLODF(sys, tol = 10e-3)
@info "calculating n-1 interface limits"
interface_lims = find_interface_limits(sys, solver, security = true);


@info "calculating n-0 interface limits for just 1 interface"
interfaces = InterfaceLimits.find_interfaces(sys)
interface_key = first(collect(keys(interfaces)))
interface = interfaces[interface_key]
@time interface_lims =
    find_interface_limits(sys, solver, interface_key, interface);
@time interface_lims = find_interface_limits(
    sys,
    solver,
    interface_key,
    interface,
    security = true,
);

# Example with custom injection limits
@time interface_lims = find_interface_limits(
    sys,
    solver,
    interface_key,
    interface,
    security = true,
    injection_limits = InjectionLimits(genbus_upper_bound=5.0, loadbus_bounds=(0.0,2.0), enforce_ldfs=true),
);
