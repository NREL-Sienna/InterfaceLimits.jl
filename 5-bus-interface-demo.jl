using InterfaceLimits
using PowerSystems
using DataFrames
# use a simple 5-bus system
sys = System(joinpath(dirname(dirname(pathof(InterfaceLimits))), "5_bus.json"));

# add areas
remove_component!(sys, first(get_components(Area, sys)))
areas = [Area(name = "one"), Area(name = "two"), Area(name = "three")]
add_components!(sys, areas)

# add bus-area memberships
for b in get_components(Bus, sys)
    if get_name(b) ∈ ["node_a", "node_c", "node_b"]
        set_area!(b, get_component(Area, sys, "one"))
    elseif get_name(b) == "node_d"
        set_area!(b, get_component(Area, sys, "two"))
    else
        set_area!(b, get_component(Area, sys, "three"))
    end
end
# set units to MW
set_units_base_system!(sys, "natural_units")

# solve problem
interface_lims = find_interface_limits(sys);

# compare solution with capacities
interfaces = find_interfaces(sys);

interface_cap = DataFrame(
    :interface => join.(keys(interfaces), "_"),
    :sum_capacity => sum.([get_rate.(br) for br in values(interfaces)]),
);
df = leftjoin(interface_lims, interface_cap, on = :interface);
@info("Interface Limits Calculated", df)
