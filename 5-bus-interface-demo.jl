using InterfaceLimits
using PowerSystems
using DataFrames
# use a simple 5-bus system
sys = System("5_bus.json")

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

# solve problem
interface_lims = find_interface_limits(sys);

# compare solution with capacities
interfaces = Dict{Set,Vector{Branch}}();
for br in get_components(Branch, sys)
    from_area = get_area(get_from(get_arc(br)))
    to_area = get_area(get_to(get_arc(br)))
    if from_area != to_area
        key = Set(get_name.([from_area, to_area]))
        if haskey(interfaces, key)
            push!(interfaces[key], br)
        else
            interfaces[key] = [br]
        end
    end
end

interface_cap = DataFrame(
    :interface => join.(keys(interfaces), "_"),
    :sum_capacity => sum.([get_rate.(br) for br in values(interfaces)]),
);
df = leftjoin(interface_lims, interface_cap, on = :interface);
@info("Interface Limits Calculated", df)
