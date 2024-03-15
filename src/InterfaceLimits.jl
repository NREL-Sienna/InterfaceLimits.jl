module InterfaceLimits

using PowerSystems
using JuMP
using DataFrames
using PowerNetworkMatrices

export find_interface_limits
export optimizer_with_attributes
export Security

const PSY = PowerSystems
const LOAD_TYPES = ElectricLoad #for NAERM analysis
const HVDC_TYPES = Union{TwoTerminalHVDCLine,TwoTerminalVSCDCLine}
# const LOAD_TYPES = Union{InterruptiblePowerLoad,StaticLoad}

"""
mutable struct Security
    contingency_branches::PowerSystems.InfrastructureSystems.FlattenIteratorWrapper{ACBranch}
    lodf::VirtualLODF
end

Struct to define security modeling data for interface transfer limit calculations.


"""
mutable struct Security
    "Iterator containing branches to be included in security constraints"
    contingency_branches::PSY.IS.FlattenIteratorWrapper{ACBranch}
    "Line outage distribution factor matrix"
    lodf::VirtualLODF
end

"""
Function to create `Security` struct with default values
"""
function Security(sys::System)
    return Security(get_available_components(ACBranch, sys), VirtualLODF(sys))
end

"""
Function to create `Security` struct with default LODF matrix
"""
function Security(
    sys::System,
    contingency_branches::PSY.IS.FlattenIteratorWrapper{ACBranch},
)
    return Security(contingency_branches, VirtualLodf(sys))
end

"""
Function to create `Security` struct with default contingency branches (all branches included)
"""
function Security(sys::System, lodf::VirtualLODF)
    return Security(get_available_components(Branch, sys), lodf)
end

function find_interfaces(sys::System, branch_filter = x -> get_available(x))
    interfaces = Dict{Set,Vector{ACBranch}}()
    for br in get_components(branch_filter, ACBranch, sys)
        from_area = get_area(get_from(get_arc(br)))
        to_area = get_area(get_to(get_arc(br)))
        if from_area != to_area
            key = Set([get_name(from_area), get_name(to_area)])
            if haskey(interfaces, key)
                push!(interfaces[key], br)
            else
                interfaces[key] = [br]
            end
        end
    end
    return Dict(zip([k[1] => k[2] for k in collect.(keys(interfaces))], values(interfaces)))
end

function find_neighbor_lines(
    sys::System,
    interface_key::Pair{String,String},
    branch_filter = x -> get_available(x),
    hops::Int64 = 3, # since only considering lines, update default hops to 3
)
    in_branches = collect(get_components(branch_filter, ACBranch, sys)) #for filtering lines at end
    all_branches = collect(get_components(x -> get_available(x), ACBranch, sys)) #for finding buses
    areas = Set(interface_key)
    # get all the buses in the two areas of the interface
    bus_community = collect(get_components(x -> get_name(get_area(x)) ∈ areas, ACBus, sys))
    line_neighbors = Set{ACBranch}
    for hop = 1:hops
        line_indexes = findall(
            x -> (
                get_from(get_arc(x)) ∈ bus_community || get_to(get_arc(x)) ∈ bus_community
            ),
            all_branches,
        )
        line_neighbors = all_branches[line_indexes]
        # collect all buses regardless of branch filter
        bus_community = Set(
            union(get_to.(get_arc.(line_neighbors)), get_from.(get_arc.(line_neighbors))),
        )
    end
    # filter out line neighbors to only those in branch filter
    line_community = intersect(in_branches, line_neighbors)
    return line_community, bus_community
end

function find_neighbor_lines(
    sys::System,
    interfaces::Dict{Pair{String,String},Vector{ACBranch}},
    branch_filter = x -> get_available(x),
    hops::Int64 = 3,
)
    line_neighbors = Dict{Pair{String,String},Vector{ACBranch}}()
    bus_neighbors = Dict{Pair{String,String},Set{}}()
    for interface in collect(keys(interfaces))
        line_neighbors[interface], bus_neighbors[interface] =
            find_neighbor_lines(sys, interface, branch_filter, hops)
    end
    return line_neighbors, bus_neighbors
end

function find_neighbor_interfaces(
    interfaces::Dict{Pair{String,String},Vector{ACBranch}},
    interface_key::Pair{String,String},
    hops::Int64 = 1,
)
    neighbors = Dict{Pair{String,String},Set{String}}()
    all_interfaces = collect(keys(interfaces))
    community = Set(interface_key)
    for hop = 1:hops
        add_neighbors =
            findall(x -> (first(x) ∈ community || last(x) ∈ community), all_interfaces)
        union!(community, union(Set.(all_interfaces[add_neighbors])...))
    end
    return community
end

function find_neighbor_interfaces(
    interfaces::Dict{Pair{String,String},Vector{ACBranch}},
    hops::Int64 = 1,
)
    neighbors = Dict{Pair{String,String},Set{String}}()
    for interface in collect(keys(interfaces))
        neighbors[interface] = find_neighbor_interfaces(interfaces, interface, hops)
    end
    return neighbors
end

function find_injector_type(
    gen_buses::Set{ACBus},
    load_buses::Set{ACBus},
    hvdc_buses::Set{ACBus},
)
    injector_type = Dict{ACBus,Type{<:StaticInjection}}()
    for bus in union(gen_buses, load_buses, hvdc_buses)
        if bus in intersect(gen_buses, load_buses)
            injector_type[bus] = StaticInjection
        elseif bus in gen_buses
            injector_type[bus] = Generator
        elseif bus in load_buses
            injector_type[bus] = ElectricLoad
        else
            injector_type[bus] = InterconnectingConverter
        end
    end
    return injector_type
end

function add_injector_constraint!( # for buses that have loads and generators
    m::Model,
    injector_type::Type{StaticInjection},
    p_var, #injection variable,
    hvdc_inj;
    ldf_lim = nothing,
    max_gen = nothing, # pass this if enforce_gen_limits = true
)
    if !isnothing(ldf_lim)
        @constraint(m, p_var - hvdc_inj >= ldf_lim)
        !isnothing(max_gen) && (max_gen += ldf_lim)
    end
    !isnothing(max_gen) && @constraint(m, p_var - hvdc_inj <= max_gen)
end

function add_injector_constraint!( # for buses that have loads only
    m::Model,
    injector_type::Type{ElectricLoad},
    p_var, #injection variable,
    hvdc_inj;
    ldf_lim = nothing,
    max_gen = nothing, # pass this if enforce_gen_limits = true
)
    if !isnothing(ldf_lim)
        @constraint(m, p_var - hvdc_inj == ldf_lim)
    else
        @constraint(m, p_var - hvdc_inj <= 0.0)
    end

end

function add_injector_constraint!( # for buses that have generators only
    m::Model,
    injector_type::Type{Generator},
    p_var, #injection variable,
    hvdc_inj;
    ldf_lim = nothing,
    max_gen = nothing, # pass this if enforce_gen_limits = true
)
    @constraint(m, p_var - hvdc_inj >= 0.0)
    !isnothing(max_gen) && @constraint(m, p_var - hvdc_inj <= max_gen)
end

function add_injector_constraint!( # for TwoTerminalHVDCLine buses with no gens or loads
    m::Model,
    injector_type::Type{InterconnectingConverter},
    p_var, #injection variable,
    hvdc_inj;
    ldf_lim = nothing,
    max_gen = nothing, # pass this if enforce_gen_limits = true
)
    @constraint(m, p_var - hvdc_inj == 0.0)
end

function add_variables!(
    m::Model,
    inames::Vector{String},
    in_branches::Vector{ACBranch},
    gen_buses::Set,#PSY.FlattenIteratorWrapper{Bus},
    load_buses::Set,# PSY.FlattenIteratorWrapper{Bus},
    hvdc_buses::Set,# PSY.FlattenIteratorWrapper{Bus},
    security::Union{Bool, Security},
    enforce_load_distribution::Bool,
)
    # create flow variables for branches
    @variable(m, F[inames, get_name.(in_branches)])
    @variable(m, I[inames])
    @variable(m, P[inames, get_name.(union(gen_buses, load_buses, hvdc_buses))])
    vars = Dict{String,Any}("flow" => F, "interface" => I, "injection" => P)
    if enforce_load_distribution
        @variable(m, L, upper_bound = 0.0)
        vars["load"] = L
    end
    if isa(security, Security)
        @variable(
            m,
            CF[
                inames,
                get_name.(security.contingency_branches),
                get_name.(security.contingency_branches),
            ]
        )
        vars["cont_flow"] = CF
    end

    return vars
end

function line_direction(br::ACBranch, ikey::Pair{String,String})
    forward = get_name(get_area(get_from(get_arc(br)))) == first(ikey) ? 1.0 : -1.0
    return forward
end

function get_flow_lims(br::ACBranch)
    return (forward = get_rate(br), reverse = -1 * get_rate(br))
end

function get_flow_lims(br::TwoTerminalHVDCLine)
    return (
        forward = get_active_power_limits_to(br).max,
        reverse = -1 * get_active_power_limits_from(br).max,
    )
end

function get_flow_lims(br::Union{Branch,TwoTerminalVSCDCLine})
    @warn "$(typeof(br)) not considered in interface limit calculations"
    return (forward = Inf, reverse = -Inf)
end

function get_directional_flow_lim(br::ACBranch, ikey::Pair{String,String})
    flow_lim = get_flow_lims(br)
    dir = line_direction(br, ikey)
    lim = dir == 1.0 ? flow_lim.forward : flow_lim.reverse * dir
    return lim
end

function find_hvdc_buses(sys::System)
    dc_br = get_available_components(TwoTerminalHVDCLine, sys)
    from_b = get_from.(get_arc.(dc_br))
    to_b = get_to.(get_arc.(dc_br))
    return Set{ACBus}(union(from_b, to_b))
end

function get_hvdc_inj(b, iname, F, sys)
    dc_brs_from = get_components( 
        x -> (
            get_from(get_arc(x)) == b
        ),
        TwoTerminalHVDCLine, sys)
    dc_brs_to = get_components( 
        x -> (
            get_to(get_arc(x)) == b
        ),
        TwoTerminalHVDCLine, sys)    
    hvdc_inj = 0.0
    length(dc_brs_from) > 0 && (hvdc_inj -= sum(F[iname,get_name.(dc_brs_from)]))
    length(dc_brs_to) > 0 && (hvdc_inj += sum(F[iname,get_name.(dc_brs_to)]))
    return hvdc_inj
end


function add_constraints!(
    m::Model,
    vars,
    interface_key,
    interface,
    gen_buses,
    load_buses,
    hvdc_buses,
    in_branches,
    ptdf,
    security,
    sys,
    enforce_gen_limits,
    enforce_load_distribution,
)
    F = vars["flow"]
    I = vars["interface"]
    P = vars["injection"]
    if enforce_load_distribution
        L = vars["load"]
        ldf = find_ldfs(sys, load_buses)
    end

    if isa(security, Security)
        CF = vars["cont_flow"]
    end

    injector_types = find_injector_type(gen_buses, load_buses, hvdc_buses)

    for ikey in [interface_key, reverse(interface_key)]
        iname = join(ikey, "_")

        for br in in_branches
            name = get_name(br)
            flow_lims = get_flow_lims(br)
            @constraint(m, F[iname, name] >= flow_lims.reverse)
            @constraint(m, F[iname, name] <= flow_lims.forward)

            br isa TwoTerminalHVDCLine && continue

            ptdf_expr = [
                ptdf[name, get_number(b)] * P[iname, get_name(b)] for
                b in union(gen_buses, load_buses, hvdc_buses)
            ]
            push!(ptdf_expr, 0.0)
            @constraint(m, F[iname, name] == sum(ptdf_expr))
            # OutageFlowX = PreOutageFlowX + LODFx,y* PreOutageFlowY
            if isa(security, Security) && br in security.contingency_branches
                c_branches = setdiff(
                    intersect(security.contingency_branches, in_branches),
                    get_components(HVDC_TYPES, sys),
                )
                for cbr in c_branches
                    cname = get_name(cbr)
                    @constraint(m, CF[iname, name, cname] >= flow_lims.reverse)
                    @constraint(m, CF[iname, name, cname] <= flow_lims.forward)
                    @constraint(
                        m,
                        CF[iname, name, cname] ==
                        F[iname, name] + security.lodf[name, cname] * F[iname, cname]
                    )
                end
            end
        end

        for b in union(gen_buses, load_buses, hvdc_buses)
            bus_name = get_name(b)
            ldf_lim =
                (enforce_load_distribution && (b in load_buses)) ? ldf[bus_name] * L :
                nothing
            if enforce_gen_limits && (b in gen_buses)
                max_gen = sum(
                    get_max_active_power.(
                        get_components(
                            x -> get_available(x) && get_bus(x) == b,
                            Generator,
                            sys,
                        )
                    ),
                )
            else
                max_gen = nothing
            end

            hvdc_inj = get_hvdc_inj(b, iname, F, sys)

            add_injector_constraint!(
                m,
                injector_types[b],
                P[iname, bus_name],
                hvdc_inj,
                ldf_lim = ldf_lim,
                max_gen = max_gen,
            )
        end

        @constraint(
            m,
            I[iname] == 
            sum(F[iname, get_name(br)] * line_direction(br, ikey) for br in interface)
        )
    end
end

function ensure_injector!(inj_buses, bus_neighbors, bustype, sys)
    !isempty(inj_buses) && return # only add an injector if set is empty
    #Line hops:
    inj_buses =
        get_components(x -> (x ∈ bus_neighbors) && (get_bustype(x) == bustype), Bus, sys)

    #Region hops:
    #inj_buses = get_components(
    #    x -> (get_name(get_area(x)) ∈ bus_neighbors) && (get_bustype(x) == bustype),
    #    Bus,
    #    sys,
    #)

    if isempty(inj_buses)
        @warn("No no neighboring $bustype buses")
    end
    return inj_buses
end

function find_gen_buses(sys, bus_neighbors)
    # Line hops:
    gen_buses = filter(
        x -> x ∈ bus_neighbors,
        Set(get_bus.(get_available_components(Generator, sys))),
    )

    # Region hops:
    #gen_buses = filter(
    #    x -> get_name(get_area(x)) ∈ bus_neighbors,
    #    Set(get_bus.(get_components(Generator, sys))),
    #)

    ensure_injector!(gen_buses, bus_neighbors, ACBusTypes.PV, sys)
    return gen_buses
end

function find_load_buses(sys, bus_neighbors)
    # Line hops:
    load_buses = Set(
        get_bus.(
            get_components(
                x -> get_available(x) && get_bus(x) ∈ bus_neighbors,
                LOAD_TYPES,
                sys,
            )
        ),
    )

    # Region hops:
    #load_buses = filter(
    #    x -> get_name(get_area(x)) ∈ bus_neighbors,
    #    Set(get_bus.(get_components(LOAD_TYPES, sys))),
    #)

    ensure_injector!(load_buses, bus_neighbors, ACBusTypes.PQ, sys)
    return load_buses
end

function get_peak_load(load::ElectricLoad)
    return get_max_active_power(load)
end

function get_peak_load(load::FixedAdmittance)
    return real(get_base_voltage(get_bus(load))^2 * get_Y(load))
end

function find_ldfs(sys, load_buses)
    bus_loads = Dict(zip(get_name.(load_buses), zeros(length(load_buses))))
    all_loads =
        get_components(x -> get_available(x) && get_bus(x) ∈ load_buses, LOAD_TYPES, sys)
    total_load = sum(get_peak_load.(all_loads))
    for ld in all_loads
        bus_loads[get_name(get_bus(ld))] += get_peak_load(ld) / total_load
    end
    return bus_loads
end

"""
Calculates the bi-directional interface transfer limits for each interface in a `System`.

# Arguments

- `sys::System` : PowerSystems System data
- `solver::JuMP.MOI.OptimizerWithAttributes` : Solver
- `interface_key::Pair{String, String}` :  key to select single interface (optional arg to run calculation for single interface)
- `interface::Vector{ACBranch}` : vector of branches included in interface (optional arg to run calculation for single interface)
- `interfaces::Dict{Pair{String, String}, Vector{ACBranch}}` : dict of all interfaces (optional arg to run calculation for single interface)

# Keyword Arguments

- `branch_filter::Function = x -> get_available(x)` : generic function to filter branches to be included in transfer limit calculation
- `ptdf::VirtualPTDF = VirtualPTDF(sys),` : power transfer distribution factor matrix
- `security::Union{Bool, Security} = false` : enforce n-1 transmission security constraints
- `enforce_gen_limits::Bool = false` : enforce generator capacity limits defined by available generators in System
- `enforce_load_distribution::Bool = false` : enforce constant load distribution factor constraints
- `hops::Int = 3` : topological distance to include neighboring (outside of interface regions) transmission lines in interface transfer limit calculations

# Examples

```julia
using InterfaceLimits
using PowerSystems
using HiGHS
solver = optimizer_with_attributes(HiGHS.Optimizer)
sys = System("matpower_file_path.m")
interface_lims = find_interface_limits(sys, solver);

# enforce generator capacity limits and load distribution factors
interface_lims = find_interface_limits(sys, solver, enforce_gen_limits = true, enforce load_distribution = true);

# n-1 interface limits
interface_lims = find_interface_limits(sys, solver, security = true);

# omit lines below 230kV from calculation
bf = x -> get_available(x) && all(get_base_voltage.([get_from(get_arc(x)), get_to(get_arc(x))]) .>= 230.0)
interface_lims = find_interface_limits(sys, solver, branch_filter = bf);

# single interface calculation
interfaces = InterfaceLimits.find_interfaces(sys)
interface_key = first(collect(keys(interfaces)))
interface = interfaces[interface_key]
interface_lims = find_interface_limits(sys, solver, interface_key, interface, interfaces);
```
"""

function find_interface_limits(
    sys::System,
    solver::JuMP.MOI.OptimizerWithAttributes;
    branch_filter::Function = x -> get_available(x),
    ptdf::VirtualPTDF = VirtualPTDF(sys),
    security::Union{Bool,Security} = false,
    enforce_gen_limits::Bool = false,
    enforce_load_distribution::Bool = false,
    hops::Int = 3,
)

    interfaces = find_interfaces(sys, branch_filter)
    results_dfs = []

    ik = 0
    for (interface_key, interface) in interfaces
        ik += 1
        @info " $ik/$(length(interfaces)) Building interface limit optimization model " interface_key
        df = find_interface_limits(
            sys,
            solver,
            interface_key,
            interface,
            branch_filter = branch_filter,
            ptdf = ptdf,
            security = security,
            enforce_gen_limits = enforce_gen_limits,
            enforce_load_distribution = enforce_load_distribution,
            hops = hops,
        )
        push!(results_dfs, df)
    end

    df = vcat(results_dfs...)

    @info "Interface limits calculated" df
    return df
end

function find_interface_limits(
    sys::System,
    solver::JuMP.MOI.OptimizerWithAttributes,
    interface_key::Pair{String,String},
    interface::Vector{ACBranch};
    branch_filter::Function = x -> get_available(x),
    ptdf::VirtualPTDF = VirtualPTDF(sys),
    security::Union{Bool,Security} = false,
    enforce_gen_limits::Bool = false,
    enforce_load_distribution::Bool = false,
    hops::Int = 3,
)
    if security == true
        security = Security(sys)
    end

    # Line hops:
    in_branches, bus_neighbors =
        find_neighbor_lines(sys, interface_key, branch_filter, hops)

    gen_buses = find_gen_buses(sys, bus_neighbors)
    load_buses = find_load_buses(sys, bus_neighbors)
    hvdc_buses = find_hvdc_buses(sys)

    roundtrip_ikey = vcat(interface_key, reverse(interface_key))
    inames = join.(roundtrip_ikey, "_")

    # Build a JuMP Model
    m = direct_model(solver)
    # set_attribute(m, "BarHomogeneous", 1)
    vars = add_variables!(
        m,
        inames,
        in_branches,
        gen_buses,
        load_buses,
        hvdc_buses,
        security,
        enforce_load_distribution,
    )
    add_constraints!(
        m,
        vars,
        interface_key,
        interface,
        gen_buses,
        load_buses,
        hvdc_buses,
        in_branches,
        ptdf,
        security,
        sys,
        enforce_gen_limits,
        enforce_load_distribution,
    )
    # make max objective
    @objective(m, Max, sum(vars["interface"]))

    # solve the problem
    @info "Solving interface limit problem with" solver
    optimize!(m)

    # return the interface values
    interface_lims = value.(vars["interface"])
    df = DataFrame(
        :interface => interface_lims.axes[1],
        :transfer_limit => interface_lims.data,
    )

    flow_lims = map(x -> sum(get_directional_flow_lim.(interface, x)), roundtrip_ikey)

    interface_cap = DataFrame(:interface => inames, :sum_capacity => flow_lims)
    df = leftjoin(df, interface_cap, on = :interface)
    return df
end
end # module
