module InterfaceLimits

using PowerSystems
using JuMP
using DataFrames
using PowerNetworkMatrices

export find_interface_limits
export optimizer_with_attributes
export Security
export InjectionLimits

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
    return Security(contingency_branches, VirtualLODF(sys))
end

"""
Function to create `Security` struct with default contingency branches (all branches included)
"""
function Security(sys::System, lodf::VirtualLODF)
    return Security(get_available_components(Branch, sys), lodf)
end

"""
mutable struct InjectionLimits
    genbus_upper_bound::Float64
    loadbus_bounds::Tuple{Float64,Float64}
    enforce_ldfs::Bool
end

Struct to define injection limits at generation and load buses for interface transfer limit calculations.


"""
mutable struct InjectionLimits
    "Constrains max generation as a multiple of existing generation capacity. Must be greater than 0"
    genbus_upper_bound::Float64
    "(Lower, Upper) injection bounds for load buses as a multiple of nominal load at each bus. 
    The upper bound must be positive (maximum power withdrawal) while the lower bound may be 
    positive (minimum power withdrawal) or negative (allows new generation as a multiple of nominal load)."
    loadbus_bounds::Tuple{Float64, Float64}
    "If true, constrains load buses to nominal load distribution factors."
    enforce_ldfs::Bool 

    function InjectionLimits(genbus_upper_bound, loadbus_bounds, enforce_ldfs)
        (genbus_upper_bound <= 0.0) && throw(ArgumentError("genbus_upper_bound must be greater than 0"))
        (loadbus_bounds[1] > loadbus_bounds[2]) && throw(ArgumentError("lower load bound must be less than upper load bound"))
        (loadbus_bounds[2] <= 0.0) && throw(ArgumentError("upper load bound must be greater than 0"))
        if (loadbus_bounds[1] < 0.0) && enforce_ldfs
            @warn("enforce_ldfs and loadbus_bounds[1] < 0 cannot exist together. loadbus_bounds[1] will be set to 0.0.")
            loadbus_bounds = (0.0, loadbus_bounds[2])
        end
        new(genbus_upper_bound, loadbus_bounds, enforce_ldfs)
    end
end

"""
Function to create `InjectionLimits` struct with default values or allow users to specify values as kwargs
"""
# TODO: Keep this distinct function to also allow the creation of InjectionLimits without providing kwargs 
# or combine this into the inner constructor and only enable calling InjectionLimits if kwargs are provided?

function InjectionLimits(;
    genbus_upper_bound::Float64=1.0,
    loadbus_bounds::Tuple{Float64, Float64}=(0.0, 1.0),
    enforce_ldfs::Bool = false
)
    return InjectionLimits(genbus_upper_bound, loadbus_bounds, enforce_ldfs)
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

function add_injector_constraint!( # for buses that have loads and generators and enforce_ldfs==true
    m::Model,
    injector_type::Type{StaticInjection},
    p_var, #injection variable,
    hvdc_inj,
    load_lim::JuMP.AffExpr,
    max_gen,
)
    @constraint(m, p_var - hvdc_inj >= load_lim)
    @constraint(m, p_var - hvdc_inj <= max_gen + load_lim)
end

function add_injector_constraint!( # for buses that have loads and generators and enforce_ldfs==false
    m::Model,
    injector_type::Type{StaticInjection},
    p_var, #injection variable,
    hvdc_inj,
    load_lim::Tuple{Float64,Float64},
    max_gen,
)
    @constraint(m, p_var - hvdc_inj >= load_lim[2])
    @constraint(m, p_var - hvdc_inj <= max_gen + load_lim[1])
end

function add_injector_constraint!( # for buses that have loads only and enforce_ldfs == true
    m::Model,
    injector_type::Type{ElectricLoad},
    p_var, #injection variable,
    hvdc_inj,
    load_lim::JuMP.AffExpr,
    max_gen,
)
    @constraint(m, p_var - hvdc_inj == load_lim)
end

function add_injector_constraint!( # for buses that have loads only and enforce_ldfs == false
    m::Model,
    injector_type::Type{ElectricLoad},
    p_var, #injection variable,
    hvdc_inj,
    load_lim::Tuple{Float64,Float64},
    max_gen,
)
    @constraint(m, p_var - hvdc_inj >= load_lim[2])
    @constraint(m, p_var - hvdc_inj <= load_lim[1])

end

function add_injector_constraint!( # for buses that have generators only
    m::Model,
    injector_type::Type{Generator},
    p_var, #injection variable,
    hvdc_inj,
    load_lim,
    max_gen, 
)
    @constraint(m, p_var - hvdc_inj >= 0.0)
    @constraint(m, p_var - hvdc_inj <= max_gen)
end

function add_injector_constraint!( # for TwoTerminalHVDCLine buses with no gens or loads
    m::Model,
    injector_type::Type{InterconnectingConverter},
    p_var, #injection variable,
    hvdc_inj,
    load_lim,
    max_gen,
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
    injection_limits::InjectionLimits,
)
    # create flow variables for branches
    @variable(m, F[inames, get_name.(in_branches)])
    @variable(m, I[inames])
    @variable(m, P[inames, get_name.(union(gen_buses, load_buses, hvdc_buses))])
    vars = Dict{String,Any}("flow" => F, "interface" => I, "injection" => P)
    if injection_limits.enforce_ldfs
        @variable(m, L[inames]) # constrain in add_constraints!
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

function find_hvdc_buses(sys::System, bus_neighbors)
    dc_br = get_available_components(TwoTerminalHVDCLine, sys)
    from_b = get_from.(get_arc.(dc_br))
    to_b = get_to.(get_arc.(dc_br))
    return Set{ACBus}(intersect(union(from_b, to_b),bus_neighbors))
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
    injection_limits,
)
    F = vars["flow"]
    I = vars["interface"]
    P = vars["injection"]

    total_load, ldf = find_ldfs(sys, load_buses)
    if injection_limits.enforce_ldfs
        L = vars["load"]
        @constraint(m, L .≥ -injection_limits.loadbus_bounds[2] * total_load)
        @constraint(m, L .≤ -injection_limits.loadbus_bounds[1] * total_load)
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
            load_lim = nothing
            if b in load_buses
                bus_ldf = ldf[bus_name]
                if injection_limits.enforce_ldfs
                    load_lim = bus_ldf * L[iname]
                else
                    bus_peak_load = -bus_ldf * total_load
                    load_lim = (injection_limits.loadbus_bounds[1]*bus_peak_load, 
                                injection_limits.loadbus_bounds[2]*bus_peak_load)
                end
            end
            max_gen = nothing
            if b in gen_buses
                max_gen = injection_limits.genbus_upper_bound * sum(
                    get_max_active_power.(
                        get_components(
                            x -> get_available(x) && get_bus(x) == b,
                            Generator,
                            sys,
                        )
                    ),
                )
            end

            hvdc_inj = get_hvdc_inj(b, iname, F, sys)

            add_injector_constraint!(
                m,
                injector_types[b],
                P[iname, bus_name],
                hvdc_inj,
                load_lim,
                max_gen,
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
    
    inj_buses =
        get_components(x -> (x ∈ bus_neighbors) && (get_bustype(x) == bustype), Bus, sys)

    if isempty(inj_buses)
        @warn("No no neighboring $bustype buses")
    end
    return inj_buses
end

function find_gen_buses(sys, bus_neighbors)
    gen_buses = filter(
        x -> x ∈ bus_neighbors,
        Set(get_bus.(get_available_components(Generator, sys))),
    )

    ensure_injector!(gen_buses, bus_neighbors, ACBusTypes.PV, sys)
    return gen_buses
end

function find_load_buses(sys, bus_neighbors)
    load_buses = Set(
        get_bus.(
            get_components(
                x -> get_available(x) && get_bus(x) ∈ bus_neighbors,
                LOAD_TYPES,
                sys,
            )
        ),
    )

    ensure_injector!(load_buses, bus_neighbors, ACBusTypes.PQ, sys)
    return load_buses
end

function get_peak_load(load::ElectricLoad)
    return max(0.0, get_max_active_power(load))
end

function get_peak_load(load::FixedAdmittance)
    return max(0.0, real(get_base_voltage(get_bus(load))^2 * get_Y(load)))
end

function find_ldfs(sys, load_buses)
    bus_loads = Dict(zip(get_name.(load_buses), zeros(length(load_buses))))
    all_loads =
        get_components(x -> get_available(x) && get_bus(x) ∈ load_buses, LOAD_TYPES, sys)
    total_load = sum(get_peak_load.(all_loads))
    for ld in all_loads
        bus_loads[get_name(get_bus(ld))] += get_peak_load(ld) / total_load
    end
    return total_load, bus_loads
end

"""
Calculates the bi-directional interface transfer limits for each interface in a `System`.

# Arguments

- `sys::System` : PowerSystems System data
- `solver::JuMP.MOI.OptimizerWithAttributes` : Solver
- `interface_key::Pair{String, String}` :  key to select single interface (optional arg to run calculation for single interface)
- `interface::Vector{ACBranch}` : vector of branches included in interface (optional arg to run calculation for single interface)

# Keyword Arguments

- `branch_filter::Function = x -> get_available(x)` : generic function to filter branches to be included in transfer limit calculation
- `ptdf::VirtualPTDF = VirtualPTDF(sys),` : power transfer distribution factor matrix
- `security::Union{Bool, Security} = false` : enforce n-1 transmission security constraints
- `injection_limits::InjectionLimits = InjectionLimits()` : constrain generation and load injections
- `hops::Int = 3` : topological distance to include neighboring (outside of interface regions) transmission lines in interface transfer limit calculations

# Examples

```julia
using InterfaceLimits
using PowerSystems
using HiGHS
solver = optimizer_with_attributes(HiGHS.Optimizer)
sys = System("matpower_file_path.m")
interface_lims = find_interface_limits(sys, solver);

# enforce load distribution factors
interface_lims = find_interface_limits(sys, solver, injection_limits=InjectionLimits(enforce_ldfs=true));

# n-1 interface limits
interface_lims = find_interface_limits(sys, solver, security = true);

# omit lines below 230kV from calculation
bf = x -> get_available(x) && all(get_base_voltage.([get_from(get_arc(x)), get_to(get_arc(x))]) .>= 230.0)
interface_lims = find_interface_limits(sys, solver, branch_filter = bf);

# single interface calculation
interfaces = InterfaceLimits.find_interfaces(sys)
interface_key = first(collect(keys(interfaces)))
interface = interfaces[interface_key]
interface_lims = find_interface_limits(sys, solver, interface_key, interface);
```
"""

function find_interface_limits(
    sys::System,
    solver::JuMP.MOI.OptimizerWithAttributes;
    branch_filter::Function = x -> get_available(x),
    ptdf::VirtualPTDF = VirtualPTDF(sys),
    security::Union{Bool,Security} = false,
    injection_limits::InjectionLimits = InjectionLimits(),
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
            injection_limits = injection_limits,
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
    injection_limits::InjectionLimits = InjectionLimits(),
    hops::Int = 3,
    verbose_output::Bool=false,
)
    in_branches, bus_neighbors =
        find_neighbor_lines(sys, interface_key, branch_filter, hops)

    if security == true
        security = Security(sys, get_components(x->(x ∈ in_branches),ACBranch,sys))
    end

    gen_buses = find_gen_buses(sys, bus_neighbors)
    load_buses = find_load_buses(sys, bus_neighbors)
    hvdc_buses = find_hvdc_buses(sys, bus_neighbors)

    roundtrip_ikey = vcat(interface_key, reverse(interface_key))
    inames = join.(roundtrip_ikey, "_")

    # Build a JuMP Model
    m = Model(solver)
    # set_attribute(m, "BarHomogeneous", 1)
    vars = add_variables!(
        m,
        inames,
        in_branches,
        gen_buses,
        load_buses,
        hvdc_buses,
        security,
        injection_limits,
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
        injection_limits,
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

    if verbose_output
        ### CODE FOR DEBUGGING TO OUTPUT MORE VARIABLES ###
        
        inj = vars["injection"]
        flows = vars["flow"]
        total_nom_load, ldfs = find_ldfs(sys, load_buses)

        L_for = 1.0
        L_rev= 1.0
        if injection_limits.enforce_ldfs
            L_for = value.(vars["load"]).data[1]
            L_rev = value.(vars["load"]).data[2]
        end

        inj_df = DataFrame(
            :bus => inj.axes[2],
            :forward_inj => value.(inj[first(inj.axes[1]),:]).data,
            :reverse_inj => value.(inj[last(inj.axes[1]),:]).data,
        )

        transform!(inj_df, :bus =>ByRow(b->get_name(get_area(get_bus(sys,b))))=> :area)

        gen_and_load_buses = intersect(load_buses, gen_buses)

        loads = filter(row -> row.bus ∈ get_name.(setdiff(load_buses, gen_and_load_buses)), inj_df)
        gens = filter(row -> row.bus ∈ get_name.(setdiff(gen_buses, gen_and_load_buses)), inj_df)
        genloads = filter(row -> row.bus ∈ get_name.(gen_and_load_buses), inj_df)

        load_summary = DataFrame(
            :Total_Nom_Load => [-total_nom_load],
            :Total_Fwd_Load => [sum(loads.forward_inj) + sum(genloads[genloads.forward_inj .<0, :forward_inj])], #sum of all buses < 0
            :Total_Rev_Load => [sum(loads.reverse_inj) + sum(genloads[genloads.reverse_inj .<0, :reverse_inj])],
            :Total_Fwd_LoadBuses => [sum(loads.forward_inj) + sum(genloads.forward_inj)], #sum of all buses with loads
            :Total_Rev_LoadBuses => [sum(loads.reverse_inj) + sum(genloads.reverse_inj)],
            :L_fwd => [L_for],
            :L_rev => [L_rev],

        )
        # Add original ldfs to loads
        transform!(loads, :bus =>ByRow(x->ldfs[x]*load_summary.Total_Nom_Load[1])=> :nominal_load)
        transform!(genloads, :bus =>ByRow(x->ldfs[x]*load_summary.Total_Nom_Load[1])=> :nominal_load)
        loads[!, :MaxLoad] = loads.nominal_load .* injection_limits.loadbus_bounds[2]
        loads[!, :MinLoad] = loads.nominal_load .* injection_limits.loadbus_bounds[1]
        genloads[!, :MaxLoad] = genloads.nominal_load .* injection_limits.loadbus_bounds[2]
        genloads[!, :MinLoad] = genloads.nominal_load .* injection_limits.loadbus_bounds[1]

        if injection_limits.enforce_ldfs
            transform!(loads, :bus =>ByRow(x->ldfs[x])=> :original_ldf)
            transform!(genloads, :bus =>ByRow(x->ldfs[x])=> :original_ldf)
            genloads[!, :Ldf_MaxLoad_Fwd] = genloads.original_ldf .* L_for
            genloads[!, :Ldf_MaxLoad_Rev] = genloads.original_ldf .* L_rev

            loads[!,:L_ldf_fwd] = loads.forward_inj ./ L_for #LDF using L variable to calculate
            loads[!,:L_ldf_rev] = loads.reverse_inj ./ L_rev
        end

        # Add gen lims to gens
        transform!(gens, :bus =>ByRow(b->sum(
            get_max_active_power.(
                get_components(
                    x -> get_available(x) && get_name(get_bus(x)) == b,
                    Generator,
                    sys,
                )
            ),
        ) * injection_limits.genbus_upper_bound)=> :gen_lims)

        transform!(genloads, :bus =>ByRow(b->sum(
            get_max_active_power.(
                get_components(
                    x -> get_available(x) && get_name(get_bus(x)) == b,
                    Generator,
                    sys,
                )
            ),
        ) * injection_limits.genbus_upper_bound)=> :gen_lims)
        if injection_limits.enforce_ldfs
            genloads[!, :MaxGen_Fwd] = genloads.Ldf_MaxLoad_Fwd + genloads.gen_lims
            genloads[!, :MaxGen_Rev] = genloads.Ldf_MaxLoad_Rev + genloads.gen_lims
        else
            genloads.MinLoad += genloads.gen_lims
        end

        ## Output Flows
        flow_df = DataFrame(
            :branch => flows.axes[2],
            :forward_flow => value.(flows[first(flows.axes[1]),:]).data,
            :reverse_flow => value.(flows[last(flows.axes[1]),:]).data,
        )
        # add to/from ratings and thermal capacity
        transform!(flow_df, :branch =>ByRow(br->get_name(get_area(get_from(get_arc(get_component(ACBranch, sys, br))))))=> :from_area)
        transform!(flow_df, :branch =>ByRow(br->get_name(get_area(get_to(get_arc(get_component(ACBranch, sys, br))))))=> :to_area)
        transform!(flow_df, :branch =>ByRow(br->get_rate(get_component(ACBranch, sys, br)))=> :flow_limit)

        ### END OF DEBUGGING CODE ###
        return df, loads, gens, genloads, load_summary, flow_df 
    end
    return df
end
end # module
