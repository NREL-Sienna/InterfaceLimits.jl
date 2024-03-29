module InterfaceLimits

using PowerSystems
using JuMP
using DataFrames
using PowerNetworkMatrices

export find_interface_limits
export find_sparse_interface_limits
export find_interfaces
export find_neighbor_interfaces
export find_monolithic_interface_limits
export optimizer_with_attributes

const LOAD_TYPES = Union{InterruptiblePowerLoad,StaticLoad}

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

function add_variables!(
    m,
    inames,
    in_branches,
    gen_buses,
    load_buses,
    security,
    enforce_load_distribution,
)
    # create flow variables for branches
    @variable(m, F[inames, get_name.(in_branches)])
    @variable(m, I[inames])
    @variable(m, P[inames, get_name.(union(gen_buses, load_buses))])
    vars = Dict{String,Any}("flow" => F, "interface" => I, "injection" => P)
    if enforce_load_distribution
        @variable(m, L, upper_bound = 0.0)
        vars["load"] = L
    end
    if security
        @variable(m, CF[inames, get_name.(in_branches), get_name.(in_branches)])
        vars["cont_flow"] = CF
    end

    return vars
end

function add_constraints!(
    m::Model,
    vars,
    interface_key,
    interface,
    gen_buses,
    load_buses,
    in_branches,
    ptdf,
    security,
    lodf,
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

    if security
        CF = vars["cont_flow"]
    end

    for ikey in [interface_key, reverse(interface_key)]
        forward = ikey == interface_key ? 1 : -1
        iname = join(ikey, "_")

        for b in gen_buses # only gens connected
            if enforce_gen_limits
                max_gen = sum(
                    get_max_active_power.(
                        get_components(
                            x -> get_available(x) && get_bus(x) == b,
                            Generator,
                            sys,
                        )
                    ),
                )
                @constraint(m, P[iname, get_name(b)] <= max_gen)
            end
            @constraint(m, P[iname, get_name(b)] >= 0.0)
        end
        for b in load_buses # only loads connected
            if enforce_load_distribution
                @constraint(m, P[iname, get_name(b)] == ldf[get_name(b)] * L)
            else
                @constraint(m, P[iname, get_name(b)] <= 0.0)
            end
        end

        for br in in_branches
            name = get_name(br)
            @constraint(m, get_rate(br) >= F[iname, name] >= get_rate(br) * -1)

            ptdf_expr = [
                ptdf[name, get_number(b)] * P[iname, get_name(b)] for
                b in union(gen_buses, load_buses)
            ]
            push!(ptdf_expr, 0.0)
            @constraint(m, F[iname, name] == sum(ptdf_expr) * forward)
            # OutageFlowX = PreOutageFlowX + LODFx,y* PreOutageFlowY
            if security
                isnothing(lodf) && error("lodf must be defined")
                for cbr in in_branches
                    cname = get_name(cbr)
                    @constraint(
                        m,
                        get_rate(br) >= CF[iname, name, cname] >= get_rate(br) * -1
                    )
                    @constraint(
                        m,
                        CF[iname, name, cname] ==
                        F[iname, name] + lodf[name, cname] * F[iname, cname]
                    )
                end
            end
        end

        @constraint(m, I[iname] == sum(F[iname, get_name(br)] for br in interface))
    end
end

function ensure_injector!(inj_buses, neighbors, bustype, sys)
    !isempty(inj_buses) && return # only add an injector if set is empty
    inj_buses = get_components(
        x -> (get_name(get_area(x)) ∈ neighbors) && (get_bustype(x) == bustype),
        Bus,
        sys,
    )
    if isempty(buses)
        @warn("No no neighboring $bustype buses")
    end
    return inj_buses
end

function find_gen_buses(sys, neighbors)
    gen_buses = filter(
        x -> get_name(get_area(x)) ∈ neighbors,
        Set(get_bus.(get_components(Generator, sys))),
    )
    ensure_injector!(gen_buses, neighbors, BusTypes.PV, sys)
    return gen_buses
end

function find_load_buses(sys, neighbors)
    load_buses = filter(
        x -> get_name(get_area(x)) ∈ neighbors,
        Set(get_bus.(get_components(LOAD_TYPES, sys))),
    )
    ensure_injector!(load_buses, neighbors, BusTypes.PQ, sys)
    return load_buses
end

function find_ldfs(sys, load_buses)
    total_load = sum(get_max_active_power.(get_components(get_available, LOAD_TYPES, sys)))
    ldf = Dict([
        get_name(l) =>
            sum(
                get_max_active_power.(
                    get_components(x -> get_bus(x) == l, LOAD_TYPES, sys)
                ),
            ) / total_load for l in load_buses
    ])
    return ldf
end

"""
find the interface limits (between `Area`s) for a specified `System`

Args:
- `sys::System`
- `solver::OptimizerWithAttributes`

Supported kwargs:
- `branch_filter::function` : function to filter branches considered for interface limit calculation
- `ptdf::VirtualPTDF` : PTDF matrix
- `lodf::VirtualLODF` : LODF matrix
- `security::Bool = False` : include n-1 transmission security constraints in calculation
- `enforce_gen_limits::Bool = False` : enforce generator max power limits in calculation
- `enforce_load_distribution::Bool = False` : enforce constant load distribution factor in calculation
- `hops::Int = 1` : number of hops to travel to include neighboring regions in calculation
"""
function find_interface_limits(
    sys::System,
    solver::JuMP.MOI.OptimizerWithAttributes,
    interface_key,
    interface,
    interfaces;
    branch_filter = x -> get_available(x),
    ptdf = VirtualPTDF(sys),
    lodf = nothing,
    security = false, # n-1 security
    enforce_gen_limits = false,
    enforce_load_distribution = false,
    hops = 1, # neighboring areas to include
)
    interface_neighbors = find_neighbor_interfaces(interfaces, hops)
    neighbors = interface_neighbors[interface_key]
    branches = get_components(branch_filter, ACBranch, sys) # could filter for monitored lines here
    in_branches = filter(
        x -> (
            get_name(get_area(get_from(get_arc(x)))) ∈ neighbors ||
            get_name(get_area(get_to(get_arc(x)))) ∈ neighbors
        ),
        collect(branches),
    )

    inames = join.(vcat(interface_key, reverse(interface_key)), "_")

    gen_buses = find_gen_buses(sys, neighbors)
    load_buses = find_load_buses(sys, neighbors)

    # Build a JuMP Model
    m = direct_model(solver)
    vars = add_variables!(
        m,
        inames,
        in_branches,
        gen_buses,
        load_buses,
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
        in_branches,
        ptdf,
        security,
        lodf,
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

    interface_cap = DataFrame(
        :interface => inames,
        :sum_capacity => ones(2) .* sum(get_rate.(interface)),
    )
    @show df = leftjoin(df, interface_cap, on = :interface)
    return df
end

function find_interface_limits(
    sys::System,
    solver::JuMP.MOI.OptimizerWithAttributes;
    branch_filter = x -> get_available(x),
    ptdf = VirtualPTDF(sys),
    lodf = nothing,
    security = false, # n-1 security
    enforce_gen_limits = false,
    enforce_load_distribution = false,
    hops = 1, # neighboring areas to include
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
            interfaces,
            branch_filter = branch_filter,
            ptdf = ptdf,
            lodf = lodf,
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

function find_monolithic_interface_limits(
    sys::System,
    solver::JuMP.MOI.OptimizerWithAttributes;
    branch_filter = x -> get_available(x),
    ptdf = VirtualPTDF(sys),
    lodf = nothing,
    security = false,
    enforce_load_distribution = false,
    enforce_gen_limits = false,
)
    # Build a JuMP Model
    @info "Building interface limit optimization model"
    m = direct_model(solver)

    interfaces = find_interfaces(sys, branch_filter)
    in_branches = get_components(branch_filter, ACBranch, sys) # could filter for monitored lines here
    gen_buses = Set(get_bus.(get_components(get_available, Generator, sys)))
    load_buses = Set(get_bus.(get_components(get_available, LOAD_TYPES, sys)))

    inames = join.(vcat(collect(keys(interfaces)), reverse.(keys(interfaces))), "_")
    vars = add_variables!(
        m,
        inames,
        in_branches,
        gen_buses,
        load_buses,
        security,
        enforce_load_distribution,
    )
    for (interface_key, interface) in interfaces
        add_constraints!(
            m,
            vars,
            interface_key,
            interface,
            gen_buses,
            load_buses,
            in_branches,
            ptdf,
            security,
            lodf,
            sys,
            enforce_gen_limits,
            enforce_load_distribution,
        )
    end

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

    # add capacities to df
    interface_cap = DataFrame(
        :interface => inames,
        :sum_capacity => repeat(sum.([get_rate.(br) for br in values(interfaces)]), 2),
    )
    df = leftjoin(df, interface_cap, on = :interface)

    @info "Interface limits calculated" df
    return df
end

end # module
