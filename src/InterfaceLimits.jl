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
export find_neighbor_lines
export optimizer_with_attributes

const LOAD_TYPES = ElectricLoad #for NAERM analysis
# const LOAD_TYPES = Union{InterruptiblePowerLoad,StaticLoad}

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
        line_indexes =
            findall(x -> (get_from(get_arc(x)) ∈ bus_community || get_to(get_arc(x)) ∈ bus_community), all_branches)
        line_neighbors = all_branches[line_indexes]
        # collect all buses regardless of branch filter
        bus_community = Set(union(get_to.(get_arc.(line_neighbors)),get_from.(get_arc.(line_neighbors))))
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
        line_neighbors[interface], bus_neighbors[interface] = find_neighbor_lines(sys, interface, branch_filter, hops)
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

function find_injector_type(gen_buses, load_buses)
    injector_type = Dict{ACBus, Type{<: StaticInjection}}()
    for bus in union(gen_buses, load_buses)
        if bus in intersect(gen_buses, load_buses)
            injector_type[bus] = StaticInjection
        elseif bus in gen_buses
            injector_type[bus] = Generator
        else
            injector_type[bus] = ElectricLoad
        end
    end
    return injector_type
end

function add_injector_constraint!( # for buses that have loads and generators
    m::Model,
    injector_type::Type{StaticInjection},
    p_var; #injection variable,
    ldf_lim = nothing,
    max_gen = nothing, # pass this if enforce_gen_limits = true
)
    if !isnothing(ldf_lim)
        @constraint(m, p_var >= ldf_lim)
        !isnothing(max_gen) && (max_gen += ldf_lim)
    end
    !isnothing(max_gen) && @constraint(m, p_var <= max_gen)
end

function add_injector_constraint!( # for buses that have loads only
    m::Model,
    injector_type::Type{ElectricLoad},
    p_var; #injection variable,
    ldf_lim = nothing,
    max_gen = nothing, # pass this if enforce_gen_limits = true
)
    if !isnothing(ldf_lim)
        @constraint(m, p_var == ldf_lim)
    else
        @constraint(m, p_var <= 0.0)
    end
end

function add_injector_constraint!( # for buses that have generators only
    m::Model,
    injector_type::Type{Generator},
    p_var; #injection variable,
    ldf_lim = nothing,
    max_gen = nothing, # pass this if enforce_gen_limits = true
)
    @constraint(m, p_var >= 0.0)
    !isnothing(max_gen) && @constraint(m, p_var <= max_gen)
end

function add_variables!(
    m,
    inames,
    in_branches,
    gen_buses,
    load_buses,
    security,
    c_branches, #must not be nothing if security = true
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
        isnothing(c_branches) && error("c_branches must be defined")
        @variable(m, CF[inames, get_name.(c_branches), get_name.(c_branches)])
        vars["cont_flow"] = CF
    end

    return vars
end

function line_direction(br, iname)
    forward = get_name(get_area(get_from(get_arc(br)))) == first(iname) ? 1.0 : -1.0
    return forward
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
    c_branches, #must not be nothing if security = true
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

    injector_types = find_injector_type(gen_buses, load_buses)

    for ikey in [interface_key, reverse(interface_key)]
        iname = join(ikey, "_")

        for b in union(gen_buses, load_buses)
            bus_name = get_name(b)
            ldf_lim = (enforce_load_distribution && (b in load_buses)) ? ldf[bus_name] * L : nothing
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

            add_injector_constraint!(
                m,
                injector_types[b],
                P[iname, bus_name],
                ldf_lim = ldf_lim,
                max_gen = max_gen,
            )
        end

        for br in in_branches
            name = get_name(br)
            @constraint(m, F[iname, name] >= get_rate(br) * -1)
            @constraint(m, F[iname, name] <= get_rate(br))

            ptdf_expr = [
                ptdf[name, get_number(b)] * P[iname, get_name(b)] for
                b in union(gen_buses, load_buses)
            ]
            push!(ptdf_expr, 0.0)
            @constraint(m, F[iname, name] == sum(ptdf_expr))
            # OutageFlowX = PreOutageFlowX + LODFx,y* PreOutageFlowY
            if security
                isnothing(lodf) && error("lodf must be defined")
                if br in c_branches
                    for cbr in c_branches
                        cname = get_name(cbr)
                        @constraint(
                            m,
                            CF[iname, name, cname] >= get_rate(br) * -1
                        )
                        @constraint(
                            m,
                            CF[iname, name, cname] <= get_rate(br)
                        )
                        @constraint(
                            m,
                            CF[iname, name, cname] ==
                            F[iname, name] + lodf[name, cname] * F[iname, cname]
                        )
                    end
                end
            end
        end

        @constraint(m, I[iname] == sum(F[iname, get_name(br)] * line_direction(br, iname) for br in interface))
    end
end

function ensure_injector!(inj_buses, bus_neighbors, bustype, sys)
    !isempty(inj_buses) && return # only add an injector if set is empty
    #Line hops:
    inj_buses = get_components(x -> (x ∈ bus_neighbors) && (get_bustype(x) == bustype), Bus, sys)

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
        Set(get_bus.(get_components(Generator, sys))),
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
    load_buses = filter(
        x -> x ∈ bus_neighbors,
        Set(get_bus.(get_components(LOAD_TYPES, sys))),
    )

    # Region hops:
    #load_buses = filter(
    #    x -> get_name(get_area(x)) ∈ bus_neighbors,
    #    Set(get_bus.(get_components(LOAD_TYPES, sys))),
    #)

    ensure_injector!(load_buses, bus_neighbors, ACBusTypes.PQ, sys)
    return load_buses
end

function find_ldfs(sys, load_buses)
    fa_loads = get_components(FixedAdmittance, sys)
    total_load = sum(get_max_active_power.(get_components(get_available, StandardLoad, sys))) +
                    sum(real.(get_base_voltage.(get_bus.(fa_loads)) .^ 2 .* get_Y.(fa_loads)))
    ldf = Dict([
        get_name(l) =>
            (sum(
                get_max_active_power.(
                    get_components(x -> get_bus(x) == l, StandardLoad, sys)
                ), init=0
            ) +
            sum(
                real.(get_base_voltage.(get_bus.(get_components(x -> get_bus(x) == l, FixedAdmittance, sys)))
                    .^ 2 .* get_Y.(get_components(x -> get_bus(x) == l, FixedAdmittance, sys))), init=0
            )) / total_load for l in load_buses
    ])
    return ldf
end

function find_interface_limits(
    sys::System,
    solver,
    interface_key,
    interface,
    interfaces;
    branch_filter = x -> get_available(x),
    ptdf = VirtualPTDF(sys),
    lodf = nothing,
    c_branches = nothing,
    security = false, # n-1 security
    enforce_gen_limits = false,
    enforce_load_distribution = false,
    hops = 3, # neighboring areas to include
)
    # Region hops:
    #interface_neighbors = find_neighbor_interfaces(interfaces, hops)
    #neighbors = interface_neighbors[interface_key]
    #branches = get_components(branch_filter, ACBranch, sys) # could filter for monitored lines here
    #in_branches = filter(
    #    x -> (
    #        get_name(get_area(get_from(get_arc(x)))) ∈ neighbors ||
    #        get_name(get_area(get_to(get_arc(x)))) ∈ neighbors
    #    ),
    #    collect(branches),
    #)
    #gen_buses = find_gen_buses(sys, neighbors)
    #load_buses = find_load_buses(sys, neighbors)

    # Line hops:
    in_branches, bus_neighbors = find_neighbor_lines(sys, interface_key, branch_filter, hops)
    gen_buses = find_gen_buses(sys, bus_neighbors)
    load_buses = find_load_buses(sys, bus_neighbors)

    inames = join.(vcat(interface_key, reverse(interface_key)), "_")

    # Build a JuMP Model
    m = direct_model(solver)
    # set_attribute(m, "BarHomogeneous", 1)
    vars = add_variables!(
        m,
        inames,
        in_branches,
        gen_buses,
        load_buses,
        security,
        c_branches,
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
        c_branches,
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
    solver;
    branch_filter = x -> get_available(x),
    ptdf = VirtualPTDF(sys),
    lodf = nothing,
    c_branches = nothing,
    security = false, # n-1 security
    enforce_gen_limits = false,
    enforce_load_distribution = false,
    hops = 3, # neighboring areas to include
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
            c_branches = c_branches,
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
    solver;
    branch_filter = x -> get_available(x),
    ptdf = VirtualPTDF(sys),
    lodf = nothing,
    c_branches = nothing,
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
        c_branches,
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
            c_branches,
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
