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

function find_interfaces(sys, branch_filter = x -> get_available(x))
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

# need to remove DC lines and inactive components to create a "basic_network" in PowerModels
function clean_sys!(sys)
    remove_components!(!get_available, sys, Device)
    remove_components!(!get_available, sys, HVDCLine)
end

function add_variables!(m, inames, in_branches, injection_buses, security)
    # create flow variables for branches
    @variable(m, F[inames, get_name.(in_branches)])
    @variable(m, I[inames])
    @variable(m, P[inames, get_name.(injection_buses)])
    vars = Dict("flow" => F, "interface" => I, "injection" => P)
    if security
        @variable(m, CF[inames, get_name.(in_branches), get_name.(in_branches)])
        vars["cont_flow"] = CF
    end

    return vars
end

function add_constraints!(m, vars, interface_key, interface, injection_buses, gen_buses, load_buses, in_branches, ptdf, ptdf_threshold, security, lodf)
    F = vars["flow"]
    I = vars["interface"]
    P = vars["injection"]
    if security
        CF = vars["cont_flow"]
    end

    for ikey in [interface_key, reverse(interface_key)]
        forward = ikey == interface_key ? 1 : -1
        iname = join(ikey, "_")
        for b in injection_buses
            if b ∈ setdiff(gen_buses, load_buses) # only gens connected
                @constraint(m, P[iname, get_name(b)] >= 0.0)
            elseif b ∈ setdiff(load_buses, gen_buses) # only loads connected
                @constraint(m, P[iname, get_name(b)] <= 0.0)
            end
        end

        for br in in_branches
            name = get_name(br)
            @constraint(m, get_rate(br) >= F[iname, name] >= get_rate(br) * -1)

            ptdf_expr = [
                ptdf[name, get_number(b)] * P[iname, get_name(b)] for
                b in injection_buses if
                abs(ptdf[name, get_number(b)]) > ptdf_threshold
            ]
            push!(ptdf_expr, 0.0)
            @constraint(
                m,
                F[iname, name] ==sum(ptdf_expr) * forward
            )
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

function ensure_injector!(inj_buses, neighbors, bustype,  sys)
    !isempty(inj_buses) && return # only add an injector if set is empty
    inj_buses = get_components(x->(get_name(get_area(x)) ∈ neighbors) && (get_bustype(x) == bustype), Bus, sys)
    if isempty(buses)
        @warn("No no neighboring $bustype buses")
    end
    return inj_buses
end

function find_interface_limits(
    sys,
    solver,
    interface_key,
    interface,
    interfaces;
    branch_filter = x -> get_available(x),
    ptdf = VirtualPTDF(sys),
    lodf = nothing,
    security = false, # n-1 security
    hops = 1, # neighboring areas to include
    ptdf_threshold = 1e-6, #rounding threshold for including ptdf
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

    interface_areas = get_components(x->get_name(x) ∈ interface_key, Area, sys)
    #inames = join.(vcat(collect(keys(interfaces)), reverse.(keys(interfaces))), "_")
    inames = join.(vcat(interface_key, reverse(interface_key)), "_")
    gen_buses = filter(
        x -> get_name(get_area(x)) ∈ neighbors,
        Set(get_bus.(get_components(Generator, sys))),
    )
    ensure_injector!(gen_buses, neighbors, BusTypes.PV, sys)
    load_buses = filter(
        x -> get_name(get_area(x)) ∈ neighbors,
        Set(get_bus.(get_components(ElectricLoad, sys))),
    )
    ensure_injector!(load_buses, neighbors, BusTypes.PQ, sys)
    injection_buses = union(gen_buses, load_buses)

    # Build a JuMP Model
    m = direct_model(solver)
    vars = add_variables!(m, inames, in_branches, injection_buses, security)
    add_constraints!(m, vars, interface_key, interface, injection_buses, gen_buses, load_buses, in_branches, ptdf, ptdf_threshold, security, lodf)
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
    sys,
    solver;
    branch_filter = x -> get_available(x),
    ptdf = VirtualPTDF(sys),
    lodf = nothing,
    security = false, # n-1 security
    hops = 1, # neighboring areas to include
    ptdf_threshold = 1e-6, #rounding threshold for including ptdf
)

    interfaces = find_interfaces(sys, branch_filter)
    results_dfs = []

    ik = 0
    for (interface_key, interface) in interfaces
        ik+=1
        @info " $ik/$(length(interfaces)) Building interface limit optimization model " interface_key
        df = find_interface_limits(sys, solver, interface_key, interface, interfaces, branch_filter = branch_filter, ptdf = ptdf, lodf = lodf, security = security, hops = hops, ptdf_threshold = ptdf_threshold)
        push!(results_dfs, df)
    end

    df = vcat(results_dfs...)

    @info "Interface limits calculated" df
    return df
end

function find_monolithic_interface_limits(
    sys,
    solver;
    branch_filter = x -> get_available(x),
    ptdf = VirtualPTDF(sys),
    lodf = nothing,
    security = false,
    ptdf_threshold = 1e-6, #rounding threshold for including ptdf
)
    # Build a JuMP Model
    @info "Building interface limit optimization model"
    m = direct_model(solver)

    interfaces = find_interfaces(sys, branch_filter)
    in_branches = get_components(branch_filter, ACBranch, sys) # could filter for monitored lines here
    gen_buses = Set(get_bus.(get_components(get_available, Generator, sys)))
    load_buses = Set(get_bus.(get_components(get_available, ElectricLoad, sys)))
    injection_buses = union(gen_buses, load_buses)

    inames = join.(vcat(collect(keys(interfaces)), reverse.(keys(interfaces))), "_")
    vars = add_variables!(m, inames, in_branches, injection_buses, security)
    for (interface_key, interface) in interfaces
        add_constraints!(m, vars, interface_key, interface, injection_buses, gen_buses, load_buses, in_branches, ptdf, ptdf_threshold, security, lodf)
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
