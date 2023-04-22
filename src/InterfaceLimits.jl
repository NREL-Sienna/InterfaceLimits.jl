module InterfaceLimits

using PowerSystems
using JuMP
using Ipopt
using DataFrames
import PowerModelsInterface
const PMI = PowerModelsInterface

export find_interface_limits
export find_sparse_interface_limits
export find_interfaces

function find_interfaces(sys, branch_filter = x -> get_available(x))
    interfaces = Dict{Set,Vector{ACBranch}}()
    for br in get_components(ACBranch, sys, branch_filter)
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
    return interfaces
end

# need to remove DC lines and inactive components to create a "basic_network" in PowerModels
function clean_sys!(sys)
    remove_components!(!get_available, sys, Device)
    remove_components!(!get_available, sys, HVDCLine)
end

function find_interface_limits(
    sys;
    solver = Ipopt.Optimizer,
    branch_filter = x -> get_available(x),
    ptdf = nothing
)
    if isnothing(ptdf)
        @info "calculating interface limits using PowerModels sparse PTDF technique"
        clean_sys!(sys)
        pm_data = PMI.PM.make_basic_network(PMI.get_pm_data(sys))
        pm_bus_map = Dict(value["name"] => parse(Int64, ix) for (ix, value) in pm_data["bus"])
        pm_branch_map = Dict(get_name(value) => parse(Int64, ix) for (ix, value) in PMI.get_pm_map(sys, ACBranch))
    elseif ptdf isa PowerSystems.PowerNetworkMatrix
        @info "calculating interface limits using PowerSystems PTDF"
        buses = get_components(Bus, sys)
        pm_bus_map = Dict(zip(get_name.(buses), 1:length(buses)))
    end

    # Build a JuMP Model
    @info "Building interface limit optimization model"
    m = Model(solver)

    interfaces = find_interfaces(sys, branch_filter)
    branches = get_components(branch_filter, ACBranch, sys) # could filter for monitored lines here

    inames = join.(keys(interfaces), "_")
    # create flow variables for branches
    @variable(m, F[inames, get_name.(branches)])
    # create variables and interfaces
    @variable(m, I[inames])

    gen_buses = get_bus.(get_components(get_available, Generator, sys))
    load_buses = get_bus.(get_components(get_available, Generator, sys))
    injection_buses = union(gen_buses, load_buses)

    @variable(m, P[inames, get_name.(injection_buses)])

    for (ikey, interface) in interfaces
        iname = join(ikey, "_")
        for b in injection_buses
            if b ∈ setdiff(gen_buses, load_buses) # only gens connected
                @constraint(m, P[iname, get_name(b)] >= 0.0)
            elseif b ∈ setdiff(load_buses, gen_buses) # only loads connected
                @constraint(m, P[iname, get_name(b)] <= 0.0)
            end
        end

        for br in branches
            name = get_name(br)
            ptdf_row = isnothing(ptdf) ? PMI.PM.calc_basic_ptdf_row(pm_data, pm_branch_map[name]) : ptdf[name,:]
            @constraint(m, get_rate(br) >= F[iname, name] >= get_rate(br) * -1)

            @constraint(
                m,
                F[iname, name] == sum([
                    ptdf_row[pm_bus_map[get_name(b)]] * P[iname, get_name(b)] for
                    b in injection_buses
                ])
            )
        end

        @constraint(m, I[iname] == sum(F[iname, get_name(br)] for br in interface))
    end

    # make max objective
    @objective(m, Max, sum(I))

    # solve the problem
    optimize!(m)

    # return the interface values
    interface_lims = value.(I)
    df = DataFrame(
        :interface => interface_lims.axes[1],
        :transfer_limit => interface_lims.data,
    )

    # add capacities to df
    interface_cap = DataFrame(
        :interface => join.(keys(interfaces), "_"),
        :sum_capacity => sum.([get_rate.(br) for br in values(interfaces)]),
    )
    df = leftjoin(df, interface_cap, on = :interface)

    @info "Interface limits calculated" df
    return df
end

end # module
