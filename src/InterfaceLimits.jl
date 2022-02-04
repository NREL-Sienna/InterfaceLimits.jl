module InterfaceLimits

using PowerSystems
using JuMP
using Ipopt
using DataFrames

export find_interface_limits
export find_interfaces

function find_interfaces(sys)
    interfaces = Dict{Set,Vector{ACBranch}}();
    for br in get_components(ACBranch, sys)
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

function find_interface_limits(sys)
    # calculate the PTDF
    @info "Building PTDF"
    ptdf = PTDF(sys)

    # Build a JuMP Model
    @info "Building interface limit optimization model"
    m = Model(Ipopt.Optimizer)

    interfaces = find_interfaces(sys)
    branches = get_components(ACBranch, sys) # could filter for monitored lines here

    inames = join.(keys(interfaces), "_")
    # free power injectin variables
    @variable(m, P[inames, get_name.(get_components(Bus, sys))])
    # create flow variables for branches
    @variable(m, F[inames, get_name.(branches)])
    # create variables and interfaces
    @variable(m, I[inames])

    for (ikey, interface) in interfaces
        iname = join(ikey, "_")
        for br in branches
            name = get_name(br)
            @constraint(m,  get_rate(br) >= F[iname, name] >= get_rate(br) * -1)

            @constraint(
                m,
                F[iname, name] == sum([
                    ptdf[name, get_number(b)] * P[iname, get_name(b)] for b in get_components(Bus, sys)
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
    );
    df = leftjoin(df, interface_cap, on = :interface);

    @info "Interface limits calculated" df
    return df
end

end # module
