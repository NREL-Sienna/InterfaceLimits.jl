module InterfaceLimits

using PowerSystems
using JuMP
using Ipopt
using DataFrames

export find_interface_limits
export find_interfaces

function find_interfaces(sys)
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
    return interfaces
end

function find_interface_limits(sys)
    # calculate the PTDF
    ptdf = PTDF(sys)

    # Build a JuMP Model
    m = Model(Ipopt.Optimizer)

    interfaces = find_interfaces(sys)
    branches = get_components(Branch, sys) # could filter for monitored lines here

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
            @constraint(m, F[iname, name] >= get_rate(br) * -1)
            @constraint(m, F[iname, name] <= get_rate(br))
            @constraint(
                m,
                F[iname, name] == sum([
                    ptdf[name, get_number(b)] * P[iname, get_name(b)] for b in get_components(Bus, sys)
                ])
            )
        end

        @constraint(m, I[join(ikey, "_")] == sum(F[iname, get_name(br)] for br in interface))
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

    return df
end

end # module
