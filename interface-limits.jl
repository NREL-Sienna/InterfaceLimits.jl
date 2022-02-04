using PowerSystems
using JuMP
using Ipopt
using DataFrames


function find_interface_limits(sys)
    # calculate the PTDF
    ptdf = PTDF(sys)

    # Build a JuMP Model
    m = Model(Ipopt.Optimizer)

    # free power injectin variables
    @variable(m, P[get_name.(get_components(Bus, sys))])

    # create variables and constraints for branches and find interfaces
    branches = get_components(Branch, sys)
    @variable(m, F[get_name.(branches)])

    interfaces = Dict{Set,Vector{Branch}}()
    for br in get_components(Branch, sys)
        name = get_name(br)
        @constraint(m, F[name] >= get_rate(br) * -1)
        @constraint(m, F[name] <= get_rate(br))
        @constraint(
            m,
            F[name] == sum([
                ptdf[name, get_number(b)] * P[get_name(b)] for b in get_components(Bus, sys)
            ])
        )
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

    # create variables and constraints for interfaces
    @variable(m, I[join.(keys(interfaces), "_")])

    for (ikey, branches) in interfaces
        @constraint(m, I[join(ikey, "_")] == sum(F[br] for br in get_name.(branches)))
    end

    # make min objective
    @objective(m, Min, sum(I))

    # solve the problem
    optimize!(m)

    # return the interface values
    interface_lims = value.(I)
    min_df = DataFrame(
        :interface => interface_lims.axes[1],
        :min_transfer_limit => interface_lims.data,
    )

    # make max objective
    @objective(m, Max, sum(I))

    # solve the problem
    optimize!(m)

    # return the interface values
    interface_lims = value.(I)
    max_df = DataFrame(
        :interface => interface_lims.axes[1],
        :max_transfer_limit => interface_lims.data,
    )

    return leftjoin(min_df, max_df, on = :interface)
end
