# Demo for TAMU interface limits

using PowerSystems
using InterfaceLimits
using Xpress #Use a good solver for this

#Download the data from https://electricgrids.engr.tamu.edu/

file_path = joinpath(homedir(), "Downloads/ACTIVSg10k/ACTIVSg10k.RAW")
sys = System(
    file_path,
    bus_name_formatter = x -> string(x["name"] * "-" * string(x["index"])),
)
set_units_base_system!(sys, "natural_units") # set units to MW

# calculate interfafce limits
v = 160.0 #set voltage threshold
# make a branch filter function to select available branches with endpoints above the v theshold
bf = x -> get_available(x) &&
        get_base_voltage(get_from(get_arc(x))) >= v &&
        get_base_voltage(get_to(get_arc(x))) >= v

interface_lims = find_interface_limits(sys, Xpress.Optimizer, bf)
