# Demo for TAMU interface limits

using PowerSystems
using InterfaceLimits

#Download the data from https://electricgrids.engr.tamu.edu/

file_path = joinpath(homedir(),"Downloads/ACTIVSg10k/ACTIVSg10k.RAW")
sys = System(file_path, bus_name_formatter = x->string(x["name"]*"-"*string(x["index"])))
set_units_base_system!(sys, "natural_units") # set units to MW

# calculate interfafce limits
interface_lims = find_interface_limits(sys)