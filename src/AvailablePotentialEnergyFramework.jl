VERSION < v"0.1" && __precompile__()

module AvailablePotentialEnergyFramework

using ImageFiltering,
    Statistics,
    NCDatasets,
    Interpolations,
    Parameters,
    ImageSegmentation

include("apehelperfunctions.jl")
include("compositehelperfunctions.jl")
include("apebudgets.jl")
include("physicalconstants.jl")
include("physicsfunctions.jl")
#include("arraytools.jl")
export
# Filters
    filter_array!,
    filter_array_2!,
    filter_array_time,
    filter_array,
    getsmoothdata!,
    getsmoothdata_nospace,
    #Data structures
    ape_budget,
    cat_ape_budget,
    cutborders!,
    surf_quantities,
    cyclone_comp_timemean,
    Composite_Cyclone,
    Composite_Cyclone_v2,
    Composite_Cyclone_v3,
    #Methods
    findlocalmaxima,
    cyclonecompositer,
    shifter,
    shifter!,
    smoothfilter,
    cyclonecompositer_v2,
    cyclonecompositer_v3,
    cyclonecompositer_v3_partial_60d,
    cyclonecompositer_v4_partial,
    timemean_nofalseframe,
    removefalseframes,
    getapebudget,
    buoyancybudget,
    ##Constants
    R,             
    heat_capacity,
    L,             
    epsilon,       
    g,
    #Physics functions
    compute_N2,
    compute_mse,
    get_mse_budget,
    get_vorticity,
    get_okubo_weiss,
    get_divergence,
    #math functions
    integrate_vertically
"""
Data filters:    `filter_array`, `filter_array_2`, `filter_array_time`, `getsmoothdata`, `getsmoothdata_nospace`
Data structures:  `ape_budget`, `cat_ape_budget`, `cutborders!`, `surf_quantities`, `cyclone_comp_timemean`, `Composite_Cyclone`, `Composite_Cyclone_v2`, `Composite_Cyclone_v3`
Methods:    `cyclonecompositer`, `shifter`, `smoothfilter`, `cyclonecompositer_v2`, `cyclonecompositer_v3`, `timemean_nofalseframe`, `removefalseframes`, `getapebudget`, `buoyancybudget`

"""
    
SAMtools

end
