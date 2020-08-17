VERSION < v"0.1" && __precompile__()

"""
Data filters:    `filter_array`, `filter_array_2`, `filter_array_time`, `getsmoothdata`, `getsmoothdata_nospace`
Data structures:  `ape_budget`, `cat_ape_budget`, `cutborders!`, `surf_quantities`, `cyclone_comp_timemean`, `Composite_Cyclone`, `Composite_Cyclone_v2`, `Composite_Cyclone_v3`
Methods:    `cyclonecompositer`, `shifter`, `smoothfilter`, `cyclonecompositer_v2`, `cyclonecompositer_v3`, `timemean_nofalseframe`, `removefalseframes`, `getapebudget`, `buoyancybudget`

"""
module AvailablePotentialEnergyFramework

using ImageFiltering,
    Statistics,
    NCDatasets,
    Interpolations,
    Parameters,
    ImageSegmentation,
    SparseArrays

include("apehelperfunctions.jl")
include("compositehelperfunctions.jl")
include("apebudgets.jl")
include("physicalconstants.jl")
include("physicsfunctions.jl")
include("useful_diagnostics.jl")
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
    add_allcyclones!,
    averageallindistance,
    detect_cyclones,
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
    integrate_vertically,
    average_precipitation_per_pw_bin
    
end
