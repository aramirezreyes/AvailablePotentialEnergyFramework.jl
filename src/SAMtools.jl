VERSION < v"0.1" && __precompile__()

module SAMtools

using ImageFiltering,
    Statistics,
    NCDatasets,
    Interpolations

include("apehelperfunctions.jl")
include("compositehelperfunctions.jl")
include("apebudgets.jl")

export
# Filters
    filter_array,
    filter_array_2!,
    filter_array_time,
    getsmoothdata,
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
    cyclonecompositer,
    shifter,
    smoothfilter,
    cyclonecompositer_v2,
    cyclonecompositer_v3,
    cyclonecompositer_v3_partial_60d,
    timemean_nofalseframe,
    removefalseframes,
    getapebudget,
    buoyancybudget
    
"""
Data filters:    `filter_array`, `filter_array_2`, `filter_array_time`, `getsmoothdata`, `getsmoothdata_nospace`
Data structures:  `ape_budget`, `cat_ape_budget`, `cutborders!`, `surf_quantities`, `cyclone_comp_timemean`, `Composite_Cyclone`, `Composite_Cyclone_v2`, `Composite_Cyclone_v3`
Methods:    `cyclonecompositer`, `shifter`, `smoothfilter`, `cyclonecompositer_v2`, `cyclonecompositer_v3`, `timemean_nofalseframe`, `removefalseframes`, `getapebudget`, `buoyancybudget`

"""
    
SAMtools

end