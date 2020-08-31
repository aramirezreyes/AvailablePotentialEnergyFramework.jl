VERSION < v"0.1" 

"""
Data filters:    `filter_array`, `filter_array_2`, `filter_array_time`, `getsmoothdata`, `getsmoothdata_nospace`
Data structures:  `ape_budget`, `cat_ape_budget`, `cutborders!`, `surf_quantities`, `cyclone_comp_timemean`, `Composite_Cyclone`, `Composite_Cyclone_v2`, `Composite_Cyclone_v3`
Methods:    `cyclonecompositer`, `shifter`, `smoothfilter`, `cyclonecompositer_v2`, `cyclonecompositer_v3`, `timemean_nofalseframe`, `removefalseframes`, `getapebudget`, `buoyancybudget`

"""
module AvailablePotentialEnergyFramework


using DataStructures: OrderedDict
using ImageFiltering: imfilter, imfilter!, centered, kernelfactors, mapwindow, Kernel, Inner
using Images: findlocalminima
using ImageSegmentation: SegmentedImage, segment_labels, region_adjacency_graph, seeded_region_growing, labels_map
using Interpolations: LinearInterpolation, interpolate, Gridded, Linear
using NCDatasets: Dataset, variable
using OffsetArrays:OffsetArray
using SparseArrays: SparseMatrixCSC
using Statistics: mean
using Parameters: @with_kw



include("apehelperfunctions.jl")
include("compositehelperfunctions.jl")
include("apebudgets.jl")
include("physicalconstants.jl")
include("physicsfunctions.jl")
include("useful_diagnostics.jl")
include("ape_computation_from_julia_output.jl") #testing purposes
include("arrayfiltering.jl")
include("datamanagement.jl")
export
# Filters
    filter_array!,
    filter_array_2!,
    filter_array_time,
    filter_array,
    getsmoothdata!,
    getsmoothdata_nospace,
    #Data structures
    cat_ape_budget,
    cutborders!,
    #Methods
    findlocalmaxima,
    cyclonecompositer,
    shifter,
    shifter!,
    smoothfilter,
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
    average_precipitation_per_pw_bin,
    #Files and datamanagement
    smooth_vars_and_write_to_netcdf!
    
end
