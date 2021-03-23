VERSION < v"0.1" 

"""
Data filters:    `filter_array`, `filter_array_2`, `filter_array_time`, `getsmoothdata`, `getsmoothdata_nospace`
Data structures:  `ape_budget`, `cat_ape_budget`, `cutborders!`, `surf_quantities`, `cyclone_comp_timemean`, `Composite_Cyclone`, `Composite_Cyclone_v2`, `Composite_Cyclone_v3`
Methods:    `cyclonecompositer`, `shifter`, `smoothfilter`, `cyclonecompositer_v2`, `cyclonecompositer_v3`, `timemean_nofalseframe`, `removefalseframes`, `getapebudget`, `buoyancybudget`

"""
module AvailablePotentialEnergyFramework


using DataStructures: OrderedDict
using ImageFiltering: imfilter, imfilter!, centered, kernelfactors, mapwindow, mapwindow!, Kernel, Inner
using Images: findlocalminima
using ImageSegmentation: SegmentedImage, segment_labels, region_adjacency_graph, seeded_region_growing, labels_map
using Interpolations: LinearInterpolation, interpolate, Gridded, Linear
using NCDatasets: Dataset, variable
using OffsetArrays: OffsetArray
using SparseArrays: SparseMatrixCSC
using Statistics: mean, median!, mean!
using Unitful: @u_str, unit, ustrip, Quantity



include("apehelperfunctions.jl")
include("compositehelperfunctions.jl")
include("apebudgets.jl")
include("physicalconstants.jl")
include("physicsfunctions.jl")
include("useful_diagnostics.jl")
#include("ape_computation_from_julia_output.jl") #testing purposes
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
    smooth_and_mask,
    smooth_and_mask!,
    getapebudget,
    buoyancybudget,
    add_allcyclones!,
    averageallindistance,
    detect_cyclones,
    detect_cyclones!,
    get_diabatic_as_residual_buoyancy,
    run_distributed_test,
    ##Constants
    R,             
    Dryair,
    Watervapor,
    Liquidwater,             
    epsilon,       
    g,
    #Physics functions
    compute_N2,
    compute_mse,
    get_mse_budget,
    get_vorticity,
    get_okubo_weiss,
    get_divergence,
    get_saturation_vapor_pressure,
    get_partial_vapor_pressure,
    get_mixing_ratio,
    get_virtual_temperature,
    get_lifted_condensation_level,
    get_specific_entropy,
    get_potential_temperature,
    get_virtual_temperature,
    mixing_ratio_to_specific_humidity,
    specific_humidity_to_mixing_ratio,
    get_buoyancy_of_lifted_parcel,
    surface_sensible_heat_flux_to_buoyancy,
    surface_latent_heat_flux_to_buoyancy,
    get_buoyancy,
    radiative_heating_rate_to_buoyancy,
    #math functions
    integrate_vertically,
    average_precipitation_per_pw_bin,
    average_precipitation_per_pw_bin_dayang,
    #Files and datamanagement
    smooth_vars_and_write_to_netcdf!
    
end
