using NCDatasets: Dataset
using DataStructures: OrderedDict

"""
    smooth_vars_and_write_to_netcdf!(output_file,input_file,vars_to_smooth,window_h,window_t)
Take a netcdf file, a list of variables and two integers and performs a moving mean smoothing of these variables. Write the smoothed fields into a netcdf file called output_file
"""
function smooth_vars_and_write_to_netcdf!(output_file,input_file,vars_to_smooth,window_h,window_t)

    ds_orig = Dataset(input_file)
    ds_new = Dataset(output_file,"c", attrib = OrderedDict(
        "history"                   => "Filtered data from file: input_file, using $window_h points in the space dimensions and $window_t points in the time dimension"
    ))
    # Dimensions
    for (dim,dim_size) in ds_orig.dim
        ds_new.dim[dim] = dim_size # unlimited dimension
        v = ds_orig[dim]
        create_var = defVar(ds_new,dim, eltype(v.var), dimnames(v),attrib = v.var.attrib)
        create_var[:] = v[:]
    end
    # Declare variables
    for current_var in vars_to_smooth
        v = ds_orig[current_var]
        create_var = defVar(ds_new,current_var, eltype(v.var), dimnames(v),attrib = v.var.attrib)
        create_var[:] = filter_array(v[:],window_h,window_t)
    end
    close(ds_new)
    close(ds_orig)
end





