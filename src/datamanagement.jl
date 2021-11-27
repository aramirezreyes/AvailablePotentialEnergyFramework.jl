"""
    smooth_vars_and_write_to_netcdf!(output_file,input_file,vars_to_smooth,window_h,window_t)
Take a netcdf file, a list of variables and two integers and performs a moving mean smoothing of these variables. Write the smoothed fields into a netcdf file called output_file
"""
function smooth_vars_and_write_to_netcdf!(output_file,input_file,vars_to_smooth,window_h,window_t)

    ds_orig = Dataset(input_file)
    ds_new = Dataset(output_file,"c", attrib = OrderedDict(
        "history"                   => "Filtered data from file: $input_file, using $window_h points in the space dimensions and $window_t points in the time dimension"
    ))
    # Dimensions
    for (dim,dim_size) in ds_orig.dim
        @info "Writing $dim coordinate in netcdf file"
        flush(stdout)
        ds_new.dim[dim] = dim_size # unlimited dimension
        v = ds_orig[dim]
        create_var = defVar(ds_new,dim, eltype(v.var), dimnames(v),attrib = v.var.attrib)
        create_var[:] = v[:]
    end
    # Declare variables
    for current_var in vars_to_smooth
        @info "Processing $current_var"
        flush(stdout)
        v = ds_orig[current_var]
        create_var = defVar(ds_new,current_var, eltype(v.var), dimnames(v),attrib = v.var.attrib)
        create_var[:] = filter_array(v[:],window_h,window_t)
        @info "Finished processing $current_var"
        flush(stdout)
    end
    close(ds_new)
    close(ds_orig)
end

const dimensions_of_APE_variables = Dict(
    "rho0"    => ("z","t"),
    "N2"      => ("z","t"),
    "B"       => ("x","y","z","t"),
    "QCONVEC" => ("x","y","z","t"),
    "QRAD"    => ("x","y","z","t"),
    "APE"     => ("t",),
    "KE"      => ("t",),
    "APERAD"  => ("t",),
    "APEDIA"  => ("t",),
    "APEWN2"  => ("t",),
    "APEUb2"  => ("t",),
    "APEVb2"  => ("t",),
    "APErate" => ("t",),
    "APESURF" => ("t",)
)

"""
Create NetCDF for storing ape budget variables
"""
function create_APE_netcdf(filename,var_size)

    Dataset(filename,"c") do ds

        # Dimensions
        
        ds.dim["x"] = var_size[1]
        ds.dim["y"] = var_size[2]
        ds.dim["z"] = var_size[3]
        ds.dim["time"] = Inf # unlimited dimension
        
        # Declare variables
        
        ncx = defVar(ds,"x", Float32, ("x",), deflatelevel=1,  attrib = OrderedDict(
            "units"                     => "m",
        ))
        
        ncy = defVar(ds,"y", Float32, ("y",), deflatelevel=1, attrib = OrderedDict(
            "units"                     => "m",
        ))
        
        ncz = defVar(ds,"z", Float32, ("z",), deflatelevel=1,  attrib = OrderedDict(
            "units"                     => "m",
            "long_name"                 => "height",
        ))
        
        nctime = defVar(ds,"time", Float32, ("time",), deflatelevel=1,  attrib = OrderedDict(
            "units"                     => "d",
            "long_name"                 => "time",
        ))
        
        ncAPE = defVar(ds,"APE", Float32, ("time",), deflatelevel=1,  attrib = OrderedDict(
            "long_name"                 => "Available Potential Energy",
            "units"                     => "J/m^2     ",
        ))
        
        ncKE = defVar(ds,"KE", Float32, ("time",), deflatelevel=1,  attrib = OrderedDict(
            "long_name"                 => "Kinetic energy",
            "units"                     => "J/m^2     ",
        ))
        ncAPERAD = defVar(ds,"APERAD", Float32, ("time",), deflatelevel=1,  attrib = OrderedDict(
            "long_name"                 => "Available Potential Energy production by radiative processes",
            "units"                     => "W/m^2     ",
        ))
        ncAPEDIA = defVar(ds,"APEDIA", Float32, ("time",), deflatelevel=1,  attrib = OrderedDict(
            "long_name"                 => "Available Potential Energy production by radiative processes other than radiative and surface fluxes",
            "units"                     => "W/m^2     ",
        ))
        ncAPEWN2 = defVar(ds,"APEWN2", Float32, ("time",), deflatelevel=1,  attrib = OrderedDict(
            "long_name"                 => "Available Potential Energy conversion into kinetic energy",
            "units"                     => "W/m^2     ",
        ))
        ncAPEUB2 = defVar(ds,"APEUb2", Float32, ("time",), deflatelevel=1,  attrib = OrderedDict(
            "long_name"                 => "Available Potential Energy advection by the U component of the wind",
            "units"                     => "W/m^2     ",
        ))
        ncAPEVB2 = defVar(ds,"APEVb2", Float32, ("time",), deflatelevel=1,  attrib = OrderedDict(
            "long_name"                 => "Available Potential Energy advection by the V component of the wind",
            "units"                     => "W/m^2     ",
        ))
        ncAPESURF = defVar(ds,"APESURF", Float32, ("time",), deflatelevel=1,  attrib = OrderedDict(
            "long_name"                 => "Available Potential Energy production by surface fluxes",
            "units"                     => "W/m^2     ",
        ))
        ncAPErate = defVar(ds,"APErate", Float32, ("time",), deflatelevel=1,  attrib = OrderedDict(
            "long_name"                 => "Available Potential Energy rate of change",
            "units"                     => "W/m^2     ",
        ))
        
        ncrho0 = defVar(ds,"rho0", Float32, ("z","time"), deflatelevel=1, shuffle=true, attrib = OrderedDict(
            "units"                     => "kg/m^3",
            "long_name"                 => "Horizontally averaged density",
        ))
        
        ncN2 = defVar(ds,"N2", Float32, ("z","time"), deflatelevel=1, shuffle=true, attrib = OrderedDict(
            "units"                     => "1/s^2",
            "long_name"                 => "Brunt-Vaisala frequency",
        ))
        
        ncQCONVEC = defVar(ds,"QCONVEC", Float32, ("x", "y", "z", "time"), deflatelevel=1, shuffle=true, attrib = OrderedDict(
            "long_name"                 => "Convective heating departure from horizontal mean",
            "units"                     => "m/s^3       ",
        ))
        
        ncQRAD = defVar(ds,"QRAD", Float32, ("x", "y", "z", "time"), deflatelevel=1, shuffle=true, attrib = OrderedDict(
            "long_name"                 => "Radiative heating rate departure from horizontal mean",
            "units"                     => "m/s^3     ",
        ))
        
        ncB = defVar(ds,"B", Float32, ("x", "y", "z", "time"), deflatelevel=1, shuffle=true, attrib = OrderedDict(
            "long_name"                 => "Buoyancy rate departure from horizontal mean",
            "units"                     => "m/s^2     ",
        ))


        for itime = 1:var_size[end]
            nctime[itime] = itime * 3600.0
        end
            
    end

end

"""
Create NetCDF for storing ape budget variables
"""
function set_netcdf_var!(filename,var,data)
    @info "Started writing: $var into file: $filename"
    n_dimensions = ndims(data)
    last_dimension_size = last(size(data))
    Dataset(filename,"a") do ds
        #for ind in 1:last_dimension_size
        #    #@show ndims(data),last_dimension_size, ind,  selectdim(data,n_dimensions,ind)
        #    @info "Writing step: " ind
        #    selectdim(ds[var],n_dimensions,ind)[:] =  selectdim(data,n_dimensions,ind)
        #end
        if n_dimensions == 1
            ds[var][:] = data[:]
        elseif n_dimensions == 2
            ds[var][:,:] = data[:,:]
        elseif n_dimensions == 3
            ds[var][:,:,:] = data[:,:,:]
        elseif n_dimensions == 4
            ds[var][:,:,:,:] = data[:,:,:,:]
        end
        @info "Finished writing: $var into file: $filename"
    end
end

