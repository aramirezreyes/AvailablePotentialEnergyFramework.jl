cd("/Users/arreyes/Documents/Research/developement/testNetCDF/")
file = "timeSlab_2d.nc"
var = "PW"
ds = NCDatasets.Dataset(file)
units = ds[var].attrib["units"]
x                   = ds3d["x"].var[:]
y                   = ds3d["y"].var[:]
z                   = ds3d["z"].var[:]
t                   = ds3d["time"].var[iterator_time_3d]

#create struct with data::AxisArray
#    and with units


function getvar(ds::Dataset,var::String,portion::tuple=(1))
    
    return ds[var].var[:]    
end


function getunits(ds::Dataset,var::String)
    return ds[var].attrib["units"]    
end

