"""
    distance(x1,x2,gridspacing :: Number)

Compute the cartesian distance between two points given their indices and the gridspacing. It asummes uniform grid.


"""
function distance(x1,x2,gridspacing :: Number,weight=1)
    return gridspacing*hypot( x2[1]-x1[1], x2[2]-x1[2] )
end


"""
    compute_N2(xBar_Tv,z)
Take a (1,1,size(z),size(t)) profile of temperature or virtual temperature and return the Brunt - Väisälä frequency at each z level and at each t.
"""
function compute_N2(xBar_Tv,z)
    T = typeof(ustrip(g)/z[1])
    N2 = zeros(T,length(z),size(xBar_Tv,4))
    factor = ustrip(g/Dryair.cp)
    @views N2[1:end-1,:]  .= ustrip(g) .* ( (xBar_Tv[1,1,2:end,:] .- xBar_Tv[1,1,1:end-1,:]) ./(z[2:end] .- z[1:end-1]) .+ factor) ./ xBar_Tv[1,1,1:end-1,:]
    @views N2[end,:]      .= N2[end-1,:]
    ind_smallN2 = findall(abs.(N2) .< 1e-6)
    one_z = CartesianIndex(1,0)
    @inbounds for ind in ind_smallN2
        if 1 < ind[1] < length(z)
           N2[ind] =  (N2[ind + one_z] + N2[ind - one_z])/2 # If N2 is small, substite by mean of neighbours
        elseif ind[1] == 1
           N2[ind] = N2[ind + one_z]
        elseif ind[1] == length(z)
            N2[ind] = N2[ind - one_z]
        end
    end
    return N2
end

function compute_N2(xBar_Tv :: Array{ <:Quantity }, z :: Array{ <:Quantity })
    T = typeof(g/z[1])
    N2 = zeros(T,length(z),size(xBar_Tv,4))
    factor = g/Dryair.cp
    @views N2[1:end-1,:]  .= g .* ( (xBar_Tv[1,1,2:end,:] .- xBar_Tv[1,1,1:end-1,:]) ./ (z[2:end] .- z[1:end-1]) .+ factor) ./ xBar_Tv[1,1,1:end-1,:]
    @views N2[end,:]      .= N2[end-1,:]
    ind_smallN2 = findall(abs.(N2) .< 1e-6u"s^-2")
    one_z = CartesianIndex(1,0)
    @inbounds for ind in ind_smallN2
        if 1 < ind[1] < length(z)
           N2[ind] = (N2[ind + one_z] + N2[ind - one_z])/2 # If N2 is small, substite by mean of neighbours
        elseif ind[1] == 1
           N2[ind] = N2[ind + one_z]
        elseif ind[1] == length(z)
            N2[ind] = N2[ind - one_z]
        end
    end
    return N2
end

#WIP
function compute_mse(T,z,qv)
    sz = size(T)
    return  Dryair.cp*T .+ g*reshape(z,(1,1,sz[3],1)) .+ liquidwater.Lv*qv

end

function get_tendency(field :: AbstractArray{T,4}; dt = error("dt is required for the budget computation")) where T
    dfield_dt = similar(field)
    @. dfield_dt[:,:,:,1:end-1] = @views (field[:,:,:,2:end] - field[:,:,:,1:end-1]) / dt
    @. dfield_dt[:,:,:,end] = dfield_dt[:,:,:,end-1]
 return dfield_dt
end


function get_tendency(field :: AbstractArray{T,1}; dt = error("dt is required for the budget computation")) where T
    dfield_dt = similar(field)
    @. dfield_dt[1:end-1] = (field[2:end] - field[1:end-1]) / dt
    dfield_dt[end] = dfield_dt[end-1]
 return dfield_dt
end

function get_advection_asresidual(tendency,sources...)
    return reduce(+,sources) .- tendency
end

# integrate_vertically(field :: AbstractArray{T,4};dz = 1, weight = 1) where {T} = reduce(+,dz*weight.*field,dims=3)

# integrate_vertically(field :: AbstractArray{T,3};dz = 1, weight = 1) where {T} = reduce(+,dz*weight.*field,dims=2)

# integrate_vertically(field :: AbstractArray{T,2};dz = 1, weight = 1) where {T} = reduce(+,dz*weight.*field,dims=2)

# integrate_vertically(field :: AbstractArray{T,1};dz = 1 , weight = 1) where {T} = reduce(+,dz*weight*field)


# function integrate_vertically(field :: AbstractArray{T,4}; coord :: AbstractArray{T,1},weight = 1) where T
#     sz = size(field)
#     sc = size(coord)
#     integral = zeros(eltype(field),(sz[1],sz[2],1,sz[4]))
#     @inbounds for ind in CartesianIndices(field)
#         if ind[3] != sz[3]
#             integral[ind[1],ind[2],1,ind[4]] = weight*(coord[ind[3]+1] - coord[ind[3]]) * field[ind]
#         end
#     end
#     return integral
# end

# function integrate_vertically(field :: AbstractArray{T,4}; coord :: AbstractArray{T,1},weight :: AbstractArray{T,1}) where T
#     sz = size(field)
#     sc = size(coord)
#     integral = zeros(eltype(field),(sz[1],sz[2],1,sz[4]))
#     @inbounds for ind in CartesianIndices(field)
#         if ind[3] != sz[3]
#             integral[ind[1],ind[2],1,ind[4]] += weight[ind[3]]*(coord[ind[3]+1] - coord[ind[3]]) * field[ind]
#         end
#     end
#     return integral
# end


function integrate_horizontally(field :: AbstractArray{T,4}; darea) where T
    return darea*reduce(+,field)

end


function integrate_vertically(field :: AbstractArray{T,4}; coord :: AbstractArray{T,1},weight :: AbstractArray{T,4}) where T
    sz = size(field)
    sc = size(coord)
    integral = zeros(eltype(field),(sz[1],sz[2],1,sz[4]))
    for ind in CartesianIndices(field)
        if ind[3] != sz[3]
            integral[ind[1],ind[2],1,ind[4]] += weight[ind]*(coord[ind[3]+1] - coord[ind[3]]) * field[ind]
        end
    end
    return integral
end

"""
    compute_virtual_temp(temperature, specific_humidity)
Compute the virtual temperature considering only water vapor mixing ratio.   
"""

function compute_virtual_temp(temperature,specific_humidity)
    return temperature*(1 + epsilon*specific_humidity)
end


function spatial_derivative!(output,field,dx,dim)
    if dim == 1
        onex = CartesianIndex((1, ntuple(_->0, ndims(u) - 1)...))

        @inbounds for ind in CartesianIndices(field)
            if ind[dim] == 1
                output[ind] = (field[ind+onex] - field[ind])/dx
            elseif ind[dim] == size(field,dim)
                output[ind] = (field[ind] - field[ind-onex])/dx
            else
                output[ind] = (field[ind+onex] - field[ind-onex])/2dx
            end
        end
    elseif dim == 2
        oney = CartesianIndex((0,1, ntuple(_->0, ndims(u) - 2)...))
        @inbounds for ind in CartesianIndices(field)
            if ind[dim] == 1
                output[ind] = (field[ind+oney] - field[ind])/dx
            elseif ind[dim] == size(field,dim)
                output[ind] = (field[ind] - field[ind-oney])/dx
            else
                output[ind] = (field[ind+oney] - field[ind-oney])/2dx
            end
        end
        
    end
    return output
end

function spatial_derivative(field, dx, dim)
    output = similar(field)
    return spatial_derivative!(output,field,dx,dim)
end

function get_vorticity!(output,u,v,dx,dy)
    ## This implementation assumes 4d input
    onex = CartesianIndex((1, ntuple(_->0, ndims(u) - 1)...))
    oney = CartesianIndex((0,1, ntuple(_->0, ndims(u) - 2)...))
    for ind in CartesianIndices(output)
        if ind[1] == 1
            output[ind] = (v[ind+onex] - v[ind])/dx
            elseif ind[1] == size(v,1)
            output[ind] = (v[ind] - v[ind-onex])/dx
        else
            output[ind] = (v[ind+onex] - v[ind-onex])/2dx
            end
    end
    for ind in CartesianIndices(output)
        if ind[2] == 1
            output[ind] -= (u[ind+oney] - u[ind])/dy
        elseif ind[2] == size(u,2)
            output[ind] -= (u[ind] - u[ind-oney])/dy
        else
            output[ind] -= (u[ind+oney] - u[ind-oney])/2dy
        end
    end
    return output
end


function get_divergence!(output,u,v,dx,dy)
    ## This implementation assumes 4d input
    onex = CartesianIndex((1, ntuple(_->0, ndims(u) - 1)...))
    oney = CartesianIndex((0,1, ntuple(_->0, ndims(u) - 2)...))
    for ind in CartesianIndices(output)
        if ind[1] == 1
            output[ind] = (v[ind+onex] - v[ind])/dy
            elseif ind[1] == size(v,1)
            output[ind] = (v[ind] - v[ind-onex])/dy
        else
            output[ind] = (v[ind+onex] - v[ind-onex])/2dy
            end
    end
    for ind in CartesianIndices(output)
        if ind[2] == 1
            output[ind] += (u[ind+oney] - u[ind])/dx
        elseif ind[2] == size(u,2)
            output[ind] += (u[ind] - u[ind-oney])/dx
        else
            output[ind] += (u[ind+oney] - u[ind-oney])/2dx
        end
    end
    return output
end

function get_divergence(u,v,dx,dy)
    output = similar(u)
    return get_divergence!(output,u,v,dx,dy)
end


function get_vorticity(u,v,dx,dy)
    output = similar(u)
    return get_vorticity!(output,u,v,dx,dy)
end

function get_okubo_weiss!(output,u,v,dx,dy)
    onex = CartesianIndex((1, ntuple(_->0, ndims(u) - 1)...))
    oney = CartesianIndex((0,1, ntuple(_->0, ndims(u) - 2)...))
    sz = size(u)
#    dvdx
    for ind in CartesianIndices(u)
        ##dudx and dvdx
        if ind[1] == 1
            dudx = (u[ind+onex] - u[ind])/dx
            dvdx = (v[ind+onex] - v[ind])/dx
        elseif ind[1] == sz[1]
            dudx = (u[ind] - u[ind-onex])/dx
            dvdx = (v[ind] - v[ind-onex])/dx
        else
            dudx = (u[ind+onex] - u[ind-onex])/2dx
            dvdx = (v[ind+onex] - v[ind-onex])/2dx
        end

        #dudy and dvdy
        if ind[2] == 1
            dvdy = (v[ind+oney] - v[ind])/dy
            dudy = (u[ind+oney] - u[ind])/dy
        elseif ind[2] == sz[2]
            dvdy = (v[ind] - v[ind-oney])/dy
            dudy = (u[ind] - u[ind-oney])/dy
        else
            dvdy = (v[ind+oney] - v[ind-oney])/2dy
            dudy = (u[ind+oney] - u[ind-oney])/2dy
        end
        output[ind] = dudx*dudx + dvdy*dvdy - 2dudx*dvdy + 2dvdx*dudy +
            2dvdx*dudy
    end
    return output
end

function get_okubo_weiss(u,v,dx,dy)
    ow = similar(u)
    get_okubo_weiss!(ow,u,v,dx,dy)
    return ow
end

"""
    get_saturation_vapor_pressure(T)
Receive temperature T in Kelvin and compute the saturation vapor pressure in hPa from the August-Roche-Magnus formula that approximates the solution to the Clausius-Clapeyron relationship (Wikipedia contributors. (2020, December 19). Clausius–Clapeyron relation. In Wikipedia, The Free Encyclopedia. Retrieved 06:57, December 20, 2020, from https://en.wikipedia.org/w/index.php?title=Clausius%E2%80%93Clapeyron_relation&oldid=995159175)
"""
function get_saturation_vapor_pressure(T)
    return 6.112*exp(17.67 * (T-273.15) / (243.5 + (T - 273.15)))
end

function get_saturation_vapor_pressure(T :: Quantity)
    return 6.112u"hPa"*exp(17.67 * (T-273.15u"K") / (243.5u"K" + (T - 273.15u"K")))
end

"""
    get_partial_vapor_pressure(mixing_ratio,pressure)
Receive a water vapor mixing ratio (unitless g/g) and environmental pressure and compute the partial pressure of water vapor in the same units as the input pressure.
"""
function get_partial_vapor_pressure(mixing_ratio,pressure)
    return mixing_ratio*pressure/(epsilon + mixing_ratio)
end

"""
    get_mixing_ratio(water_vapor_partial_pressure,env_pressure)
Receive a water vapor mixing ratio (unitless g/g) and environmental pressure and compute the partial pressure of water vapor in the same units as the incoming pressure.
"""
function get_mixing_ratio(water_vapor_partial_pressure,env_pressure)
    return epsilon*water_vapor_partial_pressure/(env_pressure - water_vapor_partial_pressure)
end

"""
    get_specific_entropy(temperature,mixing_ratio,pressure)
Receive temperature in Kelvin, water vapor mixing ratio (unitless g/g) and pressure (hPa) and compute the specific entropy of a parcel using equation in Emmanuel's (E94, EQN. 4.5.9)
"""
function get_specific_entropy(temperature,mixing_ratio,pressure)
    CPD=1005.7f0u"J/kg/K"       # [J/kg.K] Specific heat of dry air at constant pressure
    CPV=1870.0f0u"J/kg/K"       # [J/kg.K] Specific heat of water vapor at constant pressure
    CL=4190.0f0u"J/kg/K"        # [J/kg.K] Modified specific heat of liquid water
    CPVMCL=CPV-CL
    alv = Liquidwater.Lv + CPVMCL*(temperature - 273.15f0u"K")
    vapor_pressure = get_partial_vapor_pressure(mixing_ratio,pressure)
    saturation_vapor_pressure = get_saturation_vapor_pressure(temperature)
    RH = min(vapor_pressure/saturation_vapor_pressure,1.0)
    specific_entropy =  (Dryair.cp + mixing_ratio * Liquidwater.cp) *
        log(temperature/unit(temperature)) - Dryair.R * log((pressure - vapor_pressure)/unit(pressure)) +
        alv * mixing_ratio / temperature - mixing_ratio * Watervapor.R * log(RH)
end

"""
    get_specific_entropy_emanuel(temperature,mixing_ratio,pressure)
Receive temperature in Kelvin, water vapor mixing ratio (unitless g/g) and pressure (hPa) and compute the specific entropy of a parcel using equation in Emmanuel's (E94, EQN. 4.5.9)
"""
function get_specific_entropy_emanuel(temperature,mixing_ratio,pressure)
    CPD=1005.7u"J/kg/K"       # [J/kg.K] Specific heat of dry air at constant pressure
    CPV=1870.0u"J/kg/K"       # [J/kg.K] Specific heat of water vapor at constant pressure
    CL=4190.0u"J/kg/K"        # [J/kg.K] Modified specific heat of liquid water
    CPVMCL=CPV-CL
    alv = Liquidwater.Lv + CPVMCL*(temperature - 273.15u"K")
    saturation_vapor_pressure = get_saturation_vapor_pressure(temperature)
    #vapor_pressure = get_partial_vapor_pressure(saturation_vapor_pressure,pressure)
    saturation_mixing_ratio = get_mixing_ratio(saturation_vapor_pressure, pressure)
    specific_entropy =  (Dryair.cp + mixing_ratio * Liquidwater.cp) *
        log(temperature/unit(temperature)) - Dryair.R * log((pressure - saturation_vapor_pressure)/unit(pressure)) +
        alv * saturation_mixing_ratio / temperature 
end

"""

"""
function ∂specific_entropy_∂temp(temperature, mixing_ratio)
    ∂specific_entropy_∂temp = (Dryair.cp + mixing_ratio * Liquidwater.cp)/temperature - Liquidwater.Lv * mixing_ratio / temperature^2 
end

function ∂specific_entropy_∂temp_emanuel(temperature, mixing_ratio, pressure)
    CPD=1005.7u"J/kg/K"       # [J/kg.K] Specific heat of dry air at constant pressure
    CPV=1870.0u"J/kg/K"       # [J/kg.K] Specific heat of water vapor at constant pressure
    CL=4190.0u"J/kg/K"        # [J/kg.K] Modified specific heat of liquid water
    CPVMCL=CPV-CL
    alv = Liquidwater.Lv + CPVMCL*(temperature - 273.15u"K")
    saturation_vapor_pressure = get_saturation_vapor_pressure(temperature)
    saturation_mixing_ratio = get_mixing_ratio(saturation_vapor_pressure, pressure)
    ∂specific_entropy_∂temp = (Dryair.cp + mixing_ratio * Liquidwater.cp + alv^2 * saturation_mixing_ratio /(Watervapor.R*temperature^2))/temperature

end

"""
    get_lifted_condensation_level(temperature,relative_humidity,pressure)   
Receive temperature in Kelvin, relative humidity (unitless) and pressure (hPa) and compute the lifted condensation level based on Emanuel's E94 "calcsound.f" code at http://texmex.mit.edu/pub/emanuel/BOOK/
"""
function get_lifted_condensation_level(temperature,relative_humidity,pressure) 
    return pressure * (relative_humidity^(temperature/(1669.0-122.0*relative_humidity-temperature)))
end

function get_lifted_condensation_level(temperature :: Quantity ,relative_humidity ,pressure :: Quantity) 
    return pressure * (relative_humidity^(temperature/(1669.0u"K"-122.0u"K"*relative_humidity-temperature)))
end
#we need temperature to celsius
#saturation vapor pressure

"""
    specific_humidity_to_mixing_ratio(specific_humidity)
Take a specific humidity (unitless g/g) and return a mixing ratio
"""
function specific_humidity_to_mixing_ratio(specific_humidity)
return mixing_ratio = specific_humidity / (1 - specific_humidity)

end


"""
    mixing_ratio_to_specific_humidity(mixing_ratio)
Take a mixing ratio (unitless g/g) and return a specific humidity
"""
function mixing_ratio_to_specific_humidity(mixing_ratio)
    return q = mixing_ratio / (1 + mixing_ratio)
end

"""
    get_virtual_temperature(temperature,mixing_ratio_total_water,mixing_ratio_water_vapor)
Receive temperature (K) and mixing ratios of total water and water vapor (unitless g/g) and compute the virtual temperature
"""
function get_virtual_temperature(temperature,mixing_ratio_total_water,mixing_ratio_water_vapor)
    return temperature*(1 + mixing_ratio_water_vapor/epsilon)/(1 + mixing_ratio_total_water)
end

"""
    get_minimum_pressure_of_tropical_cyclone(sea_surface_temp,sea_surface_pressure, pressure, temperature, mixing_ratio; ckovercd = 0.9, reversible_ascent=1, dissipative_heating = true)
temperatures in kelvin

"""
function get_minimum_pressure_of_tropical_cyclone(sea_surface_temp,sea_surface_pressure, pressure, temperature, mixing_ratio; ckovercd = 0.9, reversible_ascent=1, dissipative_heating = true)

    initial_level_for_lifting = 1
    exponent_central_pressure = 2.0
    vreduc = 0.8

    t0 = 230.0
    vmax = 0.0
    pmin = 0.0
    ifl = 0 #what is this
    saturation_vapor_pressure0=6.112.*exp(17.67*sea_surface_temperature/(243.5+sea_surface_temperature));
    # initial values (for what?)
    ifl=1;
    NP=0;
    min_pressure=970.0;
    min_pressure_old=min_pressure;
    min_pressure_neew=0.0;
#
#   ***   Find environmental CAPE *** 
#
    tparcel=temperature[initial_level_for_lifting];
    rparcel=mixing_ratio[initial_level_for_lifting]
    pparcel=pressure[initial_level_for_lifting]
    cape_env, outflow_temp_env, index_level_of_neutral_buoyancy = get_cape_and_outflow_temp_from_sounding(tparcel,rparcel,pparcel,temperature,mixing_ratio,pressure)

    while (abs(min_pressure_new-min_pressure_old)) > 0.2
#
#   ***  Find CAPE at radius of maximum winds   ***
#
      tparcel=temperature[initial_level_for_lifting]
      pparcel_approx=min(min_pressure,1000.0) #these two are the ones we are iterating over
      rparcel_approx = epsilon*mixing_ratio[initial_level_for_lifting]*sea_surface_pressure / (pparcel*(epsilon+mixing_ratio[initial_level_for_lifting]) - mixing_ratio[initial_level_for_lifting]*sea_surface_pressure) #what in the name of god is this? it is not documented
      cape_at_rmax, outflow_temp_at_rmax, index_level_of_neutral_buoyancy = get_cape_and_outflow_temp_from_sounding(tparcel,rparcel_approx,pparcel_approx,temperature,mixing_ratio,pressure);
#
#  ***  Find saturation CAPE at radius of maximum winds   ***
#
      tparcel=sea_surface_temperature;
      pparcel_sat=min(min_pressure,1000.0)
      rparcel_sat=get_mixing_ratio(saturation_vapor_presure0,pparcel_approx)
      saturation_cape_at_rmax, temp_outflow, intdex_lnb = get_cape_and_outflow_temp_from_sounding(tparcel,rparcel_sat,pparcel_sat,temperature,mixing_ratio,pressure);
      temp_ratio=sea_surface_temperature/temp_outflow;
       
#
#  ***  Initial estimate of minimum pressure   ***
#
        virtual_temp_parcel_approx=get_virtual_temperature(temperature[initial_level_for_lifting],mixing_ratio_parcel_approx)
        virtual_temp_parcel_sst=get_virtual_temperature(sea_surface_temperature,specific_humidity_parcel_approx)
        average_virtual_temp = 0.5.*(virtual_temperature_parcel_approx + virtual_temperature_parcel_sst);
	CAT=cape_at_rmax-cape_env+0.5.*ck_over_cd.*temp_ratio.*(saturation_cape_at_rmax-cape_at_rmax);
	CAT=max(CAT,0.0);
	min_pressure_new=sea_surface_temperature*exp(-CAT./(287.04.*average_virtual_temp));

	min_pressure_old = min_pressure;
	min_pressure = min_pressure_new;
	niter = niter + 1;
    end
end


function get_cape_and_outflow_temp_from_sounding(tparcel,rparcel,pparcel,t,r,p,ptop=50*unit(pparcel))
    buoyancy_profile = get_buoyancy_of_lifted_parcel(tparcel,rparcel,pparcel,t,r,p,ptop)
#    @show buoyancy_profile
    negative_area=0.0*unit(buoyancy_profile[1]*Dryair.R)
    positive_area=0.0*unit(buoyancy_profile[1]*Dryair.R)
    level_neutral_buoyancy = searchsortedlast(buoyancy_profile, 0.0unit(buoyancy_profile[1]); rev=true)
    @info level_neutral_buoyancy
    for level in 2:(level_neutral_buoyancy + 1)
        area=Dryair.R*(buoyancy_profile[level]+buoyancy_profile[level-1])*(p[level-1]-p[level])/(p[level]+p[level-1])
        positive_area += max(area,0.0*unit(area))
        negative_area -= min(area,0.0*unit(area))
    end
    outflow_temp = t[level_neutral_buoyancy]
    return positive_area - negative_area, outflow_temp, level_neutral_buoyancy
end

"""
    get_buoyancy_of_lifted_parcel(tparcel,rparcel,pparcel,t,r,p,ptop=50)
"""
function get_buoyancy_of_lifted_parcel(tparcel,rparcel,pparcel,t,r,p,ptop=50*unit(pparcel))
    n_valid_levels = findfirst(<(ptop),p)
    p = p[begin:n_valid_levels]
    t = t[begin:n_valid_levels]
    r = r[begin:n_valid_levels]
    tvirtual_diff_parcel_env = similar(t)
    parcel_sat_vapor_pressure = get_saturation_vapor_pressure(tparcel)
    parcel_vapor_pressure = get_partial_vapor_pressure(rparcel,pparcel)
    parcel_rh = min(parcel_vapor_pressure/parcel_sat_vapor_pressure  , 1.0)
    parcel_specific_entropy = get_specific_entropy(tparcel,rparcel,pparcel)
    @info parcel_specific_entropy
    @info parcel_vapor_pressure
    @info parcel_sat_vapor_pressure
    parcel_lcl = get_lifted_condensation_level(tparcel,parcel_rh,pparcel)
    @info parcel_lcl
#    @show parcel_lcl
    levels_below_lcl = findall(>=(parcel_lcl),p)
    levels_above_lcl = findall(<(parcel_lcl),p)
    
    #These two must populate buoyancy of lifted parcel_get_vapor_pressure
    #this would be adiabatic lifting, easy enough
    for level in levels_below_lcl
        tlifted = tparcel*(p[level]/pparcel)^(Dryair.R/Dryair.cp)
        rlifted = rparcel
        tvirtual_lifted = get_virtual_temperature(tlifted,rparcel,rparcel)
        tvirtual_env = get_virtual_temperature(t[level],r[level],r[level])
        tvirtual_diff_parcel_env[level] = tvirtual_lifted - tvirtual_env
        #@info level, tvirtual_diff_parcel_env[level]
    end

    #We start with environmental values of temperature, mixing ratio, entropy etc
    #Our goal: to find the temperature Tx such that s_approx(Tx) ≈ s_parcel s being total entropy
    for level in levels_above_lcl
        initial_guess = t[level]
        target_value = parcel_specific_entropy
        tlifted = find_root_newton_raphson(temp -> get_specific_entropy_emanuel(temp,rparcel,p[level]), temp -> ∂specific_entropy_∂temp_emanuel(temp,rparcel, p[level]); target_value, initial_guess)
        #@info level, tlifted
        saturation_vapor_pressure_lifted = get_saturation_vapor_pressure(tlifted)
        mixing_ratio_lifted = get_mixing_ratio(saturation_vapor_pressure_lifted,p[level])
        tvirtual_lifted = get_virtual_temperature(tlifted, rparcel, mixing_ratio_lifted)
        tvirtual_env = get_virtual_temperature(t[level],r[level], r[level])
        tvirtual_diff_parcel_env[level] = tvirtual_lifted - tvirtual_env
        @info level, tvirtual_diff_parcel_env[level]
    end
    return tvirtual_diff_parcel_env
end

"""
    approximate_t_newton_raphson(tguess,penv,target_entropy)
Given a parcel's total entropy s_parcel and a pressure level P, use the Newton-Raphson method to find T⋆ such that s(T⋆,P) = s_parcel. It assumes reversible moist process.
Total entropy of the parcel before condensation: 
s = (cpd + rt*cl)ln(T) - Rd*ln(pd) + Lv*r/T - r*Rvln(H) : rt is total water (vapor plus liquid) 
"""
function find_root_newton_raphson(func, func_derivative; target_value = 0.0,  initial_guess = 0.0, atol = 0.001*unit(initial_guess), max_iter = 500)
    niter = 0
    initial_err = atol + 5unit(initial_guess)
    approximation = initial_guess
    err = initial_err
    step_size = 1.0
    while (abs(err) > 0.001 * unit(initial_guess))
        niter += 1
        step_size = niter < 3 ? 0.3 : 1.0
        approximation = approximation - step_size*(func(approximation) - target_value)/func_derivative(approximation)
        err = step_size*(func(approximation) - target_value)/func_derivative(approximation)
        if (niter > max_iter ) 
            error("Function didn't converge after $max_iter iterations")
        end
    end
    return approximation
end


"""
    get_potential_temperature(temperature, pressure, reference_pressure)
Compute potential temperature from temperature and pressure.
"""
function get_potential_temperature(temperature, pressure, reference_pressure)
    exponent = ustrip(Dryair.R / Dryair.cp)
    return temperature * (reference_pressure/pressure)^exponent
end

function get_potential_temperature(temperature :: Quantity, pressure :: Quantity, reference_pressure :: Quantity)
    exponent = Dryair.R / Dryair.cp
    return temperature * (reference_pressure/pressure)^exponent
end

"""
    get_virtual_temperature(temperature, specific_humidity)
Compute virtual temperature from temperature and specific humidity.
"""
function get_virtual_temperature(temperature, specific_humidity)
    return (one(temperature) + one(temperature)/1000*epsilon*specific_humidity)*temperature
end

function get_virtual_temperature(temperature :: Quantity, specific_humidity :: Quantity)
    return (one(temperature) + one(temperature)/1000*epsilon*specific_humidity)*temperature
end

"""
   function surface_sensible_heat_flux_to_buoyancy(SST , sensible_heat_flux ; rho = 1)
Convert surface energy fluxes in units of W/m^2 to units of buoyancy m^2/s^3).
"""
function surface_sensible_heat_flux_to_buoyancy(SST, sensible_heat_flux; rho = 1.0)
    return ustrip(g) /(ustrip(Dryair.cp)*SST) * sensible_heat_flux
end

function surface_sensible_heat_flux_to_buoyancy(SST :: Quantity, sensible_heat_flux :: Quantity; rho = 1u"kg/m^3")
    return g /(1u"kg/m^3"*Dryair.cp*SST) * sensible_heat_flux
end

"""
   function surface_latent_heat_flux_to_buoyancy(SST , sensible_heat_flux ; rho = 1)
Convert surface energy fluxes in units of W/m^2 to units of buoyancy m^2/s^3).
"""
function surface_latent_heat_flux_to_buoyancy(SST, latent_heat_flux; rho = 1.0)
    return ustrip(g)/(one(SST)*ustrip(Dryair.cp)*SST)*(epsilon*ustrip(Dryair.cp)*SST/ustrip(Liquidwater.Lv)*latent_heat_flux) 
end

function surface_latent_heat_flux_to_buoyancy(SST :: Quantity, latent_heat_flux :: Quantity; rho = 1u"kg/m^3")
    return g/(1u"kg/m^3"*Dryair.cp*SST)*(epsilon*Dryair.cp*SST/Liquidwater.Lv*latent_heat_flux) 
end

"""
   get_buoyancy(temperature_anomaly,mean_temperature) 
Compute buoyancy from a temperature anomaly and a temperature profile as in RamirezReyes and Yang 2021
"""
function get_buoyancy(temperature_anomaly,mean_temperature) 
    return ustrip(g) * ustrip(temperature_anomaly)/ustrip(mean_temperature)
end

function get_buoyancy(temperature_anomaly :: Quantity ,mean_temperature :: Quantity)
    return g * temperature_anomaly/mean_temperature
end

"""
   radiative_heating_rate_to_buoyancy(mean_temperature,radiative_heating_rate)
Convert radiative heating rate in K/s to buoyancy production rate in m/s^2 using the a mean sounding as reference profile
"""
function radiative_heating_rate_to_buoyancy(mean_temperature,radiative_heating_rate)
    return ustrip(g) * radiative_heating_rate / mean_temperature
end

function radiative_heating_rate_to_buoyancy(mean_temperature :: Quantity,radiative_heating_rate :: Quantity)
    return g * radiative_heating_rate / mean_temperature
end

"""
   get_density_ideal_gas(pressure, temperature)
Compute density of dry air from pressure and temperature using the Ideal gas law.
"""
function get_density_ideal_gas(pressure, temperature)
    return pressure / temperature / eltype(temperature)(ustrip(Dryair.R))
end

function get_density_ideal_gas(pressure :: Quantity, temperature :: Quantity)
    return pressure / temperature / Dryair.R
end
