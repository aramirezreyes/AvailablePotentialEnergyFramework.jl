"""
    distance(x1,x2,gridspacing :: Number)

Compute the cartesian distance between two points given their indices and the gridspacing. It asummes uniform grid.


"""
function distance(x1,x2,gridspacing :: Number,weight=1)
    return gridspacing*hypot( x2[1]-x1[1], x2[2]-x1[2] )
end

"""
    as_ints(a::AbstractArray{CartesianIndex{L}}) where L
Take an array of cartesian indices and transforms it to an array of integers
"""
as_ints(a::AbstractArray{CartesianIndex{L}}) where L = reshape(reinterpret(Int, a), (L, size(a)...))


"""
    compute_N2(xBar_Tv,z)
Take a (1,1,size(z),size(t)) profile of temperature or virtual temperature and return the Brunt - Väisälä frequency at each z level and at each t.
"""
function compute_N2(xBar_Tv,z)
    T = typeof(ustrip(g)/z[1])
    N2 = zeros(T,length(z),size(xBar_Tv,4))
    factor = ustrip(g/Dryair.cp)
    @views  N2[1:end-1,:]  .= ustrip(g) * ( (xBar_Tv[1,1,2:end,:]-xBar_Tv[1,1,1:end-1,:])./(z[2:end].-z[1:end-1]) .+ factor)./xBar_Tv[1,1,1:end-1,:]
    #@views  @. N2  = g * (N2 + factor)/xBar_Tv[1,1,:,:]
    @views  N2[end,:]      .= N2[end-1,:]
        bb1 = as_ints(findall(abs.(N2) .< 1e-6 ))
        @views bb = bb1[1,:]
        @views cc = bb1[2,:]   
    @inbounds for i in 1:length(bb)
        if 1 < bb[i] < size(z)[1]
           N2[bb[i],cc[i]] = 0.5 * (N2[bb[i]-1,cc[i]] + N2[bb[i]+1,cc[i]]) # If N2 is small, substite by mean of neighbours
        elseif bb[i] == 1
           N2[bb[i],cc[i]] = N2[bb[i]+1,cc[i]]
        elseif bb[i] == size(z)[1]
            N2[bb[i],cc[i]] = N2[bb[i] - 1,cc[i]]
        end
    end
    return N2
end

function compute_N2(xBar_Tv :: Array{ <:Quantity }, z :: Array{ <:Quantity })
    T = typeof(g/z[1])
    N2 = zeros(T,length(z),size(xBar_Tv,4))
    factor = g/Dryair.cp
    @views  N2[1:end-1,:]  .= g * ( (xBar_Tv[1,1,2:end,:]-xBar_Tv[1,1,1:end-1,:])./(z[2:end].-z[1:end-1]) .+ factor)./xBar_Tv[1,1,1:end-1,:]
    #@views  @. N2  = g * (N2 + factor)/xBar_Tv[1,1,:,:]
    @views  N2[end,:]      .= N2[end-1,:]
        bb1 = as_ints(findall(abs.(N2) .< 1e-6*unit(N2[1]) ))
        @views bb = bb1[1,:]
        @views cc = bb1[2,:]   
    @inbounds for i in 1:length(bb)
        if 1 < bb[i] < size(z)[1]
           N2[bb[i],cc[i]] = 0.5 * (N2[bb[i]-1,cc[i]] + N2[bb[i]+1,cc[i]]) # If N2 is small, substite by mean of neighbours
        elseif bb[i] == 1
           N2[bb[i],cc[i]] = N2[bb[i]+1,cc[i]]
        elseif bb[i] == size(z)[1]
            N2[bb[i],cc[i]] = N2[bb[i] - 1,cc[i]]
        end
    end
    return N2
end

"""
    compute_N2_attempt(xBar_Tv,z)
Take a (1,1,size(z),size(t)) profile of temperature or virtual temperature and return the Brunt - Väisälä frequency at each z level and at each t. Tried doing it faster that the other function but have not been succesful.
"""
function compute_N2_attempt(xBar_Tv,z)
    T = eltype(xBar_Tv)
    N2 = zeros(T,length(z),size(xBar_Tv,4))
    sz,st = size(N2)
    factor = g/Dryair.cp
    @inbounds @simd for indt in 1:st
        for indz in 1:(sz - 1)
            Tv = xBar_Tv[1,1,indz,indt]
        dtvdz = (xBar_Tv[1,1,indz+1,indt] - Tv)/(z[indz+1] - z[indz])
        N2[indz,indt] = g*(dtvdz + factor)/Tv
        end
    end
    @inbounds @simd for indt in 1:st
        N2[sz,indt] = N2[sz-1,indt]
    end
    bb1 = as_ints(findall(abs.(N2) .< 1e-6))
   
    bb = bb1[1,:]
    cc = bb1[2,:]
    @inbounds @simd for i in 1:length(bb)
        if 1 < bb[i] < sz
            N2[bb[i],cc[i]] = 0.5 * (N2[bb[i]-1,cc[i]] + N2[bb[i]+1,cc[i]]) # If N2 is small, substite by mean of neighbours
        elseif bb[i] == 1
            N2[bb[i],cc[i]] = N2[bb[i]+1,cc[i]]
        elseif bb[i] == sz      
            N2[bb[i],cc[i]] = N2[bb[i]-1,cc[i]]
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


function compute_virtual_temp(T,QV)
    return T.*(1 .+ 0.61QV)
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
    vapor_pressure = get_partial_vapor_pressure(mixing_ratio,pressure)
    saturation_vapor_pressure = get_saturation_vapor_pressure(temperature)
    RH = min(vapor_pressure/saturation_vapor_pressure,1.0)
    specific_entropy =  (Dryair.cp + mixing_ratio * Liquidwater.cp) *
        log(temperature/unit(temperature)) - Dryair.R * log((pressure - vapor_pressure)/unit(pressure)) +
        Liquidwater.Lv * mixing_ratio / temperature - mixing_ratio * Watervapor.R * log(RH)
end 

"""
    get_lifted_condensation_level(temperature,relative_humidity,pressure)   
Receive temperature in Kelvin, relative humidity (unitless) and pressure (hPa) and compute the lifted condensation level based on Emanuel's E94 "calcsound.f" code at http://texmex.mit.edu/pub/emanuel/BOOK/
"""
function get_lifted_condensation_level(temperature,relative_humidity,pressure) 
    return pressure * (relative_humidity^(temperature/(1669.0-122.0*relative_humidity-temperature)))
end

function get_lifted_condensation_level(temperature :: Quantity ,relative_humidity :: Quantity ,pressure :: Quantity) 
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
    get_buoyancy_of_lifted_parcel(tparcel,rparcel,pparcel,t,r,p,ptop=50)
"""
function get_buoyancy_of_lifted_parcel(tparcel,rparcel,pparcel,t,r,p,ptop=50u"hPa")
    n_valid_levels = findfirst(<(ptop),p)
    p = p[begin:n_valid_levels]
    t = t[begin:n_valid_levels]
    r = r[begin:n_valid_levels]
    tvirtual_diff_parcel_env = similar(t)
    parcel_sat_vapor_pressure = get_saturation_vapor_pressure(tparcel)
    parcel_vapor_pressure = get_partial_vapor_pressure(rparcel,pparcel)
    parcel_rh = min(parcel_vapor_pressure/parcel_sat_vapor_pressure  , 1.0)
    parcel_specific_entropy = get_specific_entropy(tparcel,rparcel,pparcel)
    parcel_lcl = get_lifted_condensation_level(tparcel,parcel_rh,pparcel)
    @show parcel_lcl
    below_lcl = findall(>=(parcel_lcl),p)
    above_lcl = findall(<(parcel_lcl),p)

    #These two must populate buoyancy of lifted parcel_get_vapor_pressure
    #this would be adiabatic lifting, easy enough
    for level in below_lcl
        tlifted = tparcel*(p[level]/pparcel)^(Dryair.R/Dryair.cp)
        rlifted = rparcel
        tvirtual_lifted = get_virtual_temperature(tlifted,rlifted,rlifted)
        tvirtual_env = get_virtual_temperature(t[level],r[level],r[level])
        tvirtual_diff_parcel_env[level] = tvirtual_lifted - tvirtual_env
    end

    #We start with environmental values of temperature, mixing ratio, entropy etc
    for level in above_lcl
        niter = 0
        t_previousiter = t[level]
        saturation_vapor_pressure_previousiter = get_saturation_vapor_pressure(t_previousiter)
        mixing_ratio_previousiter = get_mixing_ratio(saturation_vapor_pressure_previousiter,p[level])
        t_currentiter = 0.0u"K"
        mixing_ratio_currentiter = 0.0u"g/g"
        while (abs(t_previousiter - t_currentiter) > 0.001u"K" )
            niter += 1
            t_currentiter = t_previousiter
            saturation_vapor_pressure_currentiter = get_saturation_vapor_pressure(t_currentiter)
            mixing_ratio_currentiter = get_mixing_ratio(saturation_vapor_pressure_currentiter,p[level])
            dsdt = (Dryair.cp + rparcel*Liquidwater.cp + Liquidwater.Lv*Liquidwater.Lv*mixing_ratio_currentiter/
            (Watervapor.R*t_currentiter*t_currentiter))/t_currentiter
            vapor_pressure_currentiter = get_partial_vapor_pressure(mixing_ratio_currentiter,p[level])
            entropy_currentinter = (Dryair.cp+rparcel*Liquidwater.cp)*log(t_currentiter/unit(t_currentiter)) - 
            Dryair.R*log((p[level]-vapor_pressure_currentiter)/unit(p[level])) + Liquidwater.Lv*mixing_ratio_currentiter / t_currentiter

            if niter < 3
                temperature_step = 0.3
            else
                temperature_step = 1
            end
            t_previousiter = t_currentiter + temperature_step*(parcel_specific_entropy - entropy_currentinter)/dsdt

            if (niter > 500 ) | (vapor_pressure_currentiter > ( p[level] - 1.0u"hPa") )
                error("Temperature didn't converge during lift")
            end
        end
        tvirtual_lifted = get_virtual_temperature(t_currentiter,rparcel,mixing_ratio_currentiter)
        tvirtual_env = get_virtual_temperature(t[level],r[level],r[level])
        tvirtual_diff_parcel_env[level] = tvirtual_lifted - tvirtual_env
    end
    return tvirtual_diff_parcel_env
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
    return (one(temperature) + 1e-3*epsilon*specific_humidity)*temperature
end

function get_virtual_temperature(temperature :: Quantity, specific_humidity :: Quantity)
    return (one(temperature) + 1e-3*epsilon*specific_humidity)*temperature
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
    return ustrip(g)/(1*ustrip(Dryair.cp)*SST)*(epsilon*ustrip(Dryair.cp)*SST/ustrip(Liquidwater.Lv)*latent_heat_flux) 
end

function surface_latent_heat_flux_to_buoyancy(SST :: Quantity, latent_heat_flux :: Quantity; rho = 1u"kg/m^3")
    return g/(1u"kg/m^3"*Dryair.cp*SST)*(epsilon*Dryair.cp*SST/Liquidwater.Lv*latent_heat_flux) 
end

function get_buoyancy(temperature_anomaly,mean_temperature) 
    return ustrip(g) * ustrip(temperature_anomaly)/ustrip(mean_temperature)
end

function get_buoyancy(temperature_anomaly :: Quantity ,mean_temperature :: Quantity)
    return g * temperature_anomaly/mean_temperature
end
