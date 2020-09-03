"""
    distance(x1,x2,gridspacing :: Number)

Compute the cartesian distance between two points given their indices and the gridspacing. It asummes uniform grid.


"""
function distance(x1,x2,gridspacing :: Number,weight=1)
    return gridspacing*hypot( x2[1]-x1[1], x2[2]-x1[2] )
end

"""
    velocity_topolar(u,v,index,center)
Take a velocity vector, an origin of said vector and a center and return the tangential and azimuthal velocities with respect to that center.
"""
function velocity_topolar(u,v,index,center)
    pos = index .- center
    theta1 = atan(v,u)
    theta2 = atan(pos[2],pos[1])
    return -hypot(u,v) * cos(theta1 + theta2), hypot(u,v)*sin(theta1 + theta2)
end

function compute_N2_old_old(var_Tv,xBar_Tv,z)

    dxBar_Tv_dz                = zeros(1,length(z),size(var_Tv,4))
    dxBar_Tv_dz[1,1:end-1,:]  .= (xBar_Tv[1,1,2:end,:].-xBar_Tv[1,1,1:end-1,:])./(z[2:end]-z[1:end-1])
    dxBar_Tv_dz[1,end,:]      .= dxBar_Tv_dz[1,end-1,:]
    N2                         = g .*(dxBar_Tv_dz .+ g/Dryair.cp)[1,:,:]./xBar_Tv[1,1,:,:]  # Why is this using virtual and not potential, why g/cp check
    bb = findall(abs.(N2) .< 1e-6)
    for i=1:length(bb)
        if bb[i][1]>2
            N2[bb[i]] = 0.5 * (N2[bb[i][1]-1,bb[i][2]] + N2[bb[i][1]+1,bb[i][2]]) # If N2 is small, substite by mean of neighbours
        else
            N2[bb[i]] = N2[bb[i][1]+1,bb[i][2]]
        end
    end

    return N2
end

as_ints(a::AbstractArray{CartesianIndex{L}}) where L = reshape(reinterpret(Int, a), (L, size(a)...))



function compute_N2(xBar_Tv,z)
    T = eltype(xBar_Tv)
    N2 = zeros(T,length(z),size(xBar_Tv,4))
    @views  N2[1:end-1,:]  .= (xBar_Tv[1,1,2:end,:].-xBar_Tv[1,1,1:end-1,:])./(z[2:end].-z[1:end-1])
    @views  N2[end,:]      .= N2[end-1,:]
    @views  @. N2  = g * (N2 + g/Dryair.cp)/xBar_Tv[1,1,:,:] 
        bb1 = as_ints(findall(abs.(N2) .< 1e-6))
     
        bb = bb1[1,:]
        cc = bb1[2,:]
        
    for i in 1:length(bb)
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

function compute_N2_attempt(xBar_Tv,z)
    T = eltype(xBar_Tv)
    N2 = zeros(T,length(z),size(xBar_Tv,4))
    sz,st = size(N2)
    factor = g/Dryair.cp
    @inbounds for indt in 1:st, indz in 1:(sz - 1)
        dtvdz = (xBar_Tv[1,1,indz+1,indt] - xBar_Tv[1,1,indz,indt])/(z[indz+1] - z[indz])
        N2[indz,indt] = g*(dtvdz + factor)/xBar_Tv[1,1,indz,indt]
    end
    @inbounds for indt in 1:st
        N2[sz,indt] = N2[sz-1,indt]
    end
    bb1 = as_ints(findall(abs.(N2) .< 1e-6))
   
    bb = bb1[1,:]
    cc = bb1[2,:]
    for i in 1:length(bb)
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


