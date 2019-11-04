function compute_N2_old(var_Tv,xBar_Tv,z)

    dxBar_Tv_dz                = zeros(1,length(z),size(var_Tv,4))
    dxBar_Tv_dz[1,1:end-1,:]  .= (xBar_Tv[1,1,2:end,:].-xBar_Tv[1,1,1:end-1,:])./(z[2:end]-z[1:end-1])
    dxBar_Tv_dz[1,end,:]      .= dxBar_Tv_dz[1,end-1,:]
    N2                         = g .*(dxBar_Tv_dz .+ g/heat_capacity)[1,:,:]./xBar_Tv[1,1,:,:]  # Why is this using virtual and not potential, why g/cp check
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
    N2              = zeros(T,length(z),size(xBar_Tv,4))
    @views  N2[1:end-1,:]  .= (xBar_Tv[1,1,2:end,:].-xBar_Tv[1,1,1:end-1,:])./(z[2:end].-z[1:end-1])
    @views  N2[end,:]      .= N2[end-1,:]
    @views  @. N2  = g * (N2 + g/heat_capacity)/xBar_Tv[1,1,:,:]  # Why is this using virtual and not potential, why g/cp check
        bb1 = as_ints(findall((N2) .< 5e-5))
        bb = bb1[1,:]
        cc = bb1[2,:]
        @inbounds for i in 1:length(bb)
        if bb[i]>1
          @views  @. N2[bb[i],cc] = 0.5 * (N2[bb[i]-1,cc] + N2[bb[i]+1,cc]) # If N2 is small, substite by mean of neighbours
        else
          @views @.  N2[bb[i],cc] = N2[bb[i]+1,cc]
        end
    end
    return N2
end


#WIP
function compute_mse(T,z,qv)
    sz = size(T)
    return  dryair.cp*T .+ g*reshape(z,(1,1,sz[3],1)) .+ liquidwater.Lv*qv

end

function get_tendency(field :: AbstractArray{T,4}; dt = error("dt is required for the budget computation")) where T
    dfield_dt = similar(field)
    @. dfield_dt[:,:,:,1:end-1] = @views (field[:,:,:,2:end] - field[:,:,:,1:end-1]) / dt
    @. dfield_dt[:,:,:,end] = dfield_dt[:,:,:,end-1]
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

function integrate_vertically(field :: AbstractArray{T,4}; coord :: AbstractArray{T,1},weight :: AbstractArray{T,1}) where T
    sz = size(field)
    sc = size(coord)
    integral = zeros(eltype(field),(sz[1],sz[2],1,sz[4]))
    @inbounds for ind in CartesianIndices(field)
        if ind[3] != sz[3]
            integral[ind[1],ind[2],1,ind[4]] += weight[ind[3]]*(coord[ind[3]+1] - coord[ind[3]]) * field[ind]
        end
    end
    return integral
end


function integrate_vertically(field :: AbstractArray{T,4}; coord :: AbstractArray{T,1},weight :: AbstractArray{T,4}) where T
    sz = size(field)
    sc = size(coord)
    integral = zeros(eltype(field),(sz[1],sz[2],1,sz[4]))
    @inbounds for ind in CartesianIndices(field)
        if ind[3] != sz[3]
            integral[ind[1],ind[2],1,ind[4]] += weight[ind]*(coord[ind[3]+1] - coord[ind[3]]) * field[ind]
        end
    end
    return integral
end


function compute_virtual_temp(T,QV)
    return T.*(1 .+ 0.61QV)
end