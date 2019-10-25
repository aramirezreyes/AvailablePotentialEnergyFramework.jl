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
