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


function compute_N2(var_Tv,xBar_Tv,z)

    N2              = zeros(length(z),size(var_Tv,4))
    N2[1:end-1,:]  .= (xBar_Tv[1,1,2:end,:].-xBar_Tv[1,1,1:end-1,:])./(z[2:end].-z[1:end-1])
    N2[end,:]      .= N2[end-1,:]

     @. N2                         = g * (N2 + g/heat_capacity)/xBar_Tv[1,1,:,:]  # Why is this using virtual and not potential, why g/cp check
    @inbounds for i in findall(abs.(N2) .< 1e-6)
        if i[1]>1
            N2[i] = 0.5 * (N2[i[1]-1,i[2]] + N2[i[1]+1,i[2]]) # If N2 is small, substite by mean of neighbours
        else
            N2[i] = N2[i[1]+1,i[2]] 
        end
    end
    return N2
end
