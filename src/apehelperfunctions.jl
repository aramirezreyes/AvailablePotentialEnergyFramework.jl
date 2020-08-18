

function cutborders!(array::Array{T,4},smooth_time,position) where T <: Real

    if position == 1
        array = array[:,:,:,1+div((smooth_time-1),2):end]
    elseif position == 2
        array = array[:,:,:,1+div((smooth_time-1),2):end-div((smooth_time-1),2)]
    elseif position == 3
        array = array[:,:,:,1:end - div((smooth_time-1),2)]
    end
    return array
end

function cutborders!(array::Array{T,3},smooth_time,position) where T <: Real
 if position == 1
        array = array[:,:,1+div((smooth_time-1),2):end]
    elseif position == 2
        array = array[:,:,1+div((smooth_time-1),2):end-div((smooth_time-1),2)]
    elseif position == 3
        array = array[:,:,1:end - div((smooth_time-1),2)]
    end
    return array
end

function cutborders!(array::Array{T,1},smooth_time,position) where T <: Real

    if position == 1
        array = array[1+div((smooth_time-1),2):end]
    elseif position == 2
        array = array[1+div((smooth_time-1),2):end-div((smooth_time-1),2)]
    elseif position == 3
        array = array[1:end-div((smooth_time-1),2)]
    end
    return array

end



function getsmoothdata!(U,V, W, Tv, ThetaV, RAD, Fs, smooth_x,smooth_y,smooth_time,position)
    buf3d = similar(U)
    buf2d = similar(Fs)
    filter_array!(buf3d,U,smooth_x,smooth_time,position)
    filter_array!(buf3d,V,smooth_x,smooth_time,position)
    filter_array!(buf3d,W,smooth_x,smooth_time,position)
    filter_array!(buf3d,Tv,smooth_x,smooth_time,position)
    filter_array!(buf3d,ThetaV,smooth_x,smooth_time,position)
    filter_array!(buf3d,RAD,smooth_x,smooth_time,position)
    filter_array!(buf2d,Fs,smooth_x,smooth_time,position)
    #return U, V, W, Tv, ThetaV, RAD, Fs
end

function getsmoothdata_nospace(U::Array{T,4},V::Array{T,4}, W::Array{T,4}, Tv::Array{T,4}, ThetaV::Array{T,4}, RAD::Array{T,4}, Fs::Array{T,3}, smooth_x::Int,smooth_y::Int,smooth_time::Int,position::Int) where T <: Real
    U = filter_array_nospace(U,smooth_x,smooth_time,position)
    V = filter_array_nospace(V,smooth_x,smooth_time,position)
    W = filter_array_nospace(W,smooth_x,smooth_time,position)
    Tv = filter_array_nospace(Tv,smooth_x,smooth_time,position)
    ThetaV = filter_array_nospace(ThetaV,smooth_x,smooth_time,position)
    RAD = filter_array_nospace(RAD,smooth_x,smooth_time,position)
    Fs = filter_array_nospace(Fs,smooth_x,smooth_time,position)
    return U, V, W, Tv, ThetaV, RAD, Fs
end

