mutable struct surf_quantities
    surf_pres
    surf_speed
    precip_water
end

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

mutable struct ape_budget
    int_APE           
    int_KE            
    int_APE_RAD       
    int_APE_DIA       
    int_APE_WN2       
    int_APE_Ub2       
    int_APE_Vb2       
    int_APE_rate      
    xBar_APE_Fs       
    int_APE_BL        
    int_KE_BL         
    int_APE_RAD_BL    
    int_APE_DIA_BL    
    int_APE_WN2_BL    
    int_APE_Ub2_BL    
    int_APE_Vb2_BL    
    int_APE_rate_BL   
    xBar_APE_Fs_BL    
#    res_APE           
#    FT                
#    BL                
    dia_a             
    rad_a             
    B_a               
    dia_ape           
    rad_ape           
    ape_budget() = new()
    
ape_budget(int_APE,int_KE,int_APE_RAD,int_APE_DIA,int_APE_WN2,int_APE_Ub2,int_APE_Vb2,int_APE_rate,xBar_APE_Fs,int_APE_BL,int_KE_BL,int_APE_RAD_BL,int_APE_DIA_BL,int_APE_WN2_BL,int_APE_Ub2_BL,int_APE_Vb2_BL,int_APE_rate_BL,xBar_APE_Fs_BL,dia_a,rad_a,B_a,dia_ape,rad_ape)=new(int_APE,int_KE,int_APE_RAD,int_APE_DIA,int_APE_WN2,int_APE_Ub2,int_APE_Vb2,int_APE_rate,xBar_APE_Fs,int_APE_BL,int_KE_BL,int_APE_RAD_BL,int_APE_DIA_BL,int_APE_WN2_BL,int_APE_Ub2_BL,int_APE_Vb2_BL,int_APE_rate_BL,xBar_APE_Fs_BL,dia_a,rad_a,B_a,dia_ape,rad_ape)
end

function cat_ape_budget(a::ape_budget,b::ape_budget)
    return ape_budget(
          cat(a.int_APE             ,b.int_APE            ,dims=1),
          cat(a.int_KE              ,b.int_KE             ,dims=1),
          cat(a.int_APE_RAD         ,b.int_APE_RAD        ,dims=1),
          cat(a.int_APE_DIA         ,b.int_APE_DIA        ,dims=1),
          cat(a.int_APE_WN2         ,b.int_APE_WN2        ,dims=1),
          cat(a.int_APE_Ub2         ,b.int_APE_Ub2        ,dims=1),
          cat(a.int_APE_Vb2         ,b.int_APE_Vb2        ,dims=1),
          cat(a.int_APE_rate        ,b.int_APE_rate       ,dims=1),
          cat(a.xBar_APE_Fs         ,b.xBar_APE_Fs        ,dims=1),
          cat(a.int_APE_BL          ,b.int_APE_BL         ,dims=1),
          cat(a.int_KE_BL           ,b.int_KE_BL          ,dims=1),
          cat(a.int_APE_RAD_BL      ,b.int_APE_RAD_BL     ,dims=1),
          cat(a.int_APE_DIA_BL      ,b.int_APE_DIA_BL     ,dims=1),
          cat(a.int_APE_WN2_BL      ,b.int_APE_WN2_BL     ,dims=1),
          cat(a.int_APE_Ub2_BL      ,b.int_APE_Ub2_BL     ,dims=1),
          cat(a.int_APE_Vb2_BL      ,b.int_APE_Vb2_BL     ,dims=1),
          cat(a.int_APE_rate_BL     ,b.int_APE_rate_BL    ,dims=1),
          cat(a.xBar_APE_Fs_BL      ,b.xBar_APE_Fs_BL     ,dims=1),
#           cat(a.res_APE             ,b.res_APE            ,dims=1),
#           cat(a.FT                  ,b.FT                 ,dims=1),
#           cat(a.BL                  ,b.BL                 ,dims=1),
          cat(a.dia_a               ,b.dia_a              ,dims=4),
          cat(a.rad_a               ,b.rad_a              ,dims=4),
          cat(a.B_a                 ,b.B_a                ,dims=4),
          cat(a.dia_ape             ,b.dia_ape            ,dims=4),
          cat(a.rad_ape             ,b.rad_ape            ,dims=4)
)
end
                                
function kernelt(window_t)
   if !isodd(window_t)
        window_t += 1
    end
    if !isodd(window_t)
        window_t += 1
    end
    kernel_t = centered(ones(window_t)/(window_t))
end

function kernel2d(window_h)
    if !isodd(window_h)
        window_h += 1
    end
    kernel_h = (ones(window_h)/(window_h))
    return  kernelfactors(( centered(kernel_h), centered(kernel_h) ))
end

function kernel4d(window_h,window_t)
    if !isodd(window_h)
        window_h += 1
    end
    if !isodd(window_t)
        window_t += 1
    end

    kernel_h = ones(window_h)/(window_h)
    return kernelfactors(( centered(kernel_h), centered(kernel_h), centered([1.0]), centered([1.0]) ))
end

function kernel3d(window_h,window_t)
   if !isodd(window_h)
        window_h += 1
    end
    if !isodd(window_t)
        window_t += 1
    end
    kernel_h = ones(window_h)/(window_h)
    return  kernelfactors(( centered(kernel_h), centered(kernel_h), centered([1.0]) ))
end

function kernel4d_t(window_t)
    if !isodd(window_t)
        window_t += 1
    end
    kernel_t = ones(window_t)/(window_t)

    return kernelfactors((centered([1.0]), centered([1.0]),centered([1.0]),centered(kernel_t)))
end

function kernel3d_t(window_t)
    if !isodd(window_t)
        window_t += 1
    end
    kernel_t = ones(window_t)/(window_t)
    return  kernelfactors((centered([1.0]), centered([1.0]), centered(kernel_t)))
end


function filter_array_2!(array::Array{T,4},smooth_x,smooth_time,position) where T <: Real
    filtered = similar(array)
    if position == 2
        filtered = imfilter(imfilter(array, kernel4d(smooth_x,smooth_time),"circular"),kernel4d_t(smooth_time),"inner")
    else
    ###Filtering in space
        for t in 1:size(array,4)
            for z in 1:size(array,3)
                filtered[:,:,z,t] = imfilter(array[:,:,z,t], kernel2d(smooth_x),"circular")
            end
        end
    ### filtering in time
        for z in 1:size(array,3)
            for y in 1:size(array,2)
                for x in 1:size(array,1)
                    array[x,y,z,:] = imfilter(filtered[x,y,z,:], kernelt(smooth_time),"replicate")
                end
            end
        end
    end
       # return array
end    

function filter_array_2!(array::Array{T,3},smooth_x,smooth_time,position) where T <: Real
    filtered = similar(array)
    if position == 2
        filtered = imfilter(imfilter(array, kernel4d(smooth_x,smooth_time),"circular"),kernel4d_t(smooth_time),"inner")
    else
    ###Filtering in space

        for t in 1:size(array,3)
            filtered[:,:,t] = imfilter(array[:,:,t], kernel2d(smooth_x),"circular")
        end

    ### filtering in time
        for y in 1:size(array,2)
            for x in 1:size(array,1)
                array[x,y,:] = imfilter(filtered[x,y,:], kernelt(smooth_time),"replicate")
            end
        end
    end
        #return array
end  

function filter_array!(buf::Array{T,4},array::Array{T,4},smooth_x,smooth_time,position) where T <: Real
    if position == 2
       error("Filter_array: Inner array is not implemented yer")
    else
        imfilter!(buf,array, kernel4d(smooth_x,smooth_time),"circular")
        imfilter!(array,buf,kernel4d_t(smooth_time),"symmetric")
    end
#        return filtered
end    

function filter_array!(buf::Array{T,3},array::Array{T,3},smooth_x,smooth_time,position) where T <: Real
    if position==2
       error("Filter_array: Inner array is not implemented yer")
    else
        imfilter!(buf,array, kernel3d(smooth_x,smooth_time),"circular")
        imfilter!(array,buf,kernel3d_t(smooth_time),"symmetric")
    end
end

                               
                                        
function filter_array_nospace(array::Array{T,4},smooth_x,smooth_time,position) where T <: Real
    if position==2

        filtered = imfilter(array,kernel4d_t(smooth_time),"inner")
    else
        
        filtered = imfilter(array,kernel4d_t(smooth_time),"replicate")
    end
    return filtered
end    

function filter_array_nospace(array::Array{T,3},smooth_x,smooth_time,position) where T <: Real
    if position==2
        filtered = imfilter(array,kernel3d_t(smooth_time),"inner")
    else
        
        filtered = imfilter(array,kernel3d_t(smooth_time),"replicate")
    end
    return filtered
end 

function filter_array_time(array::Array{T,1},window,position) where T <: Real
    kern = centered(ones(window)./window)
    if position==2
        smooth_t = imfilter(array,kern,"inner")
    else
        smooth_t = imfilter(array,kern)
    end
    return smooth_t
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

