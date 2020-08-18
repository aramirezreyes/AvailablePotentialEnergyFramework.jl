                                
function kernelt(window_t)
    if !isodd(window_t)
         window_t += 1
     end
     if !isodd(window_t)
         window_t += 1
     end
     kernel_t = centered(ones(window_t)./(window_t))
 end
 
 
 
 function kernel4d(window_h,window_t,T::Type=Float64)
     if !isodd(window_h)
         window_h += 1
     end
     if !isodd(window_t)
         window_t += 1
     end
 
     kernel_h = ones(T,window_h)./(window_h)
     return kernelfactors(( centered(kernel_h), centered(kernel_h), centered([T(1.0)]), centered([T(1.0)]) ))
 end
 
 function kernel3d(window_h,window_t,T::Type=Float64)
    if !isodd(window_h)
         window_h += 1
     end
     if !isodd(window_t)
         window_t += 1
     end
     kernel_h = ones(T,window_h)./(window_h)
     return  kernelfactors(( centered(kernel_h), centered(kernel_h), centered([T(1.0)]) ))
 end
 
 
 function kernel2d(window_h,T::Type=Float64)
    if !isodd(window_h)
         window_h += 1
     end
     if !isodd(window_t)
         window_t += 1
     end
     kernel_h = ones(T,window_h)./(window_h)
     return  kernelfactors(( centered(kernel_h), centered(kernel_h) ))
 end
 
 
 # function kernel4d_t(window_t)
 #     if !isodd(window_t)
 #         window_t += 1
 #     end
 #     kernel_t = ones(window_t)/(window_t)
 
 #     return kernelfactors((centered([1.0]), centered([1.0]),centered([1.0]),centered(kernel_t)))
 # end
 
 function kernel4d_t(window_t,T::Type=Float64)
     if !isodd(window_t)
         window_t += 1
     end
     kernel_t = ones(T,window_t)./(window_t)
 
     return kernelfactors(( centered([T(1.0)]), centered([T(1.0)]), centered([T(1.0)]), centered(kernel_t) ))
 end
 
 function kernel3d_t(window_t,T::Type=Float64)
     if !isodd(window_t)
         window_t += 1
     end
     kernel_t = ones(T,window_t)./(window_t)
     return  kernelfactors(( centered([T(1.0)]), centered([T(1.0)]), centered(kernel_t) ))
 end
 """
     filter_array_2!(array,smooth_x,smooth_t,position)
 
 Filters the input array _in-place_ using a moving mean. 
 
 In the first two dimensions the window width is smooth_x and the border is treated as circular (for doubly periodic domains). In the space dimension the width is smooth_t and the border value is replicated.
 
 This function calls under the hood the imfilter function of Images.jl
 
 """
 function filter_array_2!(array::Array{T,4},smooth_x,smooth_time,position) where T <: Real
 #    filtered = similar(array)
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
                     array[x,y,z,:] = imfilter(filtered[x,y,z,:], kernelt(smooth_time),"symmetric")
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
             filtered[:,:,t] = imfilter(filtered,array[:,:,t], kernel2d(smooth_x),"circular")
         end
 
     ### filtering in time
         for y in 1:size(array,2)
             for x in 1:size(array,1)
                 array[x,y,:] = imfilter(filtered[x,y,:], kernelt(smooth_time),"symmetric")
             end
         end
     end
         #return array
 end  
 """
     filter_array!(buffer,array,smooth_x,smooth_t,position)
 
 Filters the input array _in-place_ using a moving mean. In the first two dimensions the window width is smooth_x and the border is treated as circular (for doubly periodic domains). In the space dimension the width is smooth_t and the border value is replicated.
 
 This function calls under the hood the imfilter function of Images.jl
 
 The first argument must be a buffer of the same size of array.
 
 """
 function filter_array!(buf::Array{T,4},array::Array{T,4},smooth_x,smooth_time,position) where T <: Real
     if position == 2
        error("Filter_array: Inner array is not implemented yet")
     else
         imfilter!(buf,array,kernel4d(smooth_x,smooth_time,T)[1:2],"circular")
         imfilter!(array,buf,(kernel4d_t(smooth_time,T)[4],),"symmetric")
     end
 #        return filtered
 end    
 
 
 function filter_array!(buf::Array{T,3},array::Array{T,3},smooth_x,smooth_time,position) where T <: Real
     if position==2
        error("Filter_array: Inner array is not implemented yet")
     else
         imfilter!(buf,array, kernel3d(smooth_x,smooth_time,T)[1:2],"circular")
         imfilter!(array,buf,(kernel3d_t(smooth_time,T)[3],),"symmetric")
     end
 end
 
 """
 filter_array(array,smooth_x,smooth_t,position)
 
 Filters the input arrayusing a moving mean. In the first two dimensions the window width is smooth_x and the border is treated as circular (for doubly periodic domains). In the space dimension the width is smooth_t and the border value is replicated.
 
 This function calls under the hood the imfilter function of Images.jl
 
 The first argument must be a buffer of the same size of array.
 
 """
 
 function filter_array(array::Array{T,4},smooth_x,smooth_time,position) where T <: Real
     if position == 2
        error("Filter_array: Inner array is not implemented yet")
     else
         buf = similar(array)
         imfilter!(buf,array,kernel4d(smooth_x,smooth_time,T)[1:2],"circular")
         imfilter!(array,buf,(kernel4d_t(smooth_time,T)[4],),"symmetric")
     end
         return array
 end    
 
 
 function filter_array(array::Array{T,3},smooth_x,smooth_time,position) where T <: Real
     if position==2
        error("Filter_array: Inner array is not implemented yet")
     else
         buf = similar(array)
         imfilter!(buf,array, kernel3d(smooth_x,smooth_time,T)[1:2],"circular")
         imfilter!(array,buf,(kernel3d_t(smooth_time,T)[3],),"symmetric")
     end
     return array
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
         smooth_t = imfilter(array,kern,"symmetric")
     end
     return smooth_t
 end