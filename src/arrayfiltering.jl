"""
    kernel_1d(window)    
Create a kernel to perform a moving average filtering of a 1-d array. It will use window as windows size unless the given vaule is even, in which case it will add one to the given value before creating the kernel. This kernel has to be passed into imfilter or imfilter!    
"""                                
function kernel_1d(window)
    if !isodd(window)
         window += 1
     end
     kernel_t = centered(ones(window)./(window))
 end
 
 
 """
    kernel_4d(window_h,window_t,T :: Type :: Float64)    
Create a kernel to perform a moving average filtering of a 4-d array along the dimensions 1,2 and 4. It will use window_h (dimensions 1 and 2) and window_t (dimension 4)  as windows sizes unless the given vaules are even, in which case it will add one to the given value before creating the kernel. This kernel has to be passed into imfilter or imfilter!    
""" 
 function kernel_4d(window_h,window_t,T::Type=Float64)
    kernel_h = if (window_h == 1) | (window_h == false)
        centered([T(1.0)])
    else
        !isodd(window_h) && (window_h += 1 )
        ones(T,window_h)./(window_h)
    end
     if !isodd(window_t)
         window_t += 1
     end
     
    kernel_t = if (window_t == 1) | (window_t == false)
        centered([T(1.0)])
    else
        !isodd(window_t) && (window_t += 1 )
        ones(T,window_t)./(window_t)
    end
    
     return kernelfactors(( kernel_h, kernel_h, centered([T(1.0)]) , kernel_t ))
 end

 """
kernel_3d(window_h,window_t,T :: Type :: Float64)    
Create a kernel to perform a moving average filtering of a 3-d array along the dimensions 1,2 and 3. It will use window_h (dimensions 1 and 2) and window_t (dimension 3) as windows sizes unless the given vaules are even, in which case it will add one to the given value before creating the kernel. This kernel has to be passed into imfilter or imfilter!
"""
 function kernel_3d(window_h,window_t,T::Type=Float64)
    if !isodd(window_h)
         window_h += 1
     end
     if !isodd(window_t)
         window_t += 1
     end
     kernel_h = ones(T,window_h)./(window_h)
     return  kernelfactors(( centered(kernel_h), centered(kernel_h), centered([T(1.0)]) ))
 end
 
 """
 kernel_2d(window_h,T :: Type :: Float64)    
 Create a kernel to perform a moving average filtering of a 2-d array along the dimensions 1,2. It will use window_h as windows size unless the given vaule is even, in which case it will add one to the given value before creating the kernel. This kernel has to be passed into imfilter or imfilter!
 """
 function kernel_2d(window_h,T::Type=Float64)
    if !isodd(window_h)
         window_h += 1
     end
     kernel_h = ones(T,window_h)./(window_h)
     return  kernelfactors(( centered(kernel_h), centered(kernel_h) ))
 end
  
 """
 kernel_4d_t(window_t,T :: Type :: Float64)    
 Create a kernel to perform a moving average filtering of a 4-d array along the dimension 4. It will use window_t as windows size unless the given vaule is even, in which case it will add one to the given value before creating the kernel. This kernel has to be passed into imfilter or imfilter!
 """
 function kernel_4d_t(window_t,T::Type=Float64)
     if !isodd(window_t)
         window_t += 1
     end
     kernel_t = ones(T,window_t)./(window_t)
 
     return kernelfactors(( centered([T(1.0)]), centered([T(1.0)]), centered([T(1.0)]), centered(kernel_t) ))
 end
 
 """
 kernel_3d_t(window_t,T :: Type :: Float64)    
 Create a kernel to perform a moving average filtering of a 3-d array along the dimension 3. It will use window_t as windows size unless the given vaule is even, in which case it will add one to the given value before creating the kernel. This kernel has to be passed into imfilter or imfilter!
 """
 function kernel_3d_t(window_t,T::Type=Float64)
     if !isodd(window_t)
         window_t += 1
     end
     kernel_t = ones(T,window_t)./(window_t)
     return  kernelfactors(( centered([T(1.0)]), centered([T(1.0)]), centered(kernel_t) ))
 end

 """
     filter_array_2!(array,smooth_x,smooth_t,position)
 
 Filter the input array _in-place_ using a moving mean. 
 
 In the first two dimensions the window width is smooth_x and the border is treated as circular (for doubly periodic domains). In the space dimension the width is smooth_t and the border value is replicated.
 
 This function calls under the hood the imfilter function of Images.jl 
 
 """
function filter_array_2!(buf :: Array{T,4},array::Array{T,4},smooth_x,smooth_time,position) where T <: Real
    sx,sy,sz,st = size(array)
#    filtered = similar(array)
     if position == 2
         filtered = imfilter(imfilter(array, kernel4d(smooth_x,smooth_time),"circular"),kernel_4d_t(smooth_time),"inner")
     else
         ###Filtering in space
        
         @info typeof(array) typeof(buf)
         @inbounds for t in 1:st
             @Threads.threads for z in 1:sz
                 imfilter!(view(buf,:,:,z,t),array[:,:,z,t], kernel_2d(smooth_x),"circular")
             end
         end
         @info typeof(array) typeof(buf)         
     ### filtering in time
         @inbounds for z in 1:sz, y in 1:sy
             @Threads.threads for x in 1:sx
                 array[x,y,z,:] .= imfilter(buf[x,y,z,:], kernel_1d(smooth_time),"symmetric")
             end
         end
     end
    return array
end    
 
 """
 filter_array_2!(array,smooth_x,smooth_t,position)

Filters the input array _in-place_ using a moving mean. 

In the first two dimensions the window width is smooth_x and the border is treated as circular (for doubly periodic domains). In the space dimension the width is smooth_t and the border value is replicated.

This function calls under the hood the imfilter function of Images.jl

"""
 function filter_array_2!(array::Array{T,3},smooth_x,smooth_time,position) where T <: Real
     filtered = similar(array)
     sx,sy,st = size(array)
     if position == 2
         filtered = imfilter(imfilter(array, kernel4d(smooth_x,smooth_time),"circular"),kernel4d_t(smooth_time),"inner")
     else
         ###Filtering in space
         for t in 1:st
             filtered[:,:,t] = imfilter(filtered,array[:,:,t], kernel_2d(smooth_x),"circular")
         end
         ### filtering in time
         for y in 1:sy, x in 1:sx
             array[x,y,:] = imfilter(filtered[x,y,:], kernel_1d(smooth_time),"symmetric")
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
         imfilter!(buf,array,kernel_4d(smooth_x,smooth_time,T)[1:2],"circular")
         imfilter!(array,buf,(kernel_4d_t(smooth_time,T)[4],),"symmetric")
     end
 #        return filtered
 end    
 
 
 function filter_array!(buf::Array{T,3},array::Array{T,3},smooth_x,smooth_time,position) where T <: Real
     if position==2
        error("Filter_array: Inner array is not implemented yet")
     else
         imfilter!(buf,array, kernel_3d(smooth_x,smooth_time,T)[1:2],"circular")
         imfilter!(array,buf,(kernel_3d_t(smooth_time,T)[3],),"symmetric")
     end
 end
 
 """
 filter_array(array,smooth_x,smooth_t,position)
 
 Filters the input arrayusing a moving mean. In the first two dimensions the window width is smooth_x and the border is treated as circular (for doubly periodic domains). In the space dimension the width is smooth_t and the border value is replicated.
 
 This function calls under the hood the imfilter function of Images.jl
 
 The first argument must be a buffer of the same size of array.
 
 """
 function filter_array(array::Array{T,4},smooth_x,smooth_time,position = 1) where T <: Real
     if position == 2
        asize = size(array)
        tail_size = (asize[4]-1)÷2
        buf1 = Array{T,4}(undef,asize[1],asize[2],asize[3],asize[4] - 2*tail_size)
        buf2 = Array{T,4}(undef,asize[1],asize[2],asize[3],asize[4] - 2*tail_size)
        buf_axes = axes(buf1)
        imfilter!(OffsetArray(buf1,buf_axes[1],buf_axes[2],buf_axes[3],buf_axes[4].+tail_size),array,(kernel_4d_t(smooth_time,T)[3],),Inner())
        imfilter!(buf2,buf1, kernel_4d(smooth_x,smooth_time,T)[1:2],"circular")
     else
         buf1 = similar(array)
         buf2 = similar(array)
         imfilter!(buf1,array,kernel_4d(smooth_x,smooth_time,T)[1:2],"circular")
         imfilter!(buf2,buf1,(kernel_4d_t(smooth_time,T)[4],),"symmetric")
     end
         return buf2
 end    
 
 """
 filter_array(array,smooth_x,smooth_t,position)
 
 Filters the input arrayusing a moving mean. In the first two dimensions the window width is smooth_x and the border is treated as circular (for doubly periodic domains). In the space dimension the width is smooth_t and the border value is replicated.
 
 This function calls under the hood the imfilter function of Images.jl
 
 The first argument must be a buffer of the same size of array.
 
 """
function filter_array(array::Array{T,3},smooth_x,smooth_time,position = 1) where T <: Real
    if !isodd(smooth_time)
         smooth_time += 1
    end
    if position==2
        asize = size(array)
        tail_size = (asize[3]-1)÷2
        buf1 = Array{T,3}(undef,asize[1],asize[2],asize[3] - 2*tail_size)
        buf2 = Array{T,3}(undef,asize[1],asize[2],asize[3] - 2*tail_size)
        buf_axes = axes(buf1)
        imfilter!(OffsetArray(buf1,buf_axes[1],buf_axes[2],buf_axes[3].+tail_size),array,(kernel_3d_t(smooth_time,T)[3],),Inner())
        imfilter!(buf2,buf1, kernel_3d(smooth_x,smooth_time,T)[1:2],"circular")

     else
        buf1 = similar(array)
        buf2 = similar(array)
        imfilter!(buf1,array, kernel_3d(smooth_x,smooth_time,T)[1:2],"circular")
        imfilter!(buf2,buf1,(kernel_3d_t(smooth_time,T)[3],),"symmetric")
     end
     return buf2
 end
 
 """
 filter_array_nospace(array,smooth_t,position)
 
 Filters the input array using a moving mean along the fourth dimension. In the 4th dimension the width is smooth_t and the border value is replicated except if position = 2, in which case it only takes the inner part of the smoothed array.
 
 This function calls under the hood the imfilter function of Images.jl
 """
 function filter_array_nospace(array::Array{T,4},smooth_time,position = 1) where T <: Real
     if position==2
         filtered = imfilter(array,kernel_4d_t(smooth_time),"inner")
     else     
         filtered = imfilter(array,kernel_4d_t(smooth_time),"replicate")
     end
     return filtered
 end    
 """
 filter_array_nospace(array,smooth_t,position)
 
 Filters the input array using a moving mean along the third dimension. In the 3th dimension the width is smooth_t and the border value is replicated except if position = 2, in which case it only takes the inner part of the smoothed array.
 
 This function calls under the hood the imfilter function of Images.jl
 """
 function filter_array_nospace(array::Array{T,3},smooth_time,position = 1) where T <: Real
     if position==2
         filtered = imfilter(array,kernel_3d_t(smooth_time),"inner")
     else
         
         filtered = imfilter(array,kernel_3d_t(smooth_time),"replicate")
     end
     return filtered
 end 
 """
 filter_array_nospace(array,smooth_t,position=1)
 
 Filters the input, 1-d array. The border value is replicated except if position = 2, in which case it only takes the inner part of the smoothed array.
 
 This function calls under the hood the imfilter function of Images.jl
 """
 function filter_array_time(array::Array{T,1},window,position = 1) where T <: Real
     kern = centered(ones(window)./window)
     if position==2
         smooth_t = imfilter(array,kern,"inner")
     else
         smooth_t = imfilter(array,kern,"symmetric")
     end
     return smooth_t
 end
