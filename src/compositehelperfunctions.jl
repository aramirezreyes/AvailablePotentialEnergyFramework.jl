using Images: findlocalminima
using ImageSegmentation: SegmentedImage, segment_labels, region_adjacency_graph
using SparseArrays

"""
    distance(x1,x2,gridspacing :: Number)

Compute the cartesian distance between two points given their indices and the gridspacing. It asummes uniform grid.


"""
function distance(x1,x2,gridspacing :: Number,weight=1)
    return gridspacing*hypot( x2[1]-x1[1], x2[2]-x1[2] )
end
 

"""
    isindexindistancebin(binlimits,index,center = (0,0),gridspacing=1) = (binlimits[1] < distance(index,center,gridspacing) <= binlimits[2]) ? true : false
computes the distance of one index to the origin and returns true if that distance is inside a bin
"""
isindexindistancebin(binlimits,index,center,gridspacing=1) = (binlimits[1] < distance(index,center,gridspacing) <= binlimits[2]) ? true : false

"""
    averageallindistance(radiusbin,array :: Array{T,2},mask,center,gridspacing = 1)
Create an average of the quantity in array at all the points located between radiusbin[1] and radiusbin[2] from a center.
The points should be masked by a boolean array. It assumes a uniform gridspacing.
"""
function averageallindistance(radiusbin,array :: Array{T,2},center,gridspacing = 1) where T
    average = 0.0
    count = 0    
    @inbounds for index in CartesianIndices(array)
        if isindexindistancebin(radiusbin,index,center,gridspacing)
                count += 1
                average += array[index]
        end 
    end
    if !iszero(count)
        return average./count
    else
        return average
    end
end 

"""
    averageallindistance(radiusbin,array :: Array{T,3},mask,center,gridspacing = 1)
Create an average of the quantity in array at all the points located between radiusbin[1] and radiusbin[2] from a center.
The points should be masked by a boolean array. It assumes a uniform gridspacing.
"""
function averageallindistance(radiusbin,array :: Array{T,3},center,gridspacing = 1) where T
    average = zeros(size(array,3))
    sx,sy,sz = size(array)
    count = 0
    @inbounds for index in CartesianIndices((1:sx,1:sy))
        if isindexindistancebin(radiusbin,index,center,gridspacing)
            count += 1
            @views average[:] .+= array[index,:]
        end
    end
    if !iszero(count)
        return average./count
    else
        return average
    end
end

"""
    averageallindistance!(radiusbin,array :: Array{T,3},mask,center,gridspacing = 1)
Create an average of the quantity in array at all the points located between radiusbin[1] and radiusbin[2] from a center.
The points should be masked by a boolean array. It assumes a uniform gridspacing.
"""
function averageallindistance!(average,radiusbin,array :: Array{T,3},center,gridspacing = 1) where T
    sx,sy,sz = size(array)
    count = 0
    @inbounds for index in CartesianIndices((1:sx,1:sy))
        if isindexindistancebin(radiusbin,index,center,gridspacing)
            count += 1
            @views average[:] .+= array[index,:]
        end
    end
    if !iszero(count)
        return average./count
    else
        return average
    end
end


"""

"""
function velocity_topolar(u,v,index,center)
    pos = index .- center
    theta1 = atan(v,u)
    theta2 = atan(pos[2],pos[1])
    return hypot(u,v) * cos(theta1 + theta2), hypot(u,v)*sin(theta1 + theta2)
end

"""
    add_allcyclones(buf :: Array{T,2},array :: Array{T,2},radius_bins,array :: Array{T,3},segmentedcyclones,cyclonescenters,gridspacing)

Compute the azimuthal average of some quantity around a center. Repeats the process and averages about all the tropical cyclones detected on the array.
It receives an array with the radius bins to use,the field to average, called `array`, each cyclone as a SegmentedImage,the centers of the cyclones and the gridspacing.
"""
function add_allcyclones!(addition,buf,array,segmentedcyclones,cyclonescenters,maskcyclones = true) 
    if !(size(addition) == size(buf) == size(array))
        DimensionMismatch("Addition, buffer and array must all have the same dimensions")
    end
    center = size(array)[1:2] .÷ 2
    if maskcyclones
        G, vert_map = region_adjacency_graph(segmentedcyclones, (i,j)->1)
        labelsmap = labels_map(segmentedcyclones)
        cyclonecount = 0
        for cyclone in 1:(length(segmentedcyclones.segment_labels)-1)
            if !isinteracting(segmentedcyclones,cyclone)
                cyclonecount += 1
                addition .+=  shifter!(buf,array.*(labelsmap .== cyclone),center,cyclonescenters[cyclone][1])
            end
        end
    else
        G, vert_map = region_adjacency_graph(segmentedcyclones, (i,j)->1)
        cyclonecount = 0
        for cyclone in 1:(length(segmentedcyclones.segment_labels)-1)
            if !isinteracting(segmentedcyclones,cyclone)           
                cyclonecount += 1
                addition .+=  shifter!(buf,array,center,cyclonescenters[cyclone][1])
            end
        end
    end
    return cyclonecount

end


"""
    azimuthalaverage_allcyclones(radius_bins,array :: Array{T,3},segmentedcyclones,cyclonescenters,gridspacing)

Compute the azimuthal average of some quantity around a center. Repeats the process and averages about all the tropical cyclones detected on the array.
It receives an array with the radius bins to use,the field to average, called `array`, each cyclone as a SegmentedImage,the centers of the cyclones and the gridspacing.
"""
function  azimuthalaverage_allcyclones(radius_bins,array :: Array{T,3},segmentedcyclones,cyclonescenters,gridspacing) where T
    G, vert_map = region_adjacency_graph(segmentedcyclones, (i,j)->1)
    labelsmap = labels_map(segmentedcyclones)
    adjacencymatrix = G.weights
    average = zeros(size(array,3),length(radius_bins)-1)
    cyclonecount = 0
    for cyclone in 1:(length(segmentedcyclones.segment_labels)-1)
        if !isinteracting(segmentedcyclones,cyclone)
            cyclonecount += 1
            for bin in 1:(length(radius_bins) - 1)
                average[:,bin] .+= averageallindistance((radius_bins[bin],radius_bins[bin+1]),array,labelsmap.==cyclone,cyclonescenters[cyclone][1],gridspacing)
            end
        end
    end
    if !iszero(cyclonecount)
        return    average./cyclonecount
    else
        return average
    end
end

"""
    azimuthalaverage_allcyclones(radius_bins,array :: Array{T,3},segmentedcyclones,cyclonescenters,gridspacing)

Compute the azimuthal average of some quantity around a center. Repeats the process and averages about all the tropical cyclones detected on the array.
It receives an array with the radius bins to use,the field to average, called `array`, each cyclone as a SegmentedImage,the centers of the cyclones and the gridspacing.
"""
function  azimuthalaverage_allcyclones(radius_bins,array :: Array{T,2},segmentedcyclones,cyclonescenters,gridspacing) where T
    average = zeros(length(radius_bins)-1)
    G, vert_map = region_adjacency_graph(segmentedcyclones, (i,j)->1)
    labelsmap = labels_map(segmentedcyclones)
    adjacencymatrix = G.weights
    cyclonecount = 0
    for cyclone in 1:(length(segmentedcyclones.segment_labels)-1)
        if !isinteracting(segmentedcyclones,cyclone)
            cyclonecount += 1
            for bin in 1:(length(radius_bins) - 1)
                average[bin] += averageallindistance((radius_bins[bin],radius_bins[bin+1]),array,labelsmap.==cyclone,cyclonescenters[cyclone][1],gridspacing)
            end
        end
    end
    if !iszero(cyclonecount)
        return    average./cyclonecount
    else
        return average
    end
end

"""
    neighbours_2d(arraysize,indicesCartesian)
Compute de neighboring points of an index considering periodic boundaries
"""
function neighbours_2d(arraysize,indicesCartesian)
    decrement_mod1(num, n) = num == 1 ? n : num - 1
    increment_mod1(num, n) = num == n ? one(num) : num + 1
    x_ind = indicesCartesian[1]
    y_ind = indicesCartesian[2]
    size_x = arraysize[1]
    size_y = arraysize[2]
    @inbounds neighbours = Iterators.product(
        (decrement_mod1(x_ind,size_x), x_ind, increment_mod1(x_ind,size_x)),
        (decrement_mod1(y_ind,size_y), y_ind, increment_mod1(y_ind,size_y)) )
    return CartesianIndex.(neighbours)
end

"""
    function detect_cyclones(surface_pressure,pressure_threshold,resolution)

"""
function detect_cyclones(surface_pressure,pressure_threshold,resolution)
    neighbourhood_gen(arraysize) = point -> AvailablePotentialEnergyFramework.neighbours_2d(arraysize,point)
    pres_anomaly = surface_pressure .- mean(surface_pressure,dims=(1,2))
    centers_and_labels = findcyclonecenters_aspressureminima(pres_anomaly,pressure_threshold,resolution)
    #@info centers_and_labels
    if length(centers_and_labels) == 0
        return (nothing,nothing)
    end
    centers_and_labels = [(centers_and_labels[i],i) for i in 1:length(centers_and_labels)]
    mask = pres_anomaly .<= pressure_threshold
    push!(centers_and_labels,(findfirst(!,mask),1000))
    cyclones = seeded_region_growing(pres_anomaly.*mask,centers_and_labels,neighbourhood_gen(size(surface_pressure)))
    return (centers_and_labels,cyclones)
end

"""
    isinteracting(adjacencyMatrix :: SparseMatrixCSC,regionnumber)
Computes the adjacency of identified cyclones and returns true if the cylone-number is adjacent to a region other than the background.
"""
function isinteracting(adjacencyMatrix :: SparseMatrixCSC,regionnumber)
    @inbounds for i in 1:(adjacencyMatrix.m - 1)
#        @views adjacencyMatrix[regionnumber,i] == 1 ? (return true) : nothing
        @views adjacencyMatrix[regionnumber,i] == 1  && (return true) 
    end
    return false
end

"""
    isinteracting(cyclones :: SegmentedImage,regionnumber)
Computes the adjacency of identified cyclones and returns true if the cylone-number is adjacent to a region other than the background.
"""
function isinteracting(cyclones :: SegmentedImage,regionnumber)
    adjacency, vert_map = region_adjacency_graph(cyclones,(i,j)->1)
    @inbounds for i in 1:(adjacency.weights.m - 1)
        adjacency.weights[regionnumber,i] == 1 ? (return true) : nothing
    end
    return false
end

"""
    findcyclonecenters_aspressureminima(surf_pres_anomaly,detection_threshold)

Takes a surface pressure anomaly array surf_pres_anomaly[x,y] = surf_pres[x,y] .- mean(surf_pres,dims(1,2))
a detection threshold for the anomaly, and return an array of tuples (x,y) 
where each tuple represents the centers of cyclones identified as the minima of the anomaly.
"""
function findcyclonecenters_aspressureminima(surf_pres_anomaly,detection_threshold,grid_spacing=5000)
    surf_pres_filtered = smoothfilter(surf_pres_anomaly,detection_threshold,grid_spacing);
    peaks              = findlocalminima(surf_pres_filtered,[1,2],false);
end

"""
    smoothfilter(surf_pres_anomaly, threshold = 9)

Takes a 2d array and smooths it using a rolling median filter. 
It then returns the elements of the filtered array whose values are less than the threshold.
"""
function smoothfilter(surf_pres_anomaly,treshold=-9, resolution = 2000)
    windowsize = 30000 ÷ (2*resolution) # Length order of magnitude of the eye
    windowsize % 2 == 0 ? windowsize = windowsize + 1 : windowsize
    surf_pres_median = mapwindow(median!,surf_pres_anomaly,[windowsize,windowsize],border="circular");
    surf_pres_median = surf_pres_median.*(surf_pres_median.<treshold);
    surf_pres_filtered = imfilter(surf_pres_median,Kernel.gaussian(3),"circular");
    surf_pres_filtered = surf_pres_filtered.*(surf_pres_filtered.<treshold);
end

"""
    shifter(array,domain_center,peak)

Returns an array in which a pressure perturbation center is displaced to the center of the domain using circshift.
Using it may assume periodic domain.
Receives and SAM 3D+time or 2D+time array and two tuples, the (x,y) indices of the domain center,
and the (x,y) indices of the location of the pressure perturbation peak.
"""
function shifter(array,domain_center,peak)
    if ndims(array)==2
      return  circshift(array[:,:],[domain_center[1]-peak[1],domain_center[2]-peak[2]]);
    elseif ndims(array)==3
      return  circshift(array[:,:,:],[domain_center[1]-peak[1],domain_center[2]-peak[2],0]);
    end
end
"""
    shifter!(dest,array,time,domain_center,peak)

Stores in `dest` an array in which a pressure perturbation center is displaced to the center of the domain using circshift.
Using it may assume periodic domain.
Receives and SAM 3D+time or 2D+time array and two tuples, the (x,y) indices of the domain center,
and the (x,y) indices of the location of the pressure perturbation peak.
"""
function shifter!(dest,array,domain_center,peak)
#    if 
        return  circshift!(dest,array,[domain_center[1]-peak[1],domain_center[2]-peak[2]]);
#    elseif ndims(array)==3
    #     return  circshift!(dest,array,[domain_center[1]-peak[1],domain_center[2]-peak[2]]);
    # elseif ndims(array)==4
    #   return  circshift!(dest,array,[domain_center[1]-peak[1],domain_center[2]-peak[2],0]);
#    end
end

"""
    timemean_nofalseframe(input)

Takes a 3-d or 3-d array and computes the average along the third or fourth dimension, skipping the slices input[:,:,:,i] for which the maxiumum value is 0.0

"""
function timemean_nofalseframe(input::Array{T,4}) where {T<:Real}
        sx,sy,sz,st = size(input)
        copy = zeros(T,sx,sy,sz)
        trueframes = 1
        for t in 1:size(input,4)
            if maximum(view(input,:,:,:,t)) != 0.0
               copy[:,:,:] .+= @views  input[:,:,:,t]
               trueframes = trueframes + 1
            end
        end
        return copy .= copy ./ trueframes
end
    
function timemean_nofalseframe(input::Array{T,3}) where {T<:Real}
        sx,sy,st = size(input)
        copy = zeros(T,sx,sy)
        trueframes = 1
        for t in 1:size(input,4)
            if maximum(view(input,:,:,t)) != 0.0
               copy[:,:] .+= @views  input[:,:,t]
               trueframes = trueframes + 1
            end
        end
        return copy .= copy ./ trueframes
end
    
function removefalseframes(input::Array{T,4},peaktimes) where {T<:Real}
        sizet = length(unique(peaktimes))
        copy = zeros(T,size(input,1),size(input,2),size(input,3),sizet)
        trueframes = 1
        Threads.@threads for t in 1:size(input,4)
            if maximum(view(input,:,:,:,t)) != 0
               @views copy[:,:,:,trueframes] .= input[:,:,:,t]
               trueframes = trueframes + 1
            end
        end
        return copy
end
 
   
function removefalseframes(input::Array{T,3},peaktimes) where {T<:Real}
        sizet = length(unique(peaktimes))
        copy = zeros(T,size(input,1),size(input,2),sizet)
        trueframes = 1
        Threads.@threads for t in 1:size(input,3)
            if maximum(view(input,:,:,t)) != 0
               @views copy[:,:,trueframes] .= input[:,:,t]
               trueframes = trueframes + 1
            end
        end
        return copy
end