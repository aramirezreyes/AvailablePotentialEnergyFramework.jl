using Images: findlocalminima
using ImageSegmentation: SegmentedImage, segment_labels, region_adjacency_graph
using SparseArrays

"""
    distance(x1,x2,gridspacing :: Number)

Compute the cartesian distance between two points given their indices and the gridspacing. It asummes uniform grid.


"""
function distance(x1,x2,gridspacing :: Number,weight=1)
    return gridspacing*sqrt( (x2[1]-x1[1])^2 + (x2[2]-x1[2])^2 )
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
function averageallindistance(radiusbin,array :: Array{T,2},mask,center,gridspacing = 1) where T
    average = 0.0
    count = 0    
    @inbounds for index in CartesianIndices(mask)
        if mask[index] && isindexindistancebin(radiusbin,index,center,gridspacing)
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
function averageallindistance(radiusbin,array :: Array{T,3},mask,center,gridspacing = 1) where T
    average = zeros(size(array,3))
    count = 0    
    @inbounds for index in CartesianIndices(mask)
        if mask[index] && isindexindistancebin(radiusbin,index,center,gridspacing)
            count += 1
            average .+= array[index,:]
        end 
    end
    if !iszero(count)
        return    average./count
    else
        return average
    end
end 

"""
   azimuthalaverage_allcyclones!(radiusbins,average,buf,buf2,array,segmentedcyclones,cyclonescenters,gridspacing)

Compute the azimuthal average of some quantity around a center. It receives the field to average, called `array`, each cyclone as a SegmentedImage,the centers of the cyclones and the gridspacing.
It receives 4 arrays to mutate, they are listed below.
# Arguments
- `radius :: Array{Float,1}`: contains the list of radius
- `average :: Array{Float,ndims(array) - 1}` the array to store the averaged quantity
- `buf :: Array{Float,2}` a buffer to recenter the pressure field so that the cyclone in turn appear at the center of the domain
- `buf2 :: Array{Float,ndims(array)} ` 

"""
function azimuthalaverage_allcyclones!(radiuses,averages,buf1,buf2,array :: Array{T,3},segmentedcyclones,cyclonescenters,gridspacing)  where T
    G, vert_map = region_adjacency_graph(segmentedcyclones, (i,j)->1)
    labelsmap = labels_map(segmentedcyclones)
    adjacencymatrix = G.weights
    for cyclone in 1:(length(segmentedcyclones.segment_labels)-1)
        if !isinteracting(adjacencymatrix,cyclone)
            shifter!(buf1,array,size(array).÷2,cyclonescenters[cyclone][1])
            shifter!(buf2,labelsmap,size(array).÷2,cyclonescenters[cyclone][1])
            for index in findall((==)(cyclone),buf2)
                @info index,cyclonescenters[cyclone][1]
                distfromcenter = distance(index,cyclonescenters[cyclone][1],2000)
                radiusisatindex = findall((==)(distfromcenter),radiuses)
                if (length(radiusisatindex) == 0)
                    push!(radiuses,distfromcenter)
                end 
                return nothing
            end
        end
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
function findcyclonecenters_aspressureminima(surf_pres_anomaly,detection_threshold,resolution)
    surf_pres_filtered = smoothfilter(surf_pres_anomaly,detection_threshold,resolution);
    peaks              = findlocalminima(surf_pres_filtered,[1,2],false);
end

"""
    smoothfilter(surf_pres_anomaly, threshold = 9)

Takes a 2d array and smooths it using a rolling median filter. 
It then returns the elements of the filtered array whose values are less than the threshold.
"""
function smoothfilter(surf_pres_anomaly,treshold=9, resolution = 2000)
    windowsize = 30000 ÷ (2*resolution) # Length order of magnitude of the eye
    windowsize % 2 == 0 ? windowsize = windowsize + 1 : windowsize
    surf_pres_median = mapwindow(median!,surf_pres_anomaly,[windowsize,windowsize]);
    surf_pres_median = surf_pres_median.*(surf_pres_median.<treshold);
    surf_pres_filtered = imfilter(surf_pres_median,Kernel.gaussian(3));
    surf_pres_filtered = surf_pres_filtered.*(surf_pres_filtered.<treshold);
end

"""
    shifter(array,time,domain_center,peak)

Returns an array in which a pressure perturbation center is displaced to the center of the domain using circshift.
Using it may assume periodic domain.
Receives and SAM 3D+time or 2D+time array and two tuples, the (x,y) indices of the domain center,
and the (x,y) indices of the location of the pressure perturbation peak.
"""
function shifter(array,time,domain_center,peak)
    if ndims(array)==3
      return  circshift(array[:,:,time],[domain_center[1]-peak[1],domain_center[2]-peak[2]]);
    elseif ndims(array)==4
      return  circshift(array[:,:,:,time],[domain_center[1]-peak[1],domain_center[2]-peak[2],0]);
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

mutable struct Composite_Cyclone
    surfpres
    surfu
    surfv
    pw
    precip
    qv
    tabs
    qrad
    pres   
end
    
mutable struct Composite_Cyclone_v2
    surfpres
    surfu
    surfv
    pw
    precip
    qv
    tabs
    qrad
    pres
    U
    V
    W
    LHF
    SHF    
end
   
mutable struct Composite_Cyclone_v3
    surfpres
    surfu
    surfv
    pw
    precip
    qv
    tabs
    qrad
    pres
    U
    V
    W
    LHF
    SHF
    rad_a
    dia_a
    B_a    
end
    
function cyclone_comp_timemean(composite::Composite_Cyclone)
        cyclone_comp_timemean = Composite_Cyclone(mean(composite.surfpres,dims=3),mean(composite.surfu,dims=3),mean(composite.surfv,dims=3),mean(composite.pw,dims=3),mean(composite.precip,dims=3),
        mean(composite.qv,dims=4),mean(composite.tabs,dims=4),mean(composite.qrad,dims=4),mean(composite.pres,dims=4));
end
# function smooth_time_nofalseframe(input::Array{T,3}) where {T<:Real}
#         copy = zeros(size(input))
#         trueframe = 1
#         for t in 1:size(input,3)
#             if maximum(input[:,:,t]) !== 0
#                copy[:,:,trueframe] = input[:,:,t]
#                trueframe = trueframe + 1
#             end
#         end
#         return smooth_array(copy[:,:,1:trueframe])
# end
# function smooth_time_nofalseframe(input::Array{T,4}) where {T<:Real}
#         copy = zeros(size(input))
#         trueframe = 1
#         for t in 1:size(input,4)
#             if maximum(input[:,:,:,t]) !== 0
#                copy[:,:,:,trueframe] = input[:,:,:,t]
#                trueframe = trueframe + 1
#             end
#         end
#         return smooth_array(copy[:,:,:,1:trueframe])
# end
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


function cyclonecompositer(exp_name,output_interval,final_time_days,detection_tres,domain_center)
print("\r $exp_name : Let's do this! ")
flush(stdout) 
path_to_file        = "/global/cscratch1/sd/aramreye/for_postprocessing/all200days/";
file3d              = string(path_to_file,exp_name,"_3d.nc");
file2d              = string(path_to_file,exp_name,"_2d.nc");  
dayLength           = 60*60*24/output_interval; #How many data points make one day
final_time          = floor(Int,final_time_days*dayLength);
final_time_minus_20 = floor(Int,(final_time_days - 20)*dayLength);
iterator_time_2d    = (final_time_minus_20*2):2:final_time*2;
iterator_time_3d    = final_time_minus_20:1:final_time;
print("\r $exp_name : Reading files... ")
flush(stdout)   
ds2d = Dataset(file2d);
ds3d = Dataset(file3d);
USFC = ds2d["USFC"].var[:,:,iterator_time_2d];
VSFC = ds2d["VSFC"].var[:,:,iterator_time_2d];
Precip = ds2d["Prec"].var[:,:,iterator_time_2d];
PW                = ds2d["PW"].var[:,:,iterator_time_2d];
surf_pres           = ds2d["PSFC"].var[:,:,iterator_time_2d];
QV                  = ds3d["QV"].var[:,:,:,iterator_time_3d];
TABS                = ds3d["TABS"].var[:,:,:,iterator_time_3d];
QRAD                = ds3d["QRAD"].var[:,:,:,iterator_time_3d];
t                   = ds2d["time"].var[iterator_time_2d];
PP                  = ds3d["PP"].var[:,:,:,iterator_time_3d]; #pressure perturbation
close(ds2d)
close(ds3d)
surf_pres_anomaly   = surf_pres .- mean(surf_pres,dims=[1,2]);  
composite_cyclone   = Composite_Cyclone(similar(surf_pres),similar(surf_pres),similar(surf_pres),similar(surf_pres),similar(surf_pres),
similar(QV),similar(QV),similar(QV),similar(PP));
print("\r $exp_name : Starting analysis... ")
flush(stdout) 
for time in 1:length(t)
    surf_pres_filtered = smoothfilter(-1*surf_pres_anomaly[:,:,time],detection_tres);
    peaks              = findlocalmaxima(surf_pres_filtered);
    if length(peaks) != 0
        for peak in peaks
            composite_cyclone.surfpres[:,:,time] += shifter(surf_pres,time,domain_center,peak)/length(peaks);     
            composite_cyclone.surfu[:,:,time]    += shifter(USFC,time,domain_center,peak)/length(peaks);     
            composite_cyclone.surfv[:,:,time]    += shifter(VSFC,time,domain_center,peak)/length(peaks);     
            composite_cyclone.pw[:,:,time]       += shifter(PW,time,domain_center,peak)/length(peaks);     
            composite_cyclone.precip[:,:,time]   += shifter(Precip,time,domain_center,peak)/length(peaks);     
            composite_cyclone.qv[:,:,:,time]     += shifter(QV,time,domain_center,peak)/length(peaks);     
            composite_cyclone.tabs[:,:,:,time]   += shifter(TABS,time,domain_center,peak)/length(peaks);     
            composite_cyclone.qrad[:,:,:,time]   += shifter(QRAD,time,domain_center,peak)/length(peaks);     
            composite_cyclone.pres[:,:,:,time]   += shifter(PP,time,domain_center,peak)/length(peaks);     
        end
    end
end
print("\r $exp_name : Done ")
flush(stdout) 
    return composite_cyclone
end
    
function cyclonecompositer_v2(exp_name,output_interval,final_time_days,detection_tres,domain_center)
print("\r $exp_name : Let's do this! ")
flush(stdout) 
path_to_file        = "/global/cscratch1/sd/aramreye/for_postprocessing/all200days/";
file3d              = string(path_to_file,exp_name,"_3d.nc");
file2d              = string(path_to_file,exp_name,"_2d.nc");  
dayLength           = 60*60*24/output_interval; #How many data points make one day
final_time          = floor(Int,final_time_days*dayLength);
final_time_minus_20 = floor(Int,(final_time_days - 20)*dayLength);
iterator_time_2d    = (final_time_minus_20*2):2:final_time*2;
iterator_time_3d    = final_time_minus_20:1:final_time;
print("\r $exp_name : Reading files... ")
flush(stdout)   
ds2d = Dataset(file2d);
ds3d = Dataset(file3d);
USFC = ds2d["USFC"].var[:,:,iterator_time_2d];
VSFC = ds2d["VSFC"].var[:,:,iterator_time_2d];
Precip = ds2d["Prec"].var[:,:,iterator_time_2d];
PW                = ds2d["PW"].var[:,:,iterator_time_2d];
surf_pres           = ds2d["PSFC"].var[:,:,iterator_time_2d];
QV                  = ds3d["QV"].var[:,:,:,iterator_time_3d];
TABS                = ds3d["TABS"].var[:,:,:,iterator_time_3d];
QRAD                = ds3d["QRAD"].var[:,:,:,iterator_time_3d];
t                   = ds2d["time"].var[iterator_time_2d];
PP                  = ds3d["PP"].var[:,:,:,iterator_time_3d]; #pressure perturbation
U                   = ds3d["U"].var[:,:,:,iterator_time_3d];
V                   = ds3d["V"].var[:,:,:,iterator_time_3d];
W                   = ds3d["W"].var[:,:,:,iterator_time_3d];
LHF                 = ds2d["LHF"].var[:,:,iterator_time_2d];
SHF                 = ds2d["SHF"].var[:,:,iterator_time_2d];
close(ds2d)
close(ds3d)
surf_pres_anomaly   = surf_pres .- mean(surf_pres,dims=[1,2]);  
composite_cyclone   = Composite_Cyclone_v2(zeros(size(surf_pres)),zeros(size(surf_pres)),zeros(size(surf_pres)),zeros(size(surf_pres)),zeros(size(surf_pres)),
zeros(size(QV)),zeros(size(QV)),zeros(size(QV)),zeros(size(PP)),zeros(size(U)),zeros(size(U)),zeros(size(U)),zeros(size(LHF)),zeros(size(SHF)));
print("\r $exp_name : Starting analysis... ")
flush(stdout) 
for time in 1:length(t)
    surf_pres_filtered = smoothfilter(-1*surf_pres_anomaly[:,:,time],detection_tres);
    peaks              = findlocalmaxima(surf_pres_filtered);
    if length(peaks) != 0
        for peak in peaks
            composite_cyclone.surfpres[:,:,time] += shifter(surf_pres,time,domain_center,peak)/length(peaks);     
            composite_cyclone.surfu[:,:,time]    += shifter(USFC,time,domain_center,peak)/length(peaks);     
            composite_cyclone.surfv[:,:,time]    += shifter(VSFC,time,domain_center,peak)/length(peaks);     
            composite_cyclone.pw[:,:,time]       += shifter(PW,time,domain_center,peak)/length(peaks);     
            composite_cyclone.precip[:,:,time]   += shifter(Precip,time,domain_center,peak)/length(peaks);     
            composite_cyclone.qv[:,:,:,time]     += shifter(QV,time,domain_center,peak)/length(peaks);     
            composite_cyclone.tabs[:,:,:,time]   += shifter(TABS,time,domain_center,peak)/length(peaks);     
            composite_cyclone.qrad[:,:,:,time]   += shifter(QRAD,time,domain_center,peak)/length(peaks);     
            composite_cyclone.pres[:,:,:,time]   += shifter(PP,time,domain_center,peak)/length(peaks);
            composite_cyclone.LHF[:,:,time]       += shifter(LHF,time,domain_center,peak)/length(peaks); 
            composite_cyclone.SHF[:,:,time]       += shifter(SHF,time,domain_center,peak)/length(peaks); 
            composite_cyclone.U[:,:,:,time]   += shifter(U,time,domain_center,peak)/length(peaks);     
            composite_cyclone.V[:,:,:,time]   += shifter(V,time,domain_center,peak)/length(peaks);     
            composite_cyclone.W[:,:,:,time]   += shifter(W,time,domain_center,peak)/length(peaks);
        end
    end
end
print("\r $exp_name : Done ")
flush(stdout) 
    return composite_cyclone
end
    
function cyclonecompositer_v3(exp_name,output_interval,final_time_days,detection_tres,domain_center)
print("\r $exp_name : Let's do this! ")
flush(stdout) 
path_to_file        = "/global/cscratch1/sd/aramreye/for_postprocessing/all200days/";
file3d              = string(path_to_file,exp_name,"_3d.nc");
file2d              = string(path_to_file,exp_name,"_2d.nc");  
dayLength           = 60*60*24/output_interval; #How many data points make one day
final_time          = floor(Int,final_time_days*dayLength);
final_time_minus_20 = floor(Int,(final_time_days - 60)*dayLength);
iterator_time_2d    = (final_time_minus_20*2):2:final_time*2;
iterator_time_3d    = final_time_minus_20:1:final_time;
print("\r $exp_name : Reading files... ")
flush(stdout)   
ds2d = Dataset(file2d);
ds3d = Dataset(file3d);
USFC = ds2d["USFC"].var[:,:,iterator_time_2d];
VSFC = ds2d["VSFC"].var[:,:,iterator_time_2d];
Precip = ds2d["Prec"].var[:,:,iterator_time_2d];
PW                = ds2d["PW"].var[:,:,iterator_time_2d];
surf_pres           = ds2d["PSFC"].var[:,:,iterator_time_2d];
QV                  = ds3d["QV"].var[:,:,:,iterator_time_3d];
TABS                = ds3d["TABS"].var[:,:,:,iterator_time_3d];
QRAD                = ds3d["QRAD"].var[:,:,:,iterator_time_3d];
t                   = ds2d["time"].var[iterator_time_2d];
PP                  = ds3d["PP"].var[:,:,:,iterator_time_3d]; #pressure perturbation
U                   = ds3d["U"].var[:,:,:,iterator_time_3d];
V                   = ds3d["V"].var[:,:,:,iterator_time_3d];
W                   = ds3d["W"].var[:,:,:,iterator_time_3d];
LHF                 = ds2d["LHF"].var[:,:,iterator_time_2d];
SHF                 = ds2d["SHF"].var[:,:,iterator_time_2d];
x                   = ds3d["x"].var[:]
y                   = ds3d["y"].var[:]
z                   = ds3d["z"].var[:]
t                   = ds3d["time"].var[iterator_time_3d]
P0                  = ds3d["p"].var[:]
close(ds2d)
close(ds3d)
surf_pres_anomaly   = surf_pres .- mean(surf_pres,dims=[1,2]);  
dia_a,rad_a,B_a = buoyancy_4d(exp_name,output_interval,U,V,W,QRAD,TABS,QV,PP,PW,SHF,LHF,x,y,z,t,P0)
composite_cyclone   = Composite_Cyclone_v3(zeros(size(surf_pres)),zeros(size(surf_pres)),zeros(size(surf_pres)),zeros(size(surf_pres)),zeros(size(surf_pres)),
zeros(size(QV)),zeros(size(QV)),zeros(size(QV)),zeros(size(PP)),zeros(size(U)),zeros(size(U)),zeros(size(U)),zeros(size(LHF)),zeros(size(SHF)),zeros(size(QV)),
zeros(size(QV)),zeros(size(QV)));
print("\r $exp_name : Starting analysis... ")
flush(stdout) 
for time in 1:length(t)
    surf_pres_filtered = smoothfilter(-1*surf_pres_anomaly[:,:,time],detection_tres);
    peaks              = findlocalmaxima(surf_pres_filtered);
    if length(peaks) != 0
        for peak in peaks
            composite_cyclone.surfpres[:,:,time] += shifter(surf_pres,time,domain_center,peak)/length(peaks);     
            composite_cyclone.surfu[:,:,time]    += shifter(USFC,time,domain_center,peak)/length(peaks);     
            composite_cyclone.surfv[:,:,time]    += shifter(VSFC,time,domain_center,peak)/length(peaks);     
            composite_cyclone.pw[:,:,time]       += shifter(PW,time,domain_center,peak)/length(peaks);     
            composite_cyclone.precip[:,:,time]   += shifter(Precip,time,domain_center,peak)/length(peaks);     
            composite_cyclone.qv[:,:,:,time]     += shifter(QV,time,domain_center,peak)/length(peaks);     
            composite_cyclone.tabs[:,:,:,time]   += shifter(TABS,time,domain_center,peak)/length(peaks);     
            composite_cyclone.qrad[:,:,:,time]   += shifter(QRAD,time,domain_center,peak)/length(peaks);     
            composite_cyclone.pres[:,:,:,time]   += shifter(PP,time,domain_center,peak)/length(peaks);
            composite_cyclone.LHF[:,:,time]       += shifter(LHF,time,domain_center,peak)/length(peaks); 
            composite_cyclone.SHF[:,:,time]       += shifter(SHF,time,domain_center,peak)/length(peaks); 
            composite_cyclone.U[:,:,:,time]   += shifter(U,time,domain_center,peak)/length(peaks);     
            composite_cyclone.V[:,:,:,time]   += shifter(V,time,domain_center,peak)/length(peaks);     
            composite_cyclone.W[:,:,:,time]   += shifter(W,time,domain_center,peak)/length(peaks);
            composite_cyclone.dia_a[:,:,:,time]   += shifter(dia_a,time,domain_center,peak)/length(peaks);     
            composite_cyclone.rad_a[:,:,:,time]   += shifter(rad_a,time,domain_center,peak)/length(peaks);     
            composite_cyclone.B_a[:,:,:,time]   += shifter(B_a,time,domain_center,peak)/length(peaks);
        end
    end
end
print("\r $exp_name : Done ")
flush(stdout) 
    return composite_cyclone
end

    function cyclonecompositer_v3_partial_60d(exp_name,output_interval,initial_time_days,final_time_days,detection_tres,domain_center,position)
print("\r $exp_name : Let's do this! ")
flush(stdout) 
path_to_file        = "/global/cscratch1/sd/aramreye/for_postprocessing/all200days/";
file3d              = string(path_to_file,exp_name,"_3d.nc");
file2d              = string(path_to_file,exp_name,"_2d.nc");  
dayLength           = 60*60*24/output_interval; #How many data points make one day
final_time          = floor(Int,final_time_days*dayLength);
initial_time          = floor(Int,initial_time_days*dayLength);
final_time_minus_20 = floor(Int,(final_time_days - 60)*dayLength);
iterator_time_2d    = (final_time_minus_20*2):2:final_time*2;
iterator_time_3d    = final_time_minus_20:1:final_time;
if output_interval == 7200
        smooth_x    = 11
        smooth_y    = 11
end
    smooth_time = floor(Int,dayLength*5)

if position == 1
        iterator_time_2d    = initial_time*2-1:2:(final_time+div(smooth_time-1,2))*2  
        iterator_time_3d    = initial_time:1:final_time+div(smooth_time-1,2)
    elseif position == 2
        iterator_time_2d    = (initial_time-div(smooth_time-1,2))*2-1:2:(final_time+div(smooth_time-1,2))*2  
        iterator_time_3d    = initial_time-div(smooth_time-1,2):1:final_time+div(smooth_time-1,2)
    elseif position == 3
        iterator_time_2d    = (initial_time-div(smooth_time-1,2))*2-1:2:final_time*2  
        iterator_time_3d    = initial_time-div(smooth_time-1,2):1:final_time
    elseif position == 4

    end
println(iterator_time_2d)         
println(iterator_time_3d)
print("\r $exp_name : Reading files... ")
flush(stdout)   
ds2d = Dataset(file2d);
ds3d = Dataset(file3d);
USFC = ds2d["USFC"].var[:,:,iterator_time_2d];
VSFC = ds2d["VSFC"].var[:,:,iterator_time_2d];
Precip = ds2d["Prec"].var[:,:,iterator_time_2d];
PW                = ds2d["PW"].var[:,:,iterator_time_2d];
surf_pres           = ds2d["PSFC"].var[:,:,iterator_time_2d];
QV                  = ds3d["QV"].var[:,:,:,iterator_time_3d];
TABS                = ds3d["TABS"].var[:,:,:,iterator_time_3d];
QRAD                = ds3d["QRAD"].var[:,:,:,iterator_time_3d];
t                   = ds2d["time"].var[iterator_time_2d];
PP                  = ds3d["PP"].var[:,:,:,iterator_time_3d]; #pressure perturbation
U                   = ds3d["U"].var[:,:,:,iterator_time_3d];
V                   = ds3d["V"].var[:,:,:,iterator_time_3d];
W                   = ds3d["W"].var[:,:,:,iterator_time_3d];
LHF                 = ds2d["LHF"].var[:,:,iterator_time_2d];
SHF                 = ds2d["SHF"].var[:,:,iterator_time_2d];
x                   = ds3d["x"].var[:]
y                   = ds3d["y"].var[:]
z                   = ds3d["z"].var[:]
t                   = ds3d["time"].var[iterator_time_3d]
P0                  = ds3d["p"].var[:]
close(ds2d)
close(ds3d)

surf_pres_anomaly   = surf_pres .- mean(surf_pres,dims=[1,2]);  
dia_a,rad_a,B_a = buoyancy_4d(exp_name,output_interval,U,V,W,QRAD,TABS,QV,PP,PW,SHF,LHF,x,y,z,t,P0)

USFC                = cutborders!(USFC               ,smooth_time,position)                
VSFC                = cutborders!(VSFC               ,smooth_time,position)
Precip              = cutborders!(Precip             ,smooth_time,position)
PW                  = cutborders!(PW                 ,smooth_time,position)
surf_pres           = cutborders!(surf_pres          ,smooth_time,position)
QV                  = cutborders!(QV                 ,smooth_time,position)
TABS                = cutborders!(TABS               ,smooth_time,position)
QRAD                = cutborders!(QRAD               ,smooth_time,position)
t                   = cutborders!(t                  ,smooth_time,position)
PP                  = cutborders!(PP                 ,smooth_time,position)
U                   = cutborders!(U                  ,smooth_time,position)
V                   = cutborders!(V                  ,smooth_time,position)
W                   = cutborders!(W                  ,smooth_time,position)
LHF                 = cutborders!(LHF                ,smooth_time,position)
SHF                 = cutborders!(SHF                ,smooth_time,position)
surf_pres_anomaly   = cutborders!(surf_pres_anomaly  ,smooth_time,position)
dia_a               = cutborders!(dia_a              ,smooth_time,position)
rad_a               = cutborders!(rad_a              ,smooth_time,position)
B_a                 = cutborders!(B_a                ,smooth_time,position)


composite_cyclone   = Composite_Cyclone_v3(zeros(size(surf_pres)),zeros(size(surf_pres)),zeros(size(surf_pres)),zeros(size(surf_pres)),zeros(size(surf_pres)),
zeros(size(QV)),zeros(size(QV)),zeros(size(QV)),zeros(size(PP)),zeros(size(U)),zeros(size(U)),zeros(size(U)),zeros(size(LHF)),zeros(size(SHF)),zeros(size(QV)),
zeros(size(QV)),zeros(size(QV)));
print("\r $exp_name : Starting analysis... ")
flush(stdout) 
for time in 1:length(t)
    surf_pres_filtered = smoothfilter(-1*surf_pres_anomaly[:,:,time],detection_tres);
    peaks              = findlocalmaxima(surf_pres_filtered);
    if length(peaks) != 0
        for peak in peaks
            composite_cyclone.surfpres[:,:,time] += shifter(surf_pres,time,domain_center,peak)/length(peaks);     
            composite_cyclone.surfu[:,:,time]    += shifter(USFC,time,domain_center,peak)/length(peaks);     
            composite_cyclone.surfv[:,:,time]    += shifter(VSFC,time,domain_center,peak)/length(peaks);     
            composite_cyclone.pw[:,:,time]       += shifter(PW,time,domain_center,peak)/length(peaks);     
            composite_cyclone.precip[:,:,time]   += shifter(Precip,time,domain_center,peak)/length(peaks);     
            composite_cyclone.qv[:,:,:,time]     += shifter(QV,time,domain_center,peak)/length(peaks);     
            composite_cyclone.tabs[:,:,:,time]   += shifter(TABS,time,domain_center,peak)/length(peaks);     
            composite_cyclone.qrad[:,:,:,time]   += shifter(QRAD,time,domain_center,peak)/length(peaks);     
            composite_cyclone.pres[:,:,:,time]   += shifter(PP,time,domain_center,peak)/length(peaks);
            composite_cyclone.LHF[:,:,time]       += shifter(LHF,time,domain_center,peak)/length(peaks); 
            composite_cyclone.SHF[:,:,time]       += shifter(SHF,time,domain_center,peak)/length(peaks); 
            composite_cyclone.U[:,:,:,time]   += shifter(U,time,domain_center,peak)/length(peaks);     
            composite_cyclone.V[:,:,:,time]   += shifter(V,time,domain_center,peak)/length(peaks);     
            composite_cyclone.W[:,:,:,time]   += shifter(W,time,domain_center,peak)/length(peaks);
            composite_cyclone.dia_a[:,:,:,time]   += shifter(dia_a,time,domain_center,peak)/length(peaks);     
            composite_cyclone.rad_a[:,:,:,time]   += shifter(rad_a,time,domain_center,peak)/length(peaks);     
            composite_cyclone.B_a[:,:,:,time]   += shifter(B_a,time,domain_center,peak)/length(peaks);
        end
    end
end
print("\r $exp_name : Done ")
flush(stdout) 
    return composite_cyclone
end


    function cyclonecompositer_v4_partial(pathto,exp_name,output_interval,initial_time_days,final_time_days,detection_tres,domain_center,position)
print("\r $exp_name : Let's do this! ")
flush(stdout) 
        path_to_file        = "/global/cscratch1/sd/aramreye/for_postprocessing/all200days/";
path3d = string(pathto,"OUT_3D/")
path2d = string(pathto,"OUT_2D/")
file3d = filter(x -> occursin(r".*\.nc", x), string.(path3d, readdir(path3d)))
file2d = filter(x -> occursin(r".*\.nc", x), string.(path2d, readdir(path2d)))
#file3d              = string(path_to_file,exp_name,"_3d.nc");
#file2d              = string(path_to_file,exp_name,"_2d.nc");  
dayLength           = 60*60*24/output_interval; #How many data points make one day
final_time          = floor(Int,final_time_days*dayLength);
initial_time        = floor(Int,initial_time_days*dayLength);
final_time_minus_20 = floor(Int,(final_time_days - 60)*dayLength);
iterator_time_2d    = (final_time_minus_20*2):2:final_time*2;
iterator_time_3d    = final_time_minus_20:1:final_time;
if output_interval == 7200
        smooth_x    = 11
        smooth_y    = 11
end
    smooth_time = floor(Int,dayLength*5)

if position == 1
        iterator_time_2d    = initial_time*2-1:2:(final_time+div(smooth_time-1,2))*2  
        iterator_time_3d    = initial_time:1:final_time+div(smooth_time-1,2)
    elseif position == 2
        iterator_time_2d    = (initial_time-div(smooth_time-1,2))*2-1:2:(final_time+div(smooth_time-1,2))*2  
        iterator_time_3d    = initial_time-div(smooth_time-1,2):1:final_time+div(smooth_time-1,2)
    elseif position == 3
        iterator_time_2d    = (initial_time-div(smooth_time-1,2))*2-1:2:final_time*2  
        iterator_time_3d    = initial_time-div(smooth_time-1,2):1:final_time
    elseif position == 4

    end
println(iterator_time_2d)         
println(iterator_time_3d)
print("\r $exp_name : Reading files... ")
flush(stdout)   
ds2d = Dataset(file2d,deferopen=true,aggdim="time");
ds3d = Dataset(file3d,deferopen=true,aggdim="time");
USFC = ds2d["USFC"].var[:,:,iterator_time_2d];
VSFC = ds2d["VSFC"].var[:,:,iterator_time_2d];
Precip = ds2d["Prec"].var[:,:,iterator_time_2d];
PW                = ds2d["PW"].var[:,:,iterator_time_2d];
surf_pres           = ds2d["PSFC"].var[:,:,iterator_time_2d];
QV                  = ds3d["QV"].var[:,:,:,iterator_time_3d];
TABS                = ds3d["TABS"].var[:,:,:,iterator_time_3d];
QRAD                = ds3d["QRAD"].var[:,:,:,iterator_time_3d];
t                   = ds2d["time"].var[iterator_time_2d];
PP                  = ds3d["PP"].var[:,:,:,iterator_time_3d]; #pressure perturbation
U                   = ds3d["U"].var[:,:,:,iterator_time_3d];
V                   = ds3d["V"].var[:,:,:,iterator_time_3d];
W                   = ds3d["W"].var[:,:,:,iterator_time_3d];
LHF                 = ds2d["LHF"].var[:,:,iterator_time_2d];
SHF                 = ds2d["SHF"].var[:,:,iterator_time_2d];
x                   = ds3d["x"].var[:]
y                   = ds3d["y"].var[:]
z                   = ds3d["z"].var[:]
t                   = ds3d["time"].var[iterator_time_3d]
P0                  = ds3d["p"].var[:]
close(ds2d)
close(ds3d)

surf_pres_anomaly   = surf_pres .- mean(surf_pres,dims=[1,2]);  
dia_a,rad_a,B_a = buoyancy_4d(exp_name,output_interval,U,V,W,QRAD,TABS,QV,PP,PW,SHF,LHF,x,y,z,t,P0)

USFC                = cutborders!(USFC               ,smooth_time,position)                
VSFC                = cutborders!(VSFC               ,smooth_time,position)
Precip              = cutborders!(Precip             ,smooth_time,position)
PW                  = cutborders!(PW                 ,smooth_time,position)
surf_pres           = cutborders!(surf_pres          ,smooth_time,position)
QV                  = cutborders!(QV                 ,smooth_time,position)
TABS                = cutborders!(TABS               ,smooth_time,position)
QRAD                = cutborders!(QRAD               ,smooth_time,position)
t                   = cutborders!(t                  ,smooth_time,position)
PP                  = cutborders!(PP                 ,smooth_time,position)
U                   = cutborders!(U                  ,smooth_time,position)
V                   = cutborders!(V                  ,smooth_time,position)
W                   = cutborders!(W                  ,smooth_time,position)
LHF                 = cutborders!(LHF                ,smooth_time,position)
SHF                 = cutborders!(SHF                ,smooth_time,position)
surf_pres_anomaly   = cutborders!(surf_pres_anomaly  ,smooth_time,position)
dia_a               = cutborders!(dia_a              ,smooth_time,position)
rad_a               = cutborders!(rad_a              ,smooth_time,position)
B_a                 = cutborders!(B_a                ,smooth_time,position)


composite_cyclone   = Composite_Cyclone_v3(zeros(size(surf_pres)),zeros(size(surf_pres)),zeros(size(surf_pres)),zeros(size(surf_pres)),zeros(size(surf_pres)),
zeros(size(QV)),zeros(size(QV)),zeros(size(QV)),zeros(size(PP)),zeros(size(U)),zeros(size(U)),zeros(size(U)),zeros(size(LHF)),zeros(size(SHF)),zeros(size(QV)),
zeros(size(QV)),zeros(size(QV)));
print("\r $exp_name : Starting analysis... ")
flush(stdout) 
for time in 1:length(t)
    surf_pres_filtered = smoothfilter(-1*surf_pres_anomaly[:,:,time],detection_tres);
    peaks              = findlocalmaxima(surf_pres_filtered);
    if length(peaks) != 0
        for peak in peaks
            composite_cyclone.surfpres[:,:,time] .+= shifter(surf_pres,time,domain_center,peak)./length(peaks);     
            composite_cyclone.surfu[:,:,time]    .+= shifter(USFC,time,domain_center,peak)./length(peaks);     
            composite_cyclone.surfv[:,:,time]    .+= shifter(VSFC,time,domain_center,peak)./length(peaks);     
            composite_cyclone.pw[:,:,time]       .+= shifter(PW,time,domain_center,peak)./length(peaks);     
            composite_cyclone.precip[:,:,time]   .+= shifter(Precip,time,domain_center,peak)./length(peaks);     
            composite_cyclone.qv[:,:,:,time]     .+= shifter(QV,time,domain_center,peak)./length(peaks);     
            composite_cyclone.tabs[:,:,:,time]   .+= shifter(TABS,time,domain_center,peak)./length(peaks);     
            composite_cyclone.qrad[:,:,:,time]   .+= shifter(QRAD,time,domain_center,peak)./length(peaks);     
            composite_cyclone.pres[:,:,:,time]   .+= shifter(PP,time,domain_center,peak)./length(peaks);
            composite_cyclone.LHF[:,:,time]      .+= shifter(LHF,time,domain_center,peak)./length(peaks); 
            composite_cyclone.SHF[:,:,time]      .+= shifter(SHF,time,domain_center,peak)./length(peaks); 
            composite_cyclone.U[:,:,:,time]      .+= shifter(U,time,domain_center,peak)./length(peaks);     
            composite_cyclone.V[:,:,:,time]      .+= shifter(V,time,domain_center,peak)./length(peaks);     
            composite_cyclone.W[:,:,:,time]      .+= shifter(W,time,domain_center,peak)./length(peaks);
            composite_cyclone.dia_a[:,:,:,time]  .+= shifter(dia_a,time,domain_center,peak)./length(peaks);     
            composite_cyclone.rad_a[:,:,:,time]  .+= shifter(rad_a,time,domain_center,peak)./length(peaks);     
            composite_cyclone.B_a[:,:,:,time]    .+= shifter(B_a,time,domain_center,peak)./length(peaks);
        end
    end
end
print("\r $exp_name : Done ")
flush(stdout) 
    return composite_cyclone
end

function cyclonecompositer_v4_partial2(pathto,exp_name,output_interval,initial_time_days,final_time_days,detection_tres,domain_center,position)
    print("\r $exp_name : Let's do this! ")
    flush(stdout) 
            path_to_file        = "/global/cscratch1/sd/aramreye/for_postprocessing/all200days/";
    path3d = string(pathto,"OUT_3D/")
    path2d = string(pathto,"OUT_2D/")
    file3d = filter(x -> occursin(r".*\.nc", x), string.(path3d, readdir(path3d)))
    file2d = filter(x -> occursin(r".*\.nc", x), string.(path2d, readdir(path2d)))
    #file3d              = string(path_to_file,exp_name,"_3d.nc");
    #file2d              = string(path_to_file,exp_name,"_2d.nc");  
    dayLength           = 60*60*24/output_interval; #How many data points make one day
    final_time          = floor(Int,final_time_days*dayLength);
    initial_time        = floor(Int,initial_time_days*dayLength);
    final_time_minus_20 = floor(Int,(final_time_days - 60)*dayLength);
    iterator_time_2d    = (final_time_minus_20*2):2:final_time*2;
    iterator_time_3d    = final_time_minus_20:1:final_time;
    if output_interval == 7200
            smooth_x    = 11
            smooth_y    = 11
    end
        smooth_time = floor(Int,dayLength*5)
    
    if position == 1
            iterator_time_2d    = initial_time*2-1:2:(final_time+div(smooth_time-1,2))*2  
            iterator_time_3d    = initial_time:1:final_time+div(smooth_time-1,2)
        elseif position == 2
            iterator_time_2d    = (initial_time-div(smooth_time-1,2))*2-1:2:(final_time+div(smooth_time-1,2))*2  
            iterator_time_3d    = initial_time-div(smooth_time-1,2):1:final_time+div(smooth_time-1,2)
        elseif position == 3
            iterator_time_2d    = (initial_time-div(smooth_time-1,2))*2-1:2:final_time*2  
            iterator_time_3d    = initial_time-div(smooth_time-1,2):1:final_time
        elseif position == 4
    
        end
    println(iterator_time_2d)         
    println(iterator_time_3d)
    print("\r $exp_name : Reading files... ")
    flush(stdout)   
    ds2d = Dataset(file2d,deferopen=true,aggdim="time");
    ds3d = Dataset(file3d,deferopen=true,aggdim="time");
    USFC = ds2d["USFC"].var[:,:,iterator_time_2d];
    VSFC = ds2d["VSFC"].var[:,:,iterator_time_2d];
    Precip = ds2d["Prec"].var[:,:,iterator_time_2d];
    PW                = ds2d["PW"].var[:,:,iterator_time_2d];
    surf_pres           = ds2d["PSFC"].var[:,:,iterator_time_2d];
    QV                  = ds3d["QV"].var[:,:,:,iterator_time_3d];
    TABS                = ds3d["TABS"].var[:,:,:,iterator_time_3d];
    QRAD                = ds3d["QRAD"].var[:,:,:,iterator_time_3d];
    t                   = ds2d["time"].var[iterator_time_2d];
    PP                  = ds3d["PP"].var[:,:,:,iterator_time_3d]; #pressure perturbation
    U                   = ds3d["U"].var[:,:,:,iterator_time_3d];
    V                   = ds3d["V"].var[:,:,:,iterator_time_3d];
    W                   = ds3d["W"].var[:,:,:,iterator_time_3d];
    LHF                 = ds2d["LHF"].var[:,:,iterator_time_2d];
    SHF                 = ds2d["SHF"].var[:,:,iterator_time_2d];
    x                   = ds3d["x"].var[:]
    y                   = ds3d["y"].var[:]
    z                   = ds3d["z"].var[:]
    t                   = ds3d["time"].var[iterator_time_3d]
    P0                  = ds3d["p"].var[:]
    close(ds2d)
    close(ds3d)
    
    surf_pres_anomaly   = surf_pres .- mean(surf_pres,dims=[1,2]);  
    dia_a,rad_a,B_a = buoyancy_4d(exp_name,output_interval,U,V,W,QRAD,TABS,QV,PP,PW,SHF,LHF,x,y,z,t,P0)
    
    USFC                = cutborders!(USFC               ,smooth_time,position)                
    VSFC                = cutborders!(VSFC               ,smooth_time,position)
    Precip              = cutborders!(Precip             ,smooth_time,position)
    PW                  = cutborders!(PW                 ,smooth_time,position)
    surf_pres           = cutborders!(surf_pres          ,smooth_time,position)
    QV                  = cutborders!(QV                 ,smooth_time,position)
    TABS                = cutborders!(TABS               ,smooth_time,position)
    QRAD                = cutborders!(QRAD               ,smooth_time,position)
    t                   = cutborders!(t                  ,smooth_time,position)
    PP                  = cutborders!(PP                 ,smooth_time,position)
    U                   = cutborders!(U                  ,smooth_time,position)
    V                   = cutborders!(V                  ,smooth_time,position)
    W                   = cutborders!(W                  ,smooth_time,position)
    LHF                 = cutborders!(LHF                ,smooth_time,position)
    SHF                 = cutborders!(SHF                ,smooth_time,position)
    surf_pres_anomaly   = cutborders!(surf_pres_anomaly  ,smooth_time,position)
    dia_a               = cutborders!(dia_a              ,smooth_time,position)
    rad_a               = cutborders!(rad_a              ,smooth_time,position)
    B_a                 = cutborders!(B_a                ,smooth_time,position)
    
    
    composite_cyclone   = Composite_Cyclone_v3(zeros(size(surf_pres)),zeros(size(surf_pres)),zeros(size(surf_pres)),zeros(size(surf_pres)),zeros(size(surf_pres)),
    zeros(size(QV)),zeros(size(QV)),zeros(size(QV)),zeros(size(PP)),zeros(size(U)),zeros(size(U)),zeros(size(U)),zeros(size(LHF)),zeros(size(SHF)),zeros(size(QV)),
    zeros(size(QV)),zeros(size(QV)));
    print("\r $exp_name : Starting analysis... ")
    flush(stdout)
    buf3d = Array{typeof(B[1])}(undef, length(x),length(y), length(z))
    buf2d = Array{typeof(B[1])}(undef, length(x),length(y))
    for time in 1:length(t)
        surf_pres_filtered = smoothfilter(-1*surf_pres_anomaly[:,:,time],detection_tres);
        peaks              = findlocalmaxima(surf_pres_filtered);
        if length(peaks) != 0
            for peak in peaks
                composite_cyclone.surfpres[:,:,time] .+= shifter!(buf2d,surf_pres[:,:,time],domain_center,peak)./length(peaks);     
                composite_cyclone.surfu[:,:,time]    .+= shifter!(buf2d,USFC[:,:,time],domain_center,peak)./length(peaks);     
                composite_cyclone.surfv[:,:,time]    .+= shifter!(buf2d,VSFC[:,:,time],domain_center,peak)./length(peaks);     
                composite_cyclone.pw[:,:,time]       .+= shifter!(buf2d,PW[:,:,time],domain_center,peak)./length(peaks);     
                composite_cyclone.precip[:,:,time]   .+= shifter!(buf2d,Precip[:,:,time],domain_center,peak)./length(peaks);     
                composite_cyclone.qv[:,:,:,time]     .+= shifter!(buf3d,QV[:,:,:,time],domain_center,peak)./length(peaks);     
                composite_cyclone.tabs[:,:,:,time]   .+= shifter!(buf3d,TABS[:,:,:,time],domain_center,peak)./length(peaks);     
                composite_cyclone.qrad[:,:,:,time]   .+= shifter!(buf3d,QRAD[:,:,:,time],domain_center,peak)./length(peaks);     
                composite_cyclone.pres[:,:,:,time]   .+= shifter!(buf3d,PP[:,:,:,time],domain_center,peak)./length(peaks);
                composite_cyclone.LHF[:,:,time]      .+= shifter!(buf2d,LHF[:,:,time],domain_center,peak)./length(peaks); 
                composite_cyclone.SHF[:,:,time]      .+= shifter!(buf2d,SHF[:,:,time],domain_center,peak)./length(peaks); 
                composite_cyclone.U[:,:,:,time]      .+= shifter!(buf3d,U[:,:,:,time],domain_center,peak)./length(peaks);     
                composite_cyclone.V[:,:,:,time]      .+= shifter!(buf3d,V[:,:,:,time],domain_center,peak)./length(peaks);     
                composite_cyclone.W[:,:,:,time]      .+= shifter!(buf3d,W[:,:,:,time],domain_center,peak)./length(peaks);
                composite_cyclone.dia_a[:,:,:,time]  .+= shifter!(buf3d,dia_a[:,:,:,time],domain_center,peak)./length(peaks);     
                composite_cyclone.rad_a[:,:,:,time]  .+= shifter!(buf3d,rad_a[:,:,:,time],domain_center,peak)./length(peaks);     
                composite_cyclone.B_a[:,:,:,time]    .+= shifter!(buf3d,B_a[:,:,:,time],domain_center,peak)./length(peaks);
            end
        end
    end
    print("\r $exp_name : Done ")
    flush(stdout) 
        return composite_cyclone
    end
    
