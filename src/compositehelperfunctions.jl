import Images: findlocalmaxima
function smoothfilter(surf_pres_anomaly,treshold=9)
    surf_pres_median = mapwindow(median!,surf_pres_anomaly,[3,3]);
    surf_pres_median = surf_pres_median.*(surf_pres_median.>treshold);
    surf_pres_filtered = imfilter(surf_pres_median,Kernel.gaussian(3));
    surf_pres_filtered = surf_pres_filtered.*(surf_pres_filtered.>treshold);
end
function shifter(array,time,domain_center,peak)
    if ndims(array)==3
      return  circshift(array[:,:,time],[domain_center[1]-peak[1],domain_center[2]-peak[2]]);
    elseif ndims(array)==4
      return  circshift(array[:,:,:,time],[domain_center[1]-peak[1],domain_center[2]-peak[2],0]);
    end
end

function shifter!(dest,array,time,domain_center,peak)
    if ndims(array)==3
      return  circshift(dest,array[:,:,time],[domain_center[1]-peak[1],domain_center[2]-peak[2]]);
    elseif ndims(array)==4
      return  circshift(dest,array[:,:,:,time],[domain_center[1]-peak[1],domain_center[2]-peak[2],0]);
    end
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
    
function timemean_nofalseframe(input::Array{T,4}) where {T<:Real}
        copy = zeros(size(input))
        trueframes = 1
        for t in 1:size(input,4)
            if maximum(input[:,:,:,t]) != 0
               copy[:,:,:,trueframes] = input[:,:,:,t]
               trueframes = trueframes + 1
            end
        end
        return mean(copy[:,:,:,1:trueframes-1],dims=4)
end
    
function timemean_nofalseframe(input::Array{T,3}) where {T<:Real}
        copy = zeros(size(input))
        trueframes = 1
        for t in 1:size(input,3)
            if maximum(input[:,:,t]) != 0
               copy[:,:,trueframes] = input[:,:,t]
               trueframes = trueframes + 1
            end
        end
        return mean(copy[:,:,1:trueframes-1],dims=3)
end
    
function removefalseframes(input::Array{T,4}) where {T<:Real}
        copy = zeros(size(input))
        trueframes = 1
        for t in 1:size(input,4)
            if maximum(input[:,:,:,t]) != 0
               copy[:,:,:,trueframes] = input[:,:,:,t]
               trueframes = trueframes + 1
            end
        end
        return copy[:,:,:,1:trueframes-1]
end
 
   
function removefalseframes(input::Array{T,3}) where {T<:Real}
        copy = zeros(size(input))
        trueframes = 1
        for t in 1:size(input,3)
            if maximum(input[:,:,t]) != 0
               copy[:,:,trueframes] = input[:,:,t]
               trueframes = trueframes + 1
            end
        end
        return copy[:,:,1:trueframes-1]
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
    