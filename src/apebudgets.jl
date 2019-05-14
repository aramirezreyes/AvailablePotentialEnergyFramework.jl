### Two methods to compute APE budgets and Buoyancy budgets from SAM v3 outputs, modified from Da Yang, and based on the mathematics from
# Yang, Da. “Boundary Layer Diabatic Processes, the Virtual Effect, and Convective Self-Aggregation.” Journal of Advances in Modeling Earth Systems 10, no. 9 (September 2018): 2163–76. https://doi.org/10.1029/2017MS001261.

"""
Returns an apebudget object.
Input is 
B (buoyancy)
U,V,W The three dimensional velocities U,V and W
N2 The Brunt Va"isala frequency squared
Fs the Surface fluxes
Diabatic_other other sources of diabatic heating
rho0 the mean density
x,y,z,t the coordinate vectors
dx,dy,dz,dt the steps in each coordinate
z_up the maximum height that will be used

"""

function getapebudget_old(B, U,V, W, N2, RAD_b, Fs, Diabatic_other, rho0, x,y, z, t, dx,dy, dz, dt, z_up)

#***********Empty array generation***********#
Udb2dx = similar(U)
Vdb2dy = similar(V)
#*************  APE **************
b2      = B.*B/2
APE_b2  = mean(b2,dims=(1,2))[1,1,:,:]./N2; 
# # # figure(4)
# # # contourf(APE_b2)

#************ KE ********************
#KE      = U.*U/2 + V.*V/2
xBar_KE = mean(U.*U/2 + V.*V/2,dims=(1,2))[1,1,:,:]

#************ APE rate ***************
xBar_APE_rate = Array{typeof(B[1])}(undef,length(z), length(t))
xBar_APE_rate[:,1:end-1] = (APE_b2[:,2:end] - APE_b2[:,1:end-1])/dt; 
xBar_APE_rate[:,end] = xBar_APE_rate[:,end-1]

#*************  UdxB2 **************
b2_ghost= Array{typeof(B[1])}(undef, length(x)+1,length(y)+1, length(z), length(t))
b2_ghost[1:end-1,1:end-1,:,:] = b2

b2_ghost[end,1:end-1,:,:] = b2[1,:,:,:]
b2_ghost[1:end-1,end,:,:] = b2[:,1,:,:]

@. Udb2dx = U*(b2_ghost[2:end,1:end-1,:,:]-b2_ghost[1:end-1,1:end-1,:,:])/dx
@. Vdb2dy = V*(b2_ghost[1:end-1,2:end,:,:]-b2_ghost[1:end-1,1:end-1,:,:])/dy
xBar_APE_Ub2 = mean(Udb2dx,dims=(1,2))[1,1,:,:]./N2
xBar_APE_Vb2 = mean(Vdb2dy,dims=(1,2))[1,1,:,:]./N2

# static stability WN2
#APE_WN2     = W.*B
xBar_APE_WN2= mean(W.*B,dims=(1,2))[1,1,:,:]

# RAD generation
xBar_APE_RAD     = mean(RAD_b.*B,dims=(1,2))[1,1,:,:]./N2;  


# Diabatic_other
xBar_APE_DIA     = mean(Diabatic_other.*B,dims=(1,2))[1,1,:,:]./N2;  


# interpolation 
k_up              = argmin(abs.(z.-z_up));
z1                = z[1]:dz:z[k_up];
rho01 = zeros(length(z1),length(t))
rho02 = zeros(length(z1),length(t))
xBar_APE_b21 = similar(rho01)
xBar_APE_RAD1 = similar(rho01)
xBar_APE_DIA1 = similar(rho01)
xBar_APE_WN21 = similar(rho01)
xBar_APE_Ub21 = similar(rho01)
xBar_APE_Vb21 = similar(rho01)
xBar_KE1 = similar(rho01)
xBar_APE_rate1 = similar(rho01)
for time in 1:length(t)
    rho01_itp         = LinearInterpolation(z, rho0[:,time]);
    xBar_APE_b21_itp  = LinearInterpolation(z, APE_b2[:,time]);
    xBar_APE_RAD1_itp = LinearInterpolation(z, xBar_APE_RAD[:,time]);
    xBar_APE_DIA1_itp = LinearInterpolation(z, xBar_APE_DIA[:,time]);
    xBar_APE_WN21_itp = LinearInterpolation(z, xBar_APE_WN2[:,time]);
    xBar_APE_Ub21_itp = LinearInterpolation(z, xBar_APE_Ub2[:,time]);
    xBar_APE_Vb21_itp = LinearInterpolation(z, xBar_APE_Vb2[:,time]);
    xBar_KE1_itp      = LinearInterpolation(z, xBar_KE[:,time]);
    xBar_APE_rate1_itp    = LinearInterpolation(z, xBar_APE_rate[:,time]);
                        
                        
    rho01[:,time]              = [rho01_itp(x) for x in z1]
    xBar_APE_b21[:,time]       = [xBar_APE_b21_itp(x) for x in z1]
    xBar_APE_RAD1[:,time]      = [xBar_APE_RAD1_itp(x) for x in z1]
    xBar_APE_DIA1[:,time]      = [xBar_APE_DIA1_itp(x) for x in z1]
    xBar_APE_WN21[:,time]      = [xBar_APE_WN21_itp(x) for x in z1]
    xBar_APE_Ub21[:,time]      = [xBar_APE_Ub21_itp(x) for x in z1]
    xBar_APE_Vb21[:,time]      = [xBar_APE_Vb21_itp(x) for x in z1]
    xBar_KE1[:,time]           = [xBar_KE1_itp(x) for x in z1]
    xBar_APE_rate1[:,time]     = [xBar_APE_rate1_itp(x) for x in z1]
end



# vertical integration
int_mass     = sum(rho01.*dz,dims=1)[1,:]
int_KE       = sum(rho01.*xBar_KE1.*dz,dims=1)[1,:]
int_APE      = sum(rho01.*xBar_APE_b21.*dz,dims=1)[1,:]
int_APE_RAD  = sum(rho01.*xBar_APE_RAD1.*dz,dims=1)[1,:]
int_APE_DIA  = sum(rho01.*xBar_APE_DIA1.*dz,dims=1)[1,:]
int_APE_WN2  = sum(rho01.*xBar_APE_WN21.*dz,dims=1)[1,:]
int_APE_Ub2  = sum(rho01.*xBar_APE_Ub21.*dz,dims=1)[1,:]
int_APE_Vb2  = sum(rho01.*xBar_APE_Vb21.*dz,dims=1)[1,:]
int_APE_rate = sum(rho01.*xBar_APE_rate1.*dz,dims=1)[1,:]

# surface flux contribution
#N2S         = N2[1,:]
#APE_Fs      = B[:,:,1,:].*Fs

xBar_APE_Fs = mean(B[:,:,1,:].*Fs, dims=(1,2))[1,1,1,:]./N2[1,:];

residual    = int_APE_rate .+ int_APE_Ub2 .+ int_APE_Vb2+int_APE_WN2 .- (int_APE_RAD .+ int_APE_DIA .+ xBar_APE_Fs)

return (int_mass, int_KE, int_APE, int_APE_rate, int_APE_Ub2,int_APE_Vb2, int_APE_WN2, int_APE_RAD, int_APE_DIA, xBar_APE_Fs, residual)
end

"""
-------------Computes the APE budgets budget---------------

"""

function getapebudget(B, U,V, W, N2, RAD_b, Fs, Diabatic_other, rho0, x,y, z, t, dx,dy, dz, dt, z_up)
    N2 = reshape(N2,1,1,length(z),length(t)) 
    #***********Empty array generation***********#
    Udb2dx = similar(U)
    Vdb2dy = similar(V)
    buf    = similar(U)
    #xBar_Fs = zeros(length(z), length(t))
    xBar_KE = Array{typeof(B[1])}(undef,1,1,length(z), length(t))
    APE_b2 = Array{typeof(B[1])}(undef,1,1,length(z), length(t))
    xBar_APE_rate = Array{typeof(B[1])}(undef,length(z), length(t))
    b2_ghost= Array{typeof(B[1])}(undef, length(x)+1,length(y)+1, length(z), length(t))
    xBar_APE_Ub2 = Array{typeof(B[1])}(undef,1,1,length(z), length(t))
    xBar_APE_Vb2 = Array{typeof(B[1])}(undef,1,1,length(z), length(t))
    xBar_APE_WN2 = Array{typeof(B[1])}(undef,1,1,length(z), length(t))
    xBar_APE_RAD = Array{typeof(B[1])}(undef,1,1,length(z), length(t))
    xBar_APE_DIA = Array{typeof(B[1])}(undef,1,1,length(z), length(t))
    xBar_APE_FS = Array{typeof(B[1])}(undef,1,1,length(z), length(t))
    #*************  APE **************
    @. buf      = B*B/2
    mean!(APE_b2,buf);
    @. APE_b2 = APE_b2/N2
    @. b2_ghost[1:end-1,1:end-1,:,:] = buf
    @. b2_ghost[end,1:end-1,:,:] = buf[1,:,:,:]
    @. b2_ghost[1:end-1,end,:,:] = buf[:,1,:,:]
 
    
    #************ KE ********************
    #KE      = U.*U/2 + V.*V/2
    @. buf = U*U/2 + V*V/2
    xBar_KE = mean!(xBar_KE,buf)

    #************ APE rate ***************
    
    @. xBar_APE_rate[:,1:end-1] = (APE_b2[1,1,:,2:end] - APE_b2[1,1,:,1:end-1])/dt; 
    @. xBar_APE_rate[:,end] = xBar_APE_rate[:,end-1]
    
    #*************  UdxB2 **************
   
    @. buf =    Udb2dx = U*(b2_ghost[2:end,1:end-1,:,:]-b2_ghost[1:end-1,1:end-1,:,:])/dx
    mean!(xBar_APE_Ub2,buf)

    @. buf = V*(b2_ghost[1:end-1,2:end,:,:]-b2_ghost[1:end-1,1:end-1,:,:])/dy
    
    mean!(xBar_APE_Vb2,buf)

    @. xBar_APE_Ub2 = xBar_APE_Ub2/N2
    @. xBar_APE_Vb2 = xBar_APE_Vb2/N2
    
    # static stability WN2
    #APE_WN2     = W.*B
    @. buf = W*B
    mean!(xBar_APE_WN2,buf)
    
    # RAD generation
    @. buf = RAD_b.*B
    mean!(xBar_APE_RAD,buf);  
    @. xBar_APE_RAD = xBar_APE_RAD / N2
    
    # Diabatic_other
    @. buf = Diabatic_other.*B
    mean!(xBar_APE_DIA,buf)
    @. xBar_APE_DIA = xBar_APE_DIA/N2;  
    
    # Surface fluxes contribution 
        
    xBar_APE_Fs = mean(B[:,:,1,:].*Fs, dims=(1,2))[1,1,1,:]./N2[1,1,1,:];
    # interpolation 
    k_up              = argmin(abs.(z.-z_up));
    z1                = z[1]:dz:z[k_up];
    rho01 = zeros(length(z1),length(t))
    rho02 = zeros(length(z1),length(t))
    xBar_APE_b21 = similar(rho01)
    xBar_APE_RAD1 = similar(rho01)
    xBar_APE_DIA1 = similar(rho01)
    xBar_APE_WN21 = similar(rho01)
    xBar_APE_Ub21 = similar(rho01)
    xBar_APE_Vb21 = similar(rho01)
    xBar_KE1 = similar(rho01)
    xBar_APE_rate1 = similar(rho01)
    @inbounds for time in 1:length(t)
        rho01_itp         = LinearInterpolation(z, rho0[:,time]);
        xBar_APE_b21_itp  = LinearInterpolation(z, APE_b2[1,1,:,time]);
        xBar_APE_RAD1_itp = LinearInterpolation(z, xBar_APE_RAD[1,1,:,time]);
        xBar_APE_DIA1_itp = LinearInterpolation(z, xBar_APE_DIA[1,1,:,time]);
        xBar_APE_WN21_itp = LinearInterpolation(z, xBar_APE_WN2[1,1,:,time]);
        xBar_APE_Ub21_itp = LinearInterpolation(z, xBar_APE_Ub2[1,1,:,time]);
        xBar_APE_Vb21_itp = LinearInterpolation(z, xBar_APE_Vb2[1,1,:,time]);
        xBar_KE1_itp      = LinearInterpolation(z, xBar_KE[1,1,:,time]);
        xBar_APE_rate1_itp    = LinearInterpolation(z, xBar_APE_rate[:,time]);
                            
                            
        rho01[:,time]              = [rho01_itp(x) for x in z1]
        xBar_APE_b21[:,time]       = [xBar_APE_b21_itp(x) for x in z1]
        xBar_APE_RAD1[:,time]      = [xBar_APE_RAD1_itp(x) for x in z1]
        xBar_APE_DIA1[:,time]      = [xBar_APE_DIA1_itp(x) for x in z1]
        xBar_APE_WN21[:,time]      = [xBar_APE_WN21_itp(x) for x in z1]
        xBar_APE_Ub21[:,time]      = [xBar_APE_Ub21_itp(x) for x in z1]
        xBar_APE_Vb21[:,time]      = [xBar_APE_Vb21_itp(x) for x in z1]
        xBar_KE1[:,time]           = [xBar_KE1_itp(x) for x in z1]
        xBar_APE_rate1[:,time]     = [xBar_APE_rate1_itp(x) for x in z1]
    end
    
    
    
    # vertical integration
    dmass = rho01.*dz
    int_mass     = sum(dmass,dims=1)[1,:]
    int_KE       = sum(dmass*xBar_KE1,dims=1)[1,:]
    int_APE      = sum(dmass.*xBar_APE_b21,dims=1)[1,:]
    int_APE_RAD  = sum(dmass.*xBar_APE_RAD1,dims=1)[1,:]
    int_APE_DIA  = sum(dmass.*xBar_APE_DIA1,dims=1)[1,:]
    int_APE_WN2  = sum(dmass.*xBar_APE_WN21,dims=1)[1,:]
    int_APE_Ub2  = sum(dmass.*xBar_APE_Ub21,dims=1)[1,:]
    int_APE_Vb2  = sum(dmass.*xBar_APE_Vb21,dims=1)[1,:]
    int_APE_rate = sum(dmass.*xBar_APE_rate1,dims=1)[1,:]
    
    # surface flux contribution
    #N2S         = N2[1,:]
    #APE_Fs      = B[:,:,1,:].*Fs
    
     
    residual    = int_APE_rate .+ int_APE_Ub2 .+ int_APE_Vb2+int_APE_WN2 .- (int_APE_RAD .+ int_APE_DIA .+ xBar_APE_Fs)
    
    return (int_mass, int_KE, int_APE, int_APE_rate, int_APE_Ub2,int_APE_Vb2, int_APE_WN2, int_APE_RAD, int_APE_DIA, xBar_APE_Fs, residual)
    end

"""
-------------This function computes the buoyancy budget---------------
Inputs:
B (buoyancy)
RAD_b the radiative heating, converted to units of buoyancy
Fs the Surface fluxes
U,V,W The three dimensional velocities 
N2 The Brunt Va"isala frequency squared
dx,dy,dz,dt the steps in each coordinate
x,y,z,t the coordinate vectors

"""
function buoyancybudget(B, RAD_b, Fs, U,V, W, N2, dx,dy, dz, dt, x,y, z, t)
#************ Array creation **************#
Qs      = zeros(typeof(B[1]),length(x),length(y),length(z),length(t))
B_ghost = Array{typeof(B[1])}(undef, length(x)+1,length(y)+1, length(z), length(t))
WN2     = similar(W)
dBdt    = similar(B)
UdBdx   = similar(U)
VdBdy   = similar(U)
#************ WN2 *****************

 WN2  .= reshape(N2,(1,1,size(N2,1),size(N2,2))).*W

#*************  Advection **************

@. B_ghost[1:end-1,1:end-1,:,:] = B
@. B_ghost[end,1:end-1,:,:] = B[1,:,:,:]
@. B_ghost[1:end-1,end,:,:] = B[:,1,:,:]

@. UdBdx = U*(B_ghost[2:end,1:end-1,:,:]-B_ghost[1:end-1,1:end-1,:,:])/dx
@. VdBdy = V*(B_ghost[1:end-1,2:end,:,:]-B_ghost[1:end-1,1:end-1,:,:])/dy


# @. UdBdx[1:end-1,:,:,:] .= U[1:end-1,:,:,:]*(B[2:end,:,:,:]-B[1:end-1,:,:,:])/dx
# @. UdBdx[end,:,:,:] = U[end,:,:,:]*(B[1,:,:,:]-B[end,:,:,:])/dx

# @. VdBdy[:,1:end-1,:,:] = V[:,1:end-1,:,:].*(B[:,2:end,:,:]-B[:,1:end-1,:,:])/dy
# @. VdBdy[:,end,:,:] = V[:,end,:,:]*(B[:,1,:,:]-B[:,end,:,:])/dy

Qs[:,:,1,:] .= Fs./dz
#************ Time evolution *************#

@. dBdt[:,:,:,1:end-1] = (B[:,:,:,2:end] - B[:,:,:,1:end-1])/dt; 
@. dBdt[:,:,:,end]     = dBdt[:,:,:,end-1]/dt
#*************** Return ********************#

Diabatic_other  = dBdt .+ UdBdx .+ VdBdy .+ WN2 .- RAD_b .- Qs

return dBdt, UdBdx,VdBdy, WN2, Qs, Diabatic_other
end




function buoyancy_4d(exp_name,outputInterval,U,V,W,RAD,T,qv,PP,pw,SHF,LHF,x,y,z,t,P0)
day   = (86400)
sst   = (300)
dt = outputInterval
#########################Empty array creation #############
Tv = similar(T)
ThetaV = similar(T)

dayLength = (60*60*24/outputInterval); #How many data points make one day
R       = (287); # Ideal gas constants
heat_capacity      = (1004); #Heat capacity of air
L       = (2.5*1e6); #Enthalpy of phase change
epsilon     = (29/18-1); 
g       = (10); #acceleration of gravity
#println("Daylength: $dayLength")
#Smothing is equivalent in time to picking the slow component
if outputInterval == 14400     
    smooth_x    = 5 
    smooth_y    = 5
elseif outputInterval == 7200
    smooth_x    = 11
    smooth_y    = 11
end
smooth_time = floor(Int,dayLength*5)   
dx    = x[2]-x[1] #Grid size()
dy    = dx #Grid size()
kz      = length(z) # vertical levels
kx      = length(x) # # of horizonal grid points
ky      = length(y) # # of horizonal grid points
RAD     .= RAD/day # K/s #Heating rate per second
Tv     .= (1 .+ epsilon*qv*1e-3).*T #Virtual temperature
P0      = P0*1e2
Pref    = (1000*1e2) #Pa
Ts      = sst #Sea surface temperature
qs      = (25.7*1e-3) 
Tvs     = Ts*(1+epsilon*qs)
@. SHF = g/(1*heat_capacity*Ts)*(SHF) # rho := 1
@. LHF = g/(1*heat_capacity*Ts)*(epsilon*heat_capacity*Ts/L*LHF) # rho := 1
@. SHF  = SHF + LHF

 PP .= PP .+ reshape(P0,(1,1,kz,1))
 xBar_Pt     = mean(PP,dims=(1,2)) # use first 20-day mean as the reference
    c1 = (R/heat_capacity)
ThetaV  .= Tv.*(PP./reshape(P0,(1,1,kz,1))).^c1 #Virtual potential temp


#println(" $exp_name Smoothing data (this is the longest part)... ")
#getsmoothdata!(U,V,W, Tv, ThetaV, RAD, SHF, smooth_x,smooth_y,smooth_time,1)

println(" $exp_name: Smoothing data (this is the longest part)... ")
#getsmoothdata!(U,V,W, Tv, ThetaV, RAD, SHF, smooth_x,smooth_y,smooth_time,1)


buf3d = similar(U)
buf2d = similar(SHF)
filter_array!(buf3d,U,smooth_x,smooth_time,position)
filter_array!(buf3d,V,smooth_x,smooth_time,position)
filter_array!(buf3d,W,smooth_x,smooth_time,position)
filter_array!(buf3d,Tv,smooth_x,smooth_time,position)
filter_array!(buf3d,ThetaV,smooth_x,smooth_time,position)
filter_array!(buf3d,RAD,smooth_x,smooth_time,position)
filter_array!(buf2d,SHF,smooth_x,smooth_time,position)
buf3d = []
buf2d = []

println(" $exp_name: Finished smoothing data... ")
println(size(ThetaV))



xBar_T      = mean(T,dims=(1,2)) # use first 20-day mean as the reference
xBar_Tv     = mean(Tv,dims=(1,2)) # use first 20-day mean as the reference
var_Tv = Tv .- mean(Tv,dims=(1,2))

var_ThetaV = ThetaV .- mean(ThetaV,dims=(1,2))
xBar_ThetaV = mean(ThetaV,dims=(1,2)) 

#var_ThetaV = zeros(size(ThetaV))
#var_Tv = zeros(size(ThetaV))




rho0    = dropdims(xBar_Pt/R./xBar_Tv,dims=(1,2))

# compute N2 the Brunt-Vaisala frequency squared
dxBar_Tv_dz = zeros(1,length(z),size(var_Tv,4))
@. dxBar_Tv_dz[1,1:end-1,:] = (xBar_Tv[1,1,2:end,:]-xBar_Tv[1,1,1:end-1,:])/(z[2:end]-z[1:end-1])
@. dxBar_Tv_dz[1,end,:] = dxBar_Tv_dz[1,end-1,:]


N2  = g*(dxBar_Tv_dz .+ g/heat_capacity)[1,:,:]./xBar_Tv[1,1,:,:]  #Why is this using virtual and not potential, why g/cp check

bb = findall(abs.(N2) .< 1e-6)
for i=1:length(bb)
    if bb[i][1]>2
        N2[bb[i]] = 0.5* (N2[bb[i][1]-1,bb[i][2]] + N2[bb[i][1]+1,bb[i][2]]) #If N2 is small, substite by mean of neighbours
    else
        N2[bb[i]] = N2[bb[i][1]+1,bb[i][2]] #If N2 is small, substite by mean of neighbours
    end
end

B      = g.*var_ThetaV./xBar_ThetaV
RAD   .= RAD.*(g./xBar_Tv) # convert unit to buoyancy

PP         = []
T          = []
ThetaV     = []
Tv         = []
composite  = []
qv         = []
var_ThetaV = []
var_Tv     = []

println(" $exp_name : Computing buoyancy budget... ")
# Buoyancy budget
dz          = 50
dBdt, UdBdx,VdBd, WN2, Qs, Diabatic_other = buoyancybudget(B, RAD, SHF, U,V ,W, N2, dx,dy, dz, dt, x,y, z, 1:size(B,4))
GC.gc()
dia_a = Diabatic_other .- mean(Diabatic_other,dims=(1,2))
rad_a = RAD .- mean(RAD,dims=(1,2))
B_a   =  B .- mean(B,dims=(1,2))
#dia_ape = dia_a.*B_a
#rad_ape = rad_a.*B_a
return dia_a,rad_a,B_a

end                                                            



