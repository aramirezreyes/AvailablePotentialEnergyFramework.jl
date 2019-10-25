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

function getapebudget_old_old(B, U,V, W, N2, RAD_b, Fs, Diabatic_other, rho0, x,y, z, t, dx,dy, dz, dt, z_up)

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
xBar_APE_rate = Array{eltype(B),2}(undef,length(z), length(t))
xBar_APE_rate[:,1:end-1] = (APE_b2[:,2:end] - APE_b2[:,1:end-1])/dt; 
xBar_APE_rate[:,end] = xBar_APE_rate[:,end-1]

#*************  UdxB2 **************
b2_ghost= Array{typeof(B)}(undef, length(x)+1,length(y)+1, length(z), length(t))
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
    @info z1
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

function getapebudget_old(B, U,V, W, N2, RAD_b, Fs, Diabatic_other, rho0, x,y, z, t, dx,dy, dz, dt, z_up)
    N2            = reshape(N2,1,1,length(z),length(t)) 
    #***********Empty array generation***********#
    T             = eltype(B)
    lt            = length(t)
    lz            = length(z)
    lx            = length(y)
    ly            = length(x)
    buf           = Array{T}(undef,lx,ly,lz, lt)
    xBar_KE       = Array{T}(undef,1,1,lz, lt)
    APE_b2        = Array{T}(undef,1,1,lz, lt)
    xBar_APE_rate = Array{T}(undef,lz, lt)
    b2_ghost      = Array{T}(undef, length(x)+1,length(y)+1, lz, lt)
    xBar_APE_Ub2  = Array{T}(undef,1,1,lz, lt)
    xBar_APE_Vb2  = Array{T}(undef,1,1,lz, lt)
    xBar_APE_WN2  = Array{T}(undef,1,1,lz, lt)
    xBar_APE_RAD  = Array{T}(undef,1,1,lz, lt)
    xBar_APE_DIA  = Array{T}(undef,1,1,lz, lt)
    xBar_APE_FS   = Array{T}(undef,1,1,lz, lt)
    #*************  APE **************
    @. buf                           = B*B/2
    mean!(APE_b2,buf);
    @. APE_b2                        = APE_b2/N2
    @. b2_ghost[1:end-1,1:end-1,:,:] = buf
    @. b2_ghost[end,1:end-1,:,:]     = buf[1,:,:,:]
    @. b2_ghost[1:end-1,end,:,:]     = buf[:,1,:,:]
 
    
    #************ KE ********************
    #KE      = U.*U/2 + V.*V/2
    @. buf  = U*U/2 + V*V/2
    xBar_KE = mean!(xBar_KE,buf)

    #************ APE rate ***************
    
    @. xBar_APE_rate[:,1:end-1] = (APE_b2[1,1,:,2:end] - APE_b2[1,1,:,1:end-1])/dt; 
    @. xBar_APE_rate[:,end] = xBar_APE_rate[:,end-1]
    
    #*************  UdxB2 **************
   
    @. buf = U*(b2_ghost[2:end,1:end-1,:,:]-b2_ghost[1:end-1,1:end-1,:,:])/dx
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
    @. xBar_APE_RAD = xBar_APE_RAD/N2
    
    # Diabatic_other
    @. buf = Diabatic_other.*B
    mean!(xBar_APE_DIA,buf)
    @. xBar_APE_DIA = xBar_APE_DIA/N2;  
    
    # Surface fluxes contribution 
        
    xBar_APE_Fs  = mean(B[:,:,1,:].*Fs, dims=(1,2))[1,1,1,:]./N2[1,1,1,:];
  
    # interpolation 
    k_up              = argmin(abs.(z.-z_up));
    z1                = z[1]:dz:z[k_up];
      #    int_mass      = Array{T}(undef,lt)
    int_mass      = zeros(T,lt)
    int_KE        = zeros(T,lt)
    int_APE       = zeros(T,lt)
    int_APE_RAD   = zeros(T,lt)
    int_APE_DIA   = zeros(T,lt)
    int_APE_WN2   = zeros(T,lt)
    int_APE_Ub2   = zeros(T,lt)
    int_APE_Vb2   = zeros(T,lt)
    int_APE_rate  = zeros(T,lt)
    for timeind in 1:lt
        rho01_itp         = interpolate((z,), rho0[:,timeind],Gridded(Linear()))
        xBar_APE_b21_itp  = interpolate((z,), APE_b2[1,1,:,timeind],Gridded(Linear()))
        xBar_APE_RAD1_itp = interpolate((z,), xBar_APE_RAD[1,1,:,timeind],Gridded(Linear()))
        xBar_APE_DIA1_itp = interpolate((z,), xBar_APE_DIA[1,1,:,timeind],Gridded(Linear()))
        xBar_APE_WN21_itp = interpolate((z,), xBar_APE_WN2[1,1,:,timeind],Gridded(Linear()))
        xBar_APE_Ub21_itp = interpolate((z,), xBar_APE_Ub2[1,1,:,timeind],Gridded(Linear()))
        xBar_APE_Vb21_itp = interpolate((z,), xBar_APE_Vb2[1,1,:,timeind],Gridded(Linear()))
        xBar_KE1_itp      = interpolate((z,), xBar_KE[1,1,:,timeind],Gridded(Linear()))
        xBar_APE_rate1_itp    = interpolate((z,), xBar_APE_rate[:,timeind],Gridded(Linear()))

        for zeta in z1
                mass = dz*rho01_itp(zeta)
                int_mass[timeind]         += mass
                int_APE[timeind]          += mass*xBar_APE_b21_itp(zeta) 
                int_APE_RAD[timeind]      += mass*xBar_APE_RAD1_itp(zeta) 
                int_APE_DIA[timeind]      += mass*xBar_APE_DIA1_itp(zeta) 
                int_APE_WN2[timeind]      += mass*xBar_APE_WN21_itp(zeta) 
                int_APE_Ub2[timeind]      += mass*xBar_APE_Ub21_itp(zeta) 
                int_APE_Vb2[timeind]      += mass*xBar_APE_Vb21_itp(zeta) 
                int_KE[timeind]           += mass*xBar_KE1_itp(zeta) 
                int_APE_rate[timeind]     += mass*xBar_APE_rate1_itp(zeta) 
        end
    end
    residual    = int_APE_rate .+ int_APE_Ub2 .+ int_APE_Vb2 .+ int_APE_WN2 .- (int_APE_RAD .+ int_APE_DIA .+ xBar_APE_Fs) 
    return (int_mass, int_KE, int_APE, int_APE_rate, int_APE_Ub2,int_APE_Vb2, int_APE_WN2, int_APE_RAD, int_APE_DIA, xBar_APE_Fs, residual)
    end



"""
-------------Computes the APE budgets budget---------------

"""

function getapebudget(B, U,V, W, N2, RAD_b, Fs, Diabatic_other, rho0, x,y, z, t, dx,dy, dz, dt, z_up)
    N2            = reshape(N2,1,1,length(z),length(t)) 
    #***********Empty array generation***********#
    T             = eltype(B)
    lt            = length(t)
    lz            = length(z)
    buf           = similar(U)
    xBar_KE       = Array{T}(undef,1,1,lz, lt)
    APE_b2        = Array{T}(undef,1,1,lz, lt)
    xBar_APE_rate = Array{T}(undef,lz, lt)
    b2_ghost      = Array{T}(undef, length(x)+1,length(y)+1, lz, lt)
    xBar_APE_Ub2  = Array{T}(undef,1,1,lz, lt)
    xBar_APE_Vb2  = Array{T}(undef,1,1,lz, lt)
    xBar_APE_WN2  = Array{T}(undef,1,1,lz, lt)
    xBar_APE_RAD  = Array{T}(undef,1,1,lz, lt)
    xBar_APE_DIA  = Array{T}(undef,1,1,lz, lt)
    xBar_APE_FS   = Array{T}(undef,1,1,lz, lt)
    #*************  APE **************
    @. buf                           = B*B/2
    @. b2_ghost[1:end-1,1:end-1,:,:] = buf
    @. b2_ghost[end,1:end-1,:,:]     = buf[1,:,:,:]
    @. b2_ghost[1:end-1,end,:,:]     = buf[:,1,:,:]
    @. buf = buf/N2
    mean!(APE_b2,buf);
 
    #************ KE ********************
    #KE      = U.*U/2 + V.*V/2
    @. buf  = U*U/2 + V*V/2
    xBar_KE = mean!(xBar_KE,buf)

    #************ APE rate ***************
    @.  xBar_APE_rate[:,1:end-1] = (APE_b2[1,1,:,2:end] - APE_b2[1,1,:,1:end-1])/dt; 
    @.  xBar_APE_rate[:,end] = xBar_APE_rate[:,end-1]
    
    #*************  Advection **************
   
    @. buf = @views U*(b2_ghost[2:end,1:end-1,:,:]-b2_ghost[1:end-1,1:end-1,:,:])/N2/dx
    mean!(xBar_APE_Ub2,buf)
    @. buf = @views V*(b2_ghost[1:end-1,2:end,:,:]-b2_ghost[1:end-1,1:end-1,:,:])/N2/dy
    mean!(xBar_APE_Vb2,buf)
    ################################# static stability WN2
    #APE_WN2     = W.*B
    @. buf = W*B
    mean!(xBar_APE_WN2,buf)
    # RAD generation
    @. buf = RAD_b*B/N2
    mean!(xBar_APE_RAD,buf);  
    # Diabatic_other
    @. buf = Diabatic_other*B/N2
    mean!(xBar_APE_DIA,buf)
    # Surface fluxes contribution 
    xBar_APE_Fs  = @views mean(B[:,:,1,:].*Fs, dims=(1,2))[1,1,1,:]./N2[1,1,1,:];
    # interpolation 
    k_up              = argmin(abs.(z.-z_up));
    z1                = z[1]:dz:z[k_up];
    int_mass      = zeros(T,lt)
    int_KE        = zeros(T,lt)
    int_APE       = zeros(T,lt)
    int_APE_RAD   = zeros(T,lt)
    int_APE_DIA   = zeros(T,lt)
    int_APE_WN2   = zeros(T,lt)
    int_APE_Ub2   = zeros(T,lt)
    int_APE_Vb2   = zeros(T,lt)
    int_APE_rate  = zeros(T,lt)
    @inbounds for time in 1:lt
        rho01_itp         =  interpolate((z,), rho0[:,time],Gridded(Linear()))
        xBar_APE_b21_itp  =  interpolate((z,), APE_b2[1,1,:,time],Gridded(Linear()))
        xBar_APE_RAD1_itp =  interpolate((z,), xBar_APE_RAD[1,1,:,time],Gridded(Linear()))
        xBar_APE_DIA1_itp =  interpolate((z,), xBar_APE_DIA[1,1,:,time],Gridded(Linear()))
        xBar_APE_WN21_itp =  interpolate((z,), xBar_APE_WN2[1,1,:,time],Gridded(Linear()))
        xBar_APE_Ub21_itp =  interpolate((z,), xBar_APE_Ub2[1,1,:,time],Gridded(Linear()))
        xBar_APE_Vb21_itp =  interpolate((z,), xBar_APE_Vb2[1,1,:,time],Gridded(Linear()))
        xBar_KE1_itp      =  interpolate((z,), xBar_KE[1,1,:,time],Gridded(Linear()))
        xBar_APE_rate1_itp    =  interpolate((z,), xBar_APE_rate[:,time],Gridded(Linear()))
        @inbounds for zeta in z1
            mass = dz*rho01_itp(zeta)
            int_mass[time]         += mass
            int_APE[time]          += mass*xBar_APE_b21_itp(zeta) 
            int_APE_RAD[time]      += mass*xBar_APE_RAD1_itp(zeta) 
            int_APE_DIA[time]      += mass*xBar_APE_DIA1_itp(zeta) 
            int_APE_WN2[time]      += mass*xBar_APE_WN21_itp(zeta) 
            int_APE_Ub2[time]      += mass*xBar_APE_Ub21_itp(zeta) 
            int_APE_Vb2[time]      += mass*xBar_APE_Vb21_itp(zeta) 
            int_KE[time]           += mass*xBar_KE1_itp(zeta) 
            int_APE_rate[time]     += mass*xBar_APE_rate1_itp(zeta) 

        end
    end
    residual    = int_APE_rate .+ int_APE_Ub2 .+ int_APE_Vb2 .+ int_APE_WN2 .- (int_APE_RAD .+ int_APE_DIA .+ xBar_APE_Fs) 
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
function buoyancybudget_old(B, RAD_b, Fs, U,V, W, N2, dx,dy, dz, dt, x,y, z, t)
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
    #Qs      = zeros(eltype(B),length(x),length(y),length(z),length(t))
    #B_ghost = Array{typeof(B[1])}(undef, length(x)+1,length(y)+1, length(z), length(t))
    WN2     = similar(W)
    dBdt    = similar(B)
    UdBdx   = similar(U)
    VdBdy   = similar(U)
    #************ WN2 *****************
    
    WN2  .= reshape(N2,(1,1,size(N2,1),size(N2,2))).*W
    
    #*************  Advection **************
    
    @. @views UdBdx[1:end-1,:,:,:] .= U[1:end-1,:,:,:]*(B[2:end,:,:,:]-B[1:end-1,:,:,:])/dx
    @. @views UdBdx[end,:,:,:] = U[end,:,:,:]*(B[1,:,:,:]-B[end,:,:,:])/dx
    
    @. @views VdBdy[:,1:end-1,:,:] = V[:,1:end-1,:,:].*(B[:,2:end,:,:]-B[:,1:end-1,:,:])/dy
    @. @views VdBdy[:,end,:,:] = V[:,end,:,:]*(B[:,1,:,:]-B[:,end,:,:])/dy

    #@. Qs[:,:,1,:] = Fs/dz
#************ Time evolution *************#

    @. @views dBdt[:,:,:,1:end-1] = (B[:,:,:,2:end] - B[:,:,:,1:end-1])/dt
    @. @views dBdt[:,:,:,end]     = dBdt[:,:,:,end-1]
    #*************** Return ********************#
    
    Diabatic_other  = dBdt .+ UdBdx .+ VdBdy .+ WN2 .- RAD_b
    Diabatic_other[:,:,1,:] .= @views  Diabatic_other[:,:,1,:] .+  Fs./dz
    
    #return dBdt, UdBdx,VdBdy, WN2, Qs, Diabatic_other
    return Diabatic_other
end





