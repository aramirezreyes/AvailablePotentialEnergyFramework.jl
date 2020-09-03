
using NCDatasets
using AvailablePotentialEnergyFramework
using Statistics
using Plots
using Distributed
using JLD


"""
    getapeanalysis(file2d,file3d,output_name,outputInterval,FloatType::Type=Float64)



"""
function testAPEBudget()
file2d = "/Users/arreyes/Documents/Research/Developement/testNetCDF/timeSlab_2d.nc"
file3d = "/Users/arreyes/Documents/Research/Developement/testNetCDF/timeSlab_3d.nc"
    
Float_type = Float32
   
    day           = 86400
    sst           = 300
    dt            = 14400
    Pref          = 1000*1e2                 #Pa
    Ts            = sst                      #Sea surface temperature
    qs            = 25.7*1e-3 
    Tvs           = Ts*(1+epsilon*qs)
    c1            = (R/Dryair.cp)
    


    
    ds3d                = Dataset(file3d)
    ds2d                = Dataset(file2d)
    x                   = variable(ds3d,"x")[:]    :: Array{Float_type,1}
    y                   = variable(ds3d,"y")[:]    :: Array{Float_type,1}
    z                   = variable(ds3d,"z")[:]    :: Array{Float_type,1}
    t                   = variable(ds3d,"time")[:] :: Array{Float_type,1}
    P0                  = variable(ds3d,"p")[:]    :: Array{Float_type,1}
    U                   = variable(ds3d,"U")[:,:,:,:] :: Array{Float_type,4}
    V                   = variable(ds3d,"V")[:,:,:,:] :: Array{Float_type,4}
    W                   = variable(ds3d,"W")[:,:,:,:] :: Array{Float_type,4}
    RAD                 = variable(ds3d,"QRAD")[:,:,:,:] :: Array{Float_type,4}
    T                   = variable(ds3d,"TABS")[:,:,:,:] :: Array{Float_type,4}
    Tv                  = variable(ds3d,"QV")[:,:,:,:] :: Array{Float_type,4}
    PP                  = variable(ds3d,"PP")[:,:,:,:] :: Array{Float_type,4}
    SHF                 = variable(ds2d,"SHF")[:,:,:] :: Array{Float_type,3}
    LHF                 = variable(ds2d,"LHF")[:,:,:] :: Array{Float_type,3}
    close(ds3d)   
    close(ds2d)   

    ########## FilFloat_Typeering and chopping ##########
    
    @. SHF    = g/(1*Dryair.cp*Ts)*(SHF)                            
    @. LHF    = g/(1*Dryair.cp*Ts)*(epsilon*Dryair.cp*Ts/Liquidwater.Lv*LHF) 
    @. SHF    +=  LHF     # Now it is transformed

    
    ThetaV       = similar(T)
    xBar_Pt      = Array{eltype(T),4}(undef,1,1,size(PP,3),size(PP,4))
    xBar_Tv      = Array{eltype(T),4}(undef,1,1,size(Tv,3),size(Tv,4))
    xBar_ThetaV  = Array{eltype(T),4}(undef,1,1,size(ThetaV,3),size(ThetaV,4))
    #******
    dx           = x[2]-x[1]
    dy           = y[2]-y[1]                            # Grid size
    kz           = length(z)                            # vertical levels
    kx           = length(x)                            # # of horizonal grid points
    ky           = length(y)                            # # of horizonal grid points
    @.  RAD      = RAD/day                              # K/s #Heating rate per second;
    @.  Tv       = (1 + 1e-3*epsilon*Tv)*T              # Virtual temperature
    @.  P0       = P0*1e2
    PP          .= PP .+ reshape(P0,(1,1,kz,1))
    #ThetaV      .= Tv.*(PP./reshape(P0,(1,1,kz,1))).^c1 # Virtual potential temp
    ThetaV      .= Tv.*(PP/P0[1]).^c1 # Virtual potential temp
    mean!(xBar_Pt,PP)                                     
    mean!(xBar_Tv,Tv)              
    mean!(xBar_ThetaV,ThetaV)   
    var_Tv       =  Tv     .- xBar_Tv
    var_ThetaV   =  ThetaV .- xBar_ThetaV
    rho0         = dropdims(xBar_Pt./R./xBar_Tv,dims=(1,2))
    B            = g .* var_ThetaV./xBar_ThetaV
    @. RAD       = RAD*(g/xBar_Tv)                    # convert unit to buoyancy/s
    
    
    N2           = compute_N2(xBar_Tv,z)

    # PP           = []
    # T          = []
    # ThetaV     = []
    # Tv         = []
    # var_ThetaV = []
    # var_Tv     = []
 
    ### NOTE that SHF is now the sum, saving memory
    # Buoyancy budget
    dz          = 50
    #@info size(B), size(RAD), size(SHF), size(U),size(V) ,size(W), size(N2), size(dx),size(dy), size(dz), size(dt), size(x),size(y), size(z), size(t)
    #    dBdt, UdBdx,VdBdy, WN2, Qs, Diabatic_other = buoyancybudget(B, RAD, SHF, U,V ,W, N2, dx,dy, dz, dt, x,y, z, t);
    Diabatic_other = get_diabatic_as_residual_buoyancy(B, RAD, SHF, U,V ,W, N2, dx,dy, dz, dt);
    @info "hola"

    # APE budget
    z_up        = 15000
    z_BL        = 2000

        (int_mass,
     int_KE,
     int_APE,
     int_APE_rate,
     int_APE_Ub2,
     int_APE_Vb2,
     int_APE_WN2,
     int_APE_RAD,
     int_APE_DIA,
     xBar_APE_Fs,
     residual) = getapebudget(B, U,V, W, N2, RAD, SHF, Diabatic_other, rho0, x,y, z, t, dx,dy, dz, dt, z_up)



    Diabatic_other .= Diabatic_other .- mean(Diabatic_other,dims=(1,2)) #They are now perturbations
    RAD            .= RAD .- mean(RAD,dims=(1,2))
    B              .= B .- mean(B,dims=(1,2))



    dia_ape = Diabatic_other.*B
    rad_ape = RAD.*B   

 

    # jldopen(string("outputfile.jld"), "w") do file
    #     write(file,"int_APE",int_APE)
    #     write(file,"int_KE",int_KE)
    #     write(file,"int_RAD",int_APE_RAD)
    #     write(file,"int_DIA",int_APE_DIA)
    #     write(file,"int_WN2",int_APE_WN2)
    #     write(file,"int_Ub2",int_APE_Ub2)
    #     write(file,"int_Vb2",int_APE_Vb2)
    #     write(file,"int_APE_rate",int_APE_rate)
    #     write(file,"APE_Fs",xBar_APE_Fs)
    #     write(file,"convec_heating_anomaly",Diabatic_other)
    #     write(file,"rad_heating_anomaly",RAD)
    #     write(file,"buoyancy_anomaly",B)
    #     write(file,"radiative_ape_production",dia_ape)
    #     write(file,"convective_ape_production",rad_ape)
    # end
    
    plot([int_APE_rate,
    int_APE_Ub2,
    int_APE_Vb2,
    int_APE_WN2,
    int_APE_RAD,
    int_APE_DIA,
    xBar_APE_Fs],label=["rate" "advu" "advv" "wn2" "rad" "dia" "fs" ])
plot!(residual,lw=4,label="residual")

end


"""
    getapeanalysis(file2d,file3d,output_name,outputInterval,FloatType::Type=Float64)



"""
function testAPEBudget_distributed()

    file2d = "/Users/arreyes/Documents/Research/Developement/testNetCDF/timeSlab_2d.nc"
    file3d = "/Users/arreyes/Documents/Research/Developement/testNetCDF/timeSlab_3d.nc"
    
    Float_type = Float32
   
    day           = 86400
    sst           = 300
    dt            = 14400
    Pref          = 1000*1e2                 #Pa
    Ts            = sst                      #Sea surface temperature
    qs            = 25.7*1e-3 
    Tvs           = Ts*(1+epsilon*qs)
    c1            = (R/Dryair.cp)
    


    
    # ds3d                = Dataset(file3d)
    # ds2d                = Dataset(file2d)
    # x                   = variable(ds3d,"x")[:]    :: Array{Float_type,1}
    # y                   = variable(ds3d,"y")[:]    :: Array{Float_type,1}
    # z                   = variable(ds3d,"z")[:]    :: Array{Float_type,1}
    # t                   = variable(ds3d,"time")[:] :: Array{Float_type,1}
    # P0                  = variable(ds3d,"p")[:]    :: Array{Float_type,1}
    # U                   = variable(ds3d,"U")[:,:,:,:] :: Array{Float_type,4}
    # V                   = variable(ds3d,"V")[:,:,:,:] :: Array{Float_type,4}
    # W                   = variable(ds3d,"W")[:,:,:,:] :: Array{Float_type,4}
    # RAD                 = variable(ds3d,"QRAD")[:,:,:,:] :: Array{Float_type,4}
    # T                   = variable(ds3d,"TABS")[:,:,:,:] :: Array{Float_type,4}
    # Tv                  = variable(ds3d,"QV")[:,:,:,:] :: Array{Float_type,4}
    # PP                  = variable(ds3d,"PP")[:,:,:,:] :: Array{Float_type,4}
    # SHF                 = variable(ds2d,"SHF")[:,:,:] :: Array{Float_type,3}
    # LHF                 = variable(ds2d,"LHF")[:,:,:] :: Array{Float_type,3}
    # close(ds3d)   
    # close(ds2d)

    smooth_x      = smooth_y = 15 #it was 11
    smooth_time   = floor(Int,12*5)+1 
    position = 2
    
    sx  = sy = 100
    sz  = 80
    st  = 120
    dx  = dy = 1000
    dz  = 50
    dt  = 7200

    x   = y = collect(range(0,step = dx,length=sx));
    z   = collect(range(0, step = dz,length = sz));
    t   = collect(range(0, step = dt, length = st));


    P0  = 1000*exp.(z./8000);
    U   = V = 5.0rand(sx,sy,sz,st);
    W   = 0.2rand(sx,sy,sz,st);
    RAD = rand(sx,sy,sz,st);;
    T   = 300 .+ rand(sx,sy,sz,st);
    Tv  = 0.2rand(sx,sy,sz,st);
    PP  = rand(sx,sy,sz,st);
    SHF = rand(sx,sy,st);
    LHF = rand(sx,sy,st);
    B   = similar(RAD);
    ########## smoothing  and chopping ##########

@info "Smoothing batch 1"
    fU   = @spawnat 2 cutborders!(filter_array(U,smooth_x,smooth_time,1),smooth_time,position)
    U = fetch(fU)
    fV   = @spawnat 3 cutborders!(filter_array(V,smooth_x,smooth_time,1),smooth_time,position)
    V = fetch(fV)
    fW   = @spawnat 2 cutborders!(filter_array(W,smooth_x,smooth_time,1),smooth_time,position)
    W = fetch(fW)
    fTv  = @spawnat 3 cutborders!(filter_array(Tv,smooth_x,smooth_time,1),smooth_time,position)
    Tv = fetch(fTv)
    
    fetch(@spawnat 2 GC.gc())
    fetch(@spawnat 3 GC.gc())

    @info "smoothing batch 2"
    fT   = @spawnat 2 cutborders!(filter_array(T,smooth_x,smooth_time,1),smooth_time,position)
    fPP  = @spawnat 3 cutborders!(filter_array(PP,smooth_x,smooth_time,1),smooth_time,position)
    fRAD = @spawnat 3 cutborders!(filter_array(RAD,smooth_x,smooth_time,1),smooth_time,position)
    SHF = cutborders!(filter_array(SHF,smooth_x,smooth_time,1),smooth_time,position)
    LHF = cutborders!(filter_array(LHF,smooth_x,smooth_time,1),smooth_time,position)
    t = cutborders!(t,smooth_time,position)
    
    T = fetch(fT)
    PP = fetch(fPP)
    RAD = fetch(fRAD)
    fetch(@spawnat 2 GC.gc())
    fetch(@spawnat 3 GC.gc())

    GC.gc()
    GC.gc()

    ## Processing

    @info "Beginning processing"
    
    @. SHF    = g/(1*Dryair.cp*Ts)*(SHF)                            
    @. LHF    = g/(1*Dryair.cp*Ts)*(epsilon*Dryair.cp*Ts/Liquidwater.Lv*LHF) 
    @. SHF    +=  LHF     # Now it is transformed

    
    ThetaV       = similar(T)
    xBar_Pt      = Array{eltype(T),4}(undef,1,1,size(PP,3),size(PP,4))
    xBar_Tv      = Array{eltype(T),4}(undef,1,1,size(Tv,3),size(Tv,4))
    xBar_ThetaV  = Array{eltype(T),4}(undef,1,1,size(ThetaV,3),size(ThetaV,4))
    #******
    dx           = x[2]-x[1]
    dy           = y[2]-y[1]                            # Grid size
    kz           = length(z)                            # vertical levels
    kx           = length(x)                            # # of horizonal grid points
    ky           = length(y)                            # # of horizonal grid points
    @.  RAD      = RAD/day                              # K/s #Heating rate per second;
    @.  Tv       = (1 + 1e-3*epsilon*Tv)*T              # Virtual temperature
    @.  P0       = P0*1e2
    PP          .= PP .+ reshape(P0,(1,1,kz,1))
    #ThetaV      .= Tv.*(PP./reshape(P0,(1,1,kz,1))).^c1 # Virtual potential temp
    ThetaV      .= Tv.*(PP/P0[1]).^c1 # Virtual potential temp
    mean!(xBar_Pt,PP)                                     
    mean!(xBar_Tv,Tv)              
    mean!(xBar_ThetaV,ThetaV)   
    var_Tv       =  Tv     .- xBar_Tv
    var_ThetaV   =  ThetaV .- xBar_ThetaV
    rho0         = dropdims(xBar_Pt./R./xBar_Tv,dims=(1,2))
    B            = g .* var_ThetaV./xBar_ThetaV
    @. RAD       = RAD*(g/xBar_Tv)                    # convert unit to buoyancy/s
    
    
    N2           = compute_N2(xBar_Tv,z)

    # PP           = []
    # T          = []
    # ThetaV     = []
    # Tv         = []
    # var_ThetaV = []
    # var_Tv     = []
 
    ### NOTE that SHF is now the sum, saving memory
    # Buoyancy budget
    dz          = 50
    #@info size(B), size(RAD), size(SHF), size(U),size(V) ,size(W), size(N2), size(dx),size(dy), size(dz), size(dt), size(x),size(y), size(z), size(t)
    @info "Computing buoyancy"
     Diabatic_other = fetch(@spawnat 3  get_diabatic_as_residual_buoyancy(B, RAD, SHF, U,V ,W, N2, dx,dy, dz, dt, x,y, z, t))
    fetch(@spawnat 3 GC.gc())

    # APE budget
    z_up        = 15000
    z_BL        = 2000
    @info "Computing APE budget"
        (int_mass,
     int_KE,
     int_APE,
     int_APE_rate,
     int_APE_Ub2,
     int_APE_Vb2,
     int_APE_WN2,
     int_APE_RAD,
     int_APE_DIA,
     xBar_APE_Fs,
         residual) = fetch( @spawnat 3 getapebudget(B, U,V, W, N2, RAD, SHF, Diabatic_other, rho0, x,y, z, t, dx,dy, dz, dt, z_up))
    fetch(@spawnat 3 GC.gc())



    Diabatic_other .= Diabatic_other .- mean(Diabatic_other,dims=(1,2)) #They are now perturbations
    RAD            .= RAD .- mean(RAD,dims=(1,2))
    B              .= B .- mean(B,dims=(1,2))



    dia_ape = Diabatic_other.*B
    rad_ape = RAD.*B   

 

    jldopen(string(ENV["TMPDIR"],"/outputfile.jld"), "w") do file
        write(file,"int_APE",int_APE)
        write(file,"int_KE",int_KE)
        write(file,"int_RAD",int_APE_RAD)
        write(file,"int_DIA",int_APE_DIA)
        write(file,"int_WN2",int_APE_WN2)
        write(file,"int_Ub2",int_APE_Ub2)
        write(file,"int_Vb2",int_APE_Vb2)
        write(file,"int_APE_rate",int_APE_rate)
        write(file,"APE_Fs",xBar_APE_Fs)
        write(file,"convec_heating_anomaly",Diabatic_other)
        write(file,"rad_heating_anomaly",RAD)
        write(file,"buoyancy_anomaly",B)
        write(file,"radiative_ape_production",dia_ape)
        write(file,"convective_ape_production",rad_ape)
    end
    
    # plot([int_APE_rate,
#     int_APE_Ub2,
#     int_APE_Vb2,
#     int_APE_WN2,
#     int_APE_RAD,
#     int_APE_DIA,
#     xBar_APE_Fs],label=["rate" "advu" "advv" "wn2" "rad" "dia" "fs" ])
# plot!(residual,lw=4,label="residual")

end

""" 
     run_distributed_test()
Don't forget to first run: \n
addprocs(2) \n
@everywhere using Pkg \n
@everywhere Pkg.activate("/users/arreyes/Documents/Research/Developement/AvailablePotentialEnergyFramework/") \n
@everywhere using AvailablePotentialEnergyFramework \n
"""    
function run_distributed_test()
    @info """ Don't forget to first run: \n
    addprocs(2) \n
    @everywhere using Pkg \n
    @everywhere Pkg.activate("/users/arreyes/Documents/Research/Developement/AvailablePotentialEnergyFramework/") \n
    @everywhere using AvailablePotentialEnergyFramework \n
"""
    for i in 1:10
        @info "Running ape budget, this is the $i iteration"
        testAPEBudget_distributed()
    end
end
