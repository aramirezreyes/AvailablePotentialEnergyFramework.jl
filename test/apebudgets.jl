### This will have profile tests for the ape budgets

@testset "APE budgets" begin
    sx,sy,sz,st = 100,100,100,10
    x = y = range(0.0,stop=10.0,length = sx);
    z = range(0.0,stop=500.0,length = sz);
    t = range(0.0,stop=50.0,length = st);
    dx = dy = x[2] - x[1];
    dz = z[2] - z[1];
    dt = t[2] - t[1];
    
    zero_profile = zeros(sx,sy,sz,st);
    one_profile = ones(sx,sy,sz,st);
    test_brunt = 1e-4*ones(sz,st)
    test_rho = ones(1,1,sz,st)
    test_surface = zeros(sx,sy,st)
    
    @test  getapebudget(zero_profile,zero_profile,zero_profile,zero_profile,test_brunt,
                    zero_profile,test_surface,zero_profile,test_rho,x,y,z,t,dx,dy,dz,dt,50.0)[3] ≈ zeros(st)


    ### If B is constant in space, advection should be zero
    B = 5ones(sx,sy,sz,st);
    for i in 1:st
        for j in 1:sy
            for k in 1:sx
                B[j,k,:,i] .= B[j,k,:,i] .* sin(2pi/st*i)
            end
        end
    end
    RAD_b = 2ones(sx,sy,sz,st);
    Fs = 2ones(sx,sy,st);
    U = 10ones(sx,sy,sz,st);
    V = 10ones(sx,sy,sz,st);
    W = 10ones(sx,sy,sz,st);
    N2 = 1e-4 .* ones(1,1,sz,st);
    rho0 = ones(1,1,sz,st);
    diabatic_other = get_diabatic_as_residual_buoyancy(B, RAD_b, Fs, U,V, W, N2, dx,dy, dz, dt);
    a = getapebudget(B,U,V,W,N2,RAD_b,Fs,diabatic_other, rho0, x,y,z,t,dx,dy,dz,dt,500)
    @test a[4] ≈ -(a[7] - a[8] - a[9] - a[10])

    ### Budget must be closed when we have Advection, Rad, and Surface fluxes
    B = 5ones(sx,sy,sz,st);
    for i in 1:st
        for j in 1:sy
            for k in 1:sx
                B[j,k,:,i] .= B[j,k,:,i] .* sin(2pi/st*i)*sin(2pi/sx*k)*sin(2pi/sy*j)
            end
        end
    end
    RAD_b = 2ones(sx,sy,sz,st);
    Fs = 2ones(sx,sy,st);
    U = 10rand(sx,sy,sz,st);
    V = 10rand(sx,sy,sz,st);
    W = 5rand(sx,sy,sz,st);
    N2 = 1e-4 .* ones(1,1,sz,st);
    rho0 = ones(1,1,sz,st);
    diabatic_other = get_diabatic_as_residual_buoyancy(B, RAD_b, Fs, U,V, W, N2, dx,dy, dz, dt);
    a = getapebudget(B,U,V,W,N2,RAD_b,Fs,diabatic_other, rho0, x,y,z,t,dx,dy,dz,dt,500)
    @test a[4] ≈ -(a[5] + a[6] + a[7] - a[8] - a[9] - a[10])


    ### Budget must be closed when we have Advection, Rad, and zero surface fluxes
    B = 5ones(sx,sy,sz,st);
    for i in 1:st
        for j in 1:sy
            for k in 1:sx
                B[j,k,:,i] .= B[j,k,:,i] .* sin(2pi/st*i)*sin(2pi/sx*k)*sin(2pi/sy*j)
            end
        end
    end
    RAD_b = 2ones(sx,sy,sz,st);
    Fs = zeros(sx,sy,st);
    U = 10rand(sx,sy,sz,st);
    V = 10rand(sx,sy,sz,st);
    W = 5rand(sx,sy,sz,st);
    N2 = 1e-4 .* ones(1,1,sz,st);
    rho0 = ones(1,1,sz,st);
    diabatic_other = get_diabatic_as_residual_buoyancy(B, RAD_b, Fs, U,V, W, N2, dx,dy, dz, dt);
    a = getapebudget(B,U,V,W,N2,RAD_b,Fs,diabatic_other, rho0, x,y,z,t,dx,dy,dz,dt,500)
    @test a[4] ≈ -(a[5] + a[6] + a[7] - a[8] - a[9] )


    ### Budget must be closed when we have Advection, zero Rad, and finite surface fluxes
    B = 5ones(sx,sy,sz,st);
    for i in 1:st
        for j in 1:sy
            for k in 1:sx
                B[j,k,:,i] .= B[j,k,:,i] .* sin(2pi/st*i)*sin(2pi/sx*k)*sin(2pi/sy*j)
            end
        end
    end
    RAD_b = zeros(sx,sy,sz,st);
    Fs = 2ones(sx,sy,st);
    U = 10rand(sx,sy,sz,st);
    V = 10rand(sx,sy,sz,st);
    W = 5rand(sx,sy,sz,st);
    N2 = 1e-4 .* ones(1,1,sz,st);
    rho0 = ones(1,1,sz,st);
    diabatic_other = get_diabatic_as_residual_buoyancy(B, RAD_b, Fs, U,V, W, N2, dx,dy, dz, dt);
    a = getapebudget(B,U,V,W,N2,RAD_b,Fs,diabatic_other, rho0, x,y,z,t,dx,dy,dz,dt,500)
    @test a[4] ≈ -(a[5] + a[6] + a[7] - a[9] - a[10])




    
   #  labels = ["mass" "KE" "APE" "∂tAPE" "AdvecU" "AdvecV" "Conversion" "Rad" "Convec" "Surf" "Res"]
   #  mycolors = distinguishable_colors(11) #NEEDS COLORS.JL
    

   #  fig = Figure(resolution = (1200, 700), backgroundcolor = RGBf0(0.98, 0.98, 0.98),title = "APE budget closure test")
   #  ax1 = fig[1,1] = Axis(fig)
   #  ylims!(ax1,extrema(reduce(vcat,a[4:11])))
   #  for i in 4:11
   #     lines!(ax1,a[i], label = labels[i], color = mycolors[i], linewidth = 5)
   #  end
   #  axislegend(ax1)
   # fig

    
end
