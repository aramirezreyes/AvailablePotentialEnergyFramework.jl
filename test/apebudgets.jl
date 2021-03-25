### This will have profile tests for the ape budgets

@testset "APE budgets" begin
    sx,sy,sz,st = 100,100,10,10

    zero_profile = zeros(sx,sy,sz,st)
    one_profile = ones(sx,sy,sz,st)
    x = y = range(0.0,stop=10.0,length = sx)
    z = range(0.0,stop=50.0,length = sz)
    t = range(0.0,stop=50.0,length = st)
    test_brunt = 1e-4*ones(sz,st)
    test_rho = ones(1,1,sz,st)
    test_surface = zeros(sx,sy,st)
    dx = dy = x[2] - x[1]
    dz = z[2] - z[1]
    dt = t[2] - t[1]
    @test  getapebudget(zero_profile,zero_profile,zero_profile,zero_profile,test_brunt,
                    zero_profile,test_surface,zero_profile,test_rho,x,y,z,t,dx,dy,dz,dt,50.0)[3] ≈ zeros(st)

end
