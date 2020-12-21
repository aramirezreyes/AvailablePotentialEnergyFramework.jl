@testset "Brunt-Vaisala" begin
    

    Γneutral = AvailablePotentialEnergyFramework.g/AvailablePotentialEnergyFramework.Dryair.cp
    Γstable = Γneutral - 10e-3
    Γunstable = Γneutral + 10e-3

    z = collect(0:50:10e3)
    tvprofile(Γ,z) = reshape(repeat(300 .- Γ*z,1,2),1,1,length(z),2)

    function plot_N2(Γ,z)
        tv_profile = tvprofile(Γ,z)
        p1 = plot(profile1[1,1,:,1],z,title = "Tv profile, Γ = $Γ")
        p2 = plot(AvailablePotentialEnergyFramework.compute_N2_attempt(tv_profile,z),z, title = "N2")
        plot(p1,p2)
    end

    @test isapprox( zeros(length(z),2),compute_N2(tvprofile(Γneutral,z),z), atol=1e-10)

    @test all(compute_N2(tvprofile(Γunstable,z),z) .< 0 )

    @test all(compute_N2(tvprofile(Γstable,z),z) .> 0 )

end

@testset "Thermodynamics" begin

    @test get_saturation_vapor_pressure(273.15) == 6.112
    @test get_partial_vapor_pressure(0,1000) == 0
    @test get_partial_vapor_pressure(1,1000) == 1000/(18.016/28.966 + 1.0)
    @test get_mixing_ratio(0,1000) == 0
    @test get_mixing_ratio(get_partial_vapor_pressure(0.5,1000.0), 1000.0) == 0.5
    @test get_virtual_temperature(300,0,0) == 300
    @test specific_humidity_to_mixing_ratio(mixing_ratio_to_specific_humidity(0.5)) ≈ 0.5
    @test mixing_ratio_to_specific_humidity(specific_humidity_to_mixing_ratio(0.5)) ≈ 0.5
    
    pres = Dataset(joinpath(@__DIR__,"testfiles/thermoprofile.nc")) do ds 
        variable(ds, "PRES")[:,:]
    end
    tabs = Dataset(joinpath(@__DIR__,"testfiles/thermoprofile.nc")) do ds 
        variable(ds, "TABS")[:,:]
    end
    qv = 1e-3 .* Dataset(joinpath(@__DIR__,"testfiles/thermoprofile.nc")) do ds 
        variable(ds, "QV")[:,:] #was originally in g/kg
    end
    @info size(pres) size(qv) size(tabs)
    r = specific_humidity_to_mixing_ratio.(qv)
    timeindex = 1200
    pparcel = pres[1,timeindex]
    tparcel = tabs[1,timeindex]
    rparcel = r[1,timeindex]
    @test_broken get_buoyancy_of_lifted_parcel(tparcel,rparcel,pparcel,tabs,r,pres)
end
