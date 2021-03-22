@testset "Brunt-Vaisala" begin
    

    Γneutral = AvailablePotentialEnergyFramework.g/AvailablePotentialEnergyFramework.Dryair.cp
    Γstable = Γneutral - 10e-3u"K/m"
    Γunstable = Γneutral + 10e-3u"K/m"


    z = 1u"m" .* collect(0:50:10e3)
    tvprofile(Γ,z) = reshape(repeat(300 .- Γ*z,1,2),1,1,length(z),2)
    tvprofile(Γ :: Quantity, z :: Array{ <: Quantity}) = reshape(repeat(300u"K" .- Γ*z,1,2),1,1,length(z),2)

    function plot_N2(Γ,z)
        tv_profile = tvprofile(Γ,z)
        p1 = plot(profile1[1,1,:,1],z,title = "Tv profile, Γ = $Γ")
        p2 = plot(AvailablePotentialEnergyFramework.compute_N2_attempt(tv_profile,z),z, title = "N2")
        plot(p1,p2)
    end

    @test isapprox( 1u"1/s/s" .* zeros(length(z),2),compute_N2(tvprofile(Γneutral,z),z), atol=1e-10u"1/s/s")

    @test all(compute_N2(tvprofile(Γunstable,z),z) .< 0u"1/s/s" )

    @test all(compute_N2(tvprofile(Γstable,z),z) .> 0u"1/s/s" )

    @test isapprox(  zeros(length(z),2),compute_N2(tvprofile(ustrip(Γneutral),ustrip.(z)),ustrip.(z)), atol=1e-10)

    @test all(compute_N2(tvprofile(ustrip(Γunstable),ustrip.(z)),ustrip.(z)) .< 0 )

    @test all(compute_N2(tvprofile(ustrip(Γstable),ustrip.(z)),ustrip.(z)) .> 0 )

end

@testset "Thermodynamics" begin

    @test get_saturation_vapor_pressure(273.15u"K") == 6.112u"hPa"
    @test get_partial_vapor_pressure(0,1000u"hPa") == 0u"hPa"
    @test get_partial_vapor_pressure(1,1000u"hPa") == 1000u"hPa"/(18.016/28.966 + 1.0)
    @test get_mixing_ratio(0u"hPa",1000u"hPa") == 0
    @test get_mixing_ratio(get_partial_vapor_pressure(0.5,1000.0), 1000.0) == 0.5
    @test get_virtual_temperature(300u"K",0,0) == 300u"K"
    @test specific_humidity_to_mixing_ratio(mixing_ratio_to_specific_humidity(0.5)) ≈ 0.5
    @test mixing_ratio_to_specific_humidity(specific_humidity_to_mixing_ratio(0.5)) ≈ 0.5
    @test unit(get_specific_entropy(300u"K",0.2,1000u"hPa"))== u"J/K/kg"
    @test get_potential_temperature(300u"K",1000u"hPa",1000u"hPa") == 300u"K"
    @test get_potential_temperature(300u"K",1010u"hPa",1000u"hPa") < 300u"K"
    @test get_potential_temperature(300u"K",900u"hPa",1000u"hPa") > 300u"K"
    @test get_virtual_temperature(300u"K",0u"g/kg") == 300u"K"
    @test get_virtual_temperature(300u"K",0u"g/kg") == 300u"K"
    @test get_virtual_temperature(300u"K",10u"g/kg") > 300u"K"
    @test 1unit(surface_sensible_heat_flux_to_buoyancy(300u"K", 100u"W/m^2")) == 1unit(g)/u"s"*u"m"
    @test 1unit(surface_latent_heat_flux_to_buoyancy(300u"K", 100u"W/m^2")) == 1unit(g)/u"s"*u"m"
    # No units

    @test get_saturation_vapor_pressure(273.15) == 6.112
    @test get_partial_vapor_pressure(0,1000) == 0
    @test get_partial_vapor_pressure(1,1000) == 1000/(18.016/28.966 + 1.0)
    @test get_mixing_ratio(0,1000) == 0
    @test get_potential_temperature(300,1000,1000) == 300
    @test get_potential_temperature(300,1010,1000) < 300
    @test get_potential_temperature(300,900,1000) > 300
    @test get_virtual_temperature(300,0) == 300
    @test get_virtual_temperature(300,0) == 300
    @test get_virtual_temperature(300,10) > 300
    @test ustrip(surface_sensible_heat_flux_to_buoyancy(300u"K", 100u"W/m^2")) ==  surface_sensible_heat_flux_to_buoyancy(300, 100)
    @test ustrip(surface_latent_heat_flux_to_buoyancy(300u"K", 100u"W/m^2")) ==  surface_latent_heat_flux_to_buoyancy(300, 100)

    
    
    pres = Dataset(joinpath(@__DIR__,"testfiles/thermoprofile.nc")) do ds 
        1u"hPa" .* variable(ds, "PRES")[:,:]
    end
    tabs = Dataset(joinpath(@__DIR__,"testfiles/thermoprofile.nc")) do ds 
        1u"K" .* variable(ds, "TABS")[:,:]
    end
    qv = 1e-3u"kg/g" .* 1u"g/kg" .* Dataset(joinpath(@__DIR__,"testfiles/thermoprofile.nc")) do ds 
        variable(ds, "QV")[:,:] #was originally in g/kg
    end
    @info size(pres) size(qv) size(tabs)
    r = specific_humidity_to_mixing_ratio.(qv)
    timeindex = 1200
    pparcel = pres[1,timeindex]
    tparcel = tabs[1,timeindex]
    rparcel = r[1,timeindex]

    #I will create a similar profile but with a perturbation to see what happens
    tabs_unstable = copy(tabs)
    tabs_unstable[2:40,:] .- 7.0u"K"
    tabs_unstable[41:end,:] .+ 7.0u"K"

    @test_broken unit(get_buoyancy_of_lifted_parcel(tparcel,rparcel,pparcel,tabs[:,timeindex],r[:,timeindex],pres[:,timeindex])[1]) == u"K"
    @test_broken get_buoyancy_of_lifted_parcel(tparcel,rparcel,pparcel,tabs_unstable[:,timeindex],r[:,timeindex],pres[:,timeindex])
end
