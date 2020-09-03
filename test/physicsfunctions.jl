@testset "Brunt-Vaisala" begin
    

@info Γneutral = AvailablePotentialEnergyFramework.g/AvailablePotentialEnergyFramework.Dryair.cp
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

