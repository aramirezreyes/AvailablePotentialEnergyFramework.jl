@testset "Brunt-Vaisala" begin
    

Γneutral = AvailablePotentialEnergyFramework.g/AvailablePotentialEnergyFramework.Dryair.cp
Γstable = Γneutral + 5
Γunstable = Γneutral - 50

z = collect(0:50:10e3)
tvprofile(Γ,z) = reshape(repeat(300 .- Γ*z./1000,1,2),1,1,length(z),2)

@test compute_N2(tvprofile(Γneutral,z),z) ≈ zeros(length(z),2)

@test all(compute_N2(tvprofile(Γunstable,z),z) .< 0 )

@test all(compute_N2(tvprofile(Γstable,z),z) .> 0 )

end

