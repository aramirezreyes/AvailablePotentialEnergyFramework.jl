
@testset "Diagnostic utilities" begin

pw = rand(100,100,100)
precipitation = 100pw
max_precip_bin = 200
binspacing = 0.25
pw_bins = 0:binspacing:max_precip_bin

average_precip_per_bin_da = AvailablePotentialEnergyFramework.average_precipiation_per_pw_bin_dayang(pw,precipitation,max_precip_bin,binspacing)

probabilities,average_precip_per_bin = AvailablePotentialEnergyFramework.average_precipiation_per_pw_bin(pw,precipitation,pw_bins,binspacing)

@test average_precip_per_bin ≈ average_precip_per_bin_da

end