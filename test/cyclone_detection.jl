


#Test matrix with cyclones

matrixwithcyclones = zeros(Float64,20,20)
matrixwithcyclones[2,1] = -2 #Center
matrixwithcyclones[8:10,5:7] .= -1 
matrixwithcyclones[9,6] = -2 #Center
matrixwithcyclones[17:19,17:19] .= -1 
matrixwithcyclones[18,18] = -2 

## Fake pressure perturbation

@testset "Cyclone detection tools" begin
    @testset "Detect centers" begin
        centerstest =  AvailablePotentialEnergyFramework.findcyclonecenters_aspressureminima(matrixwithcyclones,0)
        @test length(centerstest) == 2
        @test centerstest[1] == CartesianIndex(9, 6)
        @test centerstest[2] == CartesianIndex(18, 18)
    end
end