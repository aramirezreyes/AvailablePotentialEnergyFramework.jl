using AvailablePotentialEnergyFramework
using Test


# @test integrate_vertically(1:10,dz=2) == 110
# @test integrate_vertically(1:10,dz=2,weight=2) == 220
# testmat_2d= repeat(1:10,1,10)
# testweights_1d = reverse(1:10)
# testweights_2d = (1:10)'.*repeat(reverse(1:10),1,10)
# @test integrate_vertically(testmat_2d) == reshape(10*(1:10),10,1)
# @test integrate_vertically(testmat_2d,weight=testweights_1d) == reshape(10*(1:10)).*reverse((1:10),10,1)
# @test integrate_vertically(testmat_2d,weight=testweights_2d) == reshape(55*(1:10)).*reverse(1:10),10,1)

include("cyclone_detection.jl")