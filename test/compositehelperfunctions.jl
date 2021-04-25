
@testset "Composite helper functions" begin
    x1 = origin = (0,0)
    x2 = (1,1)
    x3 = (0,4)
    x4 = (0,2)
    binlimits = (0,2)
    gridspacing = 1
    @test AvailablePotentialEnergyFramework.distance(x1,x2,gridspacing) ==  AvailablePotentialEnergyFramework.distance(x2,x1,gridspacing) == sqrt(2)
    @test AvailablePotentialEnergyFramework.isindexindistancebin(binlimits,x2,x1) == true
    @test AvailablePotentialEnergyFramework.isindexindistancebin(binlimits,x3,x1) == false
    @test AvailablePotentialEnergyFramework.isindexindistancebin(binlimits,x4,x1) == true

    array_allindistance_2d = let
        array_allindistance_2d = zeros(100,100)
        for i in 1:100, j in 1:100
            array_allindistance_2d[i,j] = hypot(i,j)
        end
        array_allindistance_2d
    end
    @test isapprox(AvailablePotentialEnergyFramework.averageallindistance((99,100),array_allindistance_2d,origin), 100; rtol = 1)

    array_allindistance_3d = let
        array_allindistance_3d = zeros(100,100,10)
        for i in 1:100, j in 1:100
            array_allindistance_3d[i,j,:] .= hypot(i,j)
        end
        array_allindistance_3d
    end
    
    #@test isapprox.(AvailablePotentialEnergyFramework.averageallindistance((99,100),array_allindistance_3d,origin), 100; rtol = 1) == fill(true,10)

    @test all(isapprox.(AvailablePotentialEnergyFramework.velocity_topolar(10,0,x4,origin),(0,10),atol=0.1))

    @test all(isapprox.(AvailablePotentialEnergyFramework.velocity_topolar(-10,0,x4,origin),(0,-10),atol=0.1))

    @test all(isapprox.(AvailablePotentialEnergyFramework.velocity_topolar(0,10,x4,origin),(10,0),atol=0.1))

    @test all(isapprox.(AvailablePotentialEnergyFramework.velocity_topolar(0,-10,x4,origin),(-10,0),atol=0.1))


end


@testset "Cyclone detection tools" begin
    
    psfc = Dataset(joinpath(@__DIR__,"testfiles/test_composite_PSFC_TABS.nc")) do ds 
        variable(ds, "PSFC")[:,:,:]
    end
    psfc = cat(psfc,fill(mean(psfc) ,size(psfc)) + 5rand(size(psfc)...),dims=3)
    
    
    TABS = Dataset(joinpath(@__DIR__,"testfiles/test_composite_PSFC_TABS.nc")) do ds 
        variable(ds, "TABS")[:,:,:,:]
    end
    TABS = cat(TABS,fill(mean(TABS) ,size(TABS)) + 1rand(size(TABS)...),dims=3)
    @testset "Detect centers" begin
        pressure_anomaly = psfc .- mean(psfc,dims=(1,2))
        centerstest =  AvailablePotentialEnergyFramework.findcyclonecenters_aspressureminima(pressure_anomaly[:,:,1],-5)

        @test length(centerstest) == 5

        framewithcyclones = detect_cyclones(pressure_anomaly[:,:,1],-5,2000)
        
        @test length(framewithcyclones.labels) == 6
        @test all([framewithcyclones.segmented_frame.segment_pixel_count[key] > 1 for key in keys(framewithcyclones.segmented_frame.segment_pixel_count)])
        
        (cyclonecount_2d,addition_2d) = AvailablePotentialEnergyFramework.add_allcyclones(psfc[:,:,1],framewithcyclones; maskcyclones = false)
        @test cyclonecount_2d == 3
        
        (cyclonecount_2d,addition_2d) = AvailablePotentialEnergyFramework.add_allcyclones(psfc[:,:,1],framewithcyclones; maskcyclones = true)
        @test cyclonecount_2d == 3

        (cyclonecount_3d,addition_3d) = AvailablePotentialEnergyFramework.add_allcyclones(TABS[:,:,:,1],framewithcyclones; maskcyclones = false)
        @test cyclonecount_3d == 3
        
        (cyclonecount_3d,addition_3d) = AvailablePotentialEnergyFramework.add_allcyclones(TABS[:,:,:,1],framewithcyclones; maskcyclones = true)
        @test cyclonecount_3d == 3
        
        binlimits = 0:2000:300000
        bins = [(binlimits[ind] , binlimits[ind + 1]) for ind in 1:(length(binlimits) - 1) ]
        
        @test_nowarn [averageallindistance(bin,addition_2d,(128,128),4000) for bin in bins]

        @test_nowarn [averageallindistance(bin,addition_3d,(128,128),4000) for bin in bins]
### Try and add two frames, one with and one without TC
       frameswithcyclones = [detect_cyclones(pressure_anomaly[:,:,i],-5,2000) for i in 1:2]
        @test 3 == begin
            totalcyclonecount = 0
            for timeindex in 1:2
                framewithcyclones = frameswithcyclones[timeindex]
                if !iszero(framewithcyclones.count)
                    count, _ = AvailablePotentialEnergyFramework.add_allcyclones(TABS[:,:,timeindex],framewithcyclones;maskcyclones = false)
                    totalcyclonecount += count
                end                
            end
            totalcyclonecount
        end 
    end
end
