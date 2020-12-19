using AvailablePotentialEnergyFramework, ImageFiltering, Test
@testset "Array filtering" begin

    @testset "Kernel creation" begin

    @test AvailablePotentialEnergyFramework.kernel_4d(false,false) == 
    ImageFiltering.kernelfactors((centered([(1.0)]),centered([(1.0)]),centered([(1.0)]),centered([(1.0)]) ))

    @test AvailablePotentialEnergyFramework.kernel_4d(5,false) == 
    ImageFiltering.kernelfactors((centered(ones(5)./5.0),centered(ones(5)./5.0),centered([(1.0)]),centered([(1.0)]) ))


    end     
    @testset "Filter arrays 1 and 2" begin     
        @testset "Filter_array!" begin
            ones3d  = ones(20,20,20)
            ones4d  = ones(20,20,20,20)
            zeros3d = zeros(20,20,20)
            zeros4d = zeros(20,20,20,20)
            buf3d   = similar(ones3d)
            buf4d   = similar(ones4d)
            filter_array!(buf3d,zeros3d,5,5,1)    
            filter_array!(buf3d,ones3d,5,5,1)   
            filter_array!(buf4d,zeros4d,5,5,1) 
            filter_array!(buf4d,ones4d,5,5,1)
            @test ones3d   == ones(20,20,20)
            @test zeros3d  == zeros(20,20,20)
            @test ones4d   == ones(20,20,20,20)
            @test zeros4d  == zeros(20,20,20,20)
        end
        
        @testset "Filter_array_2" begin
            zeros3d = zeros(20,20,20);
            ones3d = ones(20,20,20)
            zeros4d = zeros(20,20,20,20);
            ones4d = ones(20,20,20,20)
            
            filter_array_2!(ones3d,5,5,1)
            filter_array_2!(zeros3d,5,5,1)
            filter_array_2!(ones4d,5,5,1) 
            filter_array_2!(zeros4d,5,5,1)

            
            @test ones3d == ones(20,20,20)
            @test zeros3d == zeros(20,20,20)
            @test ones4d == ones(20,20,20,20)
            @test zeros4d == zeros(20,20,20,20)
            
        end
        
         @testset "Both filters simultaneously" begin
            arr3d = rand(20,20,20)
            arr3d2 = copy(arr3d)
            arr3d3 = copy(arr3d)
            arr4d = rand(20,20,20,20)
            arr4d2 = copy(arr4d)
            arr4d3 = copy(arr4d)
            buf3d   = similar(arr3d)
            buf4d   = similar(arr4d)

            filter_array_2!(arr3d2,5,5,1)
            filter_array_2!(arr4d2,5,5,1)
            filter_array!(buf3d,arr3d3,5,5,1)
            filter_array!(buf4d,arr4d3,5,5,1)
            @test arr3d3 == arr3d2
            @test arr4d3 == arr4d2
            
        end
    end
end
