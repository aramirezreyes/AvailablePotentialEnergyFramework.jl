using SAMtools, ImageFiltering, Test


@testset "APE helper functions" begin

    
    @testset "Filter arrays 1 and 2" begin
        arr3d = rand(20,20,20)
        arr4d = rand(20,20,20,20)
        @testset "Filter_array" begin
            
            @test filter_array(ones(20,20,20),5,5,1) == ones(20,20,20)
            @test filter_array(zeros(20,20,20),5,5,1) == zeros(20,20,20)
            @test filter_array(ones(20,20,20,20),5,5,1) == ones(20,20,20,20)
            @test filter_array(zeros(20,20,20,20),5,5,1) == zeros(20,20,20,20)
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
            arr3d = ones(20,20,20)
            arr4d = ones(20,20,20,20)
            filt1 = filter_array(arr3d,5,5,1)
            filter_array_2!(arr3d,5,5,1)
            filt2 = filter_array(arr4d,5,5,1)
            filter_array_2!(arr4d,5,5,1)
            filt3 = filter_array(arr3d,5,5,1)
            filt4 = filter_array(arr4d,5,5,1)
            filter_array_2!(arr3d,5,5,1)
            filter_array_2!(arr4d,5,5,1)
            @test filt1 == arr3d
            @test filt2 == arr4d
            @test filt3 == arr3d
            @test filt4 == arr4d
            arr3d = rand(20,20,20)
            arr3d2 = copy(arr3d)
            arr4d = rand(20,20,20,20)
            arr4d2 = copy(arr4d)
            filter_array_2!(arr3d,5,5,1)
            filter_array_2!(arr4d,5,5,1)
            @test arr3d != arr3d2
            @test arr4d != arr4d2
            
        end
    end
    



end
