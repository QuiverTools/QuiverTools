# this script should be ran in the terminal with the command "julia test.jl"

using LinearAlgebra, Test
include("Quivers.jl")

@testset "All tests" begin

@testset "Testing basic methods" begin
    K = GeneralizedKroneckerQuiver(4)
    
    @test underlying_graph(K) == [0 4;4 0]
    @test number_of_vertices(K) == 2
    @test number_of_arrows(K) == 4
    @test is_acyclic(K) == true
    @test is_connected(K) == true
    @test is_source(K, 1) == true && is_source(K, 2) == false
    @test is_sink(K, 2) == true && is_sink(K, 1) == false
    @test euler_matrix(K) == [1 -4; 0 1]
    @test euler_form(K, [1, 1], [1, 1]) == -2
    @test opposite_quiver(K).adjacency == transpose(K.adjacency)
    @test double_quiver(K).adjacency == K.adjacency + transpose(K.adjacency)
end;

@testset "Testing all_slope_decreasing_sequences" begin
    Q = GeneralizedKroneckerQuiver(3)
    d = [2,3]
    theta = [3,-2]
    @test all_slope_decreasing_sequences(Q, d, theta) == [[[2, 3]],
    [[1, 1], [1, 2]],
    [[2, 2], [0, 1]],
    [[2, 1], [0, 2]],
    [[1, 0], [1, 3]],
    [[1, 0], [1, 2], [0, 1]],
    [[1, 0], [1, 1], [0, 2]],
    [[2, 0], [0, 3]]]
end;

@testset "Testing has_semistable_representation" begin
    Q = GeneralizedKroneckerQuiver(17)
    @test has_semistable_representation(Q,[1,13], [13,-1]) == true
    @test has_semistable_representation(Q,[1,13], [-13,1]) == false
    
    K = GeneralizedKroneckerQuiver(4)
    @test has_semistable_representation(K,[1,5], [5,-1]) == false
    @test has_semistable_representation(K,[1,5], [-5,1]) == false

    A2 = GeneralizedKroneckerQuiver(1)
    theta = [1,-1]
    d = [1,1]
    @test has_semistable_representation(A2, d, theta) == true

    d = [2,2]
    @test has_semistable_representation(A2, d, theta) == true

    d = [1,2]
    @test has_semistable_representation(A2, d, theta) == false

    d = [0,0]
    @test has_semistable_representation(A2, d, theta) == true


    K3 = GeneralizedKroneckerQuiver(3)
    theta = [3,-2]
    d = [2,3]
    @test has_semistable_representation(K3, d, theta) == true

    d = [1,4]
    @test has_semistable_representation(K3, d, theta) == false
end;

@testset "Testing all_harder_narasimhan_types" begin
    Q = GeneralizedKroneckerQuiver(3)
    d = [2,3]
    
    theta = [3,-2]
    @test all_harder_narasimhan_types(Q, d, theta) == [[[2, 3]],
    [[1, 1], [1, 2]],
    [[2, 2], [0, 1]],
    [[2, 1], [0, 2]],
    [[1, 0], [1, 3]],
    [[1, 0], [1, 2], [0, 1]],
    [[1, 0], [1, 1], [0, 2]],
    [[2, 0], [0, 3]]]

    theta = [-3,2]
    @test all_harder_narasimhan_types(Q, d, theta) == [[[0,3],[2,0]]]
end;

@testset "Testing is_amply_stable" begin
    Q = GeneralizedKroneckerQuiver(3)
    d = [2,3]

    @test is_amply_stable(Q, d, [3,-2]) == true
    @test is_amply_stable(Q, d, [-3,2]) == false
end;

end;