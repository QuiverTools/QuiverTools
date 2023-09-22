# this script should be ran in the terminal with the command "julia test.jl"

using LinearAlgebra, Test
include("Quivers.jl")


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

@testset "Testing has_semistable_representation" begin
    Q = GeneralizedKroneckerQuiver(17)
    @test has_semistable_representation(Q,[1,13], [13,-1]) == true
    @test has_semistable_representation(Q,[1,13], [-13,1]) == false
    
    K = GeneralizedKroneckerQuiver(4)
    @test has_semistable_representation(K,[1,5], [5,-1]) == false
    @test has_semistable_representation(K,[1,5], [-5,1]) == false

end;