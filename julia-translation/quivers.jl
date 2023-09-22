# This code is a translation of the Sage code in quivers.py. It is not necessarily meant to be published, but rather to help with high-performance hypothesis testing and profiling.

using LinearAlgebra, IterTools


"""
A quiver is represented by its adjacency matrix (a_ij) in M_{n x n}(N) where Q_0 = {1,...,n} and a_{ij} is the number of arrows i --> j.

Variables:
adjacency::Matrix{Int64} is the adjacency matrix of the quiver
name = None
"""
mutable struct Quiver
    adjacency::Matrix{Int64}
    name::String

    
    function Quiver(adjacency::Matrix{Int64}, name::String)
        if !(size(adjacency)[1] == size(adjacency)[1])
            throw(DomainError(adjacency, "adjacency matrix must be square"))
        elseif !(all(a >= 0 for a in adjacency))
            throw(DomainError(adjacency, "adjacency matrix must have non-negative entries"))
        else 
            new(adjacency, name)
        end
    end 
    function Quiver(adjacency::Matrix{Int64})
        if !(size(adjacency)[1] == size(adjacency)[2])
            throw(DomainError(adjacency, "adjacency matrix must be square"))
        elseif !(all(a >= 0 for a in adjacency))
            throw(DomainError(adjacency, "adjacency matrix must have non-negative entries"))
        else
            new(adjacency, "")
        end
    end
end

"""
Returns the adjacency matrix of the quiver.

OUTPUT: A square matrix M whose entry M[i,j] is the number of arrows from the vertex i to the vertex j.
"""
function adjacency_matrix(Q::Quiver)
    return Q.adjacency
end

"""
Returns the (necessarily symmetric) adjacency matrix of the underlying graph of the quiver.

OUTPUT: A square, symmetric matrix M whose entry M[i,j] = M[j,i] is the number of edges between the vertices i and j.
"""
function underlying_graph(Q::Quiver)
    return Matrix{Int64}(Q.adjacency + transpose(Q.adjacency) - diagm(diag(Q.adjacency)))
end
"""
Returns the number of vertices of the quiver.
"""
number_of_vertices(Q::Quiver) = size(Q.adjacency)[1]

"""
Returns the number of arrows of the quiver.
"""
number_of_arrows(Q::Quiver) = sum(Q.adjacency)

"""
Returns true if the quiver is acyclic, false otherwise.
"""
is_acyclic(Q::Quiver) = (Q.adjacency^number_of_vertices(Q) == zeros(Int64, number_of_vertices(Q), number_of_vertices(Q)))

"""
Returns true if the quiver is connected, false otherwise.

Examples:
    julia> Q = Quiver([0 1 0; 0 0 1; 1 0 0])
    julia> is_connected(Q)
    true

    julia> Q = Quiver([0 1 0; 1 0 0; 0 0 2])
    false

    The 4-Kronecker quiver:
    julia> Q = GeneralizedKroneckerQuiver(4)
    julia> is_connected(Q)
    true

    The 4-loop quiver:
    julia> Q = LoopQuiver(4)
    julia> is_connected(Q)
    true

    The 4-subspace quiver:
    julia> Q = SubspaceQuiver(4)
    julia> is_connected(Q)
    true

    The A10 quiver:

    julia> A10 = Quiver(   [0 1 0 0 0 0 0 0 0 0;
                            0 0 1 0 0 0 0 0 0 0;
                            0 0 0 1 0 0 0 0 0 0;
                            0 0 0 0 1 0 0 0 0 0;
                            0 0 0 0 0 1 0 0 0 0;
                            0 0 0 0 0 0 1 0 0 0;
                            0 0 0 0 0 0 0 1 0 0;
                            0 0 0 0 0 0 0 0 1 0;
                            0 0 0 0 0 0 0 0 0 1;
                            0 0 0 0 0 0 0 0 0 0] )
    julia> is_connected(A10)
    true

    The A10 quiver without one arrow:
    julia> A10 = Quiver(   [0 1 0 0 0 0 0 0 0 0;
                            0 0 1 0 0 0 0 0 0 0;
                            0 0 0 1 0 0 0 0 0 0;
                            0 0 0 0 1 0 0 0 0 0;
                            0 0 0 0 0 1 0 0 0 0;
                            0 0 0 0 0 0 0 0 0 0;
                            0 0 0 0 0 0 0 1 0 0;
                            0 0 0 0 0 0 0 0 1 0;
                            0 0 0 0 0 0 0 0 0 1;
                            0 0 0 0 0 0 0 0 0 0] )
    julia> is_connected(A10)
    false
    
"""
function is_connected(Q::Quiver)
    paths = underlying_graph(Q)
    for i in 2:number_of_vertices(Q) - 1
        paths += paths*underlying_graph(Q)
    end
    for i in 1:number_of_vertices(Q), j in 1:number_of_vertices(Q)
            if i != j && paths[i,j] == 0 && paths[j,i] == 0
                return false
            end
    end
    return true
end

# the docstrings on these functions are from the file quiver.py

"""
Returns the number of incoming arrows to the vertex j.

Examples:

julia> Q = GeneralizedKroneckerQuiver(4)
julia> indegree(Q, 1)
0
julia> indegree(Q, 2)
4
"""
indegree(Q::Quiver, j::Int64) = (1 <= j & j<= number_of_vertices(Q)) ? sum(adjacency_matrix(Q)[:,j]) : throw(DomainError(j, "vertex index out of bounds"))

"""
Returns the number of outgoing arrows from the vertex i.

Examples:

julia> Q = GeneralizedKroneckerQuiver(4)
julia> outdegree(Q, 1)
4
julia> outdegree(Q, 2)
0
"""
outdegree(Q::Quiver, i::Int64) = (1 <= i & i<= number_of_vertices(Q)) ? sum(adjacency_matrix(Q)[i,:]) : throw(DomainError(i, "vertex index out of bounds"))

"""
Returns true if the vertex i is a source, i.e. there are no incoming arrows into i, false otherwise.

Examples:

julia> Q = GeneralizedKroneckerQuiver(4)
julia> is_source(Q, 1)
true
julia> is_source(Q, 2)
false
"""
is_source(Q::Quiver, i::Int64) = (1 <= i & i<= number_of_vertices(Q)) ? indegree(Q, i) == 0 : throw(DomainError(i, "vertex index out of bounds"))

"""
Returns true if the vertex j is a sink, i.e. there are no outgoing arrows from j, false otherwise.

Examples:

julia> Q = GeneralizedKroneckerQuiver(4)
julia> is_sink(Q, 1)
false
julia> is_sink(Q, 2)
true
"""
is_sink(Q::Quiver, j::Int64) = (1 <= j & j<= number_of_vertices(Q)) ? outdegree(Q, j) == 0 : throw(DomainError(j, "vertex index out of bounds"))

"""
Returns the Euler matrix of the quiver.
"""
euler_matrix(Q::Quiver) = Matrix{Int64}(I, number_of_vertices(Q), number_of_vertices(Q)) - adjacency_matrix(Q)

"""
Returns the value of the Euler bilinear form of the quiver computed on the vectors x and y.
"""
euler_form(Q::Quiver, x::Vector{Int64}, y::Vector{Int64}) = (length(x) == number_of_vertices(Q) & length(y) == number_of_vertices(Q)) ? x'*euler_matrix(Q)*y : throw(DomainError("dimension vectors must have length equal to number of vertices"))

"""
The opposite quiver is given by the transpose of the adjacency matrix of the original quiver.

Returns a Quiver object with the same vertices and an arrow from j to i for every arrow from i to j in the original quiver.
"""
opposite_quiver(Q::Quiver) = Quiver(Matrix{Int64}(transpose(adjacency_matrix(Q))), "Opposite of "*Q.name)

"""
The adjacency matrix of the double of a quiver is the sum of the adjacency matrix of the original quiver and its transpose.
"""
double_quiver(Q::Quiver) = Quiver(adjacency_matrix(Q) + Matrix{Int64}(transpose(adjacency_matrix(Q))), "Double of "*Q.name)


## Everything that comes before this has been properly translated from the Sage code and should work.

thin_dimension_vectors(Q::Quiver) = [1 for i in 1:number_of_vertices(Q)]

"""
The canonical stability parameter is given by <d,_> - <_,d>
"""
canonical_stability_parameter(Q::Quiver, d::Vector{Int64}) = d*(-transpose(euler_matrix(Q)) + euler_matrix(Q))


"""
Returns the list of all sequences (d^1,...,d^l) which sum to d such that slope(d^1) > ... > slope(d^l)

Examples:

julia> Q = GeneralizedKroneckerQuiver(3)
julia> d = [2,3]
julia> theta = [3,-2]
julia> all_slope_decreasing_sequences(Q, d, theta)
8-element Array{Array{Vector}}:
    [[[2, 3]],
    [[1, 1], [1, 2]],
    [[2, 2], [0, 1]],
    [[2, 1], [0, 2]],
    [[1, 0], [1, 3]],
    [[1, 0], [1, 2], [0, 1]],
    [[1, 0], [1, 1], [0, 2]],
    [[2, 0], [0, 3]]]
    """
function all_slope_decreasing_sequences(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64}, denominator::Function = sum)

    # List all subdimension vectors e of bigger slope than d.
    subdimensions = filter(e -> (e != ZeroVector(number_of_vertices(Q))) && (slope(e,theta,denominator) > slope(d,theta,denominator)), all_subdimension_vectors(d))

    # We sort the subdimension vectors by slope because that will return the list of all HN types in ascending order with respect to the partial order from Def. 3.6 of https://mathscinet.ams.org/mathscinet-getitem?mr=1974891
    subdimensions = sort(subdimensions, by = e -> slope(e,theta,denominator))
    # The slope decreasing sequences which are not of the form (d) are given by (e,f^1,...,f^s) where e is a proper subdimension vector such that mu_theta(e) > mu_theta(d) and (f^1,...,f^s) is a HN type of f = d-e such that mu_theta(e) > mu_theta(f^1) holds.

    # I will rewrite this as functional programming later
    allSlopeDecreasing = []
    for e in subdimensions
        for fstar in filter(fstar -> slope(e,theta) > slope(fstar[1],theta), all_slope_decreasing_sequences(Q, d-e, theta, denominator))
        push!(allSlopeDecreasing, [e, fstar...])
        end
    end
    # Add d again, at the beginning, because it is smallest with respect to the partial order from Def. 3.6
    return [[d], allSlopeDecreasing...]
end


"""Checks if there is a theta-semistable representation of dimension vector d.

Examples:


julia> A2 = GeneralizedKroneckerQuiver(1)
julia> theta = [1,-1]
julia> d = [1,1]
julia> has_semistable_representation(A2, d, theta)
true

julia> d = [2,2]
julia> has_semistable_representation(A2, d, theta)
true

julia> d = [1,2]
julia> has_semistable_representation(A2, d, theta)
false

julia> d = [0,0]
julia> has_semistable_representation(A2, d, theta)
true

The 3-Kronecker quiver:
julia> K3 = GeneralizedKroneckerQuiver(3)
julia> theta = [3,-2]
julia> d = [2,3]
julia> has_semistable_representation(K3, d, theta)
true

julia> d = [1,4]
julia> has_semistable_representation(K3, d, theta)
false
"""
function has_semistable_representation(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64}, algorithm::String = "schofield") 
    if algorithm == "reineke"
        throw(ArgumentError("reineke algorithm not implemented"))

    elseif algorithm == "schofield"
        # collect the list of all subdimension vectors e of bigger slope than d
        subdimensionsBiggerSlope = filter(e -> e != ZeroVector(number_of_vertices(Q)) && e != d && slope(e, theta) > slope(d, theta), all_subdimension_vectors(d))
        # to have semistable representations, none of the vectors above must be generic subdimension vectors.
        return !any(e -> is_generic_subdimension_vector(Q, e, d), subdimensionsBiggerSlope)
    else
        throw(ArgumentError("algorithm not recognized"))
    end
end

function has_stable_representation(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64}, algorithm::String = "schofield")

    if algorithm == "al"
        throw(ArgumentError("al algorithm not implemented"))
    elseif algorithm == "schofield"
        if d == ZeroVector(number_of_vertices(Q))
            return false
        else
            subdimensionsSlopeNoLess = filter(e -> e != zeroVector && e != d && slope(e, theta) >= slope(d, theta), all_subdimension_vectors(d))
            return !any(e -> is_generic_subdimension_vector(Q, e, d), subdimensionsSlopeNoLess)
        end
    else
        throw(ArgumentError("algorithm not recognized"))
    end
end


is_schur_root(Q::Quiver, d::Vector{Int64}) = has_stable_representation(Q, d, canonical_stability_parameter(Q, d))

"""Checks if e is a generic subdimension vector of d.
        # using notation from Section 5 of https://arxiv.org/pdf/0802.2147.pdf
    A dimension vector e is called a generic subdimension vector of d if a generic representation of dimension vector d possesses a subrepresentation of dimension vector e.
    By a result of Schofield (see Thm. 5.3 of https://arxiv.org/pdf/0802.2147.pdf) e is a generic subdimension vector of d if and only if <e',d-e> is non-negative for all generic subdimension vectors e' of e."""
function is_generic_subdimension_vector(Q::Quiver, e::Vector{Int64}, d::Vector{Int64})
    if (e == ZeroVector(number_of_vertices(Q))) | (e == d)
        return true
    end
    # considering subdimension vectors that violate the numerical condition
    subdimensions = filter(eprime -> euler_form(Q,eprime, d-e) < 0, all_subdimension_vectors(e))
    # none of the subdimension vectors violating the condition should be generic
    return !any(eprime -> is_generic_subdimension_vector(Q,eprime,e), subdimensions)
end

all_generic_subdimension_vectors(Q::Quiver, d::Vector{Int64}) = filter(e -> is_generic_subdimension_vector(Q, e, d), all_subdimension_vectors(d))


generic_ext_vanishing(Q::Quiver, a::Vector{Int64}, b::Vector{Int64}) = is_generic_subdimension_vector(Q, a, a+b)

function canonical_decomposition(Q::Quiver, d::Vector{Int64}, algorithm::String = "derksen-weyman")
    if algorithm == "derksen-weyman"
        throw(ArgumentError("derksen-weyman algorithm not implemented"))
    elseif algorithm == "schofield-1"
        throw(ArgumentError("schofield-1 algorithm not implemented"))
    elseif algorithm == "schofield-2"
        throw(ArgumentError("schofield-2 algorithm not implemented"))
    else
        throw(ArgumentError("algorithm not recognized"))
    end
end


function is_harder_narasimhan_type(Q::Quiver, dstar::Vector{Vector{Int64}}, theta::Vector{Int64})

    if length(dstar) == 1
        return has_semistable_representation(Q, dstar[1], theta)
    else
        for i in 1:length(dstar)-1
            if slope(dstar[i]) <= slope(dstar[i+1])
                return false
            end
        end
        return all(e -> has_semistable_representation(Q, e, theta), dstar)
    end    
end

function all_harder_narasimhan_types(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64})
    throw(ArgumentError("not implemented"))
end

function is_luna_type(Q::Quiver, tau::Vector{Tuple{Vector{Int64},Int64}}, theta::Vector{Int64})
    n = number_of_vertices(Q)
    zeroVector = Vector{Int64}(zeros(Int64, n))
    d = sum(sum(tupl[1] .* tupl[2]) for tupl in tau) 
    if d == zeroVector
        return tau == [(zeroVector, 1)]
    else
        dstar = [tupl[1] for tupl in tau]
        equalSlope = all(e -> slope(e, theta, denominator=sum) == slope(d, theta, denominator=sum), dstar)
        semistable = all(e -> has_stable_representation(Q, e, theta, algorithm="schofield"), dstar)
        return equalSlope && semistable
    end
end

function all_luna_types(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64})
    throw(ArgumentError("not implemented"))
end

function semistable_equals_stable(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64}, algorithm::String = "schofield")
    throw(ArgumentError("not implemented"))
end

slope(d::Vector{Int64}, theta::Vector{Int64}, denominator::Function = sum) = (length(d) == length(theta) && denominator(d)>0) ? (theta'*d)/denominator(d) : throw(DomainError("dimension vector and stability parameter must have same length"))


function all_forbidden_subdimension_vectors(d::Vector{Int64}, theta::Vector{Int64})
    zeroVector = Vector{Int64}(zeros(Int64, length(d)))
    properSubdimensions = filter(e -> e != d && e != zeroVector, all_subdimension_vectors(d))
    return filter(e -> slope(e, theta) > slope(d, theta), properSubdimensions)
end

function is_coprime_for_stability_parameter(d::Vector{Int64}, theta::Vector{Int64})
    zeroVector = Vector{Int64}(zeros(Int64, length(d)))
    properSubdimensions = filter(e -> e != d && e != zeroVector, all_subdimension_vectors(d))
    return all(e -> slope(d, theta) != slope(e, theta), properSubdimensions)
end

function is_indivisible(d::Vector{Int64})
    return gcd(d) == 1
end

function GeneralizedKroneckerQuiver(m::Int64)
    return Quiver([0 m; 0 0], string(m)*"-Kronecker quiver")
end

KroneckerQuiver() = GeneralizedKroneckerQuiver(2)

ThreeVertexQuiver(m12::Int64, m13::Int64, m23::Int64) = Quiver([0 m12 m13; 0 0 m23; 0 0 0], "An acyclic 3-vertex quiver")

LoopQuiver(m::Int64) = Quiver([m], string(m)*"-loop quiver")

function SubspaceQuiver(m::Int64)
    A = zeros(Int64, m+1, m+1)
    for i in 1:m
        A[i, m+1] = 1
    end
    return Quiver(A, string(m)*"-subspace quiver")
end

DynkinQuiver(T::String) = throw(ArgumentError("not implemented"))
ExtendedDynkinQuiver(T::String) = throw(ArgumentError("not implemented"))
CyclicQuiver(n::Int64) = throw(ArgumentError("not implemented"))
BipartiteQuiver(m::Int64, n::Int64) = throw(ArgumentError("not implemented"))

ZeroVector(n::Int64) = Vector{Int64}(zeros(Int64, n))


all_subdimension_vectors(d::Vector{Int64}) = collect(collect.(product((0:di for di in d)...)))