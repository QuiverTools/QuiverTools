# This code is a translation of the Sage code in quivers.py. It is not meant to be published, but rather to help with high-performance hypothesis testing and profiling.

using LinearAlgebra, IterTools



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

function adjacency_matrix(Q::Quiver)
    return Q.adjacency
end

function underlying_graph(Q::Quiver)
    return Matrix{Int64}(Q.adjacency + transpose(Q.adjacency) - diagm(diag(Q.adjacency)))
end

number_of_vertices(Q::Quiver) = size(Q.adjacency)[1]

number_of_arrows(Q::Quiver) = sum(Q.adjacency)

is_acyclic(Q::Quiver) = (Q.adjacency^number_of_vertices(Q) == zeros(Int64, number_of_vertices(Q), number_of_vertices(Q)))

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

indegree(Q::Quiver, j::Int64) = (1 <= j & j<= number_of_vertices(Q)) ? sum(adjacency_matrix(Q)[:,j]) : throw(DomainError(j, "vertex index out of bounds"))

outdegree(Q::Quiver, i::Int64) = (1 <= i & i<= number_of_vertices(Q)) ? sum(adjacency_matrix(Q)[i,:]) : throw(DomainError(i, "vertex index out of bounds"))

is_source(Q::Quiver, i::Int64) = (1 <= i & i<= number_of_vertices(Q)) ? indegree(Q, i) == 0 : throw(DomainError(i, "vertex index out of bounds"))

is_sink(Q::Quiver, j::Int64) = (1 <= j & j<= number_of_vertices(Q)) ? outdegree(Q, j) == 0 : throw(DomainError(j, "vertex index out of bounds"))

euler_matrix(Q::Quiver) = Matrix{Int64}(I, number_of_vertices(Q), number_of_vertices(Q)) - adjacency_matrix(Q)

euler_form(Q::Quiver, x::Vector{Int64}, y::Vector{Int64}) = (length(x) == number_of_vertices(Q) & length(y) == number_of_vertices(Q)) ? x'*euler_matrix(Q)*y : throw(DomainError("dimension vectors must have length equal to number of vertices"))

opposite_quiver(Q::Quiver) = Quiver(Matrix{Int64}(transpose(adjacency_matrix(Q))), "Opposite of "*Q.name)

double_quiver(Q::Quiver) = Quiver(adjacency_matrix(Q) + Matrix{Int64}(transpose(adjacency_matrix(Q))), "Double of "*Q.name)


## Everything that comes before this has been properly translated from the Sage code and should work.

thin_dimension_vectors(Q::Quiver) = [1 for i in 1:number_of_vertices(Q)]

canonical_stability_parameter(Q::Quiver, d::Vector{Int64}) = d*(-transpose(euler_matrix(Q)) + euler_matrix(Q))

function has_semistable_representation(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64}, algorithm::String = "schofield") 
    if algorithm == "reineke"
        throw(ArgumentError("reineke algorithm not implemented"))
    elseif algorithm == "schofield"

        subdimensionsBiggerSlope = filter(e -> e != ZeroVector(number_of_vertices(Q)) && e != d && slope(e, theta) > slope(d, theta), all_subdimension_vectors(d))
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