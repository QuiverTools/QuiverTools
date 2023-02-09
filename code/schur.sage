"""
How to check alpha is a Schur root?

1) Schofield, General representations of quivers, Theorem 6.1
    alpha is Schur root if and only if for all beta generic subdimension vector <beta,alpha>-<alpha,beta> > 0

2) enumerating generic subdimension vectors: Theorem 5.3 of https://arxiv.org/pdf/0802.2147.pdf is the key?
    = recursive algorithm
"""

def euler(d, e):
    M = matrix([[1, -3], [0, 1]]) # 3-Kronecker
    return d * M * e

def is_generic_subdimension_vector(e, d):
    # using notation from Section 5 of https://arxiv.org/pdf/0802.2147.pdf

    #possible optimisation?
    #if e == d: return True

    # list of all dimension vectors e' which are strictly smaller than e
    subdimensions = cartesian_product([range(ei + 1) for ei in e])
    subdimensions = map(vector, subdimensions)
    # for the recursion we ignore e
    subdimensions = filter(lambda eprime: eprime != e, subdimensions)
    # check whether they are generic subdimension vectors of e
    subdimensions = filter(lambda eprime: is_generic_subdimension_vector(eprime, e), subdimensions)
    # add e back into the list
    subdimensions = list(subdimensions) + [e]

    # apply the numerical criterion
    return all(map(lambda eprime: euler(eprime, d - e) >= 0, subdimensions))

# some tests for the 3-Kronecker quiver

# for a generic representation of dimension vector (2,4) we have 3 matrices A, B and C
# such that Av, Bv and Cv will span a 3-dimensional subspace of the 4-dimensional space
# yet for a subrepresentation of dimension vector (1,2) they will live in a 2-dimensional subspace
assert is_generic_subdimension_vector(vector([1,2]), vector([2,4])) == False
# the following _is_ a generic subdimension vector by the above reasoning
assert is_generic_subdimension_vector(vector([1,3]), vector([2,4])) == True


def generic_ext_vanishing(a, b):
    return is_generic_subdimension_vector(a, a+b)
