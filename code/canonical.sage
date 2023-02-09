"""
input: Q quiver (in reality only the Euler form is needed?)
       dimension vector alpha
output: canonical decomposition
        ordered list of pairs of multiplicity r_i and dimension vector alpha_i

P1: r_i\geq 0
P2: alpha_i a Schur root
P3: alpha_i left orthogonal to alpha_j
P4: r_i=1 when alpha_i imaginary and not isotropic
P5: for all i<j: < alpha_j,alpha_i > >= 0

# how to check whether something is a Schur root
- use Proposition 11.2.7 (due to Schofield) to check for existence of a stable of that dimension vector for the canonical stability condition
    this check is something due to King's original paper
- by Lemma 11.2.8 (2): we can assume indivisibility to check Schur root?


# how to check left orthogonality
    observe that you can _start_ with the exceptional collection given by simples, and thus you have P3
    then the operations somehow preserve the left orthogonality?

    no, you need to apply Lemma 11.9.10, where left orthogonality is important!

    see Option 2 in code/generic-ext.sage


Comments:
    - original Schofield paper assumes Q acyclic in Section 4, but not in the other sections?
    - Derksen--Weyman in Compositio assume Q acyclic
"""

class Quiver:
    def __init__(self, M):
        # TODO some check on the matrix
        # for now M is the adjacency matrix
        self.__M = matrix(M)

    def number_of_vertices(self):
        return self.__M.nrows()

    def euler_matrix(self):
        return identity_matrix(ZZ, self.number_of_vertices()) - self.__M

    def euler_form(self, a, b):
        return a * self.euler_matrix() * b

    # taken from code/schur.sage
    # still need testing code from there
    def is_generic_subdimension_vector(self, e, d):
        # using notation from Section 5 of https://arxiv.org/pdf/0802.2147.pdf

        #possible optimisation?
        #if e == d: return True

        # list of all dimension vectors e' which are strictly smaller than e
        subdimensions = cartesian_product([range(ei + 1) for ei in e])
        subdimensions = map(vector, subdimensions)
        # for the recursion we ignore e
        subdimensions = filter(lambda eprime: eprime != e, subdimensions)
        # check whether they are generic subdimension vectors of e
        subdimensions = filter(lambda eprime: self.is_generic_subdimension_vector(eprime, e), subdimensions)
        # add e back into the list
        subdimensions = list(subdimensions) + [e]

        # apply the numerical criterion
        return all(map(lambda eprime: self.euler_form(eprime, d - e) >= 0, subdimensions))

    def generic_ext_vanishing(self, a, b):
        return self.is_generic_subdimension_vector(a, a+b)


def kronecker_quiver(n):
    return Quiver([[0, n], [0, 0]])


Q = kronecker_quiver(3)
assert Q.is_generic_subdimension_vector(vector([1,2]), vector([2,4])) == False
# the following _is_ a generic subdimension vector by the above reasoning
assert Q.is_generic_subdimension_vector(vector([1,3]), vector([2,4])) == True
