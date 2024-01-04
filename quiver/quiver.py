# from sage.matrix.constructor import matrix
from sage.all import *

class Quiver:

    """
    A quiver is represented by its adjacency matrix (a_ij) in M_{n x n}(N) where Q_0 = {1,...,n} and a_{ij} is the number of arrows i --> j.

    Variables:
    adjacency
    name = None
    """

    def __init__(self, M, name=None):
        # TODO should we raise an exception/error instead?
        assert M.is_square()
        assert all(a >= 0 for a in M.list())

        self._adjacency = M
        self._name = name

    def __repr__(self):
        # TODO this should be implemented following Sage's methodology
        # see https://github.com/pbelmans/hodge-diamond-cutter/issues/14 for suggestion
        # and https://github.com/pbelmans/hodge-diamond-cutter/commit/59cc6d575babe695c6e1668721e6cd5c4f17dba9 for example
        output = ""
        if (self._name == None):
            output += "A quiver with "
        else:
            output += str(self._name)+"; "
        output += "adjacency matrix:\n"+str(self._adjacency)
        return output

    """
    Basic graph-theoretic properties of the quiver
    """

    def adjacency_matrix(self):
        r"""Returns the adjacency matrix of the quiver.

        OUTPUT: A square matrix M whose entry M[i,j] is the number of arrows from the vertex i to the vertex j.
        """
        return self._adjacency

    def underlying_graph(self):
        r"""Returns the (necessarily symmetric) adjacency matrix of the underlying graph of the quiver.

        OUTPUT: A square, symmetric matrix M whose entry M[i,j] = M[j,i] is the number of edges between the vertices i and j.
        """
        return self.adjacency_matrix() + self.adjacency_matrix().transpose() - diagonal_matrix(self.adjacency_matrix().diagonal())

    def number_of_vertices(self):
        r""""Returns the amount of vertices that the quiver has.

        OUTPUT: The number of vertices as an Int.
        """
        return self.adjacency_matrix().nrows()

    def number_of_arrows(self):
        r""""Returns the number of arrows that the quiver has.

        OUTPUT: The number of arrows as an Int.

        """
        thin = self.thin_dimension_vector()
        return thin * self.adjacency_matrix() * thin

    def is_acyclic(self):
        r""""Returns the truth value of wether the quiver is acyclic.

        OUTPUT: Statement truth value as Bool.
        """
        A = self.adjacency_matrix()
        n = self.number_of_vertices()

        # a quiver is acyclic if and only if its adjacency matrix is nilpotent
        return (A**n == zero_matrix(ZZ, n))

    def is_connected(self):
        r""""Returns whether the underlying graph of the quiver is connected or not.

        OUTPUT: Statement truth value as Bool.

        EXAMPLES:

        The 4-Kronecker quiver::

            sage: from quiver import *
            sage: K = Quiver( matrix(  [[0, 4],
            ....:                       [0, 0]]))
            sage: K.is_connected()
            True

        The doubled 1-Kronecker quiver::

            sage: from quiver import *
            sage: C1 = Quiver(matrix(  [[0,1],
            ....:                       [1,0]]))
            sage: C1.is_connected()
            True

        The 3-loop point quiver::

            sage: from quiver import *
            sage: L = Quiver(matrix([[3]]))
            sage: L.is_connected()
            True

        The A_10 quiver::

            sage: from quiver import *
            sage: A10 = Quiver(matrix(     [[0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            ....:                           [0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
            ....:                           [0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
            ....:                           [0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            ....:                           [0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
            ....:                           [0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            ....:                           [0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
            ....:                           [0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
            ....:                           [0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
            ....:                           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]))
            sage: A10.is_connected()
            True

        The A_10 quiver without one arrow::

            sage: from quiver import *
            sage: discA10 = Quiver(matrix(     [[0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            ....:                               [0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
            ....:                               [0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
            ....:                               [0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            ....:                               [0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
            ....:                               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            ....:                               [0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
            ....:                               [0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
            ....:                               [0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
            ....:                               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]))
            sage: discA10.is_connected()
            False
        """
        # inefficient but functioning method. To improve?
        # more efficient algorithm: scan the graph in depth and list all reachable vertices.
        paths = self.underlying_graph()
        for i in range(2,self.number_of_vertices()): # -1 ?
            # add all paths of length i
            paths += paths*self.underlying_graph()
        # if every couple of vertices is connected (ij \neq 0 or ji \neq 0) then true, otherwise false.
        for i in range(self.number_of_vertices()):
            for j in range(self.number_of_vertices()):
                if i != j and paths[i,j] == 0 and paths[j,i] == 0:
                    return False
        return True

    """
    Some graph-theoretic properties of the quiver
    """

    def indegree(self,j):
        """The indegree of j is the number of incoming arrows into j."""

        """indeg(j) = sum_i a_{ij} where (a_{ij}) is the adjacency matrix."""
        # Question: Should we number the vertices 1,...,n or 0,...,n-1?

        """
        EXAMPLES

        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: Q.indegree(1)
        0
        sage: Q.indegree(2)
        3
        """

        assert (j > 0) and (j <= self.number_of_vertices())
        return sum(self._adjacency.column(j-1))


    def outdegree(self,i):
        """The outdegree of i is the number of incoming arrows into i."""

        """outdeg(i) = sum_j a_{ij} where (a_{ij}) is the adjacency matrix."""
        
        """
        EXAMPLES

        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: Q.outdegree(1)
        3
        sage: Q.outdegree(2)
        0
        """

        assert (i > 0) and (i <= self.number_of_vertices())
        return sum(self._adjacency.row(i-1))

    def is_source(self,i):
        """Checks if i is a source of the quiver, i.e. if there are no incoming arrows into i."""

        """
        EXAMPLES

        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: Q.is_source(1)
        True
        sage: Q.is_source(2)
        False
        """

        assert (i > 0) and (i <= self.number_of_vertices())
        return (self.indegree(i) == 0)

    def is_sink(self,j):
        """Checks if j is a sink of the quiver, i.e. if there are no outgoing arrows out of j."""

        """
        EXAMPLES

        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: Q.is_sink(1)
        False
        sage: Q.is_sink(2)
        True
        """

        assert (j > 0) and (j <= self.number_of_vertices())
        return (self.outdegree(j) == 0)


    """
    Basic representation-theoretical properties of the quiver
    """

    def euler_matrix(self):
        r"""Returns the Euler matrix of the quiver.

        OUTPUT: Sage matrix.
        """
        return matrix.identity(self.number_of_vertices()) - self.adjacency_matrix()

    def euler_form(self, x, y):
        r"""The Euler bilinear form of the quiver.

        INPUT:
        - ``x`` -- vector of integers
        - ``y`` -- vector of integers

        OUTPUT: the multiplication of ``x * self.adjacency_matrix() * y`` as an  Int.

        """
        assert (x.length() == self.number_of_vertices() and y.length() == self.number_of_vertices())
        return x * self.euler_matrix() * y

    """
    Constructing new quivers out of old
    """

    def opposite_quiver(self):
        """The opposite quiver is given by the transpose of the adjacency matrix of the original quiver.

        OUTPUT: a Quiver object the same vertices and an arrow from j to i for every arrow from i to j in the original quiver.
        """
        A = self.adjacency_matrix().transpose()

        if (self._name != None):
            name = "Opposite of " + self._name
        else:
            name = None

        return Quiver(A, name)

    def double_quiver(self):
        """The adjacency matrix of the double of a quiver is the sum of the adjacency matrix of the original quiver and its transpose."""
        A = self.adjacency_matrix() + self.adjacency_matrix().transpose()
        if (self._name != None):
            name = "Double of " + self._name
        else:
            name = None
        return Quiver(A, name)


    """
    Dimension vectors and stability conditions
    """

    def thin_dimension_vector(self):
        return vector([1 for i in range(self.number_of_vertices())])
    
    def simple_root(self, i):
        """Returns the simple root e_i = [0,...,1,...,0], i.e. the unit vector with a one in position i."""
        n = self.number_of_vertices()
        # Our convention is that vertices are numbered 1,...,n
        assert (i >= 1 and i <= n)
        ei = vector([0 for i in range(n)])
        ei[i-1] = 1
        return ei
    
    def support(self, d):
        """Returns the full subquiver supported on {i in Q_0 | d_i > 0}."""
        """
        EXAMPLES
        
        sage: from quiver import *
        Q = ThreeVertexQuiver(2,0,4)
        sage: d = vector([1,1,1])
        sage: Q.support(d)
        A quiver with adjacency matrix:
        [0 2 0]
        [0 0 4]
        [0 0 0]
        sage: d = vector([1,0,1])
        sage: Q.support(d)
        A quiver with adjacency matrix:
        [0 0]
        [0 0]
        sage: d = vector([1,1,0])
        sage: Q.support(d)
        A quiver with adjacency matrix:
        [0 2]
        [0 0]
        """
        n = self.number_of_vertices()
        A = self.adjacency_matrix()
        support = list(filter(lambda i: d[i] > 0, range(n)))
        # Submatrix (A_ij)_{i,j in supp(d)} is the adjacency matrix of the sought quiver
        ASupp = matrix([[A[i,j] for j in support] for i in support])
        return Quiver(ASupp)

    def in_fundamental_domain(self, d):
        """Checks if the dimension vector d is in the fundamental domain."""
        
        """The fundamental domain of Q is the set of dimension vectors d such that supp(d) is connected and <d,e_i> + <e_i,d> <= 0 for all simple roots e_i.
        Every d in the fundamental domain is an imaginary root and the set of imaginary roots is the Weyl group saturation of the fundamental domain. 
        If d is in the fundamental domain then it is Schurian and a general representation of dimension vector d is stable for the canonical stability parameter."""

        """
        EXAMPLES

        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: d = vector([1,1])
        sage: Q.in_fundamental_domain(d)
        True
        sage: d = vector([1,2])
        sage: Q.in_fundamental_domain(d)
        False
        sage: d = vector([2,3])
        sage: Q.in_fundamental_domain(d)
        True
        """

        n = self.number_of_vertices()
        # This is the condition <d,e_i> + <e_i,d> <= 0 for all i in Q_0
        eulerFormCondition = all([(self.euler_form(d,self.simple_root(i+1)) + self.euler_form(self.simple_root(i+1),d) <= 0) for i in range(n)])
        # Check if the support is connected
        connected = self.support(d).is_connected()
        return eulerFormCondition and connected

    def canonical_stability_parameter(self,d):
        """The canonical stability parameter is given by <d,_> - <_,d>"""
        E = self.euler_matrix()
        return d * (-self.euler_matrix().transpose() + E)

    def all_slope_decreasing_sequences(self, d, theta, denominator=sum):
        """Returns the list of all sequences (d^1,...,d^l) which sum to d such that slope(d^1) > ... > slope(d^l)"""

        """
        EXAMPLES

        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: d = vector([2,3])
        sage: theta = vector([3,-2])
        sage: Q.all_slope_decreasing_sequences(d,theta)
        [[(2, 3)],
        [(1, 1), (1, 2)],
        [(2, 2), (0, 1)],
        [(2, 1), (0, 2)],
        [(1, 0), (1, 3)],
        [(1, 0), (1, 2), (0, 1)],
        [(1, 0), (1, 1), (0, 2)],
        [(2, 0), (0, 3)]]
        """

        n = self.number_of_vertices()
        zeroVector = vector([0 for i in range(n)])    
        # List all subdimension vectors e of bigger slope than d.
        subdimensions = list(filter(lambda e: (e != zeroVector) and (slope(e,theta,denominator=denominator) > slope(d,theta,denominator=denominator)), all_subdimension_vectors(d)))
        # We sort the subdimension vectors by slope because that will return the list of all HN types in ascending order with respect to the partial order from Def. 3.6 of https://mathscinet.ams.org/mathscinet-getitem?mr=1974891
        subdimensions.sort(key=(lambda e: slope(e,theta,denominator=denominator)))
        # The slope decreasing sequences which are not of the form (d) are given by (e,f^1,...,f^s) where e is a proper subdimension vector such that mu_theta(e) > mu_theta(d) and (f^1,...,f^s) is a HN type of f = d-e such that mu_theta(e) > mu_theta(f^1) holds.
        allSlopeDecreasing =  [[e]+fstar for e in subdimensions for fstar in list(filter(lambda fstar: slope(e,theta) > slope(fstar[0],theta) ,self.all_slope_decreasing_sequences(d-e,theta,denominator=denominator)))]
        # Add d again, at the beginning, because it is smallest with respect to the partial order from Def. 3.6
        allSlopeDecreasing = [[d]] + allSlopeDecreasing
        return allSlopeDecreasing


    def has_semistable_representation(self, d, theta, algorithm="schofield"):
        """Checks if there is a theta-semistable representation of dimension vector d."""

        #assert algorithm == "reineke"

        # TODO implement this
        if algorithm == "reineke":
            raise NotImplementedError()

        """See Thm. 5.4(1) of Reineke's overview paper https://arxiv.org/pdf/0802.2147.pdf: A dimension vector d admits a theta-semi-stable representation if and only if mu_theta(e) <= mu_theta(d) for all generic subdimension vectors e of d."""

        """
        EXAMPLES:

        The A_2 quiver:
        sage: from quiver import *
        sage: A2 = GeneralizedKroneckerQuiver(1)
        sage: theta = vector([1,-1])
        sage: d = vector([1,1])
        sage: A2.has_semistable_representation(d,theta,algorithm="schofield")
        True
        sage: d = vector([2,2])
        sage: A2.has_semistable_representation(d,theta,algorithm="schofield")
        True
        sage: d = vector([1,2])
        sage: A2.has_semistable_representation(d,theta,algorithm="schofield")
        False
        sage: d = vector([0,0])
        sage: A2.has_semistable_representation(d,theta,algorithm="schofield")
        True

        The 3-Kronecker quiver:
        sage: from quiver import *
        sage: K3 = GeneralizedKroneckerQuiver(3)
        sage: theta = vector([3,-2])
        sage: d = vector([2,3])
        sage: K3.has_semistable_representation(d,theta,algorithm="schofield")
        True
        sage: d = vector([1,4])
        sage: K3.has_semistable_representation(d,theta,algorithm="schofield")
        False

        """


        if algorithm == "schofield":        
            n = self.number_of_vertices()
            zeroVector = vector([0 for i in range(n)])    
            subdimensionsBiggerSlope = list(filter(lambda e: e != zeroVector and e != d and slope(e,theta) > slope(d,theta), all_subdimension_vectors(d)))
            return not any([self.is_generic_subdimension_vector(e,d) for e in subdimensionsBiggerSlope])


    def has_stable_representation(self, d, theta, algorithm="schofield"):
        """Checks if there is a theta-stable representation of dimension vector d."""
        # TODO implement this
        # https://mathscinet.ams.org/mathscinet-getitem?mr=1315461
        """Question: What is King's algorithm for checking for existence of stable representations supposed to be? I can't find one in the paper."""
        if algorithm == "king":
            raise NotImplementedError()
        # TODO implement this
        # al stands for Adriaenssens--Le Bruyn
        # https://mathscinet.ams.org/mathscinet-getitem?mr=1972892
        if algorithm == "al":
            raise NotImplementedError()

        """See Thm. 5.4(1) of Reineke's overview paper https://arxiv.org/pdf/0802.2147.pdf: A dimension vector d admits a theta-stable representation if and only if mu_theta(e) < mu_theta(d) for all proper generic subdimension vectors e of d."""

        """
        EXAMPLES:

        The A2 quiver:
        sage: from quiver import *
        sage: A2 = GeneralizedKroneckerQuiver(1)
        sage: theta = vector([1,-1])
        sage: d = vector([1,1])
        sage: A2.has_stable_representation(d,theta,algorithm="schofield")
        True
        sage: d = vector([2,2])
        sage: A2.has_stable_representation(d,theta,algorithm="schofield")
        False
        sage: d = vector([0,0])
        sage: A2.has_stable_representation(d,theta,algorithm="schofield")
        False

        The 3-Kronecker quiver:
        sage: from quiver import *
        sage: K3 = GeneralizedKroneckerQuiver(3)
        sage: theta = vector([3,-2])
        sage: d = vector([2,3])
        sage: K3.has_stable_representation(d,theta,algorithm="schofield")
        True

        """
        if algorithm == "schofield":
            n = self.number_of_vertices()
            zeroVector = vector([0 for i in range(n)])
            if d == zeroVector:
                return False
            else:
                subdimensionsSlopeNoLess = list(filter(lambda e: e != zeroVector and e != d and slope(e,theta) >= slope(d,theta), all_subdimension_vectors(d)))
                return not any([self.is_generic_subdimension_vector(e,d) for e in subdimensionsSlopeNoLess])


    def is_schur_root(self,d):
        """Checks if d is a Schur root for the given quiver, i.e. a dimension vector which admits a Schurian representation."""

        """By a result of Schofield (https://mathscinet.ams.org/mathscinet/relay-station?mr=1162487) d is a Schur root if and only if d admits a stable representation for the canonical stability parameter."""

        """
        EXAMPLES:

        The 3-Kronecker quiver:
        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: d = vector([2,3])
        sage: Q.is_schur_root(d)
        True

        """

        theta = self.canonical_stability_parameter(d)
        return self.has_stable_representation(d,theta)


    # TODO dimension vectors should have .is_stable(), .is_amply_stable()?
    def is_amply_stable(self, d, theta):
        """Checks if d is amply stable for theta, which by definition means that the codimension of the theta-stable locus inside R(Q,d) is at least 2."""

        # By Prop. 4.1 of https://arxiv.org/pdf/1410.0466.pdf d is amply stable for theta provided that <e,d-e> <= -2 for every proper subdimension vector.
        # But can we find a necessary and sufficient condition?
        # If every theta-semi-stable representation of dimension vector d is theta-stable then theta-ample stability is equivalent to every proper HN stratum having codimension at least 2. 

        """
        EXAMPLES

        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: d = vector([2,3])
        sage: theta = vector([3,-2])
        sage: Q.is_amply_stable(d,theta)
        True
        sage: Q.is_amply_stable(d,-theta)
        False
        """

        if self.semistable_equals_stable(d,theta):
            hn = self.all_harder_narasimhan_types(d,theta)
            if [d] in hn:
                hn.remove([d])
            return all([self.codimension_of_harder_narasimhan_stratum(dstar) >= 2 for dstar in hn])
        else:
            raise NotImplementedError()

    def is_strongly_amply_stable(self, d, theta):
        """Checks if <e,d-e> <= -2 holds for all subdimension vectors e of d which satisfy slope(e) > slope(d)."""

        # All subdimension vectors of d
        es = all_subdimension_vectors(d)
        # Remove (0,...,0)
        zeroVector = vector([0 for i in range(d.length())])
        es.remove(zeroVector)
        # All of them which have bigger slope
        es = list(filter(lambda e: slope(e,theta) > slope(d,theta), es))
        return all([self.euler_form(e,d-e) <= -2 for e in es])

    # taken from code/snippets/canonical.sage
    # TODO still need testing code from there
    def is_generic_subdimension_vector(self, e, d):
        """Checks if e is a generic subdimension vector of d."""
        # using notation from Section 5 of https://arxiv.org/pdf/0802.2147.pdf
        """A dimension vector e is called a generic subdimension vector of d if a generic representation of dimension vector d possesses a subrepresentation of dimension vector e.
        By a result of Schofield (see Thm. 5.3 of https://arxiv.org/pdf/0802.2147.pdf) e is a generic subdimension vector of d if and only if <e',d-e> is non-negative for all generic subdimension vectors e' of e."""

        # Optimization: Check first if numerical condition is violated and then if any of the e' which does violate it is a generic subdimension vector.
        
        if e == d: 
            return True
        else:
            # list of all dimension vectors e' which are strictly smaller than e
            subdimensions = all_subdimension_vectors(e)
            # only those which violate the numerical condition
            subdimensions = filter(lambda eprime: self.euler_form(eprime, d-e) < 0, subdimensions)

            # check if the list contains no generic subdimension vector of e
            return not any([self.is_generic_subdimension_vector(eprime,e) for eprime in subdimensions])

    def all_generic_subdimension_vectors(self, d):
        """Returns the list of all generic subdimension vectors of d."""
        genericSubdimensions = all_subdimension_vectors(d)
        return list(filter(lambda e: self.is_generic_subdimension_vector(e,d), genericSubdimensions))

    def generic_ext_vanishing(self, a, b):
        return self.is_generic_subdimension_vector(a, a+b)


    def canonical_decomposition(self, d, algorithm="derksen-weyman"):
        # TODO implement this
        # https://mathscinet.ams.org/mathscinet-getitem?mr=1930979
        # this is implemented in code/snippets/canonical.sage, so include it here
        if algorithm == "derksen-weyman":
            raise NotImplementedError()
        # TODO implement this
        # https://mathscinet.ams.org/mathscinet-getitem?mr=1162487
        elif algorithm == "schofield-1":
            raise NotImplementedError()
        # TODO implement this
        # https://arxiv.org/pdf/math/9911014.pdf (see Section 5, and also Section 3 of https://mathscinet.ams.org/mathscinet/article?mr=1789222)
        # in Derksen--Weyman's https://mathscinet.ams.org/mathscinet-getitem?mr=1930979 it is claimed that there is a second Schofield algorithm
        # (they do cite the wrong Schofield preprint though...)
        elif algorithm == "schofield-2":
            raise NotImplementedError()

    def is_harder_narasimhan_type(self, dstar, theta, denominator=sum, algorithm="schofield"):
        """Checks if dstar is a HN type. Peforms the check of semistability according to algorithm"""

        n = self.number_of_vertices()
        zeroVector = vector([0 for i in range(n)])
        d = sum(dstar)
        if (d == zeroVector):
            return (dstar == [zeroVector])
        else:
            slopeDecreasing = all([(slope(dstar[i],theta,denominator=denominator) > slope(dstar[i+1],theta,denominator=denominator)) for i in range(len(dstar)-1)])
            semistable = all([self.has_semistable_representation(d,theta,algorithm=algorithm) for e in dstar])
            return (slopeDecreasing and semistable)

    def codimension_of_harder_narasimhan_stratum(self,dstar):
        """Computes the codimension of the HN stratum R_{d^*}^HN inside R_d."""
        """The codimension of the HN stratum of d^* = (d^1,...,d^s) is given by - sum_{k < l} <d^k,d^l>"""

        """
        EXAMPLES

        The 3-Kronecker quiver
        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: d = vector([2,3])
        sage: theta = vector([3,-2])
        sage: hn = Q.all_harder_narasimhan_types(d,theta)
        sage: hn
        [[(2, 3)],
        [(1, 1), (1, 2)],
        [(2, 2), (0, 1)],
        [(2, 1), (0, 2)],
        [(1, 0), (1, 3)],
        [(1, 0), (1, 2), (0, 1)],
        [(1, 0), (1, 1), (0, 2)],
        [(2, 0), (0, 3)]]
        sage: [Q.codimension_of_harder_narasimhan_stratum(dstar) for dstar in hn]
        [0, 3, 4, 10, 8, 9, 12, 18]
        """

        s = len(dstar)
        return -sum([self.euler_form(dstar[k],dstar[l]) for k in range(s-1) for l in range(k+1,s)])


    def all_harder_narasimhan_types(self, d, theta, denominator=sum):
        # TODO what to return?
        # list of the Harder-Narasimhan types?
        # denominator default being sum is total dimension, there are variations possible
        # and the strata will be different!
        # https://mathscinet.ams.org/mathscinet-getitem?mr=1974891
        """Returns the list of all HN types."""

        """A Harder--Narasimhan (HN) type of d with respect to theta is a sequence d^* = (d^1,...,d^s) of dimension vectors such that
        * d^1 + ... + d^s = d
        * mu_theta(d^1) > ... > mu_theta(d^s)
        * Every d^k is theta-semi-stable."""

        """
        EXAMPLES

        The 3-Kronecker quiver:
        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: theta = vector([3,-2])
        sage: d = vector([2,3])
        sage: Q.all_harder_narasimhan_types(d,theta)
        [[(2, 3)],
         [(1, 1), (1, 2)],
         [(2, 2), (0, 1)],
         [(2, 1), (0, 2)],
         [(1, 0), (1, 3)],
         [(1, 0), (1, 2), (0, 1)],
         [(1, 0), (1, 1), (0, 2)],
         [(2, 0), (0, 3)]]
        sage: Q.all_harder_narasimhan_types(d,-theta)
        [[(0, 3), (2, 0)]]

        The 5-subspace quiver:
        sage: from quiver import *
        sage: Q = SubspaceQuiver(5)
        sage: d = vector([1,1,1,1,1,2])
        sage: theta = vector([2,2,2,2,2,-5])
        sage: Q.all_harder_narasimhan_types(d,theta)
        [[(1, 1, 1, 1, 1, 2)],
         [(0, 0, 1, 1, 1, 1), (1, 1, 0, 0, 0, 1)],
         [(0, 1, 0, 1, 1, 1), (1, 0, 1, 0, 0, 1)],
         [(0, 1, 1, 0, 1, 1), (1, 0, 0, 1, 0, 1)],
         [(0, 1, 1, 1, 0, 1), (1, 0, 0, 0, 1, 1)],
         [(1, 0, 0, 1, 1, 1), (0, 1, 1, 0, 0, 1)],
         [(1, 0, 1, 0, 1, 1), (0, 1, 0, 1, 0, 1)],
         [(1, 0, 1, 1, 0, 1), (0, 1, 0, 0, 1, 1)],
         [(1, 1, 0, 0, 1, 1), (0, 0, 1, 1, 0, 1)],
         [(1, 1, 0, 1, 0, 1), (0, 0, 1, 0, 1, 1)],
         [(1, 1, 1, 0, 0, 1), (0, 0, 0, 1, 1, 1)],
         [(0, 1, 1, 1, 1, 1), (1, 0, 0, 0, 0, 1)],
         [(1, 0, 1, 1, 1, 1), (0, 1, 0, 0, 0, 1)],
         [(1, 1, 0, 1, 1, 1), (0, 0, 1, 0, 0, 1)],
         [(1, 1, 1, 0, 1, 1), (0, 0, 0, 1, 0, 1)],
         [(1, 1, 1, 1, 0, 1), (0, 0, 0, 0, 1, 1)],
         [(1, 1, 1, 1, 1, 1), (0, 0, 0, 0, 0, 1)],
         [(0, 0, 0, 0, 1, 0), (1, 1, 1, 1, 0, 2)],
         [(0, 0, 0, 0, 1, 0), (0, 1, 1, 1, 0, 1), (1, 0, 0, 0, 0, 1)],
         [(0, 0, 0, 0, 1, 0), (1, 0, 1, 1, 0, 1), (0, 1, 0, 0, 0, 1)],
         [(0, 0, 0, 0, 1, 0), (1, 1, 0, 1, 0, 1), (0, 0, 1, 0, 0, 1)],
         [(0, 0, 0, 0, 1, 0), (1, 1, 1, 0, 0, 1), (0, 0, 0, 1, 0, 1)],
         [(0, 0, 0, 0, 1, 0), (1, 1, 1, 1, 0, 1), (0, 0, 0, 0, 0, 1)],
         [(0, 0, 0, 1, 0, 0), (1, 1, 1, 0, 1, 2)],
         [(0, 0, 0, 1, 0, 0), (0, 1, 1, 0, 1, 1), (1, 0, 0, 0, 0, 1)],
         [(0, 0, 0, 1, 0, 0), (1, 0, 1, 0, 1, 1), (0, 1, 0, 0, 0, 1)],
         [(0, 0, 0, 1, 0, 0), (1, 1, 0, 0, 1, 1), (0, 0, 1, 0, 0, 1)],
         [(0, 0, 0, 1, 0, 0), (1, 1, 1, 0, 0, 1), (0, 0, 0, 0, 1, 1)],
         [(0, 0, 0, 1, 0, 0), (1, 1, 1, 0, 1, 1), (0, 0, 0, 0, 0, 1)],
         [(0, 0, 0, 1, 1, 0), (1, 1, 1, 0, 0, 2)],
         [(0, 0, 0, 1, 1, 0), (0, 1, 1, 0, 0, 1), (1, 0, 0, 0, 0, 1)],
         [(0, 0, 0, 1, 1, 0), (1, 0, 1, 0, 0, 1), (0, 1, 0, 0, 0, 1)],
         [(0, 0, 0, 1, 1, 0), (1, 1, 0, 0, 0, 1), (0, 0, 1, 0, 0, 1)],
         [(0, 0, 0, 1, 1, 0), (1, 1, 1, 0, 0, 1), (0, 0, 0, 0, 0, 1)],
         [(0, 0, 1, 0, 0, 0), (1, 1, 0, 1, 1, 2)],
         [(0, 0, 1, 0, 0, 0), (0, 1, 0, 1, 1, 1), (1, 0, 0, 0, 0, 1)],
         [(0, 0, 1, 0, 0, 0), (1, 0, 0, 1, 1, 1), (0, 1, 0, 0, 0, 1)],
         [(0, 0, 1, 0, 0, 0), (1, 1, 0, 0, 1, 1), (0, 0, 0, 1, 0, 1)],
         [(0, 0, 1, 0, 0, 0), (1, 1, 0, 1, 0, 1), (0, 0, 0, 0, 1, 1)],
         [(0, 0, 1, 0, 0, 0), (1, 1, 0, 1, 1, 1), (0, 0, 0, 0, 0, 1)],
         [(0, 0, 1, 0, 1, 0), (1, 1, 0, 1, 0, 2)],
         [(0, 0, 1, 0, 1, 0), (0, 1, 0, 1, 0, 1), (1, 0, 0, 0, 0, 1)],
         [(0, 0, 1, 0, 1, 0), (1, 0, 0, 1, 0, 1), (0, 1, 0, 0, 0, 1)],
         [(0, 0, 1, 0, 1, 0), (1, 1, 0, 0, 0, 1), (0, 0, 0, 1, 0, 1)],
         [(0, 0, 1, 0, 1, 0), (1, 1, 0, 1, 0, 1), (0, 0, 0, 0, 0, 1)],
         [(0, 0, 1, 1, 0, 0), (1, 1, 0, 0, 1, 2)],
         [(0, 0, 1, 1, 0, 0), (0, 1, 0, 0, 1, 1), (1, 0, 0, 0, 0, 1)],
         [(0, 0, 1, 1, 0, 0), (1, 0, 0, 0, 1, 1), (0, 1, 0, 0, 0, 1)],
         [(0, 0, 1, 1, 0, 0), (1, 1, 0, 0, 0, 1), (0, 0, 0, 0, 1, 1)],
         [(0, 0, 1, 1, 0, 0), (1, 1, 0, 0, 1, 1), (0, 0, 0, 0, 0, 1)],
         [(0, 0, 1, 1, 1, 0), (1, 1, 0, 0, 0, 2)],
         [(0, 0, 1, 1, 1, 0), (1, 1, 0, 0, 0, 1), (0, 0, 0, 0, 0, 1)],
         [(0, 1, 0, 0, 0, 0), (1, 0, 1, 1, 1, 2)],
         [(0, 1, 0, 0, 0, 0), (0, 0, 1, 1, 1, 1), (1, 0, 0, 0, 0, 1)],
         [(0, 1, 0, 0, 0, 0), (1, 0, 0, 1, 1, 1), (0, 0, 1, 0, 0, 1)],
         [(0, 1, 0, 0, 0, 0), (1, 0, 1, 0, 1, 1), (0, 0, 0, 1, 0, 1)],
         [(0, 1, 0, 0, 0, 0), (1, 0, 1, 1, 0, 1), (0, 0, 0, 0, 1, 1)],
         [(0, 1, 0, 0, 0, 0), (1, 0, 1, 1, 1, 1), (0, 0, 0, 0, 0, 1)],
         [(0, 1, 0, 0, 1, 0), (1, 0, 1, 1, 0, 2)],
         [(0, 1, 0, 0, 1, 0), (0, 0, 1, 1, 0, 1), (1, 0, 0, 0, 0, 1)],
         [(0, 1, 0, 0, 1, 0), (1, 0, 0, 1, 0, 1), (0, 0, 1, 0, 0, 1)],
         [(0, 1, 0, 0, 1, 0), (1, 0, 1, 0, 0, 1), (0, 0, 0, 1, 0, 1)],
         [(0, 1, 0, 0, 1, 0), (1, 0, 1, 1, 0, 1), (0, 0, 0, 0, 0, 1)],
         [(0, 1, 0, 1, 0, 0), (1, 0, 1, 0, 1, 2)],
         [(0, 1, 0, 1, 0, 0), (0, 0, 1, 0, 1, 1), (1, 0, 0, 0, 0, 1)],
         [(0, 1, 0, 1, 0, 0), (1, 0, 0, 0, 1, 1), (0, 0, 1, 0, 0, 1)],
         [(0, 1, 0, 1, 0, 0), (1, 0, 1, 0, 0, 1), (0, 0, 0, 0, 1, 1)],
         [(0, 1, 0, 1, 0, 0), (1, 0, 1, 0, 1, 1), (0, 0, 0, 0, 0, 1)],
         [(0, 1, 0, 1, 1, 0), (1, 0, 1, 0, 0, 2)],
         [(0, 1, 0, 1, 1, 0), (1, 0, 1, 0, 0, 1), (0, 0, 0, 0, 0, 1)],
         [(0, 1, 1, 0, 0, 0), (1, 0, 0, 1, 1, 2)],
         [(0, 1, 1, 0, 0, 0), (0, 0, 0, 1, 1, 1), (1, 0, 0, 0, 0, 1)],
         [(0, 1, 1, 0, 0, 0), (1, 0, 0, 0, 1, 1), (0, 0, 0, 1, 0, 1)],
         [(0, 1, 1, 0, 0, 0), (1, 0, 0, 1, 0, 1), (0, 0, 0, 0, 1, 1)],
         [(0, 1, 1, 0, 0, 0), (1, 0, 0, 1, 1, 1), (0, 0, 0, 0, 0, 1)],
         [(0, 1, 1, 0, 1, 0), (1, 0, 0, 1, 0, 2)],
         [(0, 1, 1, 0, 1, 0), (1, 0, 0, 1, 0, 1), (0, 0, 0, 0, 0, 1)],
         [(0, 1, 1, 1, 0, 0), (1, 0, 0, 0, 1, 2)],
         [(0, 1, 1, 1, 0, 0), (1, 0, 0, 0, 1, 1), (0, 0, 0, 0, 0, 1)],
         [(0, 1, 1, 1, 1, 0), (1, 0, 0, 0, 0, 1), (0, 0, 0, 0, 0, 1)],
         [(1, 0, 0, 0, 0, 0), (0, 1, 1, 1, 1, 2)],
         [(1, 0, 0, 0, 0, 0), (0, 0, 1, 1, 1, 1), (0, 1, 0, 0, 0, 1)],
         [(1, 0, 0, 0, 0, 0), (0, 1, 0, 1, 1, 1), (0, 0, 1, 0, 0, 1)],
         [(1, 0, 0, 0, 0, 0), (0, 1, 1, 0, 1, 1), (0, 0, 0, 1, 0, 1)],
         [(1, 0, 0, 0, 0, 0), (0, 1, 1, 1, 0, 1), (0, 0, 0, 0, 1, 1)],
         [(1, 0, 0, 0, 0, 0), (0, 1, 1, 1, 1, 1), (0, 0, 0, 0, 0, 1)],
         [(1, 0, 0, 0, 1, 0), (0, 1, 1, 1, 0, 2)],
         [(1, 0, 0, 0, 1, 0), (0, 0, 1, 1, 0, 1), (0, 1, 0, 0, 0, 1)],
         [(1, 0, 0, 0, 1, 0), (0, 1, 0, 1, 0, 1), (0, 0, 1, 0, 0, 1)],
         [(1, 0, 0, 0, 1, 0), (0, 1, 1, 0, 0, 1), (0, 0, 0, 1, 0, 1)],
         [(1, 0, 0, 0, 1, 0), (0, 1, 1, 1, 0, 1), (0, 0, 0, 0, 0, 1)],
         [(1, 0, 0, 1, 0, 0), (0, 1, 1, 0, 1, 2)],
         [(1, 0, 0, 1, 0, 0), (0, 0, 1, 0, 1, 1), (0, 1, 0, 0, 0, 1)],
         [(1, 0, 0, 1, 0, 0), (0, 1, 0, 0, 1, 1), (0, 0, 1, 0, 0, 1)],
         [(1, 0, 0, 1, 0, 0), (0, 1, 1, 0, 0, 1), (0, 0, 0, 0, 1, 1)],
         [(1, 0, 0, 1, 0, 0), (0, 1, 1, 0, 1, 1), (0, 0, 0, 0, 0, 1)],
         [(1, 0, 0, 1, 1, 0), (0, 1, 1, 0, 0, 2)],
         [(1, 0, 0, 1, 1, 0), (0, 1, 1, 0, 0, 1), (0, 0, 0, 0, 0, 1)],
         [(1, 0, 1, 0, 0, 0), (0, 1, 0, 1, 1, 2)],
         [(1, 0, 1, 0, 0, 0), (0, 0, 0, 1, 1, 1), (0, 1, 0, 0, 0, 1)],
         [(1, 0, 1, 0, 0, 0), (0, 1, 0, 0, 1, 1), (0, 0, 0, 1, 0, 1)],
         [(1, 0, 1, 0, 0, 0), (0, 1, 0, 1, 0, 1), (0, 0, 0, 0, 1, 1)],
         [(1, 0, 1, 0, 0, 0), (0, 1, 0, 1, 1, 1), (0, 0, 0, 0, 0, 1)],
         [(1, 0, 1, 0, 1, 0), (0, 1, 0, 1, 0, 2)],
         [(1, 0, 1, 0, 1, 0), (0, 1, 0, 1, 0, 1), (0, 0, 0, 0, 0, 1)],
         [(1, 0, 1, 1, 0, 0), (0, 1, 0, 0, 1, 2)],
         [(1, 0, 1, 1, 0, 0), (0, 1, 0, 0, 1, 1), (0, 0, 0, 0, 0, 1)],
         [(1, 0, 1, 1, 1, 0), (0, 1, 0, 0, 0, 1), (0, 0, 0, 0, 0, 1)],
         [(1, 1, 0, 0, 0, 0), (0, 0, 1, 1, 1, 2)],
         [(1, 1, 0, 0, 0, 0), (0, 0, 0, 1, 1, 1), (0, 0, 1, 0, 0, 1)],
         [(1, 1, 0, 0, 0, 0), (0, 0, 1, 0, 1, 1), (0, 0, 0, 1, 0, 1)],
         [(1, 1, 0, 0, 0, 0), (0, 0, 1, 1, 0, 1), (0, 0, 0, 0, 1, 1)],
         [(1, 1, 0, 0, 0, 0), (0, 0, 1, 1, 1, 1), (0, 0, 0, 0, 0, 1)],
         [(1, 1, 0, 0, 1, 0), (0, 0, 1, 1, 0, 2)],
         [(1, 1, 0, 0, 1, 0), (0, 0, 1, 1, 0, 1), (0, 0, 0, 0, 0, 1)],
         [(1, 1, 0, 1, 0, 0), (0, 0, 1, 0, 1, 2)],
         [(1, 1, 0, 1, 0, 0), (0, 0, 1, 0, 1, 1), (0, 0, 0, 0, 0, 1)],
         [(1, 1, 0, 1, 1, 0), (0, 0, 1, 0, 0, 1), (0, 0, 0, 0, 0, 1)],
         [(1, 1, 1, 0, 0, 0), (0, 0, 0, 1, 1, 2)],
         [(1, 1, 1, 0, 0, 0), (0, 0, 0, 1, 1, 1), (0, 0, 0, 0, 0, 1)],
         [(1, 1, 1, 0, 1, 0), (0, 0, 0, 1, 0, 1), (0, 0, 0, 0, 0, 1)],
         [(1, 1, 1, 1, 0, 0), (0, 0, 0, 0, 1, 1), (0, 0, 0, 0, 0, 1)],
         [(1, 1, 1, 1, 1, 0), (0, 0, 0, 0, 0, 2)]]
        """

        n = self.number_of_vertices()
        zeroVector = vector([0 for i in range(n)])
        if (d == zeroVector):
            return [[zeroVector]]
        else:
            subdimensions = all_subdimension_vectors(d)
            # We consider just those subdimension vectors which are not zero, whose slope is bigger than the slope of d and which admit a semi-stable representation
            # Note that we also eliminate d by the following
            subdimensions = list(filter(lambda e: (e != zeroVector) and (slope(e,theta,denominator=denominator) > slope(d,theta,denominator=denominator)) and self.has_semistable_representation(e,theta,algorithm="schofield"), subdimensions))
            # We sort the subdimension vectors by slope because that will return the list of all HN types in ascending order with respect to the partial order from Def. 3.6 of https://mathscinet.ams.org/mathscinet-getitem?mr=1974891
            subdimensions.sort(key=(lambda e: slope(e,theta,denominator=denominator)))
            # The HN types which are not of the form (d) are given by (e,f^1,...,f^s) where e is a proper subdimension vector such that mu_theta(e) > mu_theta(d) and (f^1,...,f^s) is a HN type of f = d-e such that mu_theta(e) > mu_theta(f^1) holds.
            allHNtypes =  [[e]+fstar for e in subdimensions for fstar in list(filter(lambda fstar: slope(e,theta) > slope(fstar[0],theta) ,self.all_harder_narasimhan_types(d-e,theta,denominator=denominator)))]
            # Possibly add d again, at the beginning, because it is smallest with respect to the partial order from Def. 3.6
            if self.has_semistable_representation(d,theta,algorithm="schofield"):
                allHNtypes = [[d]] + allHNtypes
            return allHNtypes

    def all_weight_bounds(self, d, theta,denominator=sum):
        """
        Returns, for a given dimension vector d and a given stability parameter theta, the list of all weights to apply Teleman quantization.
        For each HN type, the 1-PS lambda acts on det(N_{S/R}|_Z) with a certain weight. Teleman quantization gives a numerical condition involving these weights to compute cohmology on the quotient.
        """
        # TODO return the Hn type as well?

        #This is only relevant on the unstable locus
        HN = list(filter(lambda hntype: hntype != [d] ,self.all_harder_narasimhan_types(d,theta,denominator=denominator)))

        return list(map(lambda hntype: -sum([(slope(hntype[s],theta,denominator=denominator) - slope(hntype[t],theta,denominator=denominator))*self.euler_form(hntype[s],hntype[t]) for s in range(len(hntype)-1) for t in range(s+1,len(hntype))] ), HN))
        

    def does_rigidity_inequality_hold(self,d,theta,denominator=sum):
        """
        Returns True if the rigidity inequality holds for d and theta, i.e. if the weights of the 1-PS lambda on det(N_{S/R}|_Z) for each HN type are all strictly larger than the weights of the tensors of the universal bundles U_i^\vee \otimes U_j.
        """

        #This is only relevant on the unstable locus
        HN = list(filter(lambda hntype: hntype != [d] ,self.all_harder_narasimhan_types(d,theta,denominator=denominator)))
        
        # We compute the weights of the 1-PS lambda on det(N_{S/R}|_Z) for each HN type
        weights = list(map(lambda hntype: -sum([(slope(hntype[s],theta,denominator=denominator) - slope(hntype[t],theta,denominator=denominator))*self.euler_form(hntype[s],hntype[t]) for s in range(len(hntype)-1) for t in range(s+1,len(hntype))] ), HN))
        
        # We compute the maximum weight of the tensors of the universal bundles U_i^\vee \otimes U_j
        tensorWeights = list(map(lambda hntype: slope(hntype[0],theta,denominator=denominator) - slope(hntype[-1],theta,denominator=denominator), HN))

        return all([weights[i] > tensorWeights[i] for i in range(len(HN))])





    def is_luna_type(self, tau, theta):
        """Checks if tau is a Luna type for theta."""
        n = self.number_of_vertices()
        zeroVector = vector([0 for i in range(n)])
        d = sum([sum(dn[1])*dn[0] for dn in tau])
        if (d == zeroVector):
            return (tau == [tuple([zeroVector,[1]])])
        else:
            dstar = [dn[0] for dn in tau]
            equalSlope = all([slope(e,theta,denominator=sum) == slope(d,theta,denominator=sum) for e in dstar])
            semistable = all([self.has_stable_representation(e,theta,algorithm="schofield") for e in dstar])
            return (equalSlope and semistable)

    def all_luna_types(self, d, theta):
        """Returns the unordered list of all Luna types of d for theta."""

        """A Luna type of d for theta is an unordered sequence (i.e. multiset) ((d^1,m_1),...,(d^s,m_s)) of dimension vectors d^k and (positive) natural numbers m_k such that
        * m_1d^1 + ... + m_sd^s = d
        * mu_theta(d^k) = mu_theta(d)
        * All d^k admit a theta-stable representation
        """

        """Example: Suppose that d = 3e and e, 2e, d = 3e are the only stable subdimension vectors. Then the Luna types are:
        ((3e,1))
        ((2e,1),(e,1))
        ((e,3))
        ((e,2),(e,1))
        ((e,1),(e,1),(e,1))
        """

        """Therefore we implement it as follows. A Luna type for us is a list [(d^1,p^1),...,(d^s,p^s)] (actually it should be unordered, but that's difficult because vectors are mutable) of dimension vectors d^k and (non-empty) partitions p^k such that
        * |p^1|d^1 + ... + |p^s|d^s = d
        * same
        * same """

        """So in the above example, the Luna types are
        [(3e,[1])]
        [(2e,[1]),(e,[1])]
        [(e,[3])]
        [(e,[2,1])]
        [(e,[1,1,1])]
        """

        """
        EXAMPLES

        The 3-Kronecker quiver:
        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: d = vector([3,3])
        sage: theta = vector([1,-1])
        sage: Q.all_luna_types(d,theta)
        [[((1, 1), [3])],
         [((1, 1), [2, 1])],
         [((1, 1), [1, 1, 1])],
         [((1, 1), [1]), ((2, 2), [1])],
         [((3, 3), [1])]]

        The 6-subspace quiver:
        sage: from quiver import *
        sage: Q = SubspaceQuiver(6)
        sage: d = vector([1,1,1,1,1,1,2])
        sage: theta = vector([1,1,1,1,1,1,-3])
        sage: Q.all_luna_types(d,theta)
        [[((0, 0, 0, 1, 1, 1, 1), [1]), ((1, 1, 1, 0, 0, 0, 1), [1])],
         [((0, 0, 1, 0, 1, 1, 1), [1]), ((1, 1, 0, 1, 0, 0, 1), [1])],
         [((0, 0, 1, 1, 0, 1, 1), [1]), ((1, 1, 0, 0, 1, 0, 1), [1])],
         [((0, 0, 1, 1, 1, 0, 1), [1]), ((1, 1, 0, 0, 0, 1, 1), [1])],
         [((0, 1, 0, 0, 1, 1, 1), [1]), ((1, 0, 1, 1, 0, 0, 1), [1])],
         [((0, 1, 0, 1, 0, 1, 1), [1]), ((1, 0, 1, 0, 1, 0, 1), [1])],
         [((0, 1, 0, 1, 1, 0, 1), [1]), ((1, 0, 1, 0, 0, 1, 1), [1])],
         [((0, 1, 1, 0, 0, 1, 1), [1]), ((1, 0, 0, 1, 1, 0, 1), [1])],
         [((0, 1, 1, 0, 1, 0, 1), [1]), ((1, 0, 0, 1, 0, 1, 1), [1])],
         [((0, 1, 1, 1, 0, 0, 1), [1]), ((1, 0, 0, 0, 1, 1, 1), [1])],
         [((1, 1, 1, 1, 1, 1, 2), [1])]]

        """

        n = self.number_of_vertices()
        zeroVector = vector([0 for i in range(n)])

        def partial_luna_types(d):
            """Returns the list of sets of the form {(d^1,n_1),...,(d^s,n_s)} such that all d^k are distinct."""
            subdimensions = all_subdimension_vectors(d)
            # We consider just those subdimension vectors which are not zero or d, whose slope equals the slope of d and which admit a stable representation
            subdimensions = list(filter(lambda e: (e != zeroVector) and (e != d) and (slope(e,theta,denominator=sum) == slope(d,theta,denominator=sum)) and self.has_stable_representation(e,theta,algorithm="schofield"), subdimensions))
            # Create all partial Luna types
            partialLunaTypes = []
            for e in subdimensions:
                smallerPartialLunaTypes = partial_luna_types(d-e)
                for tau in smallerPartialLunaTypes:
                    # Check if e occurs as a dimension vector in tau.
                    # If so, say of the form (e,n) then remove this occurrence and add (e,n+1)
                    # If not, then add (e,1)
                    occurs = False
                    for dn in tau:
                        if (dn[0] == e):
                            # We remove dn from tau and add the tuple (e,dn[1]+1) instead
                            tau.remove(dn)
                            tau.append(tuple([e,dn[1]+1]))
                            occurs = True
                    if (not occurs):
                        tau.append(tuple([e,1]))
                    # Now tau is a Luna type of d the desired form
                    # We sort it, because it's supposed to be unordered
                    tau = sorted(tau)
                    if tau not in partialLunaTypes:
                        partialLunaTypes = partialLunaTypes + [tau]
            if self.has_stable_representation(d,theta,algorithm="schofield"):
                partialLunaTypes = partialLunaTypes + [[tuple([d,1])]]
            return partialLunaTypes

        if (d == zeroVector):
            return [tuple([zeroVector,[1]])]
        else:
            partialLunaTypes = partial_luna_types(d)
            allLunaTypes = []
            for tau in partialLunaTypes:
                listOfPartitions = [Partitions(dn[1]).list() for dn in tau]
                Prod = cartesian_product(listOfPartitions).list()
                allLunaTypes = allLunaTypes + [[tuple([tau[i][0],p[i]]) for i in range(len(tau))] for p in Prod]
            return allLunaTypes

    def semistable_equals_stable(self, d, theta, algorithm="schofield"):
        """Checks if every theta-semistable representation of dimension vector d is theta-stable"""

        """
        EXAMPLES
        The 3-Kronecker quiver:
        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: d = vector([3,3])
        sage: theta = vector([1,-1])
        sage: Q.semistable_equals_stable(d,theta)
        False
        sage: d = vector([2,3])
        sage: Q.semistable_equals_stable(d,theta)
        True
        """

        """Every theta-semistable representation is theta-stable if and only if there are no Luna types other than (possibly) (d,[1])."""

        # As the computation of all Luna types takes so much time, we should first tests if d is theta-coprime
        if is_coprime_for_stability_parameter(d,theta):
            return True
        else:
            # This is probably the fastest way as checking theta-coprimality is fast whereas checking for existence of a semi-stable representation might be a bit slower
            if not self.has_semistable_representation(d,theta,algorithm=algorithm):
                return True
            else:
                allLunaTypes = self.all_luna_types(d,theta)
                genericType = tuple([d,[1]])
                if genericType in allLunaTypes:
                    allLunaTypes.remove(genericType)
                return (not allLunaTypes)

    def in_fundamental_domain(self, d):
        # see e.g. page 3 of https://arxiv.org/pdf/2303.08522.pdf

        # there has to be a more elegant way to do this
        # oh well
        simples = [ZeroVector(self.number_of_vertices()) for i in range(self.number_of_vertices())]
        for i in range(self.number_of_vertices()):
            simples[i][i] = 1
        return all(self.euler_form(d,i) + self.euler_form(i,d) <= 0 for i in simples)

    def partial_order(self,d,e):
        """Checks if d << e, which means that d_i <= e_i for every source i, d_j >= e_j for every sink j, and d_k == e_k for every vertex k which is neither a source nor a sink."""

        """
        EXAMPLES

        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: d = vector([1,1])
        sage: e = vector([2,1])
        sage: f = vector([2,2])
        sage: Q.partial_order(d,e)
        True
        sage: Q.partial_order(e,d)
        False
        sage: Q.partial_order(d,f)
        False
        sage: Q.partial_order(f,d)
        False
        sage: Q.partial_order(e,f)
        False
        sage: Q.partial_order(f,e)
        True

        sage: Q = ThreeVertexQuiver(2,2,2)
        sage: Q
        An acyclic 3-vertex quiver; adjacency matrix:
        [0 2 2]
        [0 0 2]
        [0 0 0]
        sage: d = vector([1,1,1])
        sage: e = vector([1,2,1])
        sage: Q.partial_order(d,e)
        False
        sage: Q.partial_order(e,d)
        False
        """

        n = self.number_of_vertices()
        assert (d.length() == n) and (e.length() == n)
        less = all([d[i-1] <= e[i-1] for i in list(filter(lambda i: self.is_source(i), range(1,n+1)))])
        less = less and all([d[j-1] >= e[j-1] for j in list(filter(lambda j: self.is_sink(j), range(1,n+1)))])
        less = less and all([d[k-1] == e[k-1] for k in list(filter(lambda k: (not self.is_source(k)) and (not self.is_sink(k)), range(1,n+1)))])

        return less

    def all_minimal_forbidden_subdimension_vectors(self,d,theta):
        """Returns the list of all minimal forbidden subdimension vectors of d."""

        """Minimality is with respect to the partial order e << d which means e_i <= d_i for every source i, e_j >= d_j for every sink j, and e_k = d_k for every vertex which is neither a source nor a sink."""

        """
        EXAMPLES

        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: d = vector([2,3])
        sage: theta = vector([3,-2])
        sage: Q.all_minimal_forbidden_subdimension_vectors(d,theta)
        [(1, 1), (2, 2)]
        """

        forbidden = all_forbidden_subdimension_vectors(d,theta)
        return list(filter(lambda e: not any([self.partial_order(f,e) for f in list(filter(lambda f: f != e, forbidden))]), forbidden))


    """
    Some things that don't below anywhere else yet?
    """

    def is_cofree(self, d):
        # https://mathscinet.ams.org/mathscinet-getitem?mr=2394691
        # we don't really know what this means
        raise NotImplementedError()

    def perpendicular_category(self, d):
        # something from Schofield
        # see Theorem 11.4.6 in the Derksen--Weyman book
        raise NotImplementedError()

"""Auxiliary methods"""

def slope(d, theta, denominator=sum):
    """For denominator = sum, the slope mu_theta(d) is defined as theta*d/(sum_i d_i). We need to ensure that sum(d) is positive."""
    assert (d.length() == theta.length())
    assert (denominator(d) > 0)
    return (theta*d)/(denominator(d))

def all_subdimension_vectors(d):
    """Returns the list of all subdimension vectors of d."""
    return list(map(vector, cartesian_product([range(di + 1) for di in d])))

def all_forbidden_subdimension_vectors(d,theta):
    """Returns the list of all subdimension vectors d' of d for which mu_theta(d') > mu_theta(d)."""

    """
    EXAMPLES

    sage: from quiver import *
    sage: d = vector([2,3])
    sage: theta = vector([3,-2])
    sage: all_forbidden_subdimension_vectors(d,theta)
    [(1, 0), (1, 1), (2, 0), (2, 1), (2, 2)]
    """

    zeroVector = vector([0 for i in range(d.length())])
    properSubdimensions = list(filter(lambda e: e != d and e != zeroVector, all_subdimension_vectors(d)))
    return list(filter(lambda e: slope(e,theta) > slope(d,theta), properSubdimensions))


def is_coprime_for_stability_parameter(d,theta):
    """Checks if d is theta-coprime."""

    """A dimension vector d is theta-coprime if mu_theta(e) != mu_theta(e) for all proper subdimension vectors e of d."""

    """
    EXAMPLES
    sage: from quiver import *
    sage: d = vector([2,3])
    sage: theta = vector([3,-2])
    sage: is_coprime_for_stability_parameter(d,theta)
    True
    sage: d = vector([3,3])
    sage: theta = vector([1,-1])
    sage: is_coprime_for_stability_parameter(d,theta)
    False
    """

    assert (d.length() == theta.length())
    zeroVector = vector([0 for i in range(d.length())])
    properSubdimensions = list(filter(lambda e: e != d and e != zeroVector, all_subdimension_vectors(d)))
    return all([slope(d,theta) != slope(e,theta) for e in properSubdimensions])

def is_indivisible(d):
    """Checks if the gcd of all entries is 1 or not."""
    return (gcd(d) == 1)

"""Class methods"""

def disjoint_union(Q1,Q2):
    """Returns the disjoint union of two quivers Q1 and Q2."""

    """
    EXAMPLES

    sage: from quiver import *
    sage: Q1 = GeneralizedKroneckerQuiver(3)
    sage: Q2 = GeneralizedKroneckerQuiver(4)
    sage: Q = disjoint_union(Q1,Q2)
    sage: Q
    Disjoint union of 3-Kronecker quiver and 4-Kronecker quiver; adjacency matrix:
    [0 3 0 0]
    [0 0 0 0]
    [0 0 0 4]
    [0 0 0 0]
    """

    if (Q1._name != None) and (Q2._name != None):
        name = "Disjoint union of " + Q1._name + " and " + Q2._name
    else:
        name = None

    return Quiver(block_diagonal_matrix(Q1._adjacency,Q2._adjacency, subdivide=False), name=name)

"""Special quivers"""

# TODO convention for generator functions is capitalise them?

def GeneralizedKroneckerQuiver(m):
    """
    The generalized Kronecker quiver has two vertices and $m$ arrows from the first to the second.

    TESTS::

        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: Q.number_of_vertices()
        2
        sage: Q.number_of_arrows()
        3

    """
    Q = Quiver(matrix([[0, m], [0, 0]]), name = str(m)+"-Kronecker quiver")
    # TODO do Q.rename here
    return Q

def KroneckerQuiver():
    """The Kronecker quiver has two vertices and 2 arrows from the first vertex to the second."""
    return GeneralizedKroneckerQuiver(2)


def ThreeVertexQuiver(m12, m13, m23):
    """An acyclic quiver with 3 vertices and mij many arrows i --> j for 1 <= i < j <= 3."""
    Q = Quiver(matrix([[0, m12, m13], [0, 0, m23], [0, 0, 0]]), name = "An acyclic 3-vertex quiver")
    # TODO do Q.rename here
    return Q


def LoopQuiver( m):
    """A quiver with one vertex and m arrows."""
    Q = Quiver(matrix([[m]]), name = str(m)+"-loop quiver")
    # TODO do Q.rename here
    return Q


def JordanQuiver(cls):
    Q = loop_quiver(1)
    # TODO do Q.rename here
    return Q


def SubspaceQuiver(m):
    """A quiver with m sources 1,...,m and one sink m+1; one arrow from every source to the sink."""
    A = zero_matrix(ZZ, m + 1)
    for i in range(m):
        A[i, m] = 1

    Q = Quiver(A, name = str(m)+"-subspace quiver")
    # TODO do Q.rename here
    return Q


def GeneralizedSubspaceQuiver( m, k):
    """A quiver with m sources 1,...,m and one sink m+1; k_i many arrows from source i to the sink."""
    assert (k.length() == m)
    A = zero_matrix(ZZ, m + 1)
    # I'm sure you can do this without a for loop
    for i in range(m):
        A[i, m] = k[i]

    Q = Quiver(A, name = "A generalized "+str(m)+"-subspace quiver")
    # TODO do Q.rename here
    return Q


def DynkinQuiver(Tn):
    r"""Returns the Dynkin quiver of type Tn. Uses the standard Sagemath implementation of Dynkin diagrams."""
    # use https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/root_system/dynkin_diagram.html
    # TODO: this constructor calls the adjacency_matrix() method many times. Should we call it once and remove lower triangular entries?
    
    #parse the string Tn
    T = Tn[:-1]
    n = int(Tn[-1])

    return Quiver(matrix(n, n, lambda i, j: DynkinDiagram(Tn).adjacency_matrix()[i,j] if i < j else 0), "Dynkin quiver of type "+Tn)


def ExtendedDynkinQuiver( T):
    # TODO implement this
    # TODO orientation: have a default (so for A and D: linear, type E?) but make it possible to change the orientation
    raise NotImplementedError()


def CyclicQuiver( n):
    return ExtendedDynkinQuiver(["A", n])


def BipartiteQuiver( m, n):
    # TODO implement this
    # m is number of sources, n is number of sinks
    raise NotImplementedError()


# Sampling and testing methods

def RandomQuiver(vertices,arrow_bound=10,acyclic=False,connected=True):
    """Returns a random Quiver object.

    Input: 
        - vertices: the number of vertices of the desired quiver;
        - acyclic: If True, the quiver will not have cycles. If false, it might but it is not guaranteed. Defaults to True; and
        - arrow_bound: the maximum amount of arrows between any two vertices. Defaults to 10.
        - connected: If True, the adjacency matrix is invertible, so that the underlying graph of the quiver is connected. Defaults to True.
    """
        # This while loop just samples candidate quivers until a connected one is found.
        # There has to be a better way to do this!
        # what other features should such a function have?


    if connected:
        acceptable = False

        while not acceptable:
            adjacency = random_matrix(ZZ,vertices,vertices, x=0,y=arrow_bound)

            if acyclic:
                # upper triangular matrix
                for i in range(vertices):
                    for j in range(i,vertices):
                        adjacency[j,i]=0

            acceptable = Quiver(adjacency).is_connected() # unnecessary overhead in defining Quiver object
    elif not connected:
        adjacency = random_matrix(ZZ,vertices,vertices, x=0,y=arrow_bound)
        
        if acyclic:
            # upper triangular matrix
            for i in range(vertices):
                for j in range(i,vertices):
                    adjacency[j,i]=0

    return Quiver(adjacency)

def RandomRoot(quiver,positive=True,upper_bound=10):
    """Returns a random root, or dimension vector, for the given quiver.
        Inputs: 
            - quiver: a Quiver object;
            - positive: if True, the root will not have zero entries. Defaults to False; and
            - upper_bound: an upper bound on the entries. Defaults to 10.
    """ 
    # what other features should such a function have?
    # if given a stability condition theta, an option to generate theta-coprime roots;
    # an option to generate indivisible roots;
    # ?
    lower_bound = 0
    if positive:
        lower_bound += 1
    return vector([randint(lower_bound,upper_bound) for i in range(quiver.number_of_vertices())])
    
def RandomStability(quiver,bound=10):
    """Returns a random stability condition for the given quiver.
        Inputs: 
            - quiver: a Quiver object;
            - bound: upper and lower bound on the entries. Defaults to 10.
    """
    # what other features should this have?

    return vector([randint(-bound//2,bound) for i in range(quiver.number_of_vertices())])

def ZeroVector(dimension):
    """Returns a zero vector of the given dimension."""
    return vector([0 for i in range(dimension)])
