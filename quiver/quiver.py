# from sage.matrix.constructor import matrix
from concurrent.futures.process import _threads_wakeups
from sage.all import *
import copy

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
        r"""Returns the amount of vertices that the quiver has.

        OUTPUT: The number of vertices as an Int.
        """
        return self.adjacency_matrix().nrows()

    def number_of_arrows(self):
        r"""Returns the number of arrows that the quiver has.

        OUTPUT: The number of arrows as an Int.

        """
        thin = self.thin_dimension_vector()
        return thin * self.adjacency_matrix() * thin

    def is_acyclic(self):
        r"""Returns the truth value of wether the quiver is acyclic.

        OUTPUT: Statement truth value as Bool.
        """
        A = self.adjacency_matrix()
        n = self.number_of_vertices()

        # a quiver is acyclic if and only if its adjacency matrix is nilpotent
        return (A**n == zero_matrix(ZZ, n))

    def is_connected(self):
        r"""Returns whether the underlying graph of the quiver is connected or not.

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

    def indegree(self, j):
        r"""Returns the indegree of a vertex.

        INPUT:
        ``j`` -- An Int between 1 and self.number_of_vertices()

        OUTPUT:
        The indegree as an Int
        """

        """The indegree of j is the number of incoming arrows into j.
        indeg(j) = sum_i a_{ij} where (a_{ij}) is the adjacency matrix."""
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


    def outdegree(self, i):
        r"""Returns the outdegree of a vertex.

        INPUT:
        ``i`` -- An Int between 1 and self.number_of_vertices()

        OUTPUT:
        The outdegree as an Int
        """

        """The outdegree of i is the number of outgoing arrows from i.
        outdeg(i) = sum_j a_{ij} where (a_{ij}) is the adjacency matrix."""

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

        OUTPUT: the multiplication of ``x * self.euler_matrix() * y`` as an  Int.

        """
        assert (x.length() == self.number_of_vertices() and y.length() == self.number_of_vertices())
        return x * self.euler_matrix() * y
    
    def symmetrized_euler_form(self, x, y):
        r"""The symmetrization of the Euler bilinear form of the quiver.

        INPUT:
        - ``x`` -- vector of integers
        - ``y`` -- vector of integers

        OUTPUT: the sum ``self.euler_form(x,y) + self.euler_form(y,x)`` as an  Int.

        """
        assert (x.length() == self.number_of_vertices() and y.length() == self.number_of_vertices())
        return self.euler_form(x, y) + self.euler_form(y, x)

    def tits_form(self, x):
        r"""The Tits quadratic form of the quiver.

        INPUT:
        - ``x`` -- vector of integers

        OUTPUT: the expression ``self.euler_form(x,x)`` as an  Int.

        """
        assert x.length() == self.number_of_vertices()
        return self.euler_form(x, x)


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
    
    def framed_quiver(self, f):
        r"""Returns the framed quiver with framing vector f.
        
        INPUT:
        - ``f``: vector of Ints

        OUTPUT: Quiver object
        """

        """The framed quiver has one additional vertex 0 and f_i many arrows from 0 to i."""
        
        n = self.number_of_vertices()
        assert (f.length() == n)
        A =  self.adjacency_matrix()
        # Adjacency matrix of the framed quiver looks like this (block shape):
        # [[0 f]
        #  [0 A]]
        # Add f as a first row
        A = A.insert_row(0,f)
        # Add a zero column
        A = A.transpose().insert_row(0,ZeroVector(n+1)).transpose()
        return Quiver(A)
    
    def coframed_quiver(self, f):
        r"""Returns the coframed quiver with framing vector f.

        INPUT:
        - ``f``: vector of Ints

        OUTPUT: Quiver object
        """

        """The coframed quiver has one additional vertex oo and f_i many arrows from i to oo."""

        n = self.number_of_vertices()
        assert (f.length() == n)
        A =  self.adjacency_matrix()
        # Adjacency matrix of the coframed quiver looks like this (block shape):
        # [[A f]
        #  [0 0]]
        # Add f as a last column
        A = A.transpose().insert_row(n,f).transpose()
        # Add a zero row as last row
        A = A.insert_row(n,ZeroVector(n+1))
        return Quiver(A)
    
    def full_subquiver(self, I):
        r"""Returns the full subquiver supported on the given set of vertices.
        
        INPUT:
        - ``I``: List

        OUTPUT: Quiver object
        """

        assert all([i in range(1,self.number_of_vertices()+1) for i in I])
        A = self.adjacency_matrix()
        return Quiver(A[[i-1 for i in I], [i-1 for i in I]])

    

    """
    Dimension vectors and stability conditions
    """

    def zero_vector(self):
        r"""Returns the zero vector.
        
        OUTPUT: vector of Ints
        """
        return vector([0 for i in range(self.number_of_vertices())])

    def thin_dimension_vector(self):
        r"""Returns the thin dimension vector.
        
        OUTPUT: vector of Ints
        """

        """The thin dimension vector is [1,...,1]."""
        return vector([1 for i in range(self.number_of_vertices())])

    def simple_root(self, i):
        r"""Returns the simple root at the vertex.
        
        INPUT:
        - ``i``: Int

        OUTPUT: vector of Ints
        """

        """The simple root at i is e_i = [0,...,1,...,0], i.e. the unit vector with a one in position i."""
        n = self.number_of_vertices()
        # Our convention is that vertices are numbered 1,...,n
        assert (i >= 1 and i <= n)
        ei = vector([0 for i in range(n)])
        ei[i-1] = 1
        return ei
    
    def is_root(self, x):
        r"""Checks if x is a root of the underlying diagram of the quiver.
        
        INPUT:
        - ``x``: vector of Ints

        OUTPUT: statement truth value as bool
        """

        """A root is a non-zero vector in Z^n such that the Tits form of x is <= 1."""
        assert x.length() == self.number_of_vertices()
        return (x != self.zero_vector() and self.tits_form(x) <= 1)
    
    def is_real_root(self, x):
        r"""Checks if x is a real root of the underlying diagram of the quiver.
        
        INPUT:
        - ``x``: vector of Ints

        OUTPUT: statement truth value as bool
        """

        """A root is called real if its Tits form equals 1."""
        assert x.length() == self.number_of_vertices()
        return self.tits_form(x) == 1
    
    def is_imaginary_root(self, x):
        r"""Checks if x is an imaginary root of the underlying diagram of the quiver.
        
        INPUT:
        - ``x``: vector of Ints

        OUTPUT: statement truth value as bool
        """

        """A root is called real if its Tits form is non-positive."""
        assert x.length() == self.number_of_vertices()
        return (x != self.zero_vector() and self.tits_form(x) <= 0)

    def support(self, d):
        r"""Returns the support of the dimension vector.
        
        INPUT:
        - ``d``: vector of Ints

        OUTPUT: Quiver object
        """

        """The support is the full subquiver supported on {i in Q_0 | d_i > 0}."""
        """
        EXAMPLES

        sage: from quiver import *
        sage: Q = ThreeVertexQuiver(2,0,4)
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
        r"""Checks if a dimension vector is in the fundamental domain.
        
        INPUT:
        - ``d``: vector of Ints

        OUTPUT: statement truth value as Bool
        """

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
        assert (d.length() == n and all([di >= 0 for di in list(d)]))
        # This is the condition <d,e_i> + <e_i,d> <= 0 for all i in Q_0
        eulerFormCondition = all([(self.euler_form(d,self.simple_root(i+1)) + self.euler_form(self.simple_root(i+1),d) <= 0) for i in range(n)])
        # Check if the support is connected
        connected = self.support(d).is_connected()
        return eulerFormCondition and connected
    
    # The fundamental domain again! Which implementation should we keep?
    # The latter is lacking the connectivity condition on the support of d
    
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
        # TODO: Think of a better name.

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

    def canonical_stability_parameter(self,d):
        """The canonical stability parameter is given by <d,_> - <_,d>"""
        E = self.euler_matrix()
        return d * (-self.euler_matrix().transpose() + E)
    
    """
    Generic subdimension vectors and generic Hom and Ext
    """

    # taken from code/snippets/canonical.sage
    # TODO still need testing code from there
    def is_generic_subdimension_vector(self, e, d):
        r"""Checks if e is a generic subdimension vector of d.
        
        INPUT:
        - ``e``: vector of Ints
        - ``d``: vector of Ints

        OUTPUT: Statement truth value as Bool        
        """
        
        # Using notation from Section 5 of https://arxiv.org/pdf/0802.2147.pdf
        """A dimension vector e is called a generic subdimension vector of d if a generic representation of dimension vector d possesses a subrepresentation of dimension vector e.
        By a result of Schofield (see Thm. 5.3 of https://arxiv.org/pdf/0802.2147.pdf) e is a generic subdimension vector of d if and only if e is a subdimension vector of d (missing in Thm. 5.3!) and <f,d-e> is non-negative for all generic subdimension vectors f of e."""

        """
        EXAMPLES:
        sage: from quiver import *
        sage: Q = LoopQuiver(1)
        sage: dims = [vector([i]) for i in range(3)]
        sage: for e in dims:
        ....:     for d in dims:
        ....:         print(str(e)+" gen. subdim of "+str(d)+"?: "+str(Q.is_generic_subdimension_vector(e,d)))
        ....: 
        (0) gen. subdim of (0)?: True
        (0) gen. subdim of (1)?: True
        (0) gen. subdim of (2)?: True
        (1) gen. subdim of (0)?: False
        (1) gen. subdim of (1)?: True
        (1) gen. subdim of (2)?: True
        (2) gen. subdim of (0)?: False
        (2) gen. subdim of (1)?: False
        (2) gen. subdim of (2)?: True
        sage: Q = LoopQuiver(2)
        sage: for e in dims:
        ....:     for d in dims:
        ....:         print(str(e)+" gen. subdim of "+str(d)+"?: "+str(Q.is_generic_subdimension_vector(e,d)))
        ....: 
        (0) gen. subdim of (0)?: True
        (0) gen. subdim of (1)?: True
        (0) gen. subdim of (2)?: True
        (1) gen. subdim of (0)?: False
        (1) gen. subdim of (1)?: True
        (1) gen. subdim of (2)?: False
        (2) gen. subdim of (0)?: False
        (2) gen. subdim of (1)?: False
        (2) gen. subdim of (2)?: True
        sage: Q = GeneralizedKroneckerQuiver(1)
        sage: for e in dims:
        ....:     for d in dims:
        ....:         if is_subdimension_vector(e,d):
        ....:             print(str(e)+" gen. subdim of "+str(d)+"?: "+str(Q.is_generic_subdimension_vector(e,d)))
        ....: 
        (0, 0) gen. subdim of (0, 0)?: True
        (0, 0) gen. subdim of (0, 1)?: True
        (0, 0) gen. subdim of (0, 2)?: True
        (0, 0) gen. subdim of (1, 0)?: True
        (0, 0) gen. subdim of (1, 1)?: True
        (0, 0) gen. subdim of (1, 2)?: True
        (0, 0) gen. subdim of (2, 0)?: True
        (0, 0) gen. subdim of (2, 1)?: True
        (0, 0) gen. subdim of (2, 2)?: True
        (0, 1) gen. subdim of (0, 1)?: True
        (0, 1) gen. subdim of (0, 2)?: True
        (0, 1) gen. subdim of (1, 1)?: True
        (0, 1) gen. subdim of (1, 2)?: True
        (0, 1) gen. subdim of (2, 1)?: True
        (0, 1) gen. subdim of (2, 2)?: True
        (0, 2) gen. subdim of (0, 2)?: True
        (0, 2) gen. subdim of (1, 2)?: True
        (0, 2) gen. subdim of (2, 2)?: True
        (1, 0) gen. subdim of (1, 0)?: True
        (1, 0) gen. subdim of (1, 1)?: False
        (1, 0) gen. subdim of (1, 2)?: False
        (1, 0) gen. subdim of (2, 0)?: True
        (1, 0) gen. subdim of (2, 1)?: True
        (1, 0) gen. subdim of (2, 2)?: False
        (1, 1) gen. subdim of (1, 1)?: True
        (1, 1) gen. subdim of (1, 2)?: True
        (1, 1) gen. subdim of (2, 1)?: True
        (1, 1) gen. subdim of (2, 2)?: True
        (1, 2) gen. subdim of (1, 2)?: True
        (1, 2) gen. subdim of (2, 2)?: True
        (2, 0) gen. subdim of (2, 0)?: True
        (2, 0) gen. subdim of (2, 1)?: False
        (2, 0) gen. subdim of (2, 2)?: False
        (2, 1) gen. subdim of (2, 1)?: True
        (2, 1) gen. subdim of (2, 2)?: False
        (2, 2) gen. subdim of (2, 2)?: True
        sage: Q = GeneralizedKroneckerQuiver(2)
        sage: for e in dims:
        ....:     for d in dims:
        ....:         if is_subdimension_vector(e,d):
        ....:             print(str(e)+" gen. subdim of "+str(d)+"?: "+str(Q.is_generic_subdimension_vector(e,d)))
        ....: 
        (0, 0) gen. subdim of (0, 0)?: True
        (0, 0) gen. subdim of (0, 1)?: True
        (0, 0) gen. subdim of (0, 2)?: True
        (0, 0) gen. subdim of (1, 0)?: True
        (0, 0) gen. subdim of (1, 1)?: True
        (0, 0) gen. subdim of (1, 2)?: True
        (0, 0) gen. subdim of (2, 0)?: True
        (0, 0) gen. subdim of (2, 1)?: True
        (0, 0) gen. subdim of (2, 2)?: True
        (0, 1) gen. subdim of (0, 1)?: True
        (0, 1) gen. subdim of (0, 2)?: True
        (0, 1) gen. subdim of (1, 1)?: True
        (0, 1) gen. subdim of (1, 2)?: True
        (0, 1) gen. subdim of (2, 1)?: True
        (0, 1) gen. subdim of (2, 2)?: True
        (0, 2) gen. subdim of (0, 2)?: True
        (0, 2) gen. subdim of (1, 2)?: True
        (0, 2) gen. subdim of (2, 2)?: True
        (1, 0) gen. subdim of (1, 0)?: True
        (1, 0) gen. subdim of (1, 1)?: False
        (1, 0) gen. subdim of (1, 2)?: False
        (1, 0) gen. subdim of (2, 0)?: True
        (1, 0) gen. subdim of (2, 1)?: False
        (1, 0) gen. subdim of (2, 2)?: False
        (1, 1) gen. subdim of (1, 1)?: True
        (1, 1) gen. subdim of (1, 2)?: False
        (1, 1) gen. subdim of (2, 1)?: True
        (1, 1) gen. subdim of (2, 2)?: True
        (1, 2) gen. subdim of (1, 2)?: True
        (1, 2) gen. subdim of (2, 2)?: True
        (2, 0) gen. subdim of (2, 0)?: True
        (2, 0) gen. subdim of (2, 1)?: False
        (2, 0) gen. subdim of (2, 2)?: False
        (2, 1) gen. subdim of (2, 1)?: True
        (2, 1) gen. subdim of (2, 2)?: False
        (2, 2) gen. subdim of (2, 2)?: True

        """

        assert (self.number_of_vertices() == d.length() & self.number_of_vertices() == e.length())
        assert (all([di >= 0 for di in d.list()]) and all([ei >= 0 for ei in e.list()]))
        
        if e == d:
            return True
        else:
            if not is_subdimension_vector(e, d):
                return False
            else: # e is subdimension vector of d
                # List of all generic subdimension vectors of e
                genSubdims = self.all_generic_subdimension_vectors(e)
                return all([self.euler_form(f, d-e) >= 0 for f in genSubdims])
            
    def __all_generic_subdimension_vectors_helper(self, d):
        """Returns the list of lists of indexes of all generic subdimension vectors of e, where e ranges over all subdimension vectors of d. The index refers to the deglex order."""

        """
        EXAMPLES

        sage: from quiver import *
        sage: Q, d = GeneralizedKroneckerQuiver(3), vector([2,3])
        sage: i, s = Q._Quiver__all_generic_subdimension_vectors_helper(d)
        sage: i, s
        ([[0],
        [0, 1],
        [0, 2],
        [0, 1, 3],
        [0, 1, 4],
        [0, 2, 5],
        [0, 1, 3, 6],
        [0, 1, 3, 7],
        [0, 1, 4, 8],
        [0, 1, 3, 6, 9],
        [0, 1, 3, 7, 10],
        [0, 1, 3, 6, 7, 9, 11]],
        [[(0, 0)],
        [(0, 0), (0, 1)],
        [(0, 0), (1, 0)],
        [(0, 0), (0, 1), (0, 2)],
        [(0, 0), (0, 1), (1, 1)],
        [(0, 0), (1, 0), (2, 0)],
        [(0, 0), (0, 1), (0, 2), (0, 3)],
        [(0, 0), (0, 1), (0, 2), (1, 2)],
        [(0, 0), (0, 1), (1, 1), (2, 1)],
        [(0, 0), (0, 1), (0, 2), (0, 3), (1, 3)],
        [(0, 0), (0, 1), (0, 2), (1, 2), (2, 2)],
        [(0, 0), (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]])

        """

        subdims = all_subdimension_vectors(d)            
        subdims.sort(key=(lambda e: deglex_key(e, b=max(d)+1)))
        # We use the deglex order because it's a total order which extends the usual entry-wise partial order on dimension vectors.
        N = len(subdims)

        # genIndexes[j] will in the end be the list of indexes (in subdims) of all generic subdimension vectors of subdims[j]
        genIndexes = [list(filter(lambda i: is_subdimension_vector(subdims[i], subdims[j]), range(N))) for j in range(N)]
        
        for j in range(N):
            genIndexes[j] = list(filter(lambda i: all([self.euler_form(subdims[k], subdims[j]-subdims[i]) >= 0 for k in genIndexes[i]]), genIndexes[j]))

        genSubdims = [[subdims[i] for i in genIndexes[j]] for j in range(N)]

        return genIndexes, genSubdims

    def all_generic_subdimension_vectors(self, d):
        r"""Returns the list of all generic subdimension vectors of d.
        
        INPUT:
        - ``d``: vector of Ints

        OUTPUT: list of vectors
        """

        """
        EXAMPLES:

        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(1)
        sage: d = vector([3,3])
        sage: Q.all_generic_subdimension_vectors(d)
        [(0, 0),
        (0, 1),
        (0, 2),
        (1, 1),
        (0, 3),
        (1, 2),
        (1, 3),
        (2, 2),
        (2, 3),
        (3, 3)]
        sage: Q = GeneralizedKroneckerQuiver(2)
        sage: Q.all_generic_subdimension_vectors(d)
        [(0, 0),
        (0, 1),
        (0, 2),
        (1, 1),
        (0, 3),
        (1, 2),
        (1, 3),
        (2, 2),
        (2, 3),
        (3, 3)]
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: Q.all_generic_subdimension_vectors(d)
        [(0, 0), (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3), (3, 3)]

        """
        assert (self.number_of_vertices() == d.length())
        assert all([di >= 0 for di in d.list()]) 

        genIndexes, genSubdims = self.__all_generic_subdimension_vectors_helper(d)
        N = len(genSubdims)
        return genSubdims[N-1]
    
    def generic_ext(self, a, b):
        r"""Computes ext(a,b).
        
        INPUT:
        - ``a``: vector of Ints
        - ``b``: vector of Ints

        OUTPUT: Int
        """

        """"According to Thm. 5.4 in Schofield's 'General representations of quivers', we have ext(a,b) = max{-<c,b> | c gen. subdimension vector of a}."""

        """
        EXAMPLES:
        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(1)
        sage: R = [Q.simple_root(1), Q.simple_root(2), vector([1,1])]
        sage: for a in R:
        ....:     for b in R:
        ....:         print('ext('+str(a)+','+str(b)+') = '+str(Q.generic_ext(a,b)))
        ....: 
        ext((1, 0),(1, 0)) = 0
        ext((1, 0),(0, 1)) = 1
        ext((1, 0),(1, 1)) = 0
        ext((0, 1),(1, 0)) = 0
        ext((0, 1),(0, 1)) = 0
        ext((0, 1),(1, 1)) = 0
        ext((1, 1),(1, 0)) = 0
        ext((1, 1),(0, 1)) = 0
        ext((1, 1),(1, 1)) = 0

        """

        genSubdims = self.all_generic_subdimension_vectors(a)
        return max([-self.euler_form(c,b) for c in genSubdims])
    
    def generic_hom(self, a, b):
        r"""Computes hom(a,b).

        INPUT:
        - ``a``: vector of Ints
        - ``b``: vector of Ints

        OUTPUT: Int
        """

        """There is a non-empty open subset U of R(Q,a) x R(Q,b) such that dim Ext(M,N) = ext(a,b), i.e. is minimal, for all (M,N) in U. Therefore dim Hom(M,N) = <a,b> + dim Ext(M,N) is minimal and therefore hom(a,b) = <a,b> + ext(a,b)."""

        """
        EXAMPLES:
        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(1)
        sage: R = [Q.simple_root(1), Q.simple_root(2), vector([1,1])]
        sage: for a in R:
        ....:     for b in R:
        ....:         print('hom('+str(a)+','+str(b)+') = '+str(Q.generic_hom(a,b)))
        ....: 
        hom((1, 0),(1, 0)) = 1
        hom((1, 0),(0, 1)) = 0
        hom((1, 0),(1, 1)) = 0
        hom((0, 1),(1, 0)) = 0
        hom((0, 1),(0, 1)) = 1
        hom((0, 1),(1, 1)) = 1
        hom((1, 1),(1, 0)) = 1
        hom((1, 1),(0, 1)) = 0
        hom((1, 1),(1, 1)) = 1

        """

        return self.euler_form(a,b) + self.generic_ext(a,b)

    def generic_ext_vanishing(self, a, b):
        return self.is_generic_subdimension_vector(a, a+b)
    
    def generic_hom_vanishing(self, a, b):
        # TODO figure out a way to implement this.
        # How about this:
        return self.generic_hom(a,b) == 0

    def is_left_orthogonal(self, a, b):
        if self.generic_ext_vanishing(a, b):
            return self.euler_form(a, b) == 0
        else:
            return False
    

    """
    Semistability and HN
    """

    def has_semistable_representation(self, d, theta):
        r"""Checks if there is a theta-semistable representation of dimension vector d.
        
        INPUT:
        - ``d``: vector of Ints
        - ``theta``: vector of Ints

        OUTPUT: Statement truth value as Bool
        """

        """See Thm. 5.4(1) of Reineke's overview paper https://arxiv.org/pdf/0802.2147.pdf: A dimension vector d admits a theta-semi-stable representation if and only if mu_theta(e) <= mu_theta(d) for all generic subdimension vectors e of d."""
        # Thm. 5.4 in Markus's paper is actually a result of Schofield. So the algorithm should bear his name, if any.

        """
        EXAMPLES:

        The A_2 quiver:
        sage: from quiver import *
        sage: A2 = GeneralizedKroneckerQuiver(1)
        sage: theta = vector([1,-1])
        sage: d = vector([1,1])
        sage: A2.has_semistable_representation(d,theta)
        True
        sage: d = vector([2,2])
        sage: A2.has_semistable_representation(d,theta)
        True
        sage: d = vector([1,2])
        sage: A2.has_semistable_representation(d,theta)
        False
        sage: d = vector([0,0])
        sage: A2.has_semistable_representation(d,theta)
        True

        The 3-Kronecker quiver:
        sage: from quiver import *
        sage: K3 = GeneralizedKroneckerQuiver(3)
        sage: theta = vector([3,-2])
        sage: d = vector([2,3])
        sage: K3.has_semistable_representation(d,theta)
        True
        sage: d = vector([1,4])
        sage: K3.has_semistable_representation(d,theta)
        False

        """

        genSubdims = self.all_generic_subdimension_vectors(d)
        genSubdims = list(filter(lambda e: e != self.zero_vector(), genSubdims))
        return all([slope(e, theta) <= slope(d, theta) for e in genSubdims])
        
    def __all_semistable_subdimension_vectors_helper(self, d, theta):
        """Computes the list of indexes of all semistable subdimension vectors of d."""

        """
        EXAMPLES:
        sage: from quiver import *
        sage: Q, d, theta = GeneralizedKroneckerQuiver(1), vector([2,3]), vector([1,0])
        sage: i, s = Q._Quiver__all_semistable_subdimension_vectors_helper(d, theta); i, s
        ([1, 2, 3, 4, 5, 6, 10],
        [(0, 1), (1, 0), (0, 2), (1, 1), (2, 0), (0, 3), (2, 2)])

        """

        subdims = all_subdimension_vectors(d)            
        subdims.sort(key=(lambda e: deglex_key(e, b=max(d)+1)))
        # We use the deglex order because it's a total order which extends the usual entry-wise partial order on dimension vectors.
        N = len(subdims)
        genIndexes, genSubdims = self.__all_generic_subdimension_vectors_helper(d)
        sstIndexes =  list(filter(lambda j: all([slope(subdims[i], theta) <= slope(subdims[j], theta) for i in list(filter(lambda i: i != 0, genIndexes[j]))]), range(1,N)))
        sstSubdims = [subdims[j] for j in sstIndexes]
        return sstIndexes, sstSubdims
    
    def is_harder_narasimhan_type(self, dstar, theta, denominator=sum):
        r"""Checks if dstar is a HN type.
        
        INPUT:
        - ``dstar``: list of vectors of Ints
        - ``theta``: vector of Ints
        - ``denominator``: function which takes a vector of Ints and returns an Int

        OUTPUT: statement truth value as Bool
        """

        """
        EXAMPLES:
        sage: from quiver import *
        sage: Q, d, theta = GeneralizedKroneckerQuiver(3), vector([2,3]), vector([1,0])
        sage: hn = Q.all_harder_narasimhan_types(d, theta)
        sage: all([Q.is_harder_narasimhan_type(dstar, theta) for dstar in hn])
        True
        sage: dstar = [vector([1,0]), vector([1,0]), vector([0,3])]
        sage: Q.is_harder_narasimhan_type(dstar, theta)
        False

        """

        n = self.number_of_vertices()
        assert (all([e.length() == n for e in dstar]) and theta.length() == n)
        assert all([denominator(self.simple_root(i)) > 0 for i in range(1,n+1)])

        d = sum(dstar)
        if (d == self.zero_vector()):
            return (dstar == [self.zero_vector()])
        else:
            sstIndexes, sstSubdims = self.__all_semistable_subdimension_vectors_helper(d, theta)
            slopeDecreasing = all([(slope(dstar[i],theta,denominator=denominator) > slope(dstar[i+1],theta,denominator=denominator)) for i in range(len(dstar)-1)])
            semistable = all([e in sstSubdims for e in dstar])
            return (slopeDecreasing and semistable)
        
    def __codimension_of_harder_narasimhan_stratum_helper(self, dstar):
        """Computes the codimension of the HN stratum of dstar inside the representation variety.
        
        INPUT:
        - ``dstar``: list of vectors of Ints

        OUTPUT: codimension as Int
        """
        # This is private because it doesn't check if dstar is a HN type. This is fast but yields nonsense, if dstar is not a HN type.

        """The codimension of the HN stratum of d^* = (d^1,...,d^s) is given by - sum_{k < l} <d^k,d^l>"""

        """
        EXAMPLES

        The 3-Kronecker quiver
        sage: from quiver import *
        sage: Q, d, theta = GeneralizedKroneckerQuiver(3), vector([2,3]), vector([1,0])
        sage: hn = Q.all_harder_narasimhan_types(d, theta); hn
        [[(1, 0), (1, 1), (0, 2)],
        [(1, 0), (1, 2), (0, 1)],
        [(1, 0), (1, 3)],
        [(1, 1), (1, 2)],
        [(2, 0), (0, 3)],
        [(2, 1), (0, 2)],
        [(2, 2), (0, 1)],
        [(2, 3)]]
        sage: [Q._Quiver__codimension_of_harder_narasimhan_stratum_helper(dstar) for dstar in hn]
        [12, 9, 8, 3, 18, 10, 4, 0]

        """
        
        n = self.number_of_vertices()
        assert all([e.length() == n for e in dstar])

        s = len(dstar)
        return -sum([self.euler_form(dstar[k],dstar[l]) for k in range(s-1) for l in range(k+1,s)])

    def codimension_of_harder_narasimhan_stratum(self, dstar, theta, denominator=sum):
        r"""Computes the codimension of the HN stratum of dstar inside the representation variety, if dstar is a HN type.
        
        INPUT:
        - ``dstar``: list of vectors of Ints
        - ``theta``: vector of Ints
        - ``denominator``: function which takes a vector of Ints and returns an Int

        OUTPUT: codimension as Int
        """

        """The codimension of the HN stratum of d^* = (d^1,...,d^s) is given by - sum_{k < l} <d^k,d^l>"""

        """
        EXAMPLES

        The 3-Kronecker quiver
        sage: from quiver import *
        sage: Q, d, theta = GeneralizedKroneckerQuiver(3), vector([2,3]), vector([1,0])
        sage: hn = Q.all_harder_narasimhan_types(d, theta); hn
        [[(1, 0), (1, 1), (0, 2)],
        [(1, 0), (1, 2), (0, 1)],
        [(1, 0), (1, 3)],
        [(1, 1), (1, 2)],
        [(2, 0), (0, 3)],
        [(2, 1), (0, 2)],
        [(2, 2), (0, 1)],
        [(2, 3)]]
        sage: [Q.codimension_of_harder_narasimhan_stratum(dstar) for dstar in hn]
        [12, 9, 8, 3, 18, 10, 4, 0]

        """
        assert theta.length() == self.number_of_vertices()
        assert self.is_harder_narasimhan_type(dstar, theta, denominator=denominator)

        return self.__codimension_of_harder_narasimhan_stratum_helper(dstar)
    
    def codimension_unstable_locus(self, d, theta):
        r"""Computes the codimension of the unstable locus inside the representation variety.
        
        INPUT:
        - ``d``: vector of Ints
        - ``theta``: vector of Ints

        OUTPUT: codimension as Int
        """

        """"
        EXAMPLES:
        sage: from quiver import *
        sage: Q, d, theta = GeneralizedKroneckerQuiver(3), vector([2,3]), vector([1,0])
        sage: Q.codimension_unstable_locus(d, theta)
        3
        sage: Q, d = ThreeVertexQuiver(1,6,1), vector([1,6,6])
        sage: theta = Q.canonical_stability_parameter(d)
        sage: Q.codimension_unstable_locus(d, theta)
        1
        sage: Q, d, theta = GeneralizedKroneckerQuiver(1), vector([2,3]), vector([1,0])
        sage: Q.codimension_unstable_locus(d, theta)
        0

        """

        n = self.number_of_vertices()
        assert (d.length() == n and theta.length() == n)

        hn = list(filter(lambda dstar: dstar != [d], self.all_harder_narasimhan_types(d, theta, denominator=sum)))
        # Note that while the HN types and strata depend on the denominator, the list of all their codimensions does not.

        return min([self.__codimension_of_harder_narasimhan_stratum_helper(dstar) for dstar in hn])


    def all_harder_narasimhan_types(self, d, theta, denominator=sum):
        # TODO what to return?
        # list of the Harder-Narasimhan types?
        # denominator default being sum is total dimension, there are variations possible
        # and the strata will be different!
        # https://mathscinet.ams.org/mathscinet-getitem?mr=1974891
        # Can the above TODO go?

        r"""Returns the list of all HN types.
        
        INPUT:
        - ``d``: vector of Ints
        - ``theta``: vector of Ints
        - ``denominator``: function which takes a vector of Ints and returns an Int

        OUTPUT: list of list of vectors of Ints
        """

        """A Harder--Narasimhan (HN) type of d with respect to theta is a sequence d^* = (d^1,...,d^s) of dimension vectors such that
        * d^1 + ... + d^s = d
        * mu_theta(d^1) > ... > mu_theta(d^s)
        * Every d^k is theta-semi-stable."""

        """
        EXAMPLES:

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
        assert (d.length() == n and theta.length())
        assert all([denominator(self.simple_root(i)) > 0 for i in range(1,n+1)])
        
        subdimensions = all_subdimension_vectors(d)
        subdimensions.sort(key=(lambda e: deglex_key(e, b=max(d)+1)))
        N = len(subdimensions)

        # sstIndexes is the list of indexes of all non-zero semistable subdimension vectors in subdimensions
        sstIndexes, sstSubdims = self.__all_semistable_subdimension_vectors_helper(d, theta)

        # idx_diff(j, i) is the index of the difference subdimensions[j]-subdimensions[i] in the list subdimensions
        idx_diff = (lambda j, i: subdimensions.index(subdimensions[j]-subdimensions[i]))

        hn = [[[]] for j in range(N)]

        for j in range(1,N):
            # sstSub is the list of all indexes in subdimensions of semistable non-zero subdimension vectors of subdimensions[j]
            sstSub = list(filter(lambda i: is_subdimension_vector(subdimensions[i], subdimensions[j]), sstIndexes))
            # The HN types which are not of the form (d) are given by (e,f^1,...,f^s) where e is a proper subdimension vector such that mu_theta(e) > mu_theta(d) and (f^1,...,f^s) is a HN type of f = d-e such that mu_theta(e) > mu_theta(f^1) holds.
            hn[j] = [[i]+fstar for i in sstSub for fstar in list(filter(lambda fstar: fstar == [] or slope(subdimensions[i], theta, denominator=denominator) > slope(subdimensions[fstar[0]], theta, denominator=denominator), hn[idx_diff(j, i)]))]

        hn[0] = [[0]]

        return [[subdimensions[r] for r in fstar] for fstar in hn[N-1]]

    """
    Stability and Luna
    """

    def has_stable_representation(self, d, theta, algorithm="schofield"):
        r"""Checks if there is a stable representation of this dimension vector.
        
        INPUT:
        - ``d``: vector of Ints
        - ``theta``: vector of Ints
        - ``algorithm``: String

        OUTPUT: statement truth value as Bool
        """

        assert (algorithm in ["schofield", "king", "al"])
        n = self.number_of_vertices()
        assert (d.length() == n and theta.length() == n)

        # TODO implement this
        # https://mathscinet.ams.org/mathscinet-getitem?mr=1315461
        # Question concerning above TODO: What is King's algorithm for checking for existence of stable representations supposed to be? I can't find one in the paper.
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

        if (algorithm == "schofield"):
            if d == self.zero_vector():
                return False
            else:
                genSubdims = self.all_generic_subdimension_vectors(d)
                genSubdims = list(filter(lambda e: e != self.zero_vector() and e != d, genSubdims))
                return all([slope(e, theta) < slope(d, theta) for e in genSubdims])


    def is_schur_root(self, d):
        r"""Checks if d is a Schur root.
        
        INPUT:
        - ``d``: vector of Ints

        OUTPUT: statement truth value as Bool
        """

        """
        A Schur root is a dimension vector which admits a Schurian representation, i.e. a representation whose endomorphism ring is k. It's necessarily indecomposable.
        By a result of Schofield (https://mathscinet.ams.org/mathscinet/relay-station?mr=1162487) d is a Schur root if and only if d admits a stable representation for the canonical stability parameter."""

        """
        EXAMPLES:

        The 3-Kronecker quiver:
        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: d = vector([2,3])
        sage: Q.is_schur_root(d)
        True

        Examples from Derksen--Weyman's book (Ex. 11.1.4):
        sage: from quiver import *
        sage: Q = ThreeVertexQuiver(1,1,1)
        sage: a = vector([1,1,2])
        sage: Q.is_schur_root(a)
        True
        sage: b = vector([1,2,1])
        sage: Q.is_schur_root(b)
        False
        sage: c = vector([1,1,1])
        sage: Q.is_schur_root(c)
        True
        sage: d = vector([2,2,2])
        sage: Q.is_schur_root(d)
        False

        """

        assert d.length() == self.number_of_vertices()

        theta = self.canonical_stability_parameter(d)
        return self.has_stable_representation(d, theta)
    
    def __all_stable_subdimension_vectors_helper(self, d, theta, denominator=sum):
        """Computes the list of all stable subdimension vectors of d which have the same slope as d."""

        """
        EXAMPLES:

        sage: from quiver import *
        sage: Q, d, theta = GeneralizedKroneckerQuiver(3), vector([3,3]), vector([1,0])
        sage: all_subdimension_vectors(d)
        [(0, 0),
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 0),
        (1, 1),
        (1, 2),
        (1, 3),
        (2, 0),
        (2, 1),
        (2, 2),
        (2, 3),
        (3, 0),
        (3, 1),
        (3, 2),
        (3, 3)]
        sage: i, s = Q._Quiver__all_stable_subdimension_vectors_helper(d, theta)
        sage: i
        [4, 11, 15]
        sage: s
        [(1, 1), (2, 2), (3, 3)]

        sage: from quiver import *
        sage: Q, d, theta = GeneralizedKroneckerQuiver(2), vector([3,3]), vector([1,-1])
        sage: Q._Quiver__all_stable_subdimension_vectors_helper(d, theta)
        ([4], [(1, 1)])

        """
        
        subdims = all_subdimension_vectors(d)            
        subdims.sort(key=(lambda e: deglex_key(e, b=max(d)+1)))
        # We use the deglex order because it's a total order which extends the usual entry-wise partial order on dimension vectors.
        N = len(subdims)
        genIndexes, genSubdims = self.__all_generic_subdimension_vectors_helper(d)
        # slopeIndexes is the list of subdimension vectors of d of the same slope as d (in particular != 0)
        slopeIndexes = list(filter(lambda j: slope(subdims[j], theta, denominator=denominator) == slope(d, theta, denominator=denominator), range(1,N)))
        # stIndexes contains all j for which subdims[j] is stable
        # e = subdims[j] is stable if for all generic subdimension vectors f = subdims[i] of e, it holds that slope(f) < slope(e)
        stIndexes =  list(filter(lambda j: all([slope(subdims[i], theta, denominator=denominator) < slope(subdims[j], theta, denominator=denominator) for i in list(filter(lambda i: i != 0 and i != j, genIndexes[j]))]), slopeIndexes))
        stSubdims = [subdims[j] for j in stIndexes]
        return stIndexes, stSubdims

    
    def is_luna_type(self, tau, theta, denominator=sum):
        r"""Checks if tau is a Luna type for theta.
        
        INPUT:
        - ``tau``: list of tuples
        - ``theta``: vector of Ints
        - ``denominator``: Int-valued function

        OUTPUT: statement truth value as Bool
        """
        """
        EXAMPLES:

        """

        n = self.number_of_vertices()
        assert (theta.length() == n and all([dn[0].length() == n for dn in tau]))
        
        d = sum([sum(dn[1])*dn[0] for dn in tau])
        if (d == self.zero_vector()):
            return (tau == [tuple([self.zero_vector(),[1]])])
        else:
            dstar = [dn[0] for dn in tau]
            stIndexes, stSubdims = self.__all_stable_subdimension_vectors_helper(d, theta, denominator=denominator)
            return all([e in stSubdims for e in dstar]) # Note that in particular the zero vector must not lie in dstar

    def all_luna_types(self, d, theta, denominator=sum):
        r"""Returns the unordered list of all Luna types of d for theta.
        
        INPUT:
        - ``d``: vector of Ints
        - ``theta``: vector of Ints
        - ``denominator``: Int-valued function

        OUTPUT: list of tuples containing Int-vector and Int 
        """

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

        The Kronecker quiver:
        sage: from quiver import *
        sage: Q, d, theta = KroneckerQuiver(), vector([3,3]), vector([1,-1])
        sage: Q.all_luna_types(d, theta)
        [[((1, 1), [3])], [((1, 1), [2, 1])], [((1, 1), [1, 1, 1])]]

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
        
        if (d == self.zero_vector()):
            return [tuple([self.zero_vector(),[1]])]
        else: 
            subdims = all_subdimension_vectors(d)
            subdims.sort(key=(lambda e: deglex_key(e, b=max(d)+1)))
            N = len(subdims)
            # slopeIndexes is the list of indexes j such that the slope of e := subdims[j] equals the slope of d (this requires e != 0)
            slopeIndexes = list(filter(lambda j: slope(subdims[j], theta, denominator=denominator) == slope(d, theta, denominator=denominator), range(1,N)))
            # We consider all subdimension vectors which are not zero, whose slope equals the slope of d, and which admit a stable representation
            # They're in deglex order by the way the helper function works.
            stIndexes, stSubdims = self.__all_stable_subdimension_vectors_helper(d, theta, denominator=denominator)
            # idx_diff(j, i) is the index of the difference stSubdims[j]-stSubdims[i] in the list stSubdims
            idx_diff = (lambda j, i: subdims.index(subdims[j]-subdims[i]))

            # partialLunaTypes is going to hold all "partial Luna types" of e for every e in stSubdims; a partial luna type of e is an unordered sequence (i.e. multiset) {(e^1,n_1),...,(e^s,n_s)} such that all e^k are distinct, e^1+...+e^s = e and the slopes of all e^k are the same (and thus equal the slope of e).
            partialLunaTypes = [[] for j in range(N)]
            for j in range(N):
                stSub = list(filter(lambda i: is_subdimension_vector(subdims[i], subdims[j]) and i != j, stIndexes))
                for i in stSub:
                    smaller = partialLunaTypes[idx_diff(j,i)]
                    for tau in smaller:
                        # Check if f := stSubdims[i] occurs as a dimension vector in tau.
                        # If so, say of the form (f,n) then remove this occurrence and add (f,n+1)
                        # If not, then add (f,1)
                        tauNew = copy.deepcopy(tau)
                        occurs = False
                        for dn in tauNew:
                            if (dn[0] == i):
                                # We remove dn from tau and add the tuple (e,dn[1]+1) instead
                                tauNew.remove(dn)
                                tauNew.append(tuple([i,dn[1]+1]))
                                occurs = True
                        if (not occurs):
                            tauNew.append(tuple([i,1]))
                        # Now tauNew is a Luna type of e := subdims[j] the desired form
                        # We sort it, because it's supposed to be unordered
                        tauNew.sort()
                        # If tau isn't already contained, then we add it
                        if tauNew not in partialLunaTypes[j]:
                            partialLunaTypes[j] = partialLunaTypes[j] + [tauNew]
                if (j in stIndexes):
                    # If e = subdims[j] is stable then (e,1) is also a Luna type.
                    partialLunaTypes[j] = partialLunaTypes[j] + [[tuple([j,1])]]
        
            partial = partialLunaTypes[N-1]
            allLunaTypes = []
            for tau in partial:
                listOfPartitions = [Partitions(dn[1]).list() for dn in tau]
                Prod = cartesian_product(listOfPartitions).list()
                allLunaTypes = allLunaTypes + [[tuple([subdims[tau[i][0]],p[i]]) for i in range(len(tau))] for p in Prod]
            return allLunaTypes


    def semistable_equals_stable(self, d, theta):
        r"""Checks if every theta-semistable representation of dimension vector d is theta-stable
        
        INPUT:
        - ``d``: vector of Ints
        - ``theta``: vector of Ints

        OUTPUT: statement truth value as Bool
        """

        """
        EXAMPLES:

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
            # This is probably the fastest way as checking theta-coprimality is fast whereas checking for existence of a semi-stable representation is a bit slower
            if not self.has_semistable_representation(d,theta):
                return True
            else:
                allLunaTypes = self.all_luna_types(d,theta)
                genericType = tuple([d,[1]])
                if genericType in allLunaTypes:
                    allLunaTypes.remove(genericType)
                return (not allLunaTypes) # This checks if the list is empty

    """
    Ample stability
    """

    # TODO dimension vectors should have .is_stable(), .is_amply_stable()?
    def is_amply_stable(self, d, theta):
        r"""Checks if d is amply stable for theta.
        
        INPUT:
        - ``d``: vector of Ints
        - ``theta``: vector of Ints

        OUTPUT: statement truth value as Bool
        """
        
        """By definition, d is theta-amply stable if the codimension of the theta-stable locus inside R(Q,d) is at least 2."""

        # By Prop. 4.1 of https://arxiv.org/pdf/1410.0466.pdf d is amply stable for theta provided that <e,d-e> <= -2 for every proper subdimension vector.
        # But can we find a necessary and sufficient condition?
        # If every theta-semi-stable representation of dimension vector d is theta-stable then theta-ample stability is equivalent to every proper HN stratum having codimension at least 2.

        """
        EXAMPLES:

        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: d = vector([2,3])
        sage: theta = vector([3,-2])
        sage: Q.is_amply_stable(d,theta)
        True
        sage: Q.is_amply_stable(d,-theta)
        False

        A three vertex example from the rigidity paper:
        sage: from quiver import *
        sage: Q = ThreeVertexQuiver(1,6,1)
        sage: d = vector([1,6,6])
        sage: theta = Q.canonical_stability_parameter(d)
        sage: Q.is_amply_stable(d, theta)
        False

        """

        if self.semistable_equals_stable(d, theta):
            return self.codimension_unstable_locus(d, theta) >= 2
        else:
            raise NotImplementedError()

    def is_strongly_amply_stable(self, d, theta):
        r"""Checks if d is strongly amply stable for theta.
        
        INPUT:
        - ``d``: vector of Ints
        - ``theta``: vector of Ints

        OUTPUT: statement truth value as Bool
        """

        """We call d strongly amply stable for theta if <e,d-e> <= -2 holds for all subdimension vectors e of d which satisfy slope(e) >= slope(d)."""

        """
        EXAMPLES:

        sage: from quiver import *
        sage: Q, d, theta = GeneralizedKroneckerQuiver(3), vector([2,3]), vector([3,-2])
        sage: Q.is_strongly_amply_stable(d, theta)
        True

        sage: from quiver import *
        sage: Q = ThreeVertexQuiver(5,1,1)
        sage: d = vector([4,1,4])
        sage: theta = Q.canonical_stability_parameter(d)
        sage: Q.is_amply_stable(d, theta)
        True
        sage: Q.is_strongly_amply_stable(d, theta)
        False

        """

        # All subdimension vectors of d
        es = all_subdimension_vectors(d)
        # Remove 0 and d
        es.remove(self.zero_vector())
        es.remove(d)
        # Filter out those of bigger slope
        es = list(filter(lambda e: slope(e,theta) >= slope(d,theta), es))
        return all([self.euler_form(e,d-e) <= -2 for e in es])

    """
    Canonical decomposition
    """
    
    def rearrange_dw_decomposition(self,decomposition,i,j):
        # this is non functional, let alone tested
        # apply Lemma 11.9.10 of Derksen--Weyman to rearrange the roots from i to j as k, k+1
        S = []
        for k in range(i,j-1):
            # if k has a ``path to j'' satisfying the hypothesis of Lemma 11.9.10, we keep it
            # m is the j-k times j-k matrix with entry m[i,j] = 1 if !generic_hom_vanishing(Q,decomposition[j][2],decomposition[i][2]) and 0 otherwise
            m = matrix(j-k)
            for l in range(k,j-1):
                for s in range(l+1,j):
                    if not self.generic_hom_vanishing(decomposition[s][2],decomposition[l][2]):
                        m[l-k,s-k] = 1
            paths = matrix(j-k) # paths[i,j] is the number of paths from k + i - 1 to k + j - 1 of length at most j-k
            for l in range(j-k):
                paths = paths + m**l
            if paths[0,j-k-1] > 0:
                S.append(k)
        rearrangement = [l for l in range(i+1,j-1) if l not in S]
        final = S + [i,j] + rearrangement
        decomposition_temp = [decomposition[l] for l in final]
        for l in range(i,j):
            decomposition[l] = decomposition_temp[l-i]
        return i + len(S) - 1 #this is the index of the element i now.


    def canonical_decomposition(self, d, algorithm="derksen-weyman"):
        # TODO implement this
        # https://mathscinet.ams.org/mathscinet-getitem?mr=1930979
        # this is implemented in code/snippets/canonical.sage, so include it here
        """
        There's something wrong with this implementation. I ran it on the Kronecker quiver and it gave the following:

        sage: from quiver import *
        sage: Q = KroneckerQuiver()
        sage: d = vector([2,4])
        sage: Q.canonical_decomposition(d, algorithm="derksen-weyman")
        [[2, (1, 0)], [4, (0, 1)]]

        But clearly a general representation of dimension vector (2,4) isn't semisimple. It is instead a direct sum of two copies of the preprojective representation of dimension vector (1,2). So the canonical decomposition should be [2, (1, 2)].
        """
        if algorithm == "derksen-weyman":
            decomposition = [[d[i],self.simple_root(i+1)] for i in range(self.number_of_vertices())]
            while True:
                decomposition = list(filter(lambda root: root[0] > 0, decomposition))
                violating_pairs = [(i,j) for i in range(len(decomposition)-1) for j in range(i+1,len(decomposition)) if self.euler_form(decomposition[j][1],decomposition[i][1]) < 0]
                if not violating_pairs:
                    break
                violating_pairs = sorted(violating_pairs, key=(lambda pair: pair[1] - pair[0]))
                i,j = violating_pairs[0]

                # this modifies the decomposition in place
                i = self.rearrange_dw_decomposition(decomposition,i,j)

                p,xi = decomposition[i]
                q,eta = decomposition[i+1]
                xi_real = self.is_real_root(xi)
                eta_real = self.is_real_root(eta)
                zeta = p*xi + q*eta
                if xi_real and eta_real:
                    discriminant = self.euler_form(zeta,zeta)
                    if discriminant > 0:
                        pass # TODO figure out what this should do
                    elif discriminant == 0:
                        zeta_prime = zeta // gcd(zeta)
                        del decomposition[i+1]
                        decomposition[i] = [1,zeta_prime]
                    else:
                        del decomposition[i+1]
                        decomposition[i] = [1,zeta]
                elif xi_real and not eta_real:
                    if p + q*self.euler_form(eta,xi) >= 0:
                        del decomposition[i+1]
                        decomposition[i] = [1,eta - self.euler_form(eta,xi)*xi]
                    else:
                        del decomposition[i+1]
                        decomposition[i] = [1,zeta]
                elif not xi_real and eta_real:
                    if q + p*self.euler_form(eta,xi) >= 0:
                        decomposition[i] = [1,eta]
                        decomposition[i+1] = [1,xi - self.euler_form(eta,xi)*eta]
                    else:
                        del decomposition[i+1]
                        decomposition[i] = [1,zeta]
                elif not xi_real and not eta_real:
                    del decomposition[i+1]
                    decomposition[i] = [1,zeta]
            return decomposition

        # https://mathscinet.ams.org/mathscinet-getitem?mr=1162487
        elif algorithm == "schofield-1":
            raise NotImplementedError()
        # TODO implement this
        # https://arxiv.org/pdf/math/9911014.pdf (see Section 5, and also Section 3 of https://mathscinet.ams.org/mathscinet/article?mr=1789222)
        # in Derksen--Weyman's https://mathscinet.ams.org/mathscinet-getitem?mr=1930979 it is claimed that there is a second Schofield algorithm
        # (they do cite the wrong Schofield preprint though...)
        elif algorithm == "schofield-2":
            raise NotImplementedError()
        
        elif algorithm == "recursive":
            """I'm not sure if one of the Schofield algorithms is meant to be this one. But here's a very simple recursion which computes the canonical decomposition. It is based on Lem. 11.2.5 in Derksen--Weyman's book (a.k.a. 'the book with typos'):
        
            Lemma: Let a be a dimension vector and a = b+c a decomposition such that ext(b,c) = ext(c,b) = 0. If b = b_1 + ... + b_s and c = c_1 + ... + c_t are the canonical decompositions, then a = b_1 + ... + b_s + c_1 + ... + c_t is the canonical decomposition of a.

            If no non-trivial decomposition a = b+c as above exists, then a is a Schur root and therefore its own canonical decomposition. This is because a generic representation has no subdimension vector b which admits both a subrepresentation and a quotient representation. So a generic representation is indecomposable, which implies that a is a Schur root.
            """

            """
            EXAMPLES:

            sage: Q = KroneckerQuiver()
            sage: ds = [vector([i,j]) for i in range(7) for j in range(7)]
            sage: can = [Q.canonical_decomposition(d, algorithm="recursive") for d in ds]
            sage: for i in range(7):
            ....:     for j in range(7):
            ....:         print('Canonical decomp. of '+str(ds[7*i+j])+' is: '+str(can[7*i+j]))
            ....: 
            Canonical decomp. of (0, 0) is: [(0, 0)]
            Canonical decomp. of (0, 1) is: [(0, 1)]
            Canonical decomp. of (0, 2) is: [(0, 1), (0, 1)]
            Canonical decomp. of (0, 3) is: [(0, 1), (0, 1), (0, 1)]
            Canonical decomp. of (0, 4) is: [(0, 1), (0, 1), (0, 1), (0, 1)]
            Canonical decomp. of (0, 5) is: [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1)]
            Canonical decomp. of (0, 6) is: [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1)]
            Canonical decomp. of (1, 0) is: [(1, 0)]
            Canonical decomp. of (1, 1) is: [(1, 1)]
            Canonical decomp. of (1, 2) is: [(1, 2)]
            Canonical decomp. of (1, 3) is: [(0, 1), (1, 2)]
            Canonical decomp. of (1, 4) is: [(0, 1), (0, 1), (1, 2)]
            Canonical decomp. of (1, 5) is: [(0, 1), (0, 1), (0, 1), (1, 2)]
            Canonical decomp. of (1, 6) is: [(0, 1), (0, 1), (0, 1), (0, 1), (1, 2)]
            Canonical decomp. of (2, 0) is: [(1, 0), (1, 0)]
            Canonical decomp. of (2, 1) is: [(2, 1)]
            Canonical decomp. of (2, 2) is: [(1, 1), (1, 1)]
            Canonical decomp. of (2, 3) is: [(2, 3)]
            Canonical decomp. of (2, 4) is: [(1, 2), (1, 2)]
            Canonical decomp. of (2, 5) is: [(0, 1), (1, 2), (1, 2)]
            Canonical decomp. of (2, 6) is: [(0, 1), (0, 1), (1, 2), (1, 2)]
            Canonical decomp. of (3, 0) is: [(1, 0), (1, 0), (1, 0)]
            Canonical decomp. of (3, 1) is: [(1, 0), (2, 1)]
            Canonical decomp. of (3, 2) is: [(3, 2)]
            Canonical decomp. of (3, 3) is: [(1, 1), (1, 1), (1, 1)]
            Canonical decomp. of (3, 4) is: [(3, 4)]
            Canonical decomp. of (3, 5) is: [(1, 2), (2, 3)]
            Canonical decomp. of (3, 6) is: [(1, 2), (1, 2), (1, 2)]
            Canonical decomp. of (4, 0) is: [(1, 0), (1, 0), (1, 0), (1, 0)]
            Canonical decomp. of (4, 1) is: [(1, 0), (1, 0), (2, 1)]
            Canonical decomp. of (4, 2) is: [(2, 1), (2, 1)]
            Canonical decomp. of (4, 3) is: [(4, 3)]
            Canonical decomp. of (4, 4) is: [(1, 1), (1, 1), (1, 1), (1, 1)]
            Canonical decomp. of (4, 5) is: [(4, 5)]
            Canonical decomp. of (4, 6) is: [(2, 3), (2, 3)]
            Canonical decomp. of (5, 0) is: [(1, 0), (1, 0), (1, 0), (1, 0), (1, 0)]
            Canonical decomp. of (5, 1) is: [(1, 0), (1, 0), (1, 0), (2, 1)]
            Canonical decomp. of (5, 2) is: [(1, 0), (2, 1), (2, 1)]
            Canonical decomp. of (5, 3) is: [(2, 1), (3, 2)]
            Canonical decomp. of (5, 4) is: [(5, 4)]
            Canonical decomp. of (5, 5) is: [(1, 1), (1, 1), (1, 1), (1, 1), (1, 1)]
            Canonical decomp. of (5, 6) is: [(5, 6)]
            Canonical decomp. of (6, 0) is: [(1, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0)]
            Canonical decomp. of (6, 1) is: [(1, 0), (1, 0), (1, 0), (1, 0), (2, 1)]
            Canonical decomp. of (6, 2) is: [(1, 0), (1, 0), (2, 1), (2, 1)]
            Canonical decomp. of (6, 3) is: [(2, 1), (2, 1), (2, 1)]
            Canonical decomp. of (6, 4) is: [(3, 2), (3, 2)]
            Canonical decomp. of (6, 5) is: [(6, 5)]
            Canonical decomp. of (6, 6) is: [(1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1)]
            sage:
            sage: all([all([Q.generic_ext(s[i],s[j]) + Q.generic_ext(s[j],s[i]) == 0 for i in range(len(s)) for j in range(i)]) for s in can])
            True
            """

            genSubdims = self.all_generic_subdimension_vectors(d)
            genSubdims = list(filter(lambda e: e != self.zero_vector() and e != d, genSubdims))
            for e in genSubdims:
                if d-e in genSubdims:
                    return self.canonical_decomposition(e, algorithm="recursive") + self.canonical_decomposition(d-e, algorithm="recursive")
            return [d]
        
        elif algorithm == "recursive_new":
            """
            EXAMPLES:

            sage: from quiver import *
            sage: Q = KroneckerQuiver()
            sage: d = vector([4,6])
            sage: Q.canonical_decomposition(d, algorithm="recursive_new")
            [(2, 3), (2, 3)]
            sage: ds = [vector([i,j]) for i in range(7) for j in range(7)]
            sage: all([Q.canonical_decomposition(d, algorithm="recursive") == Q.canonical_decomposition(d, algorithm="recursive_new") for d in ds])
            True

            """
            subdims = all_subdimension_vectors(d)            
            subdims.sort(key=(lambda e: deglex_key(e, b=max(d)+1)))
            N = len(subdims)

            idx_diff = (lambda j, i: subdims.index(subdims[j]-subdims[i]))

            genIndexes, genSubdims = self.__all_generic_subdimension_vectors_helper(d)

            def canon_indexes(j):
                """Computes for j in range(N) the list of indexes in subdims for the canonical decomposition of subdims[j]"""
                for i in list(filter(lambda i: i != 0 and i != j, genIndexes[j])):
                    k = idx_diff(j,i)
                    if k in genIndexes[j]:
                        return canon_indexes(i) + canon_indexes(k)
                return [j]
            
            return [subdims[i] for i in canon_indexes(N-1)]   

    
    """
    Teleman!
    """            

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
        r"""
        Returns True if the rigidity inequality holds for d and theta, i.e. if the weights of the 1-PS lambda on $\det(N_{S/R}|_Z)$ for each HN type are all strictly larger than the weights of the tensors of the universal bundles $U_i^\vee \otimes U_j$.
        """

        #This is only relevant on the unstable locus
        HN = list(filter(lambda hntype: hntype != [d] ,self.all_harder_narasimhan_types(d,theta,denominator=denominator)))

        # We compute the weights of the 1-PS lambda on det(N_{S/R}|_Z) for each HN type
        weights = list(map(lambda hntype: -sum([(slope(hntype[s],theta,denominator=denominator) - slope(hntype[t],theta,denominator=denominator))*self.euler_form(hntype[s],hntype[t]) for s in range(len(hntype)-1) for t in range(s+1,len(hntype))] ), HN))

        # We compute the maximum weight of the tensors of the universal bundles U_i^\vee \otimes U_j
        tensorWeights = list(map(lambda hntype: slope(hntype[0],theta,denominator=denominator) - slope(hntype[-1],theta,denominator=denominator), HN))

        return all([weights[i] > tensorWeights[i] for i in range(len(HN))])
    
    """
    Hochschild cohomology
    """

    def first_hochschild_cohomology(self):
        r"""
        Compute the dimension of the first Hochschild cohomology

        This uses the formula of Happel from Proposition 1.6 in [MR1035222].

        EXAMPLES:

        The first Hochschild cohomology of the $m$th generalized Kronecker quiver
        is the dimension of $\mathrm{PGL}_{m+1}$::

            sage: from quiver import *
            sage: GeneralizedKroneckerQuiver(3).first_hochschild_cohomology()
            8

        The first Hochschild cohomology vanishes if and only if the quiver is a tree::

            sage: from quiver import *
            sage: SubspaceQuiver(7).first_hochschild_cohomology()
            0

        """
        assert self.is_acyclic()
        # TODO think about including G into the Quiver class
        G = Graph(self.adjacency_matrix(), format="adjacency_matrix")

        # see Proposition 1.6 in Happel's "Hochschild cohomology of finite-dimensional algebras"
        return 1 - self.number_of_vertices() + sum(len(G.all_paths(a[0], a[1], use_multiedges=True)) for a in G.edges())


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

def is_subdimension_vector(e,d):
    assert (e.length() == d.length())
    n = e.length()
    return all([e[i] <= d[i] for i in range(n)])

def deglex_key(e, b):
    """A function which satisfies e <_{deglex} d iff deglex_key(e) < deglex_key(d), provided that b >> 0."""
    n = e.length()
    return sum([e[i]*b**(n-i-1) for i in range(n)])+sum(e)*b**n

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

# TODO if optional parameter is given, call GeneralizedKroneckerQuiver
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

def RandomDimensionVector(quiver,positive=False,upper_bound=10):
    """Returns a random dimension vector for the given quiver.
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

# TODO I (= Pieter) would not write this function:
# it should be possible to give dimension vectors as tuples (coerce whenever needed)
# and [0]*dimension is then perfectly good shorthand for ZeroVector(dimension)
def ZeroVector(dimension):
    """Returns a zero vector of the given dimension."""
    return vector([0 for i in range(dimension)])
