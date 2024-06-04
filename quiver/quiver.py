# from sage.matrix.constructor import matrix
from concurrent.futures.process import _threads_wakeups
from sage.all import *


class Quiver:
    # TODO this explanation is about implementation details, which don't matter
    # it should be possible to also create a quiver from a DiGraph
    # if so: things like indegree should explicitly refer to the vertices of the graph (and not 1, ..., n)
    # if the adjacency matrix constructor is given: explain _which_ DiGraph would be created
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
        if self._name == None:
            output += "A quiver with "
        else:
            output += str(self._name) + "; "
        output += "adjacency matrix:\n" + str(self._adjacency)
        return output

    """
    Basic graph-theoretic properties of the quiver
    """

    def adjacency_matrix(self):
        r"""Returns the adjacency matrix of the quiver.

        OUTPUT: A square matrix M whose entry M[i,j] is the number of arrows from the vertex i to the vertex j.
        """
        return self._adjacency

    # TODO is there a good reason to call this underlying_graph, and not just graph?
    def underlying_graph(self):
        r"""Returns the (necessarily symmetric) adjacency matrix of the underlying graph of the quiver.

        OUTPUT: A square, symmetric matrix M whose entry M[i,j] = M[j,i] is the number of edges between the vertices i and j.
        """
        return (
            self.adjacency_matrix()
            + self.adjacency_matrix().transpose()
            - diagonal_matrix(self.adjacency_matrix().diagonal())
        )

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
        return A**n == zero_matrix(ZZ, n)

    # TODO some of the examples should really be tests (so they don't clutter the docs as much)
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
        for i in range(2, self.number_of_vertices()):  # -1 ?
            # add all paths of length i
            paths += paths * self.underlying_graph()
        # if every couple of vertices is connected (ij \neq 0 or ji \neq 0) then true, otherwise false.
        for i in range(self.number_of_vertices()):
            for j in range(self.number_of_vertices()):
                if i != j and paths[i, j] == 0 and paths[j, i] == 0:
                    return False
        return True

    """
    Some graph-theoretic properties of the quiver
    """

    # TODO docstrings are split into pieces
    def indegree(self, j):
        r"""Returns the indegree of a vertex.

        INPUT:
        ``j`` -- An Int between 1 and self.number_of_vertices()

        OUTPUT:
        The indegree as an Int

        The indegree of j is the number of incoming arrows into j.
        indeg(j) = sum_i a_{ij} where (a_{ij}) is the adjacency matrix.
        # TODO Question: Should we number the vertices 1,...,n or 0,...,n-1?

        EXAMPLES:

        In the 3-Kronecker quiver the indegree is either 0 or 3::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.indegree(1)
            0
            sage: Q.indegree(2)
            3

        """

        assert (j > 0) and (j <= self.number_of_vertices())
        return sum(self._adjacency.column(j - 1))

    def outdegree(self, i):
        r"""Returns the outdegree of a vertex.

        INPUT:
        ``i`` -- An Int between 1 and self.number_of_vertices()

        OUTPUT:
        The outdegree as an Int

        The outdegree of i is the number of outgoing arrows from i.
        outdeg(i) = sum_j a_{ij} where (a_{ij}) is the adjacency matrix.

        EXAMPLES:

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.outdegree(1)
            3
            sage: Q.outdegree(2)
            0

        """

        assert (i > 0) and (i <= self.number_of_vertices())
        return sum(self._adjacency.row(i - 1))

    def is_source(self, i):
        """Checks if i is a source of the quiver, i.e. if there are no incoming arrows into i.

        EXAMPLES

        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: Q.is_source(1)
        True
        sage: Q.is_source(2)
        False
        """

        assert (i > 0) and (i <= self.number_of_vertices())
        return self.indegree(i) == 0

    def is_sink(self, j):
        """Checks if j is a sink of the quiver, i.e. if there are no outgoing arrows out of j.

        EXAMPLES

        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: Q.is_sink(1)
        False
        sage: Q.is_sink(2)
        True
        """

        assert (j > 0) and (j <= self.number_of_vertices())
        return self.outdegree(j) == 0

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
        assert (
            x.length() == self.number_of_vertices()
            and y.length() == self.number_of_vertices()
        )
        return x * self.euler_matrix() * y

    def symmetrized_euler_form(self, x, y):
        r"""The symmetrization of the Euler bilinear form of the quiver.

        INPUT:
        - ``x`` -- vector of integers
        - ``y`` -- vector of integers

        OUTPUT: the sum ``self.euler_form(x,y) + self.euler_form(y,x)`` as an  Int.

        """
        assert (
            x.length() == self.number_of_vertices()
            and y.length() == self.number_of_vertices()
        )
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

        if self._name != None:
            name = "Opposite of " + self._name
        else:
            name = None

        return Quiver(A, name)

    def double_quiver(self):
        """The adjacency matrix of the double of a quiver is the sum of the adjacency matrix of the original quiver and its transpose."""
        A = self.adjacency_matrix() + self.adjacency_matrix().transpose()
        if self._name != None:
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
        assert f.length() == n
        A = self.adjacency_matrix()
        # Adjacency matrix of the framed quiver looks like this (block shape):
        # [[0 f]
        #  [0 A]]
        # Add f as a first row
        A = A.insert_row(0, f)
        # Add a zero column
        A = A.transpose().insert_row(0, ZeroVector(n + 1)).transpose()
        return Quiver(A)

    def coframed_quiver(self, f):
        r"""Returns the coframed quiver with framing vector f.

        INPUT:
        - ``f``: vector of Ints

        OUTPUT: Quiver object
        """

        """The coframed quiver has one additional vertex oo and f_i many arrows from i to oo."""

        n = self.number_of_vertices()
        assert f.length() == n
        A = self.adjacency_matrix()
        # Adjacency matrix of the coframed quiver looks like this (block shape):
        # [[A f]
        #  [0 0]]
        # Add f as a last column
        A = A.transpose().insert_row(n, f).transpose()
        # Add a zero row as last row
        A = A.insert_row(n, ZeroVector(n + 1))
        return Quiver(A)

    def full_subquiver(self, I):
        r"""Returns the full subquiver supported on the given set of vertices.

        INPUT:
        - ``I``: List

        OUTPUT: Quiver object

        EXAMPLES:

        sage: from quiver import *
        sage: Q = ThreeVertexQuiver(2,3,4); Q
        An acyclic 3-vertex quiver; adjacency matrix:
        [0 2 3]
        [0 0 4]
        [0 0 0]
        sage: Q.full_subquiver([1,2])
        A quiver with adjacency matrix:
        [0 2]
        [0 0]
        sage: Q.full_subquiver([1,3])
        A quiver with adjacency matrix:
        [0 3]
        [0 0]

        """

        assert all([i in range(1, self.number_of_vertices() + 1) for i in I])
        A = self.adjacency_matrix()
        return Quiver(A[[i - 1 for i in I], [i - 1 for i in I]])

    """
    Dimension vectors and roots
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
        assert i >= 1 and i <= n
        ei = vector([0 for i in range(n)])
        ei[i - 1] = 1
        return ei

    def is_root(self, x):
        r"""Checks if x is a root of the underlying diagram of the quiver.

        INPUT:
        - ``x``: vector of Ints

        OUTPUT: statement truth value as bool
        """

        """A root is a non-zero vector in Z^n such that the Tits form of x is <= 1."""
        assert x.length() == self.number_of_vertices()
        return x != self.zero_vector() and self.tits_form(x) <= 1

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
        return x != self.zero_vector() and self.tits_form(x) <= 0

    def is_schur_root(self, d):
        r"""Checks if d is a Schur root.

        INPUT:
        - ``d``: vector of Ints

        OUTPUT: statement truth value as Bool

        A Schur root is a dimension vector which admits a Schurian representation, i.e. a representation whose endomorphism ring is k. It's necessarily indecomposable.
        By a result of Schofield (https://mathscinet.ams.org/mathscinet/relay-station?mr=1162487) d is a Schur root if and only if d admits a stable representation for the canonical stability parameter.

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

    def support(self, d):
        r"""Returns the support of the dimension vector.

        INPUT:

        - ``d``: dimension vector

        OUTPUT: List

        The support is the set {i in Q_0 | d_i > 0}.
        # I know it doesn't actually depend on Q, but I think I want to make it a method of the object. The reason is that I might want to introduce the possibility to give names to the vertices and arrows. Then support() should return the subset {i in Q_0 | d_i > 0} (as a list) and not the list of indexes of the vector.

        EXAMPLES

        sage: from quiver import *
        sage: Q = ThreeVertexQuiver(2, 0, 4)
        sage: d = vector([1, 1, 1])
        sage: Q.support(d)
        [1, 2, 3]
        sage: d = vector([1, 0, 1])
        sage: Q.support(d)
        [1, 3]

        """

        supp = list(filter(lambda i: d[i] > 0, range(self.number_of_vertices())))
        return [i + 1 for i in supp]

    def in_fundamental_domain(self, d):
        # TODO optional parameter for strict interior? maybe even specify how deep into the strict interior?
        r"""Checks if a dimension vector is in the fundamental domain.

        INPUT:

        - ``d``: dimension vector

        OUTPUT: statement truth value as Bool

        The fundamental domain of Q is the set of dimension vectors d such that supp(d) is connected and <d,e_i> + <e_i,d> <= 0 for all simple roots e_i.
        Every d in the fundamental domain is an imaginary root and the set of imaginary roots is the Weyl group saturation of the fundamental domain.
        If d is in the fundamental domain then it is Schurian and a general representation of dimension vector d is stable for the canonical stability parameter.

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
        assert d.length() == n and all([di >= 0 for di in list(d)])
        # This is the condition <d,e_i> + <e_i,d> <= 0 for all i in Q_0
        eulerFormCondition = all(
            [
                (
                    self.euler_form(d, self.simple_root(i + 1))
                    + self.euler_form(self.simple_root(i + 1), d)
                    <= 0
                )
                for i in range(n)
            ]
        )
        # Check if the support is connected
        connected = self.support(d).is_connected()
        return eulerFormCondition and connected

    # The fundamental domain again! Which implementation should we keep?
    # The latter is lacking the connectivity condition on the support of d

    # def in_fundamental_domain(self, d):
    #     # see e.g. page 3 of https://arxiv.org/pdf/2303.08522.pdf

    #     # there has to be a more elegant way to do this
    #     # oh well
    #     simples = [ZeroVector(self.number_of_vertices()) for i in range(self.number_of_vertices())]
    #     for i in range(self.number_of_vertices()):
    #         simples[i][i] = 1
    #     return all(self.euler_form(d,i) + self.euler_form(i,d) <= 0 for i in simples)

    def division_order(self, d, e):
        """Checks if d << e, which means that d_i <= e_i for every source i, d_j >= e_j for every sink j, and d_k == e_k for every vertex k which is neither a source nor a sink.

        # TODO: Think of a better name.
        # Good name?

        EXAMPLES

        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: d = vector([1,1])
        sage: e = vector([2,1])
        sage: f = vector([2,2])
        sage: Q.division_order(d,e)
        True
        sage: Q.division_order(e,d)
        False
        sage: Q.division_order(d,f)
        False
        sage: Q.division_order(f,d)
        False
        sage: Q.division_order(e,f)
        False
        sage: Q.division_order(f,e)
        True

        sage: Q = ThreeVertexQuiver(2,2,2)
        sage: Q
        An acyclic 3-vertex quiver; adjacency matrix:
        [0 2 2]
        [0 0 2]
        [0 0 0]
        sage: d = vector([1,1,1])
        sage: e = vector([1,2,1])
        sage: Q.division_order(d,e)
        False
        sage: Q.division_order(e,d)
        False
        """

        n = self.number_of_vertices()
        assert (d.length() == n) and (e.length() == n)
        less = all(
            [
                d[i - 1] <= e[i - 1]
                for i in list(filter(lambda i: self.is_source(i), range(1, n + 1)))
            ]
        )
        less = less and all(
            [
                d[j - 1] >= e[j - 1]
                for j in list(filter(lambda j: self.is_sink(j), range(1, n + 1)))
            ]
        )
        less = less and all(
            [
                d[k - 1] == e[k - 1]
                for k in list(
                    filter(
                        lambda k: (not self.is_source(k)) and (not self.is_sink(k)),
                        range(1, n + 1),
                    )
                )
            ]
        )

        return less

    """
    Generic subdimension vectors and generic Hom and Ext
    """

    # taken from code/snippets/canonical.sage
    # TODO still need testing code from there
    def is_generic_subdimension_vector(self, e, d):
        r"""Checks if e is a generic subdimension vector of d.

        INPUT:

        - ``e``: dimension vector for the subrepresentation

        - ``d``: dimension vector for the ambient representation

        OUTPUT: Statement truth value as Bool

        # Using notation from Section 5 of https://arxiv.org/pdf/0802.2147.pdf
        A dimension vector e is called a generic subdimension vector of d if a generic representation of dimension vector d possesses a subrepresentation of dimension vector e.
        By a result of Schofield (see Thm. 5.3 of https://arxiv.org/pdf/0802.2147.pdf) e is a generic subdimension vector of d if and only if e is a subdimension vector of d (missing in Thm. 5.3!) and <f,d-e> is non-negative for all generic subdimension vectors f of e.

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

        assert (
            self.number_of_vertices()
            == d.length() & self.number_of_vertices()
            == e.length()
        )
        assert all([di >= 0 for di in d.list()]) and all([ei >= 0 for ei in e.list()])

        if e == d:
            return True
        else:
            if not is_subdimension_vector(e, d):
                return False
            else:  # e is subdimension vector of d
                # List of all generic subdimension vectors of e
                genSubdims = self.all_generic_subdimension_vectors(e)
                return all([self.euler_form(f, d - e) >= 0 for f in genSubdims])

    def __all_generic_subdimension_vectors_helper(self, d):
        """Returns the list of lists of indexes of all generic subdimension vectors of e, where e ranges over all subdimension vectors of d. The index refers to the deglex order.

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
        subdims.sort(key=(lambda e: deglex_key(e, b=max(d) + 1)))
        # We use the deglex order because it's a total order which extends the usual entry-wise partial order on dimension vectors.
        N = len(subdims)

        # genIndexes[j] will in the end be the list of indexes (in subdims) of all generic subdimension vectors of subdims[j]
        genIndexes = [
            list(
                filter(
                    lambda i: is_subdimension_vector(subdims[i], subdims[j]), range(N)
                )
            )
            for j in range(N)
        ]

        for j in range(N):
            genIndexes[j] = list(
                filter(
                    lambda i: all(
                        [
                            self.euler_form(subdims[k], subdims[j] - subdims[i]) >= 0
                            for k in genIndexes[i]
                        ]
                    ),
                    genIndexes[j],
                )
            )

        genSubdims = [[subdims[i] for i in genIndexes[j]] for j in range(N)]

        return genIndexes, genSubdims

    def all_generic_subdimension_vectors(self, d):
        r"""Returns the list of all generic subdimension vectors of d.

        INPUT:

        - ``d``: dimension vector

        OUTPUT: list of vectors

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
        assert self.number_of_vertices() == d.length()
        assert all([di >= 0 for di in d.list()])

        genIndexes, genSubdims = self.__all_generic_subdimension_vectors_helper(d)
        N = len(genSubdims)
        return genSubdims[N - 1]

    def generic_ext(self, a, b):
        r"""Computes ext(a,b).

        INPUT:

        - ``a``: dimension vector

        - ``b``: dimension vector

        OUTPUT: Int

        According to Thm. 5.4 in Schofield's 'General representations of quivers', we have ext(a,b) = max{-<c,b> | c gen. subdimension vector of a}.

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
        return max([-self.euler_form(c, b) for c in genSubdims])

    def generic_hom(self, a, b):
        r"""Computes hom(a,b).

        INPUT:

        - ``a``: dimension vector

        - ``b``: dimension vector

        OUTPUT: Int

        There is a non-empty open subset U of R(Q,a) x R(Q,b) such that dim Ext(M,N) = ext(a,b), i.e. is minimal, for all (M,N) in U. Therefore dim Hom(M,N) = <a,b> + dim Ext(M,N) is minimal and therefore hom(a,b) = <a,b> + ext(a,b).

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

        return self.euler_form(a, b) + self.generic_ext(a, b)

    def generic_ext_vanishing(self, a, b):
        return self.is_generic_subdimension_vector(a, a + b)

    def generic_hom_vanishing(self, a, b):
        # TODO figure out a way to implement this.
        # How about this:
        return self.generic_hom(a, b) == 0

    def is_left_orthogonal(self, a, b):
        if self.generic_ext_vanishing(a, b):
            return self.euler_form(a, b) == 0
        else:
            return False

    """
    (Semi-)stability
    """

    def canonical_stability_parameter(self, d):
        """The canonical stability parameter is given by <d,_> - <_,d>"""
        E = self.euler_matrix()
        return d * (-self.euler_matrix().transpose() + E)

    def has_semistable_representation(self, d, theta):
        r"""Checks if there is a `\theta`-semistable representation of dimension vector `d`

        INPUT:
        - ``d``: dimension vector

        - ``theta``: stability parameter

        OUTPUT: Statement truth value as Bool

        See Thm. 5.4(1) of Reineke's overview paper https://arxiv.org/pdf/0802.2147.pdf: A dimension vector d admits a theta-semi-stable representation if and only if mu_theta(e) <= mu_theta(d) for all generic subdimension vectors e of d.
        # Thm. 5.4 in Markus's paper is actually a result of Schofield. So the algorithm should bear his name, if any.

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
        """Computes the list of indexes of all semistable subdimension vectors of d.

        EXAMPLES:
        sage: from quiver import *
        sage: Q, d, theta = GeneralizedKroneckerQuiver(1), vector([2,3]), vector([1,0])
        sage: i, s = Q._Quiver__all_semistable_subdimension_vectors_helper(d, theta); i, s
        ([1, 2, 3, 4, 5, 6, 10],
        [(0, 1), (1, 0), (0, 2), (1, 1), (2, 0), (0, 3), (2, 2)])

        """

        subdims = all_subdimension_vectors(d)
        subdims.sort(key=(lambda e: deglex_key(e, b=max(d) + 1)))
        # We use the deglex order because it's a total order which extends the usual entry-wise partial order on dimension vectors.
        N = len(subdims)
        genIndexes, genSubdims = self.__all_generic_subdimension_vectors_helper(d)
        sstIndexes = list(
            filter(
                lambda j: all(
                    [
                        slope(subdims[i], theta) <= slope(subdims[j], theta)
                        for i in list(filter(lambda i: i != 0, genIndexes[j]))
                    ]
                ),
                range(1, N),
            )
        )
        sstSubdims = [subdims[j] for j in sstIndexes]
        return sstIndexes, sstSubdims

    # TODO need to specify what the input for "stability parameter" is
    # TODO always have canonical stability as default?
    def has_stable_representation(self, d, theta, algorithm="schofield"):
        r"""Checks if there is a `\theta`-stable representation of this dimension vector.

        INPUT:
        - ``d``: dimension vector

        - ``theta``: stability parameter

        - ``algorithm``: String

        OUTPUT: statement truth value as Bool

        See Thm. 5.4(1) of Reineke's overview paper https://arxiv.org/pdf/0802.2147.pdf: A dimension vector d admits a theta-stable representation if and only if mu_theta(e) < mu_theta(d) for all proper generic subdimension vectors e of d.

        EXAMPLES:

        The `\mathrm{A}_2` quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(1)
            sage: theta = (1, -1)
            sage: Q.has_stable_representation([1, 1], theta, algorithm="schofield")
            True
            sage: Q.has_stable_representation([2, 2], theta, algorithm="schofield")
            False
            sage: Q.has_stable_representation([0, 0], theta, algorithm="schofield")
            False

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuivrr(3)
            sage: d = (2, 3)
            sage: theta = Q.canonical_stability(d)
            sage: Q.has_stable_representation(d, theta, algorithm="schofield")
            True

        """

        assert algorithm in ["schofield", "king", "al"]

        # coerce dimension vector and stability parameter
        d = vector(d)
        theta = vector(theta)

        assert d.length() == self.number_of_vertices() and theta.length() == self.number_of_vertices()

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

        if algorithm == "schofield":
            if d == self.zero_vector():
                return False
            else:
                genSubdims = self.all_generic_subdimension_vectors(d)
                genSubdims = list(
                    filter(lambda e: e != self.zero_vector() and e != d, genSubdims)
                )
                return all([slope(e, theta) < slope(d, theta) for e in genSubdims])

    def __all_stable_subdimension_vectors_helper(self, d, theta, denominator=sum):
        """Computes the list of all stable subdimension vectors of d which have the same slope as d.

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
        subdims.sort(key=(lambda e: deglex_key(e, b=max(d) + 1)))
        # We use the deglex order because it's a total order which extends the usual entry-wise partial order on dimension vectors.
        N = len(subdims)
        genIndexes, genSubdims = self.__all_generic_subdimension_vectors_helper(d)
        # slopeIndexes is the list of subdimension vectors of d of the same slope as d (in particular != 0)
        slopeIndexes = list(
            filter(
                lambda j: slope(subdims[j], theta, denominator=denominator)
                == slope(d, theta, denominator=denominator),
                range(1, N),
            )
        )
        # stIndexes contains all j for which subdims[j] is stable
        # e = subdims[j] is stable if for all generic subdimension vectors f = subdims[i] of e, it holds that slope(f) < slope(e)
        stIndexes = list(
            filter(
                lambda j: all(
                    [
                        slope(subdims[i], theta, denominator=denominator)
                        < slope(subdims[j], theta, denominator=denominator)
                        for i in list(
                            filter(lambda i: i != 0 and i != j, genIndexes[j])
                        )
                    ]
                ),
                slopeIndexes,
            )
        )
        stSubdims = [subdims[j] for j in stIndexes]
        return stIndexes, stSubdims

    """
    Canonical decomposition
    """

    # TODO this is an implementation detail, not something for the public interface
    def rearrange_dw_decomposition(self, decomposition, i, j):
        # this is non functional, let alone tested
        # apply Lemma 11.9.10 of Derksen--Weyman to rearrange the roots from i to j as k, k+1
        S = []
        for k in range(i, j - 1):
            # if k has a ``path to j'' satisfying the hypothesis of Lemma 11.9.10, we keep it
            # m is the j-k times j-k matrix with entry m[i,j] = 1 if !generic_hom_vanishing(Q,decomposition[j][2],decomposition[i][2]) and 0 otherwise
            m = matrix(j - k)
            for l in range(k, j - 1):
                for s in range(l + 1, j):
                    if not self.generic_hom_vanishing(
                        decomposition[s][2], decomposition[l][2]
                    ):
                        m[l - k, s - k] = 1
            paths = matrix(
                j - k
            )  # paths[i,j] is the number of paths from k + i - 1 to k + j - 1 of length at most j-k
            for l in range(j - k):
                paths = paths + m**l
            if paths[0, j - k - 1] > 0:
                S.append(k)
        rearrangement = [l for l in range(i + 1, j - 1) if l not in S]
        final = S + [i, j] + rearrangement
        decomposition_temp = [decomposition[l] for l in final]
        for l in range(i, j):
            decomposition[l] = decomposition_temp[l - i]
        return i + len(S) - 1  # this is the index of the element i now.

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
        # TODO this implementation needs to be documented (much better)
        if algorithm == "derksen-weyman":
            decomposition = [
                [d[i], self.simple_root(i + 1)]
                for i in range(self.number_of_vertices())
            ]
            while True:
                decomposition = list(filter(lambda root: root[0] > 0, decomposition))
                violating_pairs = [
                    (i, j)
                    for i in range(len(decomposition) - 1)
                    for j in range(i + 1, len(decomposition))
                    if self.euler_form(decomposition[j][1], decomposition[i][1]) < 0
                ]
                if not violating_pairs:
                    break
                violating_pairs = sorted(
                    violating_pairs, key=(lambda pair: pair[1] - pair[0])
                )
                i, j = violating_pairs[0]

                # this modifies the decomposition in place
                i = self.rearrange_dw_decomposition(decomposition, i, j)

                p, xi = decomposition[i]
                q, eta = decomposition[i + 1]
                xi_real = self.is_real_root(xi)
                eta_real = self.is_real_root(eta)
                zeta = p * xi + q * eta
                if xi_real and eta_real:
                    discriminant = self.euler_form(zeta, zeta)
                    if discriminant > 0:
                        pass  # TODO figure out what this should do
                    elif discriminant == 0:
                        zeta_prime = zeta // gcd(zeta)
                        del decomposition[i + 1]
                        decomposition[i] = [1, zeta_prime]
                    else:
                        del decomposition[i + 1]
                        decomposition[i] = [1, zeta]
                elif xi_real and not eta_real:
                    if p + q * self.euler_form(eta, xi) >= 0:
                        del decomposition[i + 1]
                        decomposition[i] = [1, eta - self.euler_form(eta, xi) * xi]
                    else:
                        del decomposition[i + 1]
                        decomposition[i] = [1, zeta]
                elif not xi_real and eta_real:
                    if q + p * self.euler_form(eta, xi) >= 0:
                        decomposition[i] = [1, eta]
                        decomposition[i + 1] = [1, xi - self.euler_form(eta, xi) * eta]
                    else:
                        del decomposition[i + 1]
                        decomposition[i] = [1, zeta]
                elif not xi_real and not eta_real:
                    del decomposition[i + 1]
                    decomposition[i] = [1, zeta]
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
            # TODO move this to docstring
            """I'm not sure if one of the Schofield algorithms is meant to be this one. But here's a very simple recursion which computes the canonical decomposition. It is based on Lem. 11.2.5 in Derksen--Weyman's book (a.k.a. 'the book with typos'):

            Lemma: Let a be a dimension vector and a = b+c a decomposition such that ext(b,c) = ext(c,b) = 0. If b = b_1 + ... + b_s and c = c_1 + ... + c_t are the canonical decompositions, then a = b_1 + ... + b_s + c_1 + ... + c_t is the canonical decomposition of a.

            If no non-trivial decomposition a = b+c as above exists, then a is a Schur root and therefore its own canonical decomposition. This is because a generic representation has no subdimension vector b which admits both a subrepresentation and a quotient representation. So a generic representation is indecomposable, which implies that a is a Schur root.

            EXAMPLES:

            # TODO make a smaller example?

            sage: Q = KroneckerQuiver()
            sage: ds = [vector([i,j]) for i in range(7) for j in range(7)]
            sage: can = [Q.canonical_decomposition(d, algorithm="recursive") for d in ds]
            sage: for i in range(7):
            ....:     for j in range(7):
            ....:         print('decomposition of '+str(ds[7*i+j])+' is: '+str(can[7*i+j]))
            ....:
            decomposition of (0, 0) is: [(0, 0)]
            decomposition of (0, 1) is: [(0, 1)]
            decomposition of (0, 2) is: [(0, 1), (0, 1)]
            decomposition of (0, 3) is: [(0, 1), (0, 1), (0, 1)]
            decomposition of (0, 4) is: [(0, 1), (0, 1), (0, 1), (0, 1)]
            decomposition of (0, 5) is: [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1)]
            decomposition of (0, 6) is: [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1)]
            decomposition of (1, 0) is: [(1, 0)]
            decomposition of (1, 1) is: [(1, 1)]
            decomposition of (1, 2) is: [(1, 2)]
            decomposition of (1, 3) is: [(0, 1), (1, 2)]
            decomposition of (1, 4) is: [(0, 1), (0, 1), (1, 2)]
            decomposition of (1, 5) is: [(0, 1), (0, 1), (0, 1), (1, 2)]
            decomposition of (1, 6) is: [(0, 1), (0, 1), (0, 1), (0, 1), (1, 2)]
            decomposition of (2, 0) is: [(1, 0), (1, 0)]
            decomposition of (2, 1) is: [(2, 1)]
            decomposition of (2, 2) is: [(1, 1), (1, 1)]
            decomposition of (2, 3) is: [(2, 3)]
            decomposition of (2, 4) is: [(1, 2), (1, 2)]
            decomposition of (2, 5) is: [(0, 1), (1, 2), (1, 2)]
            decomposition of (2, 6) is: [(0, 1), (0, 1), (1, 2), (1, 2)]
            decomposition of (3, 0) is: [(1, 0), (1, 0), (1, 0)]
            decomposition of (3, 1) is: [(1, 0), (2, 1)]
            decomposition of (3, 2) is: [(3, 2)]
            decomposition of (3, 3) is: [(1, 1), (1, 1), (1, 1)]
            decomposition of (3, 4) is: [(3, 4)]
            decomposition of (3, 5) is: [(1, 2), (2, 3)]
            decomposition of (3, 6) is: [(1, 2), (1, 2), (1, 2)]
            decomposition of (4, 0) is: [(1, 0), (1, 0), (1, 0), (1, 0)]
            decomposition of (4, 1) is: [(1, 0), (1, 0), (2, 1)]
            decomposition of (4, 2) is: [(2, 1), (2, 1)]
            decomposition of (4, 3) is: [(4, 3)]
            decomposition of (4, 4) is: [(1, 1), (1, 1), (1, 1), (1, 1)]
            decomposition of (4, 5) is: [(4, 5)]
            decomposition of (4, 6) is: [(2, 3), (2, 3)]
            decomposition of (5, 0) is: [(1, 0), (1, 0), (1, 0), (1, 0), (1, 0)]
            decomposition of (5, 1) is: [(1, 0), (1, 0), (1, 0), (2, 1)]
            decomposition of (5, 2) is: [(1, 0), (2, 1), (2, 1)]
            decomposition of (5, 3) is: [(2, 1), (3, 2)]
            decomposition of (5, 4) is: [(5, 4)]
            decomposition of (5, 5) is: [(1, 1), (1, 1), (1, 1), (1, 1), (1, 1)]
            decomposition of (5, 6) is: [(5, 6)]
            decomposition of (6, 0) is: [(1, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0)]
            decomposition of (6, 1) is: [(1, 0), (1, 0), (1, 0), (1, 0), (2, 1)]
            decomposition of (6, 2) is: [(1, 0), (1, 0), (2, 1), (2, 1)]
            decomposition of (6, 3) is: [(2, 1), (2, 1), (2, 1)]
            decomposition of (6, 4) is: [(3, 2), (3, 2)]
            decomposition of (6, 5) is: [(6, 5)]
            decomposition of (6, 6) is: [(1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1)]
            sage:
            sage: all([all([Q.generic_ext(s[i],s[j]) + Q.generic_ext(s[j],s[i]) == 0 for i in range(len(s)) for j in range(i)]) for s in can])
            True
            """

            genSubdims = self.all_generic_subdimension_vectors(d)
            genSubdims = list(
                filter(lambda e: e != self.zero_vector() and e != d, genSubdims)
            )
            for e in genSubdims:
                if d - e in genSubdims:
                    return self.canonical_decomposition(
                        e, algorithm="recursive"
                    ) + self.canonical_decomposition(d - e, algorithm="recursive")
            return [d]

        elif algorithm == "recursive_new":
            # TODO move this to docstring
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
            subdims.sort(key=(lambda e: deglex_key(e, b=max(d) + 1)))
            N = len(subdims)

            idx_diff = lambda j, i: subdims.index(subdims[j] - subdims[i])

            genIndexes, genSubdims = self.__all_generic_subdimension_vectors_helper(d)

            def canon_indexes(j):
                """Computes for j in range(N) the list of indexes in subdims for the canonical decomposition of subdims[j]"""
                for i in list(filter(lambda i: i != 0 and i != j, genIndexes[j])):
                    k = idx_diff(j, i)
                    if k in genIndexes[j]:
                        return canon_indexes(i) + canon_indexes(k)
                return [j]

            return [subdims[i] for i in canon_indexes(N - 1)]

    """
    Nilpotent cone and Hesselink
    """

    def dimension_nullcone(self, d):
        r"""Returns the dimension of the nullcone which is the set of all nilpotent representations.

        INPUT:
        - ``d``: vector of Ints

        OUTPUT: dimension as Int
        """

        if self.is_acyclic():
            return d.transpose() * self.adjacency_matrix() * d
        else:
            # TODO where is the algorithm?
            raise NotImplementedError()

    """
    Teleman!
    """
    # TODO: This section should go into QuiverModuliSpace, I think.

    def all_weight_bounds(self, d, theta, denominator=sum):
        """
        Returns, for a given dimension vector d and a given stability parameter theta, the list of all weights to apply Teleman quantization.
        For each HN type, the 1-PS lambda acts on det(N_{S/R}|_Z) with a certain weight. Teleman quantization gives a numerical condition involving these weights to compute cohmology on the quotient.
        """
        # TODO return the Hn type as well?

        # This is only relevant on the unstable locus
        HN = list(
            filter(
                lambda hntype: hntype != [d],
                self.all_harder_narasimhan_types(d, theta, denominator=denominator),
            )
        )

        return list(
            map(
                lambda hntype: -sum(
                    [
                        (
                            slope(hntype[s], theta, denominator=denominator)
                            - slope(hntype[t], theta, denominator=denominator)
                        )
                        * self.euler_form(hntype[s], hntype[t])
                        for s in range(len(hntype) - 1)
                        for t in range(s + 1, len(hntype))
                    ]
                ),
                HN,
            )
        )

    def does_rigidity_inequality_hold(self, d, theta, denominator=sum):
        r"""
        Returns True if the rigidity inequality holds for d and theta, i.e. if the weights of the 1-PS lambda on $\det(N_{S/R}|_Z)$ for each HN type are all strictly larger than the weights of the tensors of the universal bundles $U_i^\vee \otimes U_j$.
        """

        # This is only relevant on the unstable locus
        HN = list(
            filter(
                lambda hntype: hntype != [d],
                self.all_harder_narasimhan_types(d, theta, denominator=denominator),
            )
        )

        # We compute the weights of the 1-PS lambda on det(N_{S/R}|_Z) for each HN type
        weights = list(
            map(
                lambda hntype: -sum(
                    [
                        (
                            slope(hntype[s], theta, denominator=denominator)
                            - slope(hntype[t], theta, denominator=denominator)
                        )
                        * self.euler_form(hntype[s], hntype[t])
                        for s in range(len(hntype) - 1)
                        for t in range(s + 1, len(hntype))
                    ]
                ),
                HN,
            )
        )

        # We compute the maximum weight of the tensors of the universal bundles U_i^\vee \otimes U_j
        tensorWeights = list(
            map(
                lambda hntype: slope(hntype[0], theta, denominator=denominator)
                - slope(hntype[-1], theta, denominator=denominator),
                HN,
            )
        )

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
        return (
            1
            - self.number_of_vertices()
            + sum(len(G.all_paths(a[0], a[1], use_multiedges=True)) for a in G.edges())
        )

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
    assert d.length() == theta.length()
    assert denominator(d) > 0
    return (theta * d) / (denominator(d))


def is_subdimension_vector(e, d):
    assert e.length() == d.length()
    n = e.length()
    return all([e[i] <= d[i] for i in range(n)])


def deglex_key(e, b):
    """A function which satisfies e <_{deglex} d iff deglex_key(e) < deglex_key(d), provided that b >> 0."""
    n = e.length()
    return sum([e[i] * b ** (n - i - 1) for i in range(n)]) + sum(e) * b**n


def all_subdimension_vectors(d):
    """Returns the list of all subdimension vectors of d."""
    return list(map(vector, cartesian_product([range(di + 1) for di in d])))


# TODO: This method has a stupid name (my own fault). Think of a better one.
def is_coprime_for_stability_parameter(d, theta):
    """Checks if d is theta-coprime.

    A dimension vector d is theta-coprime if mu_theta(e) != mu_theta(e) for all proper subdimension vectors e of d.

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

    assert d.length() == theta.length()
    zeroVector = vector([0 for i in range(d.length())])
    properSubdimensions = list(
        filter(lambda e: e != d and e != zeroVector, all_subdimension_vectors(d))
    )
    return all([slope(d, theta) != slope(e, theta) for e in properSubdimensions])


def is_indivisible(d):
    """Checks if the gcd of all entries is 1 or not."""
    return gcd(d) == 1


"""Class methods"""
