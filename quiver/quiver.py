from sage.arith.misc import gcd
from sage.categories.cartesian_product import cartesian_product
from sage.graphs.digraph import DiGraph
from sage.graphs.graph import Graph
from sage.matrix.constructor import matrix
from sage.matrix.special import zero_matrix
from sage.modules.free_module_element import vector, zero_vector
from sage.rings.integer_ring import ZZ
from sage.structure.element import Element


class Quiver(Element):
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

    def __init__(self, G, name=None):
        r"""Constructor for a quiver.

        This takes a directed graph as input. If it is not a DiGraph instance,
        we interpret it as an adjacency matrix.
        For other constructions, see

        - meth:`Quiver.from_digraph`
        - meth:`Quiver.from_matrix`
        - meth:`Quiver.from_string`

        INPUT:

        - ``G`` -- directed graph

        - ``name`` -- optional name for the quiver

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = Quiver([[0, 3], [0, 0]]); Q
            a quiver with 2 vertices and 3 arrows

        A Dynkin quiver of type A_3::

            sage: Q = Quiver.from_string("1-2-3"); Q.adjacency_matrix()
            [0 1 0]
            [0 0 1]
            [0 0 0]

        A triangle-shaped quiver::

            sage: Q = Quiver.from_string("1-2-3, 1-3"); Q.adjacency_matrix()
            [0 1 1]
            [0 0 1]
            [0 0 0]

        """

        if isinstance(G, DiGraph):
            self.__G = G
        else:
            self.__G = DiGraph(matrix(G))

        # if name is None this doesn't do anything
        self.rename(name)

        # for caching purposes
        self.__M = self.__G.adjacency_matrix()

    @classmethod
    def from_digraph(cls, G, name=None):
        r"""Construct a quiver a graph

        INPUT:

        - ``G`` -- directed graph

        - ``name`` -- optional name for the quiver

        OUTPUT: the quiver

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: M = [[0, 3], [0, 0]]
            sage: Quiver.from_digraph(DiGraph(matrix(M))) == Quiver.from_matrix(M)
            True

        """
        return cls(G.adjacency_matrix(), name)

    @classmethod
    def from_matrix(cls, M, name=None):
        r"""Construct a quiver from its adjacency matrix

        INPUT:

        - ``M`` -- adjacency matrix of the quiver

        - ``name`` -- optional name for the quiver

        OUTPUT: the quiver

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = Quiver.from_matrix([[0, 3], [0, 0]]); Q.adjacency_matrix()
            [0 3]
            [0 0]

        """
        return cls(DiGraph(matrix(M)), name)

    @classmethod
    def from_string(cls, Q: str, name=None):
        r"""Construct a quiver from a comma-separated list of chains of the form "i-j-k-...

        You specify an arrow from `i` to `j` by writing `i-j`.
        A multiple arrow is specified by repeating the `i`, so that `1--2` is the Kronecker quiver.
        If you write `i-j-k` then you have 1 arrow from `i` to `j` and one from `j` to `k`.
        The full quiver is specified by concatenating (multiple) arrows by commas.

        OUTPUT: the quiver

        INPUT:

        - ``Q`` -- string; a string of the format described above giving a quiver

        - ``name`` -- optional name for the quiver

        EXAMPLES:

        The 3-Kronecker quiver defined in two different ways::

            sage: from quiver import *
            sage: Quiver.from_matrix([[0, 3], [0, 0]]) == Quiver.from_string("1---2")
            True

        A more complicated example::

            sage: Quiver.from_string("1--2-3,1---3,3-1").adjacency_matrix()
            [0 2 3]
            [0 0 1]
            [1 0 0]

        The actual numbers we use don't matter::

            sage: from quiver import *
            sage: Quiver.from_matrix([[0, 3], [0, 0]]) == Quiver.from_string("12---23")
            True

        """
        # remove all whitespace from the string
        Q = "".join(Q.split())

        # determine the vertices used in the string
        vertices = list(
            set(
                [
                    int(vertex)
                    for chain in Q.split(",")
                    for vertex in chain.split("-")
                    if vertex
                ]
            )
        )

        # adjacency matrix
        M = zero_matrix(len(vertices))

        for chain in Q.split(","):
            pieces = chain.split("-")

            source = vertices.index(int(pieces[0]))
            number = 1

            for piece in pieces[1:]:
                # if the string is empty we increase the number of arrows counter
                if not piece:
                    number += 1
                # if the string is non-empty we treat it as an integer
                # this means we add the appropriate number of arrows
                # and make the target the source to start the process again
                if piece:
                    target = vertices.index(int(piece))
                    M[source, target] += number

                    number = 1
                    source, target = target, None

        return cls.from_matrix(M, name)

    def __repr__(self) -> str:
        if self.get_custom_name():
            return self.get_custom_name()
        else:
            return "a quiver with {} vertices and {} arrows".format(
                self.adjacency_matrix().nrows(), sum(sum(self.adjacency_matrix()))
            )

    def __str__(self) -> str:
        return "{}\nadjacency matrix:\n{}".format(self.repr(), self.adjacency_matrix())

    def repr(self) -> str:
        r"""
        Basic description of the quiver

        To override the output, one uses `Quiver.rename` from the `Element` class.
        The output of `Quiver.repr` is that of `Quiver.get_custom_name` if it is set,
        else it is the default specifying the number of vertices and arrows.

        OUTPUT: a basic description of the quiver

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = Quiver.from_string("1---2"); Q
            a quiver with 2 vertices and 3 arrows
            sage: Q.rename("3-Kronecker quiver"); Q
            3-Kronecker quiver

        Renaming and resetting the name::

            sage: Q = Quiver.from_string("1---2")
            sage: Q.get_custom_name()

            sage: Q.rename("3-Kronecker quiver")
            sage: Q.get_custom_name()
            '3-Kronecker quiver'
            sage: Q.reset_name()
            sage: Q.get_custom_name()

            sage: Q
            a quiver with 2 vertices and 3 arrows

        """
        return self.__repr__()

    def str(self) -> str:
        r"""
        Full description of the quiver

        This combines the output of :meth:`Quiver.repr` with the adjacency matrix.

        OUTPUT: a complete description of the quiver

        EXAMPLES::

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = Quiver.from_string("1---2"); print(Q)
            a quiver with 2 vertices and 3 arrows
            adjacency matrix:
            [0 3]
            [0 0]
            sage: Q.rename("3-Kronecker quiver"); print(Q)
            3-Kronecker quiver
            adjacency matrix:
            [0 3]
            [0 0]

        """
        return self.__str__()

    def __eq__(self, other) -> bool:
        r"""
        Checks for equality of quivers.

        Equality here refers to equality of adjacency matrices,
        but disregarding the name of the quiver.

        INPUT:

        - ``other`` -- Quiver; the quiver to compare against

        OUTPUT: whether the adjacency matrices are the same

        EXAMPLES:

        The 2-Kronecker quiver and the generalized Kronecker quiver are the same::

            sage: from quiver import *
            sage: KroneckerQuiver() == GeneralizedKroneckerQuiver(2)
            True

        """
        return self.adjacency_matrix() == other.adjacency_matrix()

    def _coerce_dimension_vector(self, d):
        r"""
        Coerces `d` to be a dimension vector of the quiver

        It must be a data structure that is indexed by the vertices of the quiver,
        so most likely a dict, list, or vector.
        If it is a list it is coerced to a vector.

        As a consistency check we verify whether the length is the number of vertices.

        INPUT:

        - ``d``: a candidate dimension vector

        OUTPUT: either a dict or vector

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q._coerce_dimension_vector([1, 2])
            (1, 2)
            sage: Q._coerce_dimension_vector([1, 2, 3, 4])
            Traceback (most recent call last):
            ...
            ValueError: The input is not an element of `\mathbb{Z}Q_0`.
            sage: Q._coerce_dimension_vector([1, -3])
            Traceback (most recent call last):
            ...
            ValueError: The input is not a dimension vector of the quiver.

        """
        d = self._coerce_vector(d)
        if all(di >= 0 for di in d):
            return d
        else:
            raise ValueError("The input is not a dimension vector of the quiver.")

    def _coerce_vector(self, x):
        r"""
        Coerces `x` to be a vector in `\mathbb{Z}Q_0`.

        It raises a `ValueError` if it is not a list of integers of length the number
        of vertices in the quiver.

        INPUT:

        - ``x``: a list of integers

        OUTPUT: a Sage vector if `x` is an element of `\mathbb{Z}Q_0`

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q._coerce_vector([-1, 2])
            (-1, 2)
            sage: Q._coerce_vector([1, 2, 3, 4])
            Traceback (most recent call last):
            ...
            ValueError: The input is not an element of `\mathbb{Z}Q_0`.

        """
        if isinstance(x, list) or isinstance(x, tuple):
            x = vector(x)

        if len(x) == self.number_of_vertices():
            return x
        else:
            raise ValueError("The input is not an element of `\mathbb{Z}Q_0`.")

    """
    Basic graph-theoretic properties of the quiver
    """

    def adjacency_matrix(self):
        r"""
        Returns the adjacency matrix of the quiver.

        OUTPUT: The square matrix `M` whose entry `M[i,j]` is the number of arrows
        from the vertex `i` to the vertex `j`

        EXAMPLES::

        The adjacency matrix of a quiver construct from an adjacency matrix::

            sage: from quiver import *
            sage: M = matrix([[0, 3], [0, 0]])
            sage: M == Quiver(M).adjacency_matrix()
            True

        """
        return self.__M

    def graph(self):
        r"""
        Return the underlying graph of the quiver

        OUTPUT: the underlying quiver as a DiGraph object

        EXAMPLES:

        The underlying graph of the quiver from a directed graph is that graph::

            sage: from quiver import *
            sage: G = DiGraph(matrix([[0, 3], [0, 0]]))
            sage: G == Quiver.from_digraph(G).graph()
            True

        """
        return self.__G

    def vertices(self):
        r"""
        Return the vertices of the quiver

        If the quiver is created from a DiGraph or string, the vertices are labelled
        using the data in the DiGraph of string, as explained in :meth:`Quiver.from_graph`
        or :meth:`Quiver.from_string`.
        If the quiver is created from a matrix, the vertices are labelled from `0`
        to `n-1`, where `n` is the number of rows or columns in the matrix.

        OUTPUT: the vertices in the underlying graph

        EXAMPLES:

        Usually the vertices will be just integers::

            sage: from quiver import *
            sage: Quiver([[0, 3], [0, 0]]).vertices()
            [0, 1]

        We can have non-trivial labels for a quiver::

            sage: Quiver.from_string("foo---bar", forget_labels=False).vertices()
            ['foo', 'bar']

        """
        return self.graph().vertices(sort=False)

    def number_of_vertices(self) -> int:
        r"""Returns the number of vertices

        OUTPUT: the number of vertices

        EXAMPLES:

        There are 3 vertices in a 3-vertex quiver::

            sage: from quiver import *
            sage: ThreeVertexQuiver(1, 2, 4).number_of_vertices()
            3

        """
        return self.graph().order()

    def number_of_arrows(self) -> int:
        r"""Returns the number of arrows

        OUTPUT: the number of arrows

        EXAMPLES:

        There are 7 vertices in this 3-vertex quiver::

            sage: from quiver import *
            sage: ThreeVertexQuiver(1, 2, 4).number_of_arrows()
            7

        """
        return self.graph().size()

    def is_acyclic(self) -> bool:
        r"""Returns whether the quiver is acyclic.

        OUTPUT: True if the quiver is acyclic, False otherwise.

        EXAMPLES:

        An acyclic graph::

            sage: from quiver import *
            sage: KroneckerQuiver(3).is_acyclic()
            True

        A non-acyclic graph::

            sage: GeneralizedJordanQuiver(5).is_acyclic()
            False

        """
        return self.graph().is_directed_acyclic()

    def is_connected(self) -> bool:
        r"""Returns whether the underlying graph of the quiver is connected or not.

        OUTPUT: True if the quiver is connected, False otherwise.

        EXAMPLES:

        The n-Kronecker quivers are connected::

            sage: from quiver import *
            sage: KroneckerQuiver(4).is_connected()
            True

        The loop quivers are connected::

            sage: GeneralizedJordanQuiver(3).is_connected()
            True
        """
        return self.graph().is_connected()

    """
    Some graph-theoretic properties of the quiver
    """

    def in_degree(self, i):
        r"""Returns the in-degree of a vertex.

        The parameter `i` must be an element of the vertices of the underlying graph.
        If constructed from a matrix or string, one has that `i` can go from `0` to `n-1`
        where `n` is the number of vertices in the graph.

        The indegree of `i` is the number of incoming arrows at `i`.

        INPUT:

        ``i`` -- a vertex of the underlying graph

        OUTPUT: The indegree of the vertex `i`

        EXAMPLES:

        In the 3-Kronecker quiver the indegree is either 0 or 3::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.in_degree(0)
            0
            sage: Q.in_degree(1)
            3

        """
        return self.graph().in_degree(i)

    def out_degree(self, i):
        r"""Returns the out-degree of a vertex.

        The parameter `i` must be an element of the vertices of the underlying graph.
        If constructed from a matrix or string, one has that `i` can go from `0` to `n-1`
        where `n` is the number of vertices in the graph.

        The out-degree of `i` is the number of outgoing arrows at `i`.

        ``i`` -- a vertex of the underlying graph

        OUTPUT: The out-degree of the vertex `i`

        EXAMPLES:

        In the 3-Kronecker quiver the outdegree is either 3 or 0::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.out_degree(0)
            3
            sage: Q.out_degree(1)
            0

        """
        return self.graph().out_degree(i)

    def is_source(self, i) -> bool:
        """Checks if `i` is a source of the quiver, i.e. if there are no incoming arrows at `i`.

        EXAMPLES:

        The 3-Kronecker quiver has one source::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.is_source(0)
            True
            sage: Q.is_source(1)
            False

        """
        return self.in_degree(i) == 0

    def is_sink(self, i) -> bool:
        """Checks if `i` is a sink of the quiver, i.e. if there are no outgoing arrows out of `i`.

        EXAMPLES

        The 3-Kronecker quiver has one sink::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.is_sink(0)
            False
            sage: Q.is_sink(1)
            True

        """
        return self.out_degree(i) == 0

    """
    Basic representation-theoretical properties of the quiver
    """

    def euler_matrix(self):
        r"""Returns the Euler matrix of the quiver.

        OUTPUT: Sage matrix.
        """
        return matrix.identity(self.number_of_vertices()) - self.adjacency_matrix()

    def euler_form(self, x, y) -> int:
        r"""The Euler bilinear form of the quiver.

        INPUT:
        - ``x`` -- vector of integers
        - ``y`` -- vector of integers

        OUTPUT: the multiplication of ``x * self.euler_matrix() * y`` as an  Int.

        """
        # TODO coercion
        assert (
            x.length() == self.number_of_vertices()
            and y.length() == self.number_of_vertices()
        )
        return x * self.euler_matrix() * y

    def symmetrized_euler_form(self, x, y) -> int:
        r"""The symmetrization of the Euler bilinear form of the quiver.

        INPUT:
        - ``x`` -- vector of integers
        - ``y`` -- vector of integers

        OUTPUT: the sum ``self.euler_form(x,y) + self.euler_form(y,x)`` as an  Int.

        """
        # TODO coercion
        assert (
            x.length() == self.number_of_vertices()
            and y.length() == self.number_of_vertices()
        )
        return self.euler_form(x, y) + self.euler_form(y, x)

    def tits_form(self, x) -> int:
        r"""The Tits quadratic form of the quiver.

        INPUT:
        - ``x`` -- vector of integers

        OUTPUT: the expression ``self.euler_form(x,x)`` as an  Int.

        """
        # TODO coercion
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

        if self.get_custom_name():
            name = "opposite of " + self.get_custom_name()
        else:
            name = None

        return Quiver(A, name)

    def double_quiver(self):
        """The adjacency matrix of the double of a quiver is the sum of the adjacency matrix of the original quiver and its transpose."""
        A = self.adjacency_matrix() + self.adjacency_matrix().transpose()
        if self.get_custom_name():
            name = "double of " + self.get_custom_name()
        else:
            name = None
        return Quiver(A, name)

    # TODO optional parameter for name of the vertex
    def framed_quiver(self, f):
        r"""Returns the framed quiver with framing vector f.

        INPUT:
        - ``f``: vector of Ints

        OUTPUT: Quiver object

        The framed quiver has one additional vertex 0 and f_i many arrows from 0 to i.
        """
        # TODO tests are missing

        n = self.number_of_vertices()
        assert f.length() == n
        A = self.adjacency_matrix()
        # Adjacency matrix of the framed quiver looks like this (block shape):
        # [[0 f]
        #  [0 A]]
        # Add f as a first row
        A = A.insert_row(0, f)
        # Add a zero column
        A = A.transpose().insert_row(0, zero_vector(n + 1)).transpose()
        return Quiver(A)

    def coframed_quiver(self, f):
        r"""Returns the coframed quiver with framing vector f.

        INPUT:
        - ``f``: vector of Ints

        OUTPUT: Quiver object

        The coframed quiver has one additional vertex oo and f_i many arrows from i to oo.
        """
        n = self.number_of_vertices()
        assert f.length() == n
        A = self.adjacency_matrix()
        # Adjacency matrix of the coframed quiver looks like this (block shape):
        # [[A f]
        #  [0 0]]
        # Add f as a last column
        A = A.transpose().insert_row(n, f).transpose()
        # Add a zero row as last row
        A = A.insert_row(n, zero_vector(n + 1))
        return Quiver(A)

    def full_subquiver(self, I):
        r"""Returns the full subquiver supported on the given set of vertices.

        INPUT:
        - ``I``: List

        OUTPUT: Quiver object

        EXAMPLES:

        The support is the set {i in Q_0 | d_i > 0}.::

            sage: from quiver import *
            sage: Q = ThreeVertexQuiver(2, 3, 4); Q.adjacency_matrix()
            [0 2 3]
            [0 0 4]
            [0 0 0]
            sage: Q.full_subquiver([1, 2]).adjacency_matrix()
            [0 2]
            [0 0]
            sage: Q.full_subquiver([1, 3]).adjacency_matrix()
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

        The thin dimension vector is [1,...,1].

        OUTPUT: vector of Ints
        """
        return vector([1 for i in range(self.number_of_vertices())])

    def simple_root(self, i):
        r"""Returns the simple root at the vertex.

        The simple root at i is e_i = [0,...,1,...,0], i.e. the unit vector with a one in position i.

        INPUT:
        - ``i``: Int

        OUTPUT: vector of Ints
        """
        n = self.number_of_vertices()
        # Our convention is that vertices are numbered 1,...,n
        # TODO no it shouldn't
        assert i >= 1 and i <= n
        ei = vector([0 for i in range(n)])
        ei[i - 1] = 1
        return ei

    def is_root(self, x):
        r"""Checks if x is a root of the underlying diagram of the quiver.

        A root is a non-zero vector in Z^n such that the Tits form of x is <= 1.

        INPUT:
        - ``x``: vector of Ints

        OUTPUT: statement truth value as bool
        """
        x = self._coerce_vector(x)

        # TODO any check like e == zero_vector should really be replaced by a all(ei == 0 for ei in e) for performance
        return x != self.zero_vector() and self.tits_form(x) <= 1

    def is_real_root(self, x):
        r"""Checks if x is a real root of the underlying diagram of the quiver.

        A root is called real if its Tits form equals 1.

        INPUT:
        - ``x``: vector of Ints

        OUTPUT: statement truth value as bool
        """
        x = self._coerce_vector(x)
        return self.tits_form(x) == 1

    def is_imaginary_root(self, x):
        r"""Checks if x is an imaginary root of the underlying diagram of the quiver.

        A root is called imaginary if its Tits form is non-positive.

        INPUT:
        - ``x``: vector of Ints

        OUTPUT: statement truth value as bool
        """
        x = self._coerce_vector(x)
        return x != self.zero_vector() and self.tits_form(x) <= 0

    def is_schur_root(self, d):
        r"""Checks if d is a Schur root.

        INPUT:
        - ``d``: vector of Ints

        OUTPUT: statement truth value as Bool

        A Schur root is a dimension vector which admits a Schurian representation, i.e. a representation whose endomorphism ring is k. It's necessarily indecomposable.
        By a result of Schofield (https://mathscinet.ams.org/mathscinet/relay-station?mr=1162487) d is a Schur root if and only if d admits a stable representation for the canonical stability parameter.

        EXAMPLES:

        The root (2, 3) is Schurian for the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.is_schur_root((2, 3))
            True

        Examples from Derksen--Weyman's book (Example 11.1.4)::

            sage: from quiver import *
            sage: Q = ThreeVertexQuiver(1, 1, 1)
            sage: Q.is_schur_root((1, 1, 2))
            True
            sage: Q.is_schur_root((1, 2, 1))
            False
            sage: Q.is_schur_root((1, 1, 1))
            True
            sage: Q.is_schur_root((2, 2, 2))
            False

        """
        d = self._coerce_dimension_vector(d)
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

        The support is the set of vertices for which the dimension vector is nonzero::

            sage: from quiver import *
            sage: Q = ThreeVertexQuiver(2, 0, 4)
            sage: d = vector([1, 1, 1])
            sage: Q.support(d)
            [1, 2, 3]
            sage: d = vector([1, 0, 1])
            sage: Q.support(d)
            [1, 3]

        """
        d = self._coerce_dimension_vector(d)
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

        EXAMPLES:

        The fundamental domain of the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.in_fundamental_domain((1, 1))
            True
            sage: Q.in_fundamental_domain((1, 2))
            False
            sage: Q.in_fundamental_domain((2, 3))
            True

        """
        d = self._coerce_dimension_vector(d)

        # check if `\langle d,e_i\rangle + \langle e_i,d\rangle \leq 0`
        # for all vertices `i\in Q_0`
        inequality = all(
            [
                (
                    self.euler_form(d, self.simple_root(i + 1))
                    + self.euler_form(self.simple_root(i + 1), d)
                    <= 0
                )
                for i in range(self.number_of_vertices())
            ]
        )

        # check if the support is connected
        connected = self.full_subquiver(self.support(d)).is_connected()

        return inequality and connected

    def division_order(self, d, e):
        """Checks if d << e, which means that d_i <= e_i for every source i, d_j >= e_j for every sink j, and d_k == e_k for every vertex k which is neither a source nor a sink.

        # TODO: Think of a better name.
        # Good name?

        EXAMPLES

        Order on some dimension vectors for the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: d = [1, 1]
            sage: e = [2, 1]
            sage: f = [2, 2]
            sage: Q.division_order(d, e)
            True
            sage: Q.division_order(e, d)
            False
            sage: Q.division_order(d, f)
            False
            sage: Q.division_order(f, d)
            False
            sage: Q.division_order(e, f)
            False
            sage: Q.division_order(f, e)
            True

        Order on some dimension vectors for a 3-vertex quiver::

            sage: Q = ThreeVertexQuiver(2, 2, 2)
            sage: d = [1, 1, 1]
            sage: e = [1, 2, 1]
            sage: Q.division_order(d, e)
            False
            sage: Q.division_order(e, d)
            False
        """
        d = self._coerce_dimension_vector(d)
        e = self._coerce_dimension_vector(e)
        n = self.number_of_vertices()

        # TODO isn't this redundant? coerce dimension vector already checks this
        assert (d.length() == n) and (e.length() == n)

        # TODO instead of range(n) we need to iterate over the vertices!
        # TODO indexing dimension vectors by vertices requires some additional care...
        # TODO make this cleaner
        less = all(
            [d[i] <= e[i] for i in list(filter(lambda i: self.is_source(i), range(n)))]
        )
        less = less and all(
            [d[j] >= e[j] for j in list(filter(lambda j: self.is_sink(j), range(n)))]
        )
        less = less and all(
            [
                d[k] == e[k]
                for k in list(
                    filter(
                        lambda k: (not self.is_source(k)) and (not self.is_sink(k)),
                        range(n),
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

        OUTPUT: True if e is a generic subdimension vector of d, False otherwise.

        # Using notation from Section 5 of https://arxiv.org/pdf/0802.2147.pdf
        A dimension vector e is called a generic subdimension vector of d if a generic representation of dimension vector d possesses a subrepresentation of dimension vector e.
        By a result of Schofield (see Thm. 5.3 of https://arxiv.org/pdf/0802.2147.pdf) e is a generic subdimension vector of d if and only if e is a subdimension vector of d (missing in Thm. 5.3!) and <f,d-e> is non-negative for all generic subdimension vectors f of e.

        EXAMPLES:

        Some examples on loop quivers::

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


        Some n-Kronecker quivers::

            sage: Q = GeneralizedKroneckerQuiver(1)
            sage: dims = Tuples(range(3), 2)
            sage: for e in dims:
            ....:     for d in dims:
            ....:         if is_subdimension_vector(e,d):
            ....:             print(str(e)+" gen. subdim of "+str(d)+"?: "+str(Q.is_generic_subdimension_vector(e,d)))
            (0, 0) gen. subdim of (0, 0)?: True
            (0, 0) gen. subdim of (1, 0)?: True
            (0, 0) gen. subdim of (2, 0)?: True
            (0, 0) gen. subdim of (0, 1)?: True
            (0, 0) gen. subdim of (1, 1)?: True
            (0, 0) gen. subdim of (2, 1)?: True
            (0, 0) gen. subdim of (0, 2)?: True
            (0, 0) gen. subdim of (1, 2)?: True
            (0, 0) gen. subdim of (2, 2)?: True
            (1, 0) gen. subdim of (1, 0)?: True
            (1, 0) gen. subdim of (2, 0)?: True
            (1, 0) gen. subdim of (1, 1)?: False
            (1, 0) gen. subdim of (2, 1)?: True
            (1, 0) gen. subdim of (1, 2)?: False
            (1, 0) gen. subdim of (2, 2)?: False
            (2, 0) gen. subdim of (2, 0)?: True
            (2, 0) gen. subdim of (2, 1)?: False
            (2, 0) gen. subdim of (2, 2)?: False
            (0, 1) gen. subdim of (0, 1)?: True
            (0, 1) gen. subdim of (1, 1)?: True
            (0, 1) gen. subdim of (2, 1)?: True
            (0, 1) gen. subdim of (0, 2)?: True
            (0, 1) gen. subdim of (1, 2)?: True
            (0, 1) gen. subdim of (2, 2)?: True
            (1, 1) gen. subdim of (1, 1)?: True
            (1, 1) gen. subdim of (2, 1)?: True
            (1, 1) gen. subdim of (1, 2)?: True
            (1, 1) gen. subdim of (2, 2)?: True
            (2, 1) gen. subdim of (2, 1)?: True
            (2, 1) gen. subdim of (2, 2)?: False
            (0, 2) gen. subdim of (0, 2)?: True
            (0, 2) gen. subdim of (1, 2)?: True
            (0, 2) gen. subdim of (2, 2)?: True
            (1, 2) gen. subdim of (1, 2)?: True
            (1, 2) gen. subdim of (2, 2)?: True
            (2, 2) gen. subdim of (2, 2)?: True
            sage: Q = GeneralizedKroneckerQuiver(2)
            sage: for e in dims:
            ....:     for d in dims:
            ....:         if is_subdimension_vector(e,d):
            ....:             print(str(e)+" gen. subdim of "+str(d)+"?: "+str(Q.is_generic_subdimension_vector(e,d)))
            ....:
            (0, 0) gen. subdim of (0, 0)?: True
            (0, 0) gen. subdim of (1, 0)?: True
            (0, 0) gen. subdim of (2, 0)?: True
            (0, 0) gen. subdim of (0, 1)?: True
            (0, 0) gen. subdim of (1, 1)?: True
            (0, 0) gen. subdim of (2, 1)?: True
            (0, 0) gen. subdim of (0, 2)?: True
            (0, 0) gen. subdim of (1, 2)?: True
            (0, 0) gen. subdim of (2, 2)?: True
            (1, 0) gen. subdim of (1, 0)?: True
            (1, 0) gen. subdim of (2, 0)?: True
            (1, 0) gen. subdim of (1, 1)?: False
            (1, 0) gen. subdim of (2, 1)?: False
            (1, 0) gen. subdim of (1, 2)?: False
            (1, 0) gen. subdim of (2, 2)?: False
            (2, 0) gen. subdim of (2, 0)?: True
            (2, 0) gen. subdim of (2, 1)?: False
            (2, 0) gen. subdim of (2, 2)?: False
            (0, 1) gen. subdim of (0, 1)?: True
            (0, 1) gen. subdim of (1, 1)?: True
            (0, 1) gen. subdim of (2, 1)?: True
            (0, 1) gen. subdim of (0, 2)?: True
            (0, 1) gen. subdim of (1, 2)?: True
            (0, 1) gen. subdim of (2, 2)?: True
            (1, 1) gen. subdim of (1, 1)?: True
            (1, 1) gen. subdim of (2, 1)?: True
            (1, 1) gen. subdim of (1, 2)?: False
            (1, 1) gen. subdim of (2, 2)?: True
            (2, 1) gen. subdim of (2, 1)?: True
            (2, 1) gen. subdim of (2, 2)?: False
            (0, 2) gen. subdim of (0, 2)?: True
            (0, 2) gen. subdim of (1, 2)?: True
            (0, 2) gen. subdim of (2, 2)?: True
            (1, 2) gen. subdim of (1, 2)?: True
            (1, 2) gen. subdim of (2, 2)?: True
            (2, 2) gen. subdim of (2, 2)?: True

        """
        d = self._coerce_dimension_vector(d)
        e = self._coerce_dimension_vector(e)

        if e == d:
            return True
        else:
            if not is_subdimension_vector(e, d):
                return False
            else:  # e is subdimension vector of d
                # List of all generic subdimension vectors of e
                genSubdims = self.all_generic_subdimension_vectors(e)
                return all([self.euler_form(f, d - e) >= 0 for f in genSubdims])

    # TODO remove this and cache the recursive one instead
    # This method computes a list of all generic subdimension vectors of e, for all e which are subdimension vectors of d.
    def __all_generic_subdimension_vectors_helper(self, d):
        """Returns the list of lists of indexes of all generic subdimension vectors of e, where e ranges over all subdimension vectors of d. The index refers to the deglex order.

        EXAMPLES:

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

        Some n-Kronecker quivers::

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
        d = self._coerce_dimension_vector(d)

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


        Generic ext on the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(1)
            sage: R = [Q.simple_root(1), Q.simple_root(2), vector([1, 1])]
            sage: for a in R:
            ....:     for b in R:
            ....:         print('ext(' + str(a) + ',' + str(b) + ') = ' + str(Q.generic_ext(a, b)))
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
        a = self._coerce_dimension_vector(a)
        b = self._coerce_dimension_vector(b)
        genSubdims = self.all_generic_subdimension_vectors(a)
        return max([-self.euler_form(c, b) for c in genSubdims])

    def generic_hom(self, a, b):
        r"""Computes hom(a, b).

        INPUT:

        - ``a``: dimension vector

        - ``b``: dimension vector

        OUTPUT: Int

        There is a non-empty open subset U of R(Q,a) x R(Q,b) such that dim Ext(M,N) = ext(a,b), i.e. is minimal, for all (M,N) in U. Therefore dim Hom(M,N) = <a,b> + dim Ext(M,N) is minimal and therefore hom(a,b) = <a,b> + ext(a,b).

        EXAMPLES:

        Generic hom for the Kronecker quiver::

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
        a = self._coerce_dimension_vector(a)
        b = self._coerce_dimension_vector(b)

        return self.euler_form(a, b) + self.generic_ext(a, b)

    # TODO remove
    def generic_ext_vanishing(self, a, b):
        return self.is_generic_subdimension_vector(a, a + b)

    # TODO remove
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
        d = self._coerce_dimension_vector(d)

        return vector(d) * (-self.euler_matrix().transpose() + self.euler_matrix())

    def has_semistable_representation(self, d, theta):
        r"""Checks if there is a `\theta`-semistable representation of dimension vector `d`

        INPUT:
        - ``d``: dimension vector

        - ``theta``: stability parameter

        OUTPUT: Statement truth value as Bool

        See Thm. 5.4(1) of Reineke's overview paper https://arxiv.org/pdf/0802.2147.pdf: A dimension vector d admits a theta-semi-stable representation if and only if mu_theta(e) <= mu_theta(d) for all generic subdimension vectors e of d.
        # Thm. 5.4 in Markus's paper is actually a result of Schofield. So the algorithm should bear his name, if any.

        EXAMPLES:

        Semistables for the A_2 quiver::

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

        Semistables for the 3-Kronecker quiver::

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
        d = self._coerce_dimension_vector(d)

        genSubdims = self.all_generic_subdimension_vectors(d)
        genSubdims = list(filter(lambda e: e != self.zero_vector(), genSubdims))
        return all([slope(e, theta) <= slope(d, theta) for e in genSubdims])

    # TODO remove and cache the recursive one instead
    def __all_semistable_subdimension_vectors_helper(self, d, theta):
        """Computes the list of indexes of all semistable subdimension vectors of d.

        EXAMPLES:

        Manual caching for all semistable subdimension vectors of (2,3) for the Kronecker quiver::

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

        OUTPUT: True if there is a `\theta`-stable representation of `d`, False otherwise.

        See Thm. 5.4(1) of Reineke's overview paper https://arxiv.org/pdf/0802.2147.pdf: A dimension vector d admits a theta-stable representation if and only if mu_theta(e) < mu_theta(d) for all proper generic subdimension vectors e of d.

        EXAMPLES:

        Stables for the `\mathrm{A}_2` quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(1)
            sage: theta = (1, -1)
            sage: Q.has_stable_representation([1, 1], theta, algorithm="schofield")
            True
            sage: Q.has_stable_representation([2, 2], theta, algorithm="schofield")
            False
            sage: Q.has_stable_representation([0, 0], theta, algorithm="schofield")
            False

        Stables for the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: d = (2, 3)
            sage: theta = Q.canonical_stability_parameter(d)
            sage: Q.has_stable_representation(d, theta, algorithm="schofield")
            True

        """
        assert algorithm in ["schofield", "king", "al"]

        # coerce stability parameter
        theta = vector(theta)
        assert theta.length() == self.number_of_vertices()

        d = self._coerce_dimension_vector(d)

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

    # TODO remove and cache the recursive one instead
    def __all_stable_subdimension_vectors_helper(self, d, theta, denominator=sum):
        """Computes the list of all stable subdimension vectors of d which have the same slope as d.

        EXAMPLES:

        Manual caching of all stable subdimension vectors of (3,3) for the 3-Kronecker quiver::

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

        Now for the 2-Kronecker quiver::

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

            Canonical decompositions of small dimension vectors for the Kronecker quiver::

                # TODO make a smaller example? don't print _everything_?
                sage: from quiver import *
                sage: Q = KroneckerQuiver()
                sage: ds = Tuples(range(7), 2)
                sage: decompositions = {d: Q.canonical_decomposition(d, algorithm="recursive") for d in ds}
                sage: for d in ds:
                ....:     print("decomposition of {} is {}".format(d, decompositions[d]))
                decomposition of (0, 0) is [(0, 0)]
                decomposition of (1, 0) is [(1, 0)]
                decomposition of (2, 0) is [(1, 0), (1, 0)]
                decomposition of (3, 0) is [(1, 0), (1, 0), (1, 0)]
                decomposition of (4, 0) is [(1, 0), (1, 0), (1, 0), (1, 0)]
                decomposition of (5, 0) is [(1, 0), (1, 0), (1, 0), (1, 0), (1, 0)]
                decomposition of (6, 0) is [(1, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0)]
                decomposition of (0, 1) is [(0, 1)]
                decomposition of (1, 1) is [(1, 1)]
                decomposition of (2, 1) is [(2, 1)]
                decomposition of (3, 1) is [(1, 0), (2, 1)]
                decomposition of (4, 1) is [(1, 0), (1, 0), (2, 1)]
                decomposition of (5, 1) is [(1, 0), (1, 0), (1, 0), (2, 1)]
                decomposition of (6, 1) is [(1, 0), (1, 0), (1, 0), (1, 0), (2, 1)]
                decomposition of (0, 2) is [(0, 1), (0, 1)]
                decomposition of (1, 2) is [(1, 2)]
                decomposition of (2, 2) is [(1, 1), (1, 1)]
                decomposition of (3, 2) is [(3, 2)]
                decomposition of (4, 2) is [(2, 1), (2, 1)]
                decomposition of (5, 2) is [(1, 0), (2, 1), (2, 1)]
                decomposition of (6, 2) is [(1, 0), (1, 0), (2, 1), (2, 1)]
                decomposition of (0, 3) is [(0, 1), (0, 1), (0, 1)]
                decomposition of (1, 3) is [(0, 1), (1, 2)]
                decomposition of (2, 3) is [(2, 3)]
                decomposition of (3, 3) is [(1, 1), (1, 1), (1, 1)]
                decomposition of (4, 3) is [(4, 3)]
                decomposition of (5, 3) is [(2, 1), (3, 2)]
                decomposition of (6, 3) is [(2, 1), (2, 1), (2, 1)]
                decomposition of (0, 4) is [(0, 1), (0, 1), (0, 1), (0, 1)]
                decomposition of (1, 4) is [(0, 1), (0, 1), (1, 2)]
                decomposition of (2, 4) is [(1, 2), (1, 2)]
                decomposition of (3, 4) is [(3, 4)]
                decomposition of (4, 4) is [(1, 1), (1, 1), (1, 1), (1, 1)]
                decomposition of (5, 4) is [(5, 4)]
                decomposition of (6, 4) is [(3, 2), (3, 2)]
                decomposition of (0, 5) is [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1)]
                decomposition of (1, 5) is [(0, 1), (0, 1), (0, 1), (1, 2)]
                decomposition of (2, 5) is [(0, 1), (1, 2), (1, 2)]
                decomposition of (3, 5) is [(1, 2), (2, 3)]
                decomposition of (4, 5) is [(4, 5)]
                decomposition of (5, 5) is [(1, 1), (1, 1), (1, 1), (1, 1), (1, 1)]
                decomposition of (6, 5) is [(6, 5)]
                decomposition of (0, 6) is [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1)]
                decomposition of (1, 6) is [(0, 1), (0, 1), (0, 1), (0, 1), (1, 2)]
                decomposition of (2, 6) is [(0, 1), (0, 1), (1, 2), (1, 2)]
                decomposition of (3, 6) is [(1, 2), (1, 2), (1, 2)]
                decomposition of (4, 6) is [(2, 3), (2, 3)]
                decomposition of (5, 6) is [(5, 6)]
                decomposition of (6, 6) is [(1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1)]

            We verify the vanishing of generic Ext::

                sage: all(all(Q.generic_ext(di, dj) + Q.generic_ext(dj, di) == 0
                ....:         for (di, dj) in Combinations(s, 2))
                ....:     for s in decompositions.values())
                True

            """
            d = self._coerce_dimension_vector(d)

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

            def idx_diff(j, i):
                return subdims.index(subdims[j] - subdims[i])

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
        d = self._coerce_dimension_vector(d)

        if self.is_acyclic():
            return d.transpose() * self.adjacency_matrix() * d
        else:
            # TODO where is the algorithm?
            raise NotImplementedError()

    """
    Teleman!
    """

    # TODO: This section should go into QuiverModuliSpace, I think.
    # TODO return weights as dictionaries with HN types as keys.
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

        INPUT:

        - ``d`` -- dimension vector

        - ``theta`` -- stability parameter

        - ``denominator`` -- function to compute the denominator of the slope. Default is sum.

        OUTPUT: True if the rigidity inequality holds for d and theta, False otherwise.

        If the weights of the 1-PS lambda on $\det(N_{S/R}|_Z)$ for each HN type
        are all strictly larger than the weights of the tensors of the universal bundles $U_i^\vee \otimes U_j$,
        then the resulting moduli space is infinitesimally rigid.
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
    # coerce dimension vectors
    d = vector(d)
    e = vector(e)

    assert e.length() == d.length()

    return all(0 <= ei and ei <= di for (ei, di) in zip(e, d))


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

    Examples of coprimality::

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
