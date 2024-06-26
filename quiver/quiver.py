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

        # for caching purposes: order along the specified order of vertices
        self.__M = self.__G.adjacency_matrix(vertices=self.__G.vertices(sort=False))

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
        return cls(G, name)

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
    def from_string(cls, Q: str, forget_labels=True, name=None):
        r"""Construct a quiver from a comma-separated list of chains of the form "i-j-k-...

        You specify an arrow from `i` to `j` by writing `i-j`.
        A multiple arrow is specified by repeating the `i`, so that `1--2` is the Kronecker
        quiver. If you write `i-j-k` then you have 1 arrow from `i` to `j` and one from `j`
        to `k`. The full quiver is specified by concatenating (multiple) arrows by commas.

        The values for a vertex can be anything, and the chosen names will be used for
        the vertices in the underlying graph. Labels are cast to an integer, if possible,
        and otherwise to strings.

        OUTPUT: the quiver

        INPUT:

        - ``Q`` -- string; a string of the format described above giving a quiver

        - ``forget_labels`` -- boolean (default: True): specifies whether to use the
          names for vertices, or whether to use the integers `0,...,n-1` where `n` is
          the number of vertices

        - ``name`` -- optional name for the quiver

        EXAMPLES:

        The 3-Kronecker quiver defined in two different ways::

            sage: from quiver import *
            sage: Quiver.from_matrix([[0, 3], [0, 0]]) == Quiver.from_string("a---b")
            True

        A more complicated example::

            sage: Q = Quiver.from_string("a--b-3,a---3,3-a")
            sage: Q.adjacency_matrix()
            [0 2 3]
            [0 0 1]
            [1 0 0]
            sage: Q.vertices()
            [0, 1, 2]

        The actual labeling we use doesn't matter for the isomorphism type of the quiver::

            sage: from quiver import *
            sage: Quiver.from_matrix([[0, 3], [0, 0]]) == Quiver.from_string("12---b")
            True

        However, it does influence the labels of the vertex if we choose so::

            sage: Quiver.from_string("12---b", forget_labels=False).vertices()
            [12, 'b']
            sage: Quiver.from_string("foo---bar", forget_labels=False).vertices()
            ['foo', 'bar']

        """
        # remove all whitespace from the string
        Q = "".join(Q.split())

        # determine the vertices used in the string, preserving the order
        vertices = list(
            dict.fromkeys(
                vertex
                for chain in Q.split(",")
                for vertex in chain.split("-")
                if vertex  # this filters out "" from chained hyphens
            )
        )

        # adjacency matrix to be built
        M = zero_matrix(len(vertices))

        for chain in Q.split(","):
            pieces = chain.split("-")

            source = vertices.index(pieces[0])
            number = 1

            for piece in pieces[1:]:
                # if the string is empty we increase the number of arrows counter
                if not piece:
                    number += 1
                # if the string is non-empty we treat it as a label
                # this means we add the appropriate number of arrows
                # and make the target the source to start the process again
                if piece:
                    target = vertices.index(piece)
                    M[source, target] += number

                    number = 1
                    source, target = target, None

        G = DiGraph(M)

        # attempt to cast vertex labels to integers, otherwise to strings
        if not forget_labels:
            labels = []
            for vertex in vertices:
                try:
                    labels.append(int(vertex))
                except ValueError:
                    labels.append(str(vertex))
            G.relabel(perm=labels, inplace=True)

        return cls.from_digraph(G, name)

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

    def _is_vector(self, x):
        r"""
        Checks whether `x` is an element of `\mathbb{Z}Q_0`

        If the quiver doesn't use vertex labels we check that it has the right length.
        If the quiver uses vertex labels, we check that `d` is a dict with the right
        set of keys.

        We actually do not care whether the values are in `\mathbb{Z}`.

        INPUT:

        - `x` -- vector

        OUTPUT: whether `x` can be used as a dimension vector for the quiver

        EXAMPLES:

        Some basic examples::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: Q._is_vector([2, 3])
            True
            sage: Q._is_vector([0, 0])
            True
            sage: Q._is_vector([-2, -2])
            True
            sage: Q._is_vector([1, 2, 3])
            False

        We allow non-integral values, because this can be useful for stability::

            sage: Q._is_vector([1/2, 3])
            True

        An example with vertex labels::

            sage: Q = Quiver.from_string("foo---bar", forget_labels=False)
            sage: Q._is_vector({"foo" : 0, "bar" : 0})
            True
            sage: Q._is_vector({"bar" : 0, "foo" : 0})
            True
            sage: Q._is_vector({"baz" : 0, "ofo" : 0})
            False

        """
        if isinstance(x, dict):
            return set(x.keys()) == set(self.vertices())

        return len(x) == self.number_of_vertices()

    def _is_dimension_vector(self, d):
        r"""
        Checks whether `d` is a dimension vector of the quiver

        If the quiver doesn't use vertex labels we check that it has the right length
        and has positive entries.
        If the quiver uses vertex labels, we check that `d` is a dict with the right
        set of keys and positive entries.

        We only check for non-negativity, not for integrality.

        INPUT:

        - `d` -- dimension vector

        OUTPUT: whether `d` can be used as a dimension vector for the quiver

        EXAMPLES:

        Some basic examples::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: Q._is_dimension_vector([2, 3])
            True
            sage: Q._is_dimension_vector([0, 0])
            True
            sage: Q._is_dimension_vector([-2, -2])
            False
            sage: Q._is_dimension_vector([1, 2, 3])
            False

        An example with vertex labels::

            sage: Q = Quiver.from_string("foo---bar", forget_labels=False)
            sage: Q._is_dimension_vector({"foo" : 2, "bar" : 3})
            True
            sage: Q._is_dimension_vector({"bar" : 0, "foo" : -1})
            False
            sage: Q._is_dimension_vector({"baz" : 0, "ofo" : 0})
            False

        """
        if not self._is_vector(d):
            return False

        if isinstance(d, dict):
            return all(di >= 0 for di in d.values())

        return all(di >= 0 for di in d)

    def _coerce_dimension_vector(self, d):
        r"""
        Coerces `d` to be a dimension vector of the quiver

        It must be a data structure that is indexed by the vertices of the quiver,
        so most likely a dict, list, or vector.
        It is coerced to a vector, see :meth:`Quiver._coerce_vector`.

        As a consistency check we verify that all entries are non-negative,
        raising a `ValueError` if it isn't.

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

        It must be a data structure that is indexed by the vertices of the quiver,
        so most likely a dict, list, or vector.

        It raises a `ValueError` if it is not a data structure of length the number
        of vertices in the quiver.

        INPUT:

        - ``x``: a list, tuple, or dict of integers

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
        if len(x) != self.number_of_vertices():
            raise ValueError("The input is not an element of `\mathbb{Z}Q_0`.")

        if isinstance(x, list) or isinstance(x, tuple):
            x = vector(ZZ, x)
        elif isinstance(x, dict):
            x = vector(ZZ, [x[i] for i in self.vertices()])

        return x

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

    def __has_vertex_labels(self) -> bool:
        r"""Check whether vertex labels are used

        EXAMPLES:

        With vertex labels::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(2)
            sage: Q._Quiver__has_vertex_labels()
            False

        With vertex labels::

            sage: Q = Quiver.from_string("foo---bar", forget_labels=False)
            sage: Q._Quiver__has_vertex_labels()
            True

        """
        return self.vertices() != list(range(self.number_of_vertices()))

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

        If we specified a non-standard labeling on the vertices we must use it::

            sage: Q = Quiver.from_string("a---b", forget_labels=False)
            sage: Q.in_degree("a")
            0
            sage: Q.in_degree("b")
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

        If we specified a non-standard labeling on the vertices we must use it::

            sage: Q = Quiver.from_string("a---b", forget_labels=False)
            sage: Q.out_degree("a")
            3
            sage: Q.out_degree("b")
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

        If we specified a non-standard labeling on the vertices we must use it::

            sage: Q = Quiver.from_string("a---b", forget_labels=False)
            sage: Q.is_source("a")
            True
            sage: Q.is_source("b")
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

        If we specified a non-standard labeling on the vertices we must use it::

            sage: Q = Quiver.from_string("a---b", forget_labels=False)
            sage: Q.is_sink("a")
            False
            sage: Q.is_sink("b")
            True

        """
        return self.out_degree(i) == 0

    def sources(self):
        r"""Return the vertices which are sources in the quiver

        OUTPUT: the list of vertices without incoming edges

        EXAMPLES:

        The 3-Kronecker quiver has one source::

            sage: from quiver import *
            sage: GeneralizedKroneckerQuiver(3).sources()
            [0]

        It is possible that a quiver has no sources::

            sage: JordanQuiver().sources()
            []

        """
        return list(filter(lambda i: self.is_source(i), self.vertices()))

    def sinks(self):
        r"""Return the vertices which are sinks in the quiver

        OUTPUT: the list of vertices without incoming edges

        EXAMPLES:

        The 3-Kronecker quiver has one source::

            sage: from quiver import *
            sage: GeneralizedKroneckerQuiver(3).sources()
            [0]

        It is possible that a quiver has no sinks::

            sage: JordanQuiver().sinks()
            []

        """
        return list(filter(lambda i: self.is_sink(i), self.vertices()))

    """
    Basic representation-theoretical properties of the quiver
    """

    def euler_matrix(self):
        r"""Returns the Euler matrix of the quiver

        This is the matrix representing the Euler form, defined by

        \begin{equation}
            \langle\mathbf{d},\mathbf{e}\rangle=
            \sum_{i\in Q_0}d_ie_i-\sum_{\alpha\in Q_1}d_{s(\alpha)}e_{t(\alpha)}
        \end{equation}

        In the basis given by the vertices, it can be written as the difference
        of the identity matrix and the adjacency matrix.

        OUTPUT: the Euler matrix of the quiver

        EXAMPLES:

        The Kronecker 3-quiver::

            sage: from quiver import *
            sage: GeneralizedKroneckerQuiver(3).euler_matrix()
            [ 1 -3]
            [ 0  1]

        It uses the basis of the vertices, so it agrees with this alternative definition::

            sage: Quiver.from_string("foo---bar", forget_labels=False).euler_matrix()
            [ 1 -3]
            [ 0  1]

        """
        return matrix.identity(self.number_of_vertices()) - self.adjacency_matrix()

    def euler_form(self, x, y) -> int:
        r"""The value `\langle x,y\rangle` of the Euler form

        INPUT:

        - ``x`` -- an element of `\mathbb{Z}Q_0`

        - ``y`` -- an element of `\mathbb{Z}Q_0`

        OUTPUT: the value of the Euler form, i.e., `x * self.euler_matrix() * y`

        EXAMPLES:

        An example using the Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.euler_form([1, 3], [2, -2])
            2

        It uses the basis of the vertices, so we specify the entries of elements of
        `\mathbb{Z}Q_0` in this order, thus the same example as before::

            sage: Q = Quiver.from_string("foo---bar", forget_labels=False)
            sage: Q.euler_form([1, 3], [2, -2])
            2

        """
        x = self._coerce_vector(x)
        y = self._coerce_vector(y)

        return x * self.euler_matrix() * y

    def cartan_matrix(self):
        r"""Returns the Cartan matrix of the quiver

        This is the matrix representing the symmetrization of the Euler form,
        see :meth:`Quiver.euler_matrix`

        OUTPUT: the Cartan matrix of the quiver

        EXAMPLES::

        The Kronecker 3-quiver::

            sage: from quiver import *
            sage: GeneralizedKroneckerQuiver(3).cartan_matrix()
            [ 2 -3]
            [-3  2]

        """
        return self.euler_matrix() + self.euler_matrix().transpose()

    def symmetrized_euler_form(self, x, y) -> int:
        r"""The value `(x,y)` of the Euler form

        INPUT:

        - ``x`` -- an element of `\mathbb{Z}Q_0`

        - ``y`` -- an element of `\mathbb{Z}Q_0`

        OUTPUT: the value of the symmetrized Euler form applied to `x` and `y`

        EXAMPLES:

        An example using the Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.symmetrized_euler_form([1, 3], [2, -2])
            -20

        It uses the basis of the vertices, so we specify the entries of elements of
        `\mathbb{Z}Q_0` in this order, thus the same example as before::

            sage: Q = Quiver.from_string("foo---bar", forget_labels=False)
            sage: Q.symmetrized_euler_form([1, 3], [2, -2])
            -20

        """
        x = self._coerce_vector(x)
        y = self._coerce_vector(y)

        return self.euler_form(x, y) + self.euler_form(y, x)

    def tits_form(self, x) -> int:
        r"""The value of the Tits quadratic form of the quiver at `x`

        This is really just the value `\langle x,x\rangle` of the Euler form,
        or half of the value `(x,x)` of the symmetrized Euler form.

        INPUT:

        - ``x`` -- an element of `\mathbb{Z}Q_0`

        OUTPUT: the value of the Tits form applied to `x`
        EXAMPLES:

        An example using the Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.tits_form([2, 3])
            -5

        It uses the basis of the vertices, so we specify the entries of elements of
        `\mathbb{Z}Q_0` in this order, thus the same example as before::

            sage: Q = Quiver.from_string("foo---bar", forget_labels=False)
            sage: Q.tits_form([2, 3])
            -5

        """
        return self.euler_form(x, x)

    """
    Constructing new quivers out of old
    """

    def opposite_quiver(self):
        r"""
        Returns the opposite quiver

        The opposite quiver is the quiver with all arrows reversed.
        Its adjacency matrix is given by the transpose of the adjacency matrix.

        OUTPUT: the opposite quiver

        EXAMPLES:

        The opposite of the 3-Kronecker quiver::

            sage: from quiver import *
            sage: print(GeneralizedKroneckerQuiver(3).opposite_quiver())
            opposite of 3-Kronecker quiver
            adjacency matrix:
            [0 0]
            [3 0]

        It preserves the labelling of the vertices::

            sage: Q = Quiver.from_string("foo---bar", forget_labels=False).opposite_quiver()
            sage: Q.vertices()
            ['foo', 'bar']
            sage: Q.adjacency_matrix()
            [0 0]
            [3 0]

        """
        name = None
        if self.get_custom_name():
            name = "opposite of " + self.get_custom_name()

        return Quiver.from_digraph(self.graph().reverse(), name)

    def doubled_quiver(self):
        r"""
        Returns the doubled quiver

        The double of a quiver is the quiver where for each arrow we add an arrow in
        the opposite direction.

        Its adjacency matrix is the sum of the adjacency matrix of the original quiver
        and its transpose.

        OUTPUT: the doubled quiver

        EXAMPLES:

        The double of the 3-Kronecker quiver::

            sage: from quiver import *
            sage: print(GeneralizedKroneckerQuiver(3).doubled_quiver())
            double of 3-Kronecker quiver
            adjacency matrix:
            [0 3]
            [3 0]

        It preserves the labelling of the vertices::

            sage: Q = Quiver.from_string("foo---bar", forget_labels=False).doubled_quiver()
            sage: Q.vertices()
            ['foo', 'bar']
            sage: Q.adjacency_matrix()
            [0 3]
            [3 0]

        """
        G = DiGraph(self.graph())
        G.add_edges(self.opposite_quiver().graph().edges())

        name = None
        if self.get_custom_name():
            name = "double of " + self.get_custom_name()

        return Quiver.from_digraph(G, name)

    def framed_quiver(self, framing, vertex="-oo"):
        r"""
        Returns the framed quiver with framing vector `framing`

        The optional parameter `vertex` determines the name of the framing vertex,
        which defaults to `-oo`.

        The framed quiver has one additional vertex, and `f_i` many arrows from
        the framing vertex to `i`, for every `i\in Q_0`.

        INPUT:

        - ``framing`` -- list of non-negative integers saying how many arrows from the
                         framed vertex to `i`

        - ``vertex`` (default: "-oo") -- name of the framing vertex

        OUTPUT: the framed quiver

        EXAMPLES:

        Framing the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3).framed_quiver([1, 0])
            sage: print(Q)
            framing of 3-Kronecker quiver
            adjacency matrix:
            [0 1 0]
            [0 0 3]
            [0 0 0]
            sage: Q.vertices()
            ['-oo', 0, 1]
            sage: Q = GeneralizedKroneckerQuiver(3).framed_quiver([2, 2], vertex="a")
            sage: print(Q)
            framing of 3-Kronecker quiver
            adjacency matrix:
            [0 2 2]
            [0 0 3]
            [0 0 0]
            sage: Q.vertices()
            ['a', 0, 1]

        If you frame twice it will have to use a different vertex label::

            sage: Q = GeneralizedKroneckerQuiver(3).framed_quiver([2, 2])
            sage: Q.framed_quiver([1, 1, 1]).vertices()
            Traceback (most recent call last):
            ...
            ValueError: -oo is already a vertex

        """
        if vertex in self.vertices():
            raise ValueError("{} is already a vertex".format(vertex))

        framing = self._coerce_dimension_vector(framing)

        G = DiGraph(self.graph())

        # adding framing the vertex (as the _first_ vertex by default)
        G.add_vertex(vertex)

        # adding the arrows according to the framing vector
        for i, v in enumerate(self.vertices()):
            G.add_edges([(vertex, v)] * framing[i])

        name = None
        if self.get_custom_name():
            name = "framing of " + self.get_custom_name()

        return Quiver.from_digraph(G, name)

    def coframed_quiver(self, coframing, vertex="+oo"):
        r"""
        Returns the coframed quiver with coframing vector `coframing`

        The optional parameter `vertex` determines the name of the coframing vertex,
        which defaults to `+oo`.

        The coframed quiver has one additional vertex, and `f_i` many arrows from
        the vertex `i` to the coframed vertex, for every `i\in Q_0`.

        INPUT:

        - ``coframing`` -- list of non-negative integers saying how many arrows from the
                         framed vertex to `i`

        - ``vertex`` (default: None) -- name of the framing vertex

        OUTPUT: the framed quiver

        EXAMPLES:

        Coframing the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3).coframed_quiver([1, 0])
            sage: print(Q)
            coframing of 3-Kronecker quiver
            adjacency matrix:
            [0 3 1]
            [0 0 0]
            [0 0 0]
            sage: Q.vertices()
            [0, 1, '+oo']
            sage: Q = GeneralizedKroneckerQuiver(3).coframed_quiver([2, 2], vertex="a")
            sage: print(Q)
            coframing of 3-Kronecker quiver
            adjacency matrix:
            [0 3 2]
            [0 0 2]
            [0 0 0]
            sage: Q.vertices()
            [0, 1, 'a']

        If you coframe twice it will have to use a different vertex label::

            sage: Q = GeneralizedKroneckerQuiver(3).coframed_quiver([2, 2])
            sage: Q.coframed_quiver([1, 1, 1]).vertices()
            Traceback (most recent call last):
            ...
            ValueError: +oo is already a vertex

        """
        if vertex in self.vertices():
            raise ValueError("{} is already a vertex".format(vertex))

        coframing = self._coerce_dimension_vector(coframing)

        # build the graph from the ground up
        G = DiGraph([self.vertices() + [vertex], []], multiedges=True, loops=True)

        # there doesn't seem to be a way to save the order of vertices otherwise
        # so make sure the coframing vertex appears last
        permutation = dict(zip(G.vertices(), self.vertices() + [vertex]))
        G.relabel(perm=permutation, inplace=True)

        # adding the existing arrows
        G.add_edges(self.graph().edges())

        # adding the arrows according to the framing vector
        for i, v in enumerate(self.vertices()):
            G.add_edges([(v, vertex)] * coframing[i])

        name = None
        if self.get_custom_name():
            name = "coframing of " + self.get_custom_name()

        return Quiver.from_digraph(G, name)

    def full_subquiver(self, vertices):
        r"""Returns the full subquiver supported on the given set of vertices

        INPUT:

        - ``vertices``: list of vertices for the subquiver

        OUTPUT: the full subquiver on the specified vertices

        EXAMPLES:

        Some basic examples::

            sage: from quiver import *
            sage: Q = ThreeVertexQuiver(2, 3, 4)
            sage: print(Q.full_subquiver([0, 1]))
            full subquiver of an acyclic 3-vertex quiver of type (2, 3, 4)
            adjacency matrix:
            [0 2]
            [0 0]
            sage: print(Q.full_subquiver([0, 2]))
            full subquiver of an acyclic 3-vertex quiver of type (2, 3, 4)
            adjacency matrix:
            [0 3]
            [0 0]

        If we specified a non-standard labeling on the vertices we must use it::

            sage: Q = Quiver.from_string("a--b----c,a---c", forget_labels=False)
            sage: Q == ThreeVertexQuiver(2, 3, 4)
            True
            sage: print(Q.full_subquiver(["a", "b"]))
            a quiver with 2 vertices and 2 arrows
            adjacency matrix:
            [0 2]
            [0 0]
            sage: print(Q.full_subquiver(["a", "c"]))
            a quiver with 2 vertices and 3 arrows
            adjacency matrix:
            [0 3]
            [0 0]

        """
        name = None
        if self.get_custom_name():
            name = "full subquiver of " + self.get_custom_name()

        return Quiver.from_digraph(self.graph().subgraph(vertices=vertices), name)

    """
    Dimension vectors and roots
    """

    def zero_vector(self):
        r"""
        Returns the zero dimension vector

        The output is adapted to the vertices.

        OUTPUT: the zero dimension vector

        EXAMPLES:

        Usually it is an actual vector::

            sage: from quiver import *
            sage: KroneckerQuiver(3).zero_vector()
            (0, 0)
            sage: type(KroneckerQuiver(3).zero_vector())
            <class 'sage.modules.vector_integer_dense.Vector_integer_dense'>

        But if the quiver has custom vertex labels it is a dict::

            sage: Q = Quiver.from_string("a--b----c,a---c", forget_labels=False)
            sage: Q.zero_vector()
            {'a': 0, 'b': 0, 'c': 0}

        """
        if self.__has_vertex_labels():
            return {i: 0 for i in self.vertices()}
        else:
            return vector([0] * self.number_of_vertices())

    def thin_dimension_vector(self):
        r"""
        Returns the thin dimension vector, i.e., all ones

        The output is adapted to the vertices.

        OUTPUT: the thin dimension vector

        EXAMPLES:

        Usually it is an actual vector::

            sage: from quiver import *
            sage: KroneckerQuiver(3).thin_dimension_vector()
            (1, 1)
            sage: type(KroneckerQuiver(3).thin_dimension_vector())
            <class 'sage.modules.vector_integer_dense.Vector_integer_dense'>

        But if the quiver has custom vertex labels it is a dict::

            sage: Q = Quiver.from_string("a--b----c,a---c", forget_labels=False)
            sage: Q.thin_dimension_vector()
            {'a': 1, 'b': 1, 'c': 1}

        """
        if self.__has_vertex_labels():
            return {i: 1 for i in self.vertices()}
        else:
            return vector([1] * self.number_of_vertices())

    def simple_root(self, i):
        r"""
        Returns the simple root at the vertex `i`

        The output is adapted to the vertices.

        OUTPUT: the simple root at the vertex `i`

        EXAMPLES:

        Usually it is an actual vector::

            sage: from quiver import *
            sage: KroneckerQuiver(3).simple_root(1)
            (0, 1)
            sage: type(KroneckerQuiver(3).simple_root(1))
            <class 'sage.modules.vector_integer_dense.Vector_integer_dense'>

        But if the quiver has custom vertex labels it is a dict::

            sage: Q = Quiver.from_string("a--b----c,a---c", forget_labels=False)
            sage: Q.simple_root("b")
            {'a': 0, 'b': 1, 'c': 0}

        """
        root = self.zero_vector()
        root[i] = 1

        return root

    def is_root(self, x) -> bool:
        r"""Checks whether `x` is a root of the underlying diagram of the quiver.

        A root is a non-zero vector in `\mathbb{Z}Q_0` such that
        the Tits form of `x` is bounded above by 1..

        INPUT:
        - ``x``: integer vector

        OUTPUT: whether `x` is a root

        EXAMPLES:

        Some roots and non-roots for the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: Q.is_root([2, 3])
            True
            sage: Q.is_root(Q.zero_vector())
            False
            sage: Q.is_root([4, 1])
            False

        """
        x = self._coerce_vector(x)

        return any(x) and self.tits_form(x) <= 1

    def is_real_root(self, x) -> bool:
        r"""Checks whether `x` is a real root of the underlying diagram of the quiver.

        A root is called real if its Tits form equals 1.

        INPUT:
        - ``x``: integer vector

        OUTPUT: whether `x` is a real root

        EXAMPLES:

        Some real and non-real for the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: Q.is_real_root([2, 3])
            False
            sage: Q.is_real_root(Q.zero_vector())
            False
            sage: Q.is_real_root([3, 1])
            True

        """
        x = self._coerce_vector(x)

        return self.tits_form(x) == 1

    def is_imaginary_root(self, x) -> bool:
        r"""Checks whether `x` is a imaginary root of the underlying diagram of the quiver.

        A root is called imaginary if its Tits form is non-positive.

        INPUT:
        - ``x``: integer vector

        OUTPUT: whether `x` is an imaginary root

        EXAMPLES:

        Some imaginary roots and non imaginary roots for the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: Q.is_imaginary_root([2, 3])
            True
            sage: Q.is_imaginary_root(Q.zero_vector())
            False
            sage: Q.is_imaginary_root([4, 1])
            False

        """
        x = self._coerce_vector(x)

        return any(x) and self.tits_form(x) <= 0

    def is_schur_root(self, d) -> bool:
        r"""Checks if `d` is a Schur root.

        INPUT:
        - ``d``: dimension vector

        OUTPUT: whether `d` is an imaginary root

        A Schur root is a dimension vector which admits a Schurian representation,
        i.e., a representation whose endomorphism ring is the field itself.
        It is necessarily indecomposable.

        By a result of Schofield (https://mathscinet.ams.org/mathscinet/relay-station?mr=1162487)
        `d` is a Schur root if and only if `d` admits a stable representation for
        the canonical stability parameter.

        EXAMPLES:

        The dimension vector `(2, 3)` is Schurian for the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.is_schur_root([2, 3])
            True

        Examples from Derksen--Weyman's book (Example 11.1.4)::

            sage: from quiver import *
            sage: Q = ThreeVertexQuiver(1, 1, 1)
            sage: Q.is_schur_root([1, 1, 2])
            True
            sage: Q.is_schur_root([1, 2, 1])
            False
            sage: Q.is_schur_root([1, 1, 1])
            True
            sage: Q.is_schur_root([2, 2, 2])
            False

        """
        d = self._coerce_dimension_vector(d)
        theta = self.canonical_stability_parameter(d)

        return self.has_stable_representation(d, theta)

    def slope(self, d, theta, denominator=sum):
        r"""
        Returns the slope of `d` with respect to `theta`

        The slope is defined as the value of `theta(d)` divided by the total dimension
        of `d`. It is possible to vary the denominator, to use a function more general
        than the sum.

        INPUT:

        - `d` -- dimension vector

        - `theta` -- stability parameter

        OUTPUT: the slope of `d` with respect to `theta` and optional `denominator`

        EXAMPLES:

        Some slopes for the Kronecker quiver, first for the canonical stability
        parameter, then for some other::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: d = [2, 3]
            sage: Q.slope(d, [9, -6])
            0
            sage: Q.slope(d, [2, -2])
            -2/5

        We can use for instance a constant denominator::

            sage: constant = lambda di: 1
            sage: Q.slope(d, Q.canonical_stability_parameter(d), denominator=constant)
            0

        The only dependence on the quiver is the set of vertices, so if we don't
        use vertex labels, the choice of quiver doesn't matter::

            sage: d, theta = [2, 3], [9, -6]
            sage: KroneckerQuiver(2).slope(d, theta) == KroneckerQuiver(3).slope(d, theta)
            True

        """

        assert denominator(d) > 0

        d = self._coerce_dimension_vector(d)
        theta = self._coerce_vector(theta)

        return (theta * d) / denominator(d)

    def is_subdimension_vector(self, e, d):
        r"""
        Determine whether `e` is a subdimension vector of `d`

        INPUT:

        -- `e` -- dimension vector

        -- `d` -- dimension vector

        OUTPUT: whether `e` is a subdimension vector of `d`

        EXAMPLES:

        Some basic examples::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: Q.is_subdimension_vector([1, 2], [2, 3])
            True
            sage: Q.is_subdimension_vector([2, 3], [2, 3])
            True
            sage: Q.is_subdimension_vector([6, 6], [2, 3])
            False

        We can also work with vertex labels::

            sage: Q = Quiver.from_string("a--b----c,a---c", forget_labels=False)
            sage: d = {"a" : 3, "b" : 3, "c" : 3}
            sage: e = {"a" : 1, "b" : 2, "c" : 3}
            sage: Q.is_subdimension_vector(e, d)
            True
            sage: Q.is_subdimension_vector(d, e)
            False

        """
        d = self._coerce_dimension_vector(d)
        e = self._coerce_dimension_vector(e)

        return all(ei <= di for (ei, di) in zip(e, d))

    def deglex_key(self, e, b):
        """A function which satisfies e <_{deglex} d iff deglex_key(e) < deglex_key(d), provided that b >> 0."""
        n = len(e)
        return sum([e[i] * b ** (n - i - 1) for i in range(n)]) + sum(e) * b**n

    def all_subdimension_vectors(self, d, proper=False, nonzero=False):
        """Returns the list of all subdimension vectors of d."""
        # TODO if Q has vertex labels, then returning vector isn't good enough!
        vectors = list(map(vector, cartesian_product([range(di + 1) for di in d])))
        # TODO clean this up: what if d is zero?
        if proper:
            vectors = vectors[:-1]
        if nonzero:
            vectors = vectors[1:]
        return vectors

    # TODO: This method has a stupid name (my own fault). Think of a better one.
    # TODO whenever a theta needs to be provided, we should default to the canonical one?
    def is_theta_coprime(self, d, theta) -> bool:
        """Checks if d is theta-coprime.

        A dimension vector d is theta-coprime if mu_theta(e) != mu_theta(e) for all proper subdimension vectors e of d.

        EXAMPLES

        Examples of coprimality::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: d = [2, 3]
            sage: Q.is_theta_coprime(d, Q.canonical_stability_parameter(d))
            True
            sage: Q.is_theta_coprime([3, 3], [1, -1])
            False

        """
        assert self._is_dimension_vector(d)
        assert self._is_vector(theta)

        vectors = self.all_subdimension_vectors(d, proper=True, nonzero=True)

        return all([self.slope(d, theta) != self.slope(e, theta) for e in vectors])

    def is_indivisible(self, d) -> bool:
        """
        Checks if the gcd of all entries of `d` is 1

        INPUT:

        -- `d` -- dimension vector

        OUTPUT: whether the dimension vector is indivisible

        EXAMPLES:

        Two examples with the Kronecker quiver::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: Q.is_indivisible([2, 3])
            True
            sage: Q.is_indivisible([2, 2])
            False

        """
        return gcd(self._coerce_dimension_vector(d)) == 1

    def support(self, d):
        r"""Returns the support of the dimension vector.

        INPUT:

        - ``d``: dimension vector

        OUTPUT: subset of vertices in the underlying graph in the support

        The support is the set `\{ i \in Q_0 \mid d_i > 0 \}`.

        EXAMPLES

        The support is the set of vertices for which the value of the dimension
        vector is nonzero::

            sage: from quiver import *
            sage: Q = ThreeVertexQuiver(2, 0, 4)
            sage: d = vector([1, 1, 1])
            sage: Q.support(d)
            [0, 1, 2]
            sage: d = vector([1, 0, 1])
            sage: Q.support(d)
            [0, 2]

        It takes into account vertex labels::

            sage: Q = Quiver.from_string("a--b----c,a---c", forget_labels=False)
            sage: d = {"a": 2, "b": 3, "c": 0}
            sage: Q.support(d)
            ['a', 'b']

        """
        # TODO we need also a Quiver._is_dimension_vector(d) method?
        # TODO if d is a dict, should we have defaultdict behavior for zeroes?

        return [i for i in self.vertices() if d[i] > 0]

    def in_fundamental_domain(self, d, depth=0):
        r"""Checks if a dimension vector is in the fundamental domain.

        The fundamental domain of `Q` is the set of dimension vectors `d` such that

        * `\operatorname{supp}(\mathbf{d})` is connected
        * `\langle d,e_i\rangle + \langle e_i,d\rangle\leq 0` for every simple root `e_i`.

        Every `d` in the fundamental domain is an imaginary root and the set of
        imaginary roots is the Weyl group saturation of the fundamental domain.
        If `d` is in the fundamental domain then it is Schurian and a general
        representation of dimension vector `d` is stable for the canonical stability parameter.

        The optional parameter `depth` allows to make the inequality stricter.

        INPUT:

        - ``d``: dimension vector

        - ``depth`` (default: 0) -- how deep the vector should be in the domain

        OUTPUT: whether `d` is in the (interior of) the fundamental domain

        EXAMPLES:

        The fundamental domain of the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.in_fundamental_domain([1, 1])
            True
            sage: Q.in_fundamental_domain([1, 2])
            False
            sage: Q.in_fundamental_domain([2, 3])
            True

        The same calculation now with vertex labels::

            sage: Q = Quiver.from_string("a---b", forget_labels=False)
            sage: Q.in_fundamental_domain({"a" : 1, "b" : 1})
            True
            sage: Q.in_fundamental_domain({"a" : 1, "b" : 2})
            False
            sage: Q.in_fundamental_domain({"a" : 2, "b" : 3})
            True

        We test for dimension vectors in the strict interior, where the depth is
        equal to 1::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.in_fundamental_domain([1, 1], depth=1)
            True
            sage: Q.in_fundamental_domain([2, 3], depth=1)
            False

        """
        # TODO we need also a Quiver._is_dimension_vector(d) method?
        # TODO we don't want to coerce here!

        # check if `\langle d,e_i\rangle + \langle e_i,d\rangle \leq 0`
        # for all vertices `i\in Q_0`
        inequality = all(
            self.symmetrized_euler_form(d, self.simple_root(i)) <= -depth
            for i in self.vertices()
        )

        # check if the support is connected
        connected = self.full_subquiver(self.support(d)).is_connected()

        return inequality and connected

    # TODO what is the use case?
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

        return (
            all([d[i] <= e[i] for i in self.sources()])
            and all([d[i] >= e[i] for i in self.sinks()])
            and all(
                [
                    d[i] == e[i]
                    for i in self.vertices()
                    if i not in self.sources() and i not in self.sinks()
                ]
            )
        )

    """
    Generic subdimension vectors and generic Hom and Ext
    """

    # taken from code/snippets/canonical.sage
    # TODO still need testing code from there
    def is_generic_subdimension_vector(self, e, d) -> bool:
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
            ....:         if Q.is_subdimension_vector(e,d):
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
            ....:         if Q.is_subdimension_vector(e,d):
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
            if not self.is_subdimension_vector(e, d):
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
        subdims = self.all_subdimension_vectors(d)
        subdims.sort(key=(lambda e: self.deglex_key(e, b=max(d) + 1)))
        # We use the deglex order because it's a total order which extends the usual entry-wise partial order on dimension vectors.
        N = len(subdims)

        # genIndexes[j] will in the end be the list of indexes (in subdims) of all generic subdimension vectors of subdims[j]
        genIndexes = [
            list(
                filter(
                    lambda i: self.is_subdimension_vector(subdims[i], subdims[j]),
                    range(N),
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
            sage: R = [Q.simple_root(0), Q.simple_root(1), vector([1, 1])]
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
            sage: R = [Q.simple_root(0), Q.simple_root(1), vector([1,1])]
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

    def is_left_orthogonal(self, a, b) -> bool:
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
        zero_vector = self._coerce_dimension_vector(self.zero_vector())

        # TODO exclude_zero parameter in all_generic_subdimension_vectors?
        genSubdims = self.all_generic_subdimension_vectors(d)
        genSubdims = list(filter(lambda e: e != zero_vector, genSubdims))
        return all([self.slope(e, theta) <= self.slope(d, theta) for e in genSubdims])

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
        subdims = self.all_subdimension_vectors(d)
        subdims.sort(key=(lambda e: self.deglex_key(e, b=max(d) + 1)))
        # We use the deglex order because it's a total order which extends the usual entry-wise partial order on dimension vectors.
        N = len(subdims)
        genIndexes, genSubdims = self.__all_generic_subdimension_vectors_helper(d)
        sstIndexes = list(
            filter(
                lambda j: all(
                    [
                        self.slope(subdims[i], theta) <= self.slope(subdims[j], theta)
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
                return all(
                    [self.slope(e, theta) < self.slope(d, theta) for e in genSubdims]
                )

    # TODO remove and cache the recursive one instead
    def __all_stable_subdimension_vectors_helper(self, d, theta, denominator=sum):
        """Computes the list of all stable subdimension vectors of d which have the same slope as d.

        EXAMPLES:

        Manual caching of all stable subdimension vectors of (3,3) for the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q, d, theta = GeneralizedKroneckerQuiver(3), vector([3,3]), vector([1,0])
            sage: Q.all_subdimension_vectors(d)
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
        subdims = self.all_subdimension_vectors(d)
        subdims.sort(key=(lambda e: self.deglex_key(e, b=max(d) + 1)))
        # We use the deglex order because it's a total order which extends the usual entry-wise partial order on dimension vectors.
        N = len(subdims)
        genIndexes, genSubdims = self.__all_generic_subdimension_vectors_helper(d)
        # slopeIndexes is the list of subdimension vectors of d of the same slope as d (in particular != 0)
        slopeIndexes = list(
            filter(
                lambda j: self.slope(subdims[j], theta, denominator=denominator)
                == self.slope(d, theta, denominator=denominator),
                range(1, N),
            )
        )
        # stIndexes contains all j for which subdims[j] is stable
        # e = subdims[j] is stable if for all generic subdimension vectors f = subdims[i] of e, it holds that slope(f) < slope(e)
        stIndexes = list(
            filter(
                lambda j: all(
                    [
                        self.slope(subdims[i], theta, denominator=denominator)
                        < self.slope(subdims[j], theta, denominator=denominator)
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
                [d[i], self.simple_root(i)]
                # TODO range over self.vertices() instead
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
            subdims = self.all_subdimension_vectors(d)
            subdims.sort(key=(lambda e: self.deglex_key(e, b=max(d) + 1)))
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
                            self.slope(hntype[s], theta, denominator=denominator)
                            - self.slope(hntype[t], theta, denominator=denominator)
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
        # TODO Quiver.all_harder_narasimhan_types needs way to filter out [d]
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
                            self.slope(hntype[s], theta, denominator=denominator)
                            - self.slope(hntype[t], theta, denominator=denominator)
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
                lambda hntype: self.slope(hntype[0], theta, denominator=denominator)
                - self.slope(hntype[-1], theta, denominator=denominator),
                HN,
            )
        )

        return all([weights[i] > tensorWeights[i] for i in range(len(HN))])

    def first_hochschild_cohomology(self):
        r"""
        Compute the dimension of the first Hochschild cohomology

        This uses the formula of Happel from Proposition 1.6 in [MR1035222].
        One needs the quiver to be acyclic for this.

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

        return (
            1
            - self.number_of_vertices()
            + sum(
                len(self.graph().all_paths(a[0], a[1], use_multiedges=True))
                for a in self.graph().edges()
            )
        )

    """
    Some things that don't below anywhere else yet?
    """

    def is_cofree(self, d) -> bool:
        # https://mathscinet.ams.org/mathscinet-getitem?mr=2394691
        # we don't really know what this means
        raise NotImplementedError()

    def perpendicular_category(self, d):
        # something from Schofield
        # see Theorem 11.4.6 in the Derksen--Weyman book
        raise NotImplementedError()


"""Auxiliary methods"""


"""Class methods"""
