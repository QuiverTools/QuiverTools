from sage.arith.misc import gcd
from sage.categories.cartesian_product import cartesian_product
from sage.combinat.root_system.cartan_matrix import CartanMatrix
from sage.graphs.digraph import DiGraph
from sage.matrix.constructor import matrix
from sage.matrix.special import zero_matrix
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ
from sage.structure.element import Element


class Quiver(Element):
    r"""
    A quiver is a (finite) directed multigraph. It is an important tool in the
    representation theory of (finite-dimensional) algebras, because it allows one to
    construct the path algebra, whose modules are equivalently described as
    representations of the quiver. These in turn can be classified using moduli spaces
    of quiver representations.

    For an introduction to the subject one is referred to

    * Harm Derksen and Jerzy Weyman: An introduction to quiver representations
    * Markus Reineke: Moduli of representations of quivers

    or one of the many other resources that exist.
    """

    def __init__(self, G, name=None):
        r"""Constructor for a quiver.

        This takes a directed graph as input. If it is not a `DiGraph` instance,
        it is interpreted it as an adjacency matrix.
        For other constructions, see

        - :meth:`Quiver.from_digraph`
        - :meth:`Quiver.from_matrix`
        - :meth:`Quiver.from_string`

        INPUT:

        - ``G`` -- directed graph

        - ``name`` -- optional name for the quiver

        EXAMPLES:

        The 3-Kronecker quiver from an adjacency matrix::

            sage: from quiver import *
            sage: Q = Quiver([[0, 3], [0, 0]]); Q
            a quiver with 2 vertices and 3 arrows

        """

        if isinstance(G, DiGraph):
            # it is the user's responsibility to not change the graph afterwards
            self.__G = G
        else:
            self.__G = DiGraph(matrix(G), immutable=True)

        # if name is None this doesn't do anything
        self.rename(name)

        # for caching purposes: order along the specified order of vertices
        self.__M = self.__G.adjacency_matrix(vertices=self.__G.vertices(sort=False))

    @classmethod
    def from_digraph(cls, G, name=None):
        r"""Construct a quiver from a `DiGraph` object.

        INPUT:

        - ``G`` -- directed graph as a `DiGraph` object

        - ``name`` -- optional name for the quiver

        OUTPUT: the quiver.

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
        r"""Construct a quiver from its adjacency matrix.

        INPUT:

        - ``M`` -- adjacency matrix of the quiver

        - ``name`` -- optional name for the quiver

        OUTPUT: the quiver.

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = Quiver.from_matrix([[0, 3], [0, 0]]); Q.adjacency_matrix()
            [0 3]
            [0 0]

        """
        return cls(matrix(M), name)

    @classmethod
    def from_string(cls, Q: str, forget_labels=True, name=None):
        r"""Construct a quiver from a comma-separated list of chains like ``i-j-k-...``

        You specify an arrow from ``i`` to ``j`` by writing ``i-j``.
        Multiple arrows are specified by repeating the hyphen, so that ``1--2`` is the
        Kronecker quiver. If you write ``i-j-k`` then you have 1 arrow from ``i`` to
        ``j`` and one from ``j`` to ``k``. The full quiver is specified by concatenating
        (multiple) arrows by commas.

        The values for a vertex can be anything, and the chosen names will be used for
        the vertices in the underlying graph. Labels are cast to an integer, if
        possible, and otherwise to strings.

        INPUT:

        - ``Q`` -- a string of the format described above giving a quiver

        - ``forget_labels`` -- (default: True): whether to use labels for vertices or
          to number them ``0,...,n-1``

        - ``name`` -- optional name for the quiver

        OUTPUT: the quiver

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

        The actual labeling we use doesn't matter for the isomorphism type of the
        quiver::

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

    def _repr_(self) -> str:
        r"""
        Give a shorthand string presentation for the quiver

        If a name is set, use that, if not, just give information on number of vertices
        and arrows.

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
        if self.get_custom_name():
            return self.get_custom_name()

        return "a quiver with {} vertices and {} arrows".format(
            self.adjacency_matrix().nrows(), sum(sum(self.adjacency_matrix()))
        )

    def __str__(self) -> str:
        r"""
        Detailed description of the quiver

        Everything you get from :meth:`Quiver.repr()` together with the adjacency
        matrix.

        EXAMPLES:

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
        return "{}\nadjacency matrix:\n{}".format(self.repr(), self.adjacency_matrix())

    def repr(self) -> str:
        r"""
        Basic description of the quiver

        To override the output, one uses :meth:`Quiver.rename` from the `Element`
        class. The output of :meth:`Quiver.repr` is that of
        :meth:`Quiver.get_custom_name` if it is set, else it is the default specifying
        the number of vertices and arrows.

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
            sage: Q.get_custom_name() is None
            True
            sage: Q.rename("3-Kronecker quiver")
            sage: Q.get_custom_name()
            '3-Kronecker quiver'
            sage: Q.reset_name()
            sage: Q.get_custom_name() is None
            True
            sage: Q
            a quiver with 2 vertices and 3 arrows

        """
        return self._repr_()

    def str(self) -> str:
        r"""
        Full description of the quiver

        This combines the output of :meth:`Quiver.repr` with the adjacency matrix.

        OUTPUT: a complete description of the quiver

        EXAMPLES:

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
        Checks whether ``x`` is an element of :math:`\mathbb{Z}Q_0`

        If the quiver doesn't use vertex labels we check that it has the right length.
        If the quiver uses vertex labels, we check that ``x`` is a dict with the right
        set of keys.

        We actually do not care whether the values are in :math:`\mathbb{Z}`.

        INPUT:

        - ``x`` -- vector

        OUTPUT: whether ``x`` can be used as a dimension vector for the quiver

        EXAMPLES:

        Some basic examples::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: Q._is_vector((2, 3))
            True
            sage: Q._is_vector((0, 0))
            True
            sage: Q._is_vector((-2, -2))
            True
            sage: Q._is_vector((1, 2, 3))
            False

        We allow non-integral values, because this can be useful for stability::

            sage: Q._is_vector((1/2, 3))
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
        Checks whether ``d`` is a dimension vector of the quiver

        If the quiver doesn't use vertex labels we check that it has the right length
        and has positive entries.
        If the quiver uses vertex labels, we check that ``d`` is a dict with the right
        set of keys and positive entries.

        We only check for non-negativity, not for integrality.

        INPUT:

        - ``d`` -- dimension vector

        OUTPUT: whether ``d`` can be used as a dimension vector for the quiver

        EXAMPLES:

        Some basic examples::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: Q._is_dimension_vector((2, 3))
            True
            sage: Q._is_dimension_vector((0, 0))
            True
            sage: Q._is_dimension_vector((-2, -2))
            False
            sage: Q._is_dimension_vector((1, 2, 3))
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
        Coerces ``d`` to be a dimension vector of the quiver

        The input ``d`` must be a data structure that is indexed
        by the vertices of the quiver, so most likely a dict, list, or vector.
        It is coerced to a vector, see :meth:`Quiver._coerce_vector`.

        As a consistency check we verify that all entries are non-negative,
        raising a `ValueError` if it isn't the case.

        INPUT:

        - ``d``: a candidate dimension vector

        OUTPUT: either a dict or vector

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q._coerce_dimension_vector((1, 2))
            (1, 2)
            sage: Q._coerce_dimension_vector((1, 2, 3, 4))
            Traceback (most recent call last):
            ...
            ValueError: The input is not an element of :math:`\mathbb{Z}Q_0`.
            sage: Q._coerce_dimension_vector((1, -3))
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
        Coerces ``x`` to be a vector in :math:`\mathbb{Z}Q_0`.

        The input ``x`` must be a data structure that is indexed by
        the vertices of the quiver,
        so most likely a dict, list, tuple, or vector.

        It raises a `ValueError` if it is not a data structure of length the number
        of vertices in the quiver.

        INPUT:

        - ``x``: a list, tuple, or dict of integers

        OUTPUT: a Sage vector if ``x`` is an element of :math:`\mathbb{Z}Q_0`

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q._coerce_vector((-1, 2))
            (-1, 2)
            sage: Q._coerce_vector((1, 2, 3, 4))
            Traceback (most recent call last):
            ...
            ValueError: The input is not an element of :math:`\mathbb{Z}Q_0`.

        """
        if len(x) != self.number_of_vertices():
            raise ValueError(r"The input is not an element of :math:`\mathbb{Z}Q_0`.")

        if isinstance(x, list) or isinstance(x, tuple):
            x = vector(ZZ, x)
        elif isinstance(x, dict):
            x = vector(ZZ, [x[i] for i in self.vertices()])

        # so that it can be used for hashing
        x.set_immutable()

        return x

    """
    Basic graph-theoretic properties of the quiver
    """

    def adjacency_matrix(self):
        r"""
        Returns the adjacency matrix of the quiver.

        OUTPUT: The square matrix ``M`` whose entry ``M[i,j]`` is the number of arrows
        from the vertex ``i`` to the vertex ``j``

        EXAMPLES:

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
        using the data in the DiGraph or string, as explained in
        :meth:`Quiver.from_digraph` or :meth:`Quiver.from_string`.
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

        There are 7 arrows in this 3-vertex quiver::

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

    def is_finite_type(self) -> bool:
        r"""
        Returns whether the quiver is of finite type representation type.

        This is the case if and only the connected components of the underlying
        undirected graph are isomorphic to Dynkin diagrams.

        EXAMPLES:

        The generalized Kronecker quiver is finite only for :math:`m=1`::

            sage: from quiver import *
            sage: GeneralizedKroneckerQuiver(1).is_finite_type()
            True
            sage: GeneralizedKroneckerQuiver(2).is_finite_type()
            False
            sage: GeneralizedKroneckerQuiver(3).is_finite_type()
            False

        """
        return CartanMatrix(self.cartan_matrix()).is_finite()

    def is_tame_type(self) -> bool:
        r"""
        Returns whether the quiver is of tame type representation type.

        This is the case if and only the connected components of the underlying
        undirected graph are isomorphic to (extended) Dynkin diagrams, with at least one
        being extended Dynkin.

        EXAMPLES:

        The generalized Kronecker quiver is tame only for :math:`m=2`::

            sage: from quiver import *
            sage: GeneralizedKroneckerQuiver(1).is_tame_type()
            False
            sage: GeneralizedKroneckerQuiver(2).is_tame_type()
            True
            sage: GeneralizedKroneckerQuiver(3).is_tame_type()
            False

        """
        M = CartanMatrix(self.cartan_matrix())
        return M.is_affine() and not M.is_finite()

    def is_wild_type(self) -> bool:
        r"""
        Returns whether the quiver is of wild type representation type.

        This is the case if and only the connected components of the underlying
        undirected graph are not all isomorphic to (extended) Dynkin diagrams.

        EXAMPLES:

        The generalized Kronecker quiver is wild for all :math:`m\geq 3`::

            sage: from quiver import *
            sage: GeneralizedKroneckerQuiver(1).is_wild_type()
            False
            sage: GeneralizedKroneckerQuiver(2).is_wild_type()
            False
            sage: GeneralizedKroneckerQuiver(3).is_wild_type()
            True

        """
        return not self.is_finite_type() and not self.is_tame_type()

    """
    Some graph-theoretic properties of the quiver
    """

    def in_degree(self, i):
        r"""Returns the in-degree of a vertex.

        The in-degree of ``i`` is the number of incoming arrows at ``i``.

        The parameter ``i`` must be an element of the vertices of the underlying graph.
        If constructed from a matrix or string, ``i`` can go from `0` to
        `n-1` where `n` is the number of vertices in the graph.

        INPUT:

        - ``i`` -- a vertex of the underlying graph

        OUTPUT: The in-degree of the vertex ``i``

        EXAMPLES:

        In the 3-Kronecker quiver the in-degree is either 0 or 3::

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

        The parameter ``i`` must be an element of the vertices of the underlying graph.
        If constructed from a matrix or string, ``i`` can go from `0` to
        `n-1` where `n` is the number of vertices in the graph.

        The out-degree of ``i`` is the number of outgoing arrows at ``i``.

        INPUT:

        - ``i`` -- a vertex of the underlying graph

        OUTPUT: The out-degree of the vertex ``i``

        EXAMPLES:

        In the 3-Kronecker quiver the out-degree is either 3 or 0::

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
        """Checks if ``i`` is a source of the quiver

        The vertex ``i`` is a source if there are no incoming arrows at ``i``.

        INPUT:

        - ``i`` -- a vertex of the quiver

        OUTPUT: whether ``i`` is a source of the quiver

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
        """Checks if ``i`` is a sink of the quiver

        The vertex ``i`` is a sink if there are no outgoing arrows out of ``i``.

        INPUT:

        - ``i`` -- a vertex of the quiver

        OUTPUT: whether ``i`` is a sink of the quiver

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

        .. MATH::

            \langle\mathbf{d},\mathbf{e}\rangle=
            \sum_{i\in Q_0}d_i e_i-\sum_{\alpha\in Q_1}d_{s(\alpha)}e_{t(\alpha)}

        In the basis given by the vertices, it can be written as the difference
        of the identity matrix and the adjacency matrix.

        OUTPUT: the Euler matrix of the quiver

        EXAMPLES:

        The Kronecker 3-quiver::

            sage: from quiver import *
            sage: GeneralizedKroneckerQuiver(3).euler_matrix()
            [ 1 -3]
            [ 0  1]

        It uses the basis of the vertices, so it agrees with this alternative
        definition::

            sage: Quiver.from_string("foo---bar", forget_labels=False).euler_matrix()
            [ 1 -3]
            [ 0  1]

        """
        return matrix.identity(self.number_of_vertices()) - self.adjacency_matrix()

    def euler_form(self, x, y) -> int:
        r"""The value :math:`\langle x,y\rangle` of the Euler form

        INPUT:

        - ``x`` -- an element of :math:`\mathbb{Z}Q_0`

        - ``y`` -- an element of :math:`\mathbb{Z}Q_0`

        OUTPUT: the value of the Euler form, i.e., ``x * self.euler_matrix() * y``

        EXAMPLES:

        An example using the Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.euler_form((1, 3), (2, -2))
            2

        It uses the basis of the vertices, so we specify the entries of elements of
        :math:`\mathbb{Z}Q_0` in this order, thus the same example as before::

            sage: Q = Quiver.from_string("foo---bar", forget_labels=False)
            sage: Q.euler_form((1, 3), (2, -2))
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

        EXAMPLES:

        The Kronecker 3-quiver::

            sage: from quiver import *
            sage: GeneralizedKroneckerQuiver(3).cartan_matrix()
            [ 2 -3]
            [-3  2]

        """
        return self.euler_matrix() + self.euler_matrix().transpose()

    def symmetrized_euler_form(self, x, y) -> int:
        r"""The value :math:`(x,y)` of the Euler form

        INPUT:

        - ``x`` -- an element of :math:`\mathbb{Z}Q_0`

        - ``y`` -- an element of :math:`\mathbb{Z}Q_0`

        OUTPUT: the value of the symmetrized Euler form applied to ``x`` and ``y``

        EXAMPLES:

        An example using the Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.symmetrized_euler_form((1, 3), (2, -2))
            -20

        It uses the basis of the vertices, so we specify the entries of elements of
        :math:`\mathbb{Z}Q_0` in this order, thus the same example as before::

            sage: Q = Quiver.from_string("foo---bar", forget_labels=False)
            sage: Q.symmetrized_euler_form((1, 3), (2, -2))
            -20

        """
        x = self._coerce_vector(x)
        y = self._coerce_vector(y)

        return self.euler_form(x, y) + self.euler_form(y, x)

    def tits_form(self, x) -> int:
        r"""The value of the Tits quadratic form of the quiver at ``x``

        This is really just the value :math:`\langle x,x\rangle` of the Euler form,
        or half of the value :math:`(x,x)` of the symmetrized Euler form.

        INPUT:

        - ``x`` -- an element of :math:`\mathbb{Z}Q_0`

        OUTPUT: the value of the Tits form applied to ``x``

        EXAMPLES:

        An example using the Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.tits_form((2, 3))
            -5

        It uses the basis of the vertices, so we specify the entries of elements of
        :math:`\mathbb{Z}Q_0` in this order, thus the same example as before::

            sage: Q = Quiver.from_string("foo---bar", forget_labels=False)
            sage: Q.tits_form((2, 3))
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

            sage: Q = Quiver.from_string("foo---bar", forget_labels=False)
            sage: Qopp = Q.opposite_quiver()
            sage: Qopp.vertices()
            ['foo', 'bar']
            sage: Qopp.adjacency_matrix()
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

        The double of a quiver is the quiver where for each arrow
        we add an arrow in the opposite direction.

        Its adjacency matrix is the sum of the adjacency matrix
        of the original quiver and its transpose.

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

            sage: Q = Quiver.from_string("foo---bar", forget_labels=False)
            sage: Qbar = Q.doubled_quiver()
            sage: Qbar.vertices()
            ['foo', 'bar']
            sage: Qbar.adjacency_matrix()
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
        Returns the framed quiver with framing vector ``framing``

        The optional parameter ``vertex`` determines the name of the framing vertex,
        which defaults to `-oo`.

        The framed quiver has one additional vertex, and :math:`f_i` many arrows from
        the framing vertex to :math:`i`, for every :math:`i\in Q_0`.

        INPUT:

        - ``framing`` -- list of non-negative integers saying how many arrows from the
          framed vertex to ``i``

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
        Returns the coframed quiver with coframing vector ``coframing``

        The optional parameter ``vertex`` determines the name of the coframing vertex,
        which defaults to `+oo`.

        The coframed quiver has one additional vertex, and :math:`f_i` many arrows from
        the vertex `i` to the coframed vertex, for every :math:`i\in Q_0`.

        INPUT:

        - ``coframing`` -- list of non-negative integers saying how many arrows go from
          the framed vertex to `i`

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
        Returns the zero dimension vector.

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
        return self.__zero_vector()

    @cached_method
    def __zero_vector(self):
        r"""The cacheable implementation of :meth:`Quiver.zero_vector`

        EXAMPLES:

        The zero vector of the 3-Kronecker quiver::

            sage: from quiver import *
            sage: KroneckerQuiver(3).zero_vector()
            (0, 0)

        If we specified a non-standard labeling on the vertices it is used::

            sage: Q = Quiver.from_string("a--b----c,a---c", forget_labels=False)
            sage: Q._Quiver__zero_vector()
            {'a': 0, 'b': 0, 'c': 0}
        """
        if self.__has_vertex_labels():
            return {i: 0 for i in self.vertices()}

        return vector([0] * self.number_of_vertices(), immutable=True)

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
        return self.__thin_dimension_vector()

    @cached_method
    def __thin_dimension_vector(self):
        r"""The cacheable implementation of :meth:`Quiver.thin_dimension_vector`

        EXAMPLES:

        The thin dimension vector of the 3-Kronecker quiver::

            sage: from quiver import *
            sage: KroneckerQuiver(3).thin_dimension_vector()
            (1, 1)

        If we specified a non-standard labeling on the vertices it is used::

            sage: Q = Quiver.from_string("a--b----c,a---c", forget_labels=False)
            sage: Q._Quiver__thin_dimension_vector()
            {'a': 1, 'b': 1, 'c': 1}
        """
        if self.__has_vertex_labels():
            return {i: 1 for i in self.vertices()}

        return vector([1] * self.number_of_vertices(), immutable=True)

    def simple_root(self, i):
        r"""
        Returns the simple root at the vertex ``i``

        The output is adapted to the vertices.

        OUTPUT: the simple root at the vertex ``i``

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
        return self.__simple_root(i)

    @cached_method
    def __simple_root(self, i):
        r"""The cacheable implementation of :meth:`Quiver.simple_root`

        EXAMPLES:

        The simple root at the source of the 3-Kronecker quiver::

            sage: from quiver import *
            sage: KroneckerQuiver(3).simple_root(0)
            (1, 0)

        If we specified a non-standard labeling on the vertices it is used::

            sage: Q = Quiver.from_string("a--b----c,a---c", forget_labels=False)
            sage: Q._Quiver__simple_root("a")
            {'a': 1, 'b': 0, 'c': 0}
        """
        if self.__has_vertex_labels():
            root = {i: 0 for i in self.vertices()}
            root[i] = 1

            return root

        root = vector([0] * self.number_of_vertices())
        root[i] = 1
        root.set_immutable()

        return root

    def is_root(self, x) -> bool:
        r"""Checks whether ``x`` is a root of the underlying diagram of the quiver.

        A root is a non-zero vector `x` in :math:`\mathbb{Z}Q_0` such that
        the Tits form of `x` is at most 1.

        INPUT:

        - ``x``: integer vector

        OUTPUT: whether ``x`` is a root

        EXAMPLES:

        Some roots and non-roots for the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: Q.is_root((2, 3))
            True
            sage: Q.is_root(Q.zero_vector())
            False
            sage: Q.is_root((4, 1))
            False

        """
        x = self._coerce_vector(x)

        return any(x) and self.tits_form(x) <= 1

    def is_real_root(self, x) -> bool:
        r"""Checks whether ``x`` is a real root of the underlying diagram of the quiver.

        A root is called real if its Tits form equals 1.

        INPUT:

        - ``x``: integer vector

        OUTPUT: whether ``x`` is a real root

        EXAMPLES:

        Some real and non-real for the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: Q.is_real_root((2, 3))
            False
            sage: Q.is_real_root(Q.zero_vector())
            False
            sage: Q.is_real_root((3, 1))
            True

        """
        x = self._coerce_vector(x)

        return self.tits_form(x) == 1

    def is_imaginary_root(self, x) -> bool:
        r"""Checks whether ``x`` is a imaginary root of the quiver.

        A root is called imaginary if its Tits form is non-positive.

        INPUT:

        - ``x``: integer vector

        OUTPUT: whether ``x`` is an imaginary root

        EXAMPLES:

        Some imaginary roots and non imaginary roots for the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: Q.is_imaginary_root((2, 3))
            True
            sage: Q.is_imaginary_root(Q.zero_vector())
            False
            sage: Q.is_imaginary_root((4, 1))
            False

        """
        x = self._coerce_vector(x)

        return any(x) and self.tits_form(x) <= 0

    def is_schur_root(self, d) -> bool:
        r"""Checks if ``d`` is a Schur root.

        INPUT:

        - ``d``: dimension vector

        OUTPUT: whether ``d`` is an imaginary root

        A Schur root is a dimension vector which admits a Schurian representation,
        i.e., a representation whose endomorphism ring is the field itself.
        It is necessarily indecomposable.

        By MR1162487_ :math:`{\bf d}` is a Schur root if and only if it admits a stable
        representation for the canonical stability parameter.

        .. _MR1162487: https://mathscinet.ams.org/mathscinet/relay-station?mr=1162487

        EXAMPLES:

        The dimension vector `(2, 3)` is Schurian for the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.is_schur_root([2, 3])
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

    def slope(self, d, theta=None, denom=sum):
        r"""
        Returns the slope of ``d`` with respect to ``theta``

        The slope is defined as the value of ``theta(d)`` divided by the total dimension
        of `d` ``sum(d)``. It is possible to vary the denominator, to use a function more general
        than the sum.

        INPUT:

        - ``d`` -- dimension vector

        - ``theta`` -- (default: canonical stability parameter) stability parameter

        - ``denom`` -- (default: sum) the denominator function

        OUTPUT: the slope of ``d`` with respect to ``theta`` and optional ``denom``

        EXAMPLES:

        Some slopes for the Kronecker quiver, first for the canonical stability
        parameter, then for some other::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: d = (2, 3)
            sage: Q.slope(d, (9, -6))
            0
            sage: Q.slope(d)
            0
            sage: Q.slope(d, (2, -2))
            -2/5

        We can use for instance a constant denominator::

            sage: constant = lambda di: 1
            sage: Q.slope(d, Q.canonical_stability_parameter(d), denom=constant)
            0

        The only dependence on the quiver is the set of vertices, so if we don't
        use vertex labels, the choice of quiver doesn't matter::

            sage: d, theta = (2, 3), (9, -6)
            sage: KroneckerQuiver(3).slope(d, theta)
            0

        """
        d = self._coerce_dimension_vector(d)
        assert denom(d) > 0, "denominator needs to be strictly positive on ``d``"

        if theta is None:
            theta = self.canonical_stability_parameter(d)
        theta = self._coerce_vector(theta)

        return (theta * d) / denom(d)

    def is_subdimension_vector(self, e, d):
        r"""
        Determine whether ``e`` is a subdimension vector of ``d``

        INPUT:

        -- ``e`` -- dimension vector

        -- ``d`` -- dimension vector

        OUTPUT: whether ``e`` is a subdimension vector of ``d``

        EXAMPLES:

        Some basic examples::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: Q.is_subdimension_vector((1, 2), (2, 3))
            True
            sage: Q.is_subdimension_vector((2, 3), (2, 3))
            True
            sage: Q.is_subdimension_vector((6, 6), (2, 3))
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

    def _deglex_key(self, e, b=None) -> int:
        r"""
        An integer representation of a dimension vector

        This is the base-b expansion of a dimension vector.

        This is a function which satisfies

            e <_{deglex} d iff deglex_key(e) < deglex_key(d),

        provided that b >> 0.

        For b >> 0 the deglex order is a _total_ order which extends the usual
        entry-wise partial order on dimension vectors.

        INPUT:

        - ``e`` -- dimension vector

        - ``b`` -- the "base" of the key (default: `max(e)+1`)

        OUTPUT: the base-`b` expansion of the dimension vector

        EXAMPLES:

        If we let `b` be the largest entry plus one we get a good key, at least for
        subdimension vectors of the original one::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: d = (2, 3)
            sage: Q._deglex_key(d, max(d) + 1)
            91
            sage: d = (3, 3)
            sage: Q._deglex_key(d)
            111

        """
        e = self._coerce_dimension_vector(e)
        if b is None:
            b = max(e) + 1

        n = self.number_of_vertices()

        return (
            sum(ei * b ** (n - i - 1) for (i, ei) in enumerate(e))
            + sum(self._coerce_dimension_vector(e)) * b**n
        )

    def all_subdimension_vectors(
        self, d, proper=False, nonzero=False, forget_labels=False
    ):
        r"""
        Returns the list of all subdimension vectors of ``d``.

        INPUT:

        - ``d`` -- dimension vector

        - ``proper`` (default: False) -- whether to exclude ``d``

        - ``nonzero`` (default: False) -- whether to exclude the zero vector

        - ``forget_labels`` (default: False) -- whether to forget the vertex labels

        OUTPUT: all subdimension vectors of ``d`` (maybe excluding zero and/or ``d``)

        EXAMPLES:

        The usual use cases::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: Q.all_subdimension_vectors((2, 3))
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
                 (2, 3)]
            sage: Q.all_subdimension_vectors((2, 3), proper=True)
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
                 (2, 2)]
            sage: Q.all_subdimension_vectors((2, 3), nonzero=True)
                [(0, 1),
                 (0, 2),
                 (0, 3),
                 (1, 0),
                 (1, 1),
                 (1, 2),
                 (1, 3),
                 (2, 0),
                 (2, 1),
                 (2, 2),
                 (2, 3)]
            sage: Q.all_subdimension_vectors((2, 3), proper=True, nonzero=True)
                [(0, 1),
                 (0, 2),
                 (0, 3),
                 (1, 0),
                 (1, 1),
                 (1, 2),
                 (1, 3),
                 (2, 0),
                 (2, 1),
                 (2, 2)]

        Some exceptional cases::

            sage: Q.all_subdimension_vectors(Q.zero_vector())
            [(0, 0)]
            sage: Q.all_subdimension_vectors(Q.zero_vector(), proper=True)
            []

        If we work with labeled vertices, then we get a list of dicts::

            sage: Q = Quiver.from_string("a---b", forget_labels=False)
            sage: Q.all_subdimension_vectors((1, 2))
            [{'a': 0, 'b': 0},
             {'a': 0, 'b': 1},
             {'a': 0, 'b': 2},
             {'a': 1, 'b': 0},
             {'a': 1, 'b': 1},
             {'a': 1, 'b': 2}]
            sage: Q.all_subdimension_vectors((1, 2), forget_labels=True)
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]

        """
        assert self._is_dimension_vector(d), "``d`` needs to be a dimension vector"

        # if zero dimension vector we deal with it separately
        if sum(self._coerce_dimension_vector(d)) == 0:
            if proper or nonzero:
                return []
            return [d]

        vectors = list(cartesian_product([range(di + 1) for di in d]))

        if proper:
            vectors = vectors[:-1]
        if nonzero:
            vectors = vectors[1:]

        if self.__has_vertex_labels() and not forget_labels:
            return list(map(lambda e: dict(zip(self.vertices(), e)), vectors))

        return list(map(vector, vectors))

    def is_theta_coprime(self, d, theta=None) -> bool:
        r"""Checks if ``d`` is ``theta``-coprime.

        A dimension vector `d` is :math:`\theta`-coprime if
        :math:`\mu_{\theta}(e)\neq \mu_{\theta}(e)`
        for all proper non-zero subdimension vectors e of d.

        The default value for ``theta`` is the canonical stability parameter,
        see :meth:`canonical_stability_parameter`.

        INPUT:

        - ``d`` -- dimension vector

        - ``theta`` -- (default: canonical stability paramter) stability parameter

        EXAMPLES:

        Examples of coprimality::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: d = (2, 3)
            sage: Q.is_theta_coprime(d, Q.canonical_stability_parameter(d))
            True
            sage: Q.is_theta_coprime(d)
            True
            sage: Q.is_theta_coprime((3, 3), (1, -1))
            False

        """
        if theta is None:
            theta = self.canonical_stability_parameter(d)

        assert self._is_dimension_vector(d), "``d`` needs to be a dimension vector"
        assert self._is_vector(theta), "`theta` needs to be a stability parameter"

        vectors = self.all_subdimension_vectors(d, proper=True, nonzero=True)

        return all(self.slope(d, theta) != self.slope(e, theta) for e in vectors)

    def is_indivisible(self, d) -> bool:
        """
        Checks if the gcd of all the entries of ``d`` is 1

        INPUT:

        -- ``d`` -- dimension vector

        OUTPUT: whether the dimension vector is indivisible

        EXAMPLES:

        Two examples with the Kronecker quiver::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: Q.is_indivisible((2, 3))
            True
            sage: Q.is_indivisible((2, 2))
            False

        """
        return gcd(self._coerce_dimension_vector(d)) == 1

    def support(self, d):
        r"""Returns the support of the dimension vector.

        INPUT:

        - ``d``: dimension vector

        OUTPUT: subset of vertices in the underlying graph in the support

        The support is the set :math:`\{ i \in Q_0 \mid d_i > 0 \}`.

        EXAMPLES:

        The support is the set of vertices for which the value of the dimension
        vector is nonzero::

            sage: from quiver import *
            sage: Q = ThreeVertexQuiver(2, 0, 4)
            sage: d = (1, 1, 1)
            sage: Q.support(d)
            [0, 1, 2]
            sage: d = (1, 0, 1)
            sage: Q.support(d)
            [0, 2]

        It takes into account vertex labels::

            sage: Q = Quiver.from_string("a--b----c,a---c", forget_labels=False)
            sage: d = {"a": 2, "b": 3, "c": 0}
            sage: Q.support(d)
            ['a', 'b']

        """
        assert self._is_dimension_vector(d), "``d`` needs to be a dimension vector"

        return [i for i in self.vertices() if d[i] > 0]

    def in_fundamental_domain(self, d, depth=0):
        r"""Checks if a dimension vector is in the fundamental domain.

        The fundamental domain of :math:`Q` is the set of dimension vectors :math:`d`
        such that

        - :math:`\operatorname{supp}(\mathbf{d})` is connected
        - :math:`\langle d,e_i\rangle + \langle e_i,d\rangle\leq 0` for every simple
          root

        Every :math:`d` in the fundamental domain is an imaginary root and the set of
        imaginary roots is the Weyl group saturation of the fundamental domain.
        If :math:`d` is in the fundamental domain then it is Schurian and a general
        representation of dimension vector :math:`d` is stable for the canonical
        stability parameter.

        The optional parameter ``depth`` allows to make the inequality stricter.

        INPUT:

        - ``d``: dimension vector

        - ``depth`` (default: 0) -- how deep the vector should be in the domain

        OUTPUT: whether ``d`` is in the (interior of) the fundamental domain

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
            sage: Q.in_fundamental_domain((1, 1), depth=1)
            True
            sage: Q.in_fundamental_domain((2, 3), depth=1)
            False

        """
        assert self._is_dimension_vector(d), "``d`` needs to be a dimension vector"

        # check if `\langle d,e_i\rangle + \langle e_i,d\rangle \leq 0`
        # for all vertices `i\in Q_0`
        inequality = all(
            self.symmetrized_euler_form(d, self.simple_root(i)) <= -depth
            for i in self.vertices()
        )

        # check if the support is connected
        connected = self.full_subquiver(self.support(d)).is_connected()

        return inequality and connected

    def division_order(self, d, e):
        r"""
        Checks if :math:`d\ll e`

        This means that

        - :math:`d_i \leq e_i` for every source `i`,
        - :math:`d_j \geq e_j` for every sink `j`, and
        - :math:`d_k = e_k` for every vertex `k` which is neither a source nor a sink.

        This is used when dealing with Chow rings of quiver moduli, see also
        :meth:`QuiverModuli.chow_ring` and
        :meth:`QuiverModuli._all_minimal_forbidden_subdimension_vectors`.

        EXAMPLES:

        The division order on some dimension vectors for the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: d = (1, 1)
            sage: e = (2, 1)
            sage: f = (2, 2)
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

        The division order on some dimension vectors for a 3-vertex quiver::

            sage: Q = ThreeVertexQuiver(2, 2, 2)
            sage: d = (1, 1, 1)
            sage: e = (1, 2, 1)
            sage: Q.division_order(d, e)
            False
            sage: Q.division_order(e, d)
            False

        """
        d = self._coerce_dimension_vector(d)
        e = self._coerce_dimension_vector(e)

        return (
            all(d[i] <= e[i] for i in self.sources())
            and all(d[i] >= e[i] for i in self.sinks())
            and all(
                d[i] == e[i]
                for i in self.vertices()
                if i not in self.sources() and i not in self.sinks()
            )
        )

    """
    Generic subdimension vectors and generic Hom and Ext
    """

    def is_generic_subdimension_vector(self, e, d) -> bool:
        r"""Checks if e is a generic subdimension vector of d.

        INPUT:

        - ``e``: dimension vector for the subrepresentation

        - ``d``: dimension vector for the ambient representation

        OUTPUT: whether e is a generic subdimension vector of d

        A dimension vector `e` is a generic subdimension vector of `d`
        if a generic representation of dimension vector `d` possesses
        a subrepresentation of dimension vector `e`.
        By MR1162487_ `e` is a generic subdimension vector of `d` if and only if `e` is
        a subdimension vector of `d` and :math:`\langle f,d-e\rangle` is non-negative
        for all generic subdimension vectors `f` of `e`.

        .. _MR1162487: https://mathscinet.ams.org/mathscinet/relay-station?mr=1162487

        EXAMPLES:

        Some examples on loop quivers::

            sage: from quiver import *
            sage: Q = LoopQuiver(1)
            sage: ds = [vector([i]) for i in range(3)]
            sage: for (e, d) in cartesian_product([ds, ds]):
            ....:     if not Q.is_subdimension_vector(e, d): continue
            ....:     print("{} is generic subdimension vector of {}: {}".format(
            ....:         e, d, Q.is_generic_subdimension_vector(e,d))
            ....:     )
            (0) is generic subdimension vector of (0): True
            (0) is generic subdimension vector of (1): True
            (0) is generic subdimension vector of (2): True
            (1) is generic subdimension vector of (1): True
            (1) is generic subdimension vector of (2): True
            (2) is generic subdimension vector of (2): True
            sage: Q = LoopQuiver(2)
            sage: for (e, d) in cartesian_product([ds]*2):
            ....:     if not Q.is_subdimension_vector(e, d): continue
            ....:     print("{} is generic subdimension vector of {}: {}".format(
            ....:         e, d, Q.is_generic_subdimension_vector(e,d))
            ....:     )
            (0) is generic subdimension vector of (0): True
            (0) is generic subdimension vector of (1): True
            (0) is generic subdimension vector of (2): True
            (1) is generic subdimension vector of (1): True
            (1) is generic subdimension vector of (2): False
            (2) is generic subdimension vector of (2): True

        Some examples on generalized Kronecker quivers::

            sage: Q = GeneralizedKroneckerQuiver(1)
            sage: ds = Tuples(range(3), 2)
            sage: for (e, d) in cartesian_product([ds]*2):
            ....:     if not Q.is_subdimension_vector(e, d): continue
            ....:     print("{} is generic subdimension vector of {}: {}".format(
            ....:         e, d, Q.is_generic_subdimension_vector(e,d))
            ....:     )
            (0, 0) is generic subdimension vector of (0, 0): True
            (0, 0) is generic subdimension vector of (1, 0): True
            (0, 0) is generic subdimension vector of (2, 0): True
            (0, 0) is generic subdimension vector of (0, 1): True
            (0, 0) is generic subdimension vector of (1, 1): True
            (0, 0) is generic subdimension vector of (2, 1): True
            (0, 0) is generic subdimension vector of (0, 2): True
            (0, 0) is generic subdimension vector of (1, 2): True
            (0, 0) is generic subdimension vector of (2, 2): True
            (1, 0) is generic subdimension vector of (1, 0): True
            (1, 0) is generic subdimension vector of (2, 0): True
            (1, 0) is generic subdimension vector of (1, 1): False
            (1, 0) is generic subdimension vector of (2, 1): True
            (1, 0) is generic subdimension vector of (1, 2): False
            (1, 0) is generic subdimension vector of (2, 2): False
            (2, 0) is generic subdimension vector of (2, 0): True
            (2, 0) is generic subdimension vector of (2, 1): False
            (2, 0) is generic subdimension vector of (2, 2): False
            (0, 1) is generic subdimension vector of (0, 1): True
            (0, 1) is generic subdimension vector of (1, 1): True
            (0, 1) is generic subdimension vector of (2, 1): True
            (0, 1) is generic subdimension vector of (0, 2): True
            (0, 1) is generic subdimension vector of (1, 2): True
            (0, 1) is generic subdimension vector of (2, 2): True
            (1, 1) is generic subdimension vector of (1, 1): True
            (1, 1) is generic subdimension vector of (2, 1): True
            (1, 1) is generic subdimension vector of (1, 2): True
            (1, 1) is generic subdimension vector of (2, 2): True
            (2, 1) is generic subdimension vector of (2, 1): True
            (2, 1) is generic subdimension vector of (2, 2): False
            (0, 2) is generic subdimension vector of (0, 2): True
            (0, 2) is generic subdimension vector of (1, 2): True
            (0, 2) is generic subdimension vector of (2, 2): True
            (1, 2) is generic subdimension vector of (1, 2): True
            (1, 2) is generic subdimension vector of (2, 2): True
            (2, 2) is generic subdimension vector of (2, 2): True
            sage: Q = GeneralizedKroneckerQuiver(2)
            sage: for (e, d) in cartesian_product([ds]*2):
            ....:     if not Q.is_subdimension_vector(e, d): continue
            ....:     print("{} is generic subdimension vector of {}: {}".format(
            ....:         e, d, Q.is_generic_subdimension_vector(e,d))
            ....:     )
            (0, 0) is generic subdimension vector of (0, 0): True
            (0, 0) is generic subdimension vector of (1, 0): True
            (0, 0) is generic subdimension vector of (2, 0): True
            (0, 0) is generic subdimension vector of (0, 1): True
            (0, 0) is generic subdimension vector of (1, 1): True
            (0, 0) is generic subdimension vector of (2, 1): True
            (0, 0) is generic subdimension vector of (0, 2): True
            (0, 0) is generic subdimension vector of (1, 2): True
            (0, 0) is generic subdimension vector of (2, 2): True
            (1, 0) is generic subdimension vector of (1, 0): True
            (1, 0) is generic subdimension vector of (2, 0): True
            (1, 0) is generic subdimension vector of (1, 1): False
            (1, 0) is generic subdimension vector of (2, 1): False
            (1, 0) is generic subdimension vector of (1, 2): False
            (1, 0) is generic subdimension vector of (2, 2): False
            (2, 0) is generic subdimension vector of (2, 0): True
            (2, 0) is generic subdimension vector of (2, 1): False
            (2, 0) is generic subdimension vector of (2, 2): False
            (0, 1) is generic subdimension vector of (0, 1): True
            (0, 1) is generic subdimension vector of (1, 1): True
            (0, 1) is generic subdimension vector of (2, 1): True
            (0, 1) is generic subdimension vector of (0, 2): True
            (0, 1) is generic subdimension vector of (1, 2): True
            (0, 1) is generic subdimension vector of (2, 2): True
            (1, 1) is generic subdimension vector of (1, 1): True
            (1, 1) is generic subdimension vector of (2, 1): True
            (1, 1) is generic subdimension vector of (1, 2): False
            (1, 1) is generic subdimension vector of (2, 2): True
            (2, 1) is generic subdimension vector of (2, 1): True
            (2, 1) is generic subdimension vector of (2, 2): False
            (0, 2) is generic subdimension vector of (0, 2): True
            (0, 2) is generic subdimension vector of (1, 2): True
            (0, 2) is generic subdimension vector of (2, 2): True
            (1, 2) is generic subdimension vector of (1, 2): True
            (1, 2) is generic subdimension vector of (2, 2): True
            (2, 2) is generic subdimension vector of (2, 2): True

        """
        return self.__is_generic_subdimension_vector(e, d)

    @cached_method(
        key=lambda self, e, d: (self._coerce_vector(e), self._coerce_vector(d))
    )
    def __is_generic_subdimension_vector(self, e, d) -> bool:
        r"""
        The cacheable implementation of :meth:`Quiver.is_generic_subdimension_vector`.

        EXAMPLES:

        Generic subdimension vectors for the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q, d, theta = GeneralizedKroneckerQuiver(3), (2, 3), (3, -2)
            sage: for e in Q.all_subdimension_vectors(d):
            ....:     print("{} is generic subdimension vector of {}: {}".format(
            ....:         e, d, Q._Quiver__is_generic_subdimension_vector(e, d))
            ....:     )
            (0, 0) is generic subdimension vector of (2, 3): True
            (0, 1) is generic subdimension vector of (2, 3): True
            (0, 2) is generic subdimension vector of (2, 3): True
            (0, 3) is generic subdimension vector of (2, 3): True
            (1, 0) is generic subdimension vector of (2, 3): False
            (1, 1) is generic subdimension vector of (2, 3): False
            (1, 2) is generic subdimension vector of (2, 3): True
            (1, 3) is generic subdimension vector of (2, 3): True
            (2, 0) is generic subdimension vector of (2, 3): False
            (2, 1) is generic subdimension vector of (2, 3): False
            (2, 2) is generic subdimension vector of (2, 3): False
            (2, 3) is generic subdimension vector of (2, 3): True
        """
        d = self._coerce_dimension_vector(d)
        e = self._coerce_dimension_vector(e)

        if e == d or all(ei == 0 for ei in e):
            return True

        if not self.is_subdimension_vector(e, d):
            return False

        ds = filter(
            lambda eprime: self.euler_form(eprime, d - e) < 0,
            self.all_subdimension_vectors(e),
        )

        return not any(self.is_generic_subdimension_vector(eprime, e) for eprime in ds)

    def all_generic_subdimension_vectors(self, d, proper=False, nonzero=False):
        r"""Returns the list of all generic subdimension vectors of ``d``.

        INPUT:

        - ``d``: dimension vector

        OUTPUT: list of vectors

        EXAMPLES:

        Some n-Kronecker quivers::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(1)
            sage: d = (3, 3)
            sage: Q.all_generic_subdimension_vectors(d)
            [(0, 0),
             (0, 1),
             (0, 2),
             (0, 3),
             (1, 1),
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
             (0, 3),
             (1, 1),
             (1, 2),
             (1, 3),
             (2, 2),
             (2, 3),
             (3, 3)]
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.all_generic_subdimension_vectors(d)
            [(0, 0), (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3), (3, 3)]
            sage: Q.all_generic_subdimension_vectors(d, nonzero=True)
            [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3), (3, 3)]
            sage: Q.all_generic_subdimension_vectors(d, proper=True)
            [(0, 0), (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

        """
        d = self._coerce_dimension_vector(d)

        return list(
            filter(
                lambda e: self.is_generic_subdimension_vector(e, d),
                self.all_subdimension_vectors(
                    d, proper=proper, nonzero=nonzero, forget_labels=True
                ),
            )
        )

    def generic_ext(self, d, e):
        r"""
        Computes :math:`\operatorname{ext}(d, e)`.

        INPUT:

        - ``d``: dimension vector

        - ``e``: dimension vector

        OUTPUT: dimension of the generic ext

        According to Theorem 5.4 in Schofield's 'General representations of quivers',
        we have

        .. MATH::

            \operatorname{ext}(a,b) =
            \operatorname{max}\{-\langle c,b\rangle\},

        where :math:`c` runs over the generic subdimension vectors of :math:`a`.

        EXAMPLES:

        Generic ext on the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: ds = [Q.simple_root(0), Q.simple_root(1), Q.thin_dimension_vector()]
            sage: for (d, e) in cartesian_product([ds]*2):
            ....:     print("ext({}, {}) = {}".format(d, e, Q.generic_ext(d, e)))
            ext((1, 0), (1, 0)) = 0
            ext((1, 0), (0, 1)) = 3
            ext((1, 0), (1, 1)) = 2
            ext((0, 1), (1, 0)) = 0
            ext((0, 1), (0, 1)) = 0
            ext((0, 1), (1, 1)) = 0
            ext((1, 1), (1, 0)) = 0
            ext((1, 1), (0, 1)) = 2
            ext((1, 1), (1, 1)) = 1

        """
        d = self._coerce_dimension_vector(d)
        e = self._coerce_dimension_vector(e)

        return max(
            -self.euler_form(f, e) for f in self.all_generic_subdimension_vectors(d)
        )

    def generic_hom(self, d, e):
        r"""
        Computes :math:`\operatorname{hom}(d, e)`.

        INPUT:

        - ``d``: dimension vector

        - ``e``: dimension vector

        OUTPUT: dimension of the generic hom

        There is a non-empty open subset `U` of :math:`R(Q,d) \times R(Q,e)` such that

        .. MATH::

            \operatorname{dim}(\operatorname{Ext}(M,N)) = \operatorname{ext}(d,e),

        i.e., :math:`\operatorname{dim}(\operatorname{Ext}(M,N))` is minimal for all
        `(M,N)` in `U`.

        Therefore, :math:`\operatorname{dim}(\operatorname{Hom}(M,N)) =
        \langle a,b\rangle + \operatorname{dim}(\operatorname{Ext}(M,N))`
        is minimal, and
        :math:`\operatorname{hom}(a,b) = \langle a,b\rangle + \operatorname{ext}(a,b)`.

        EXAMPLES:

        Generic hom on the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: ds = [Q.simple_root(0), Q.simple_root(1), Q.thin_dimension_vector()]
            sage: for (d, e) in cartesian_product([ds]*2):
            ....:     print("hom({}, {}) = {}".format(d, e, Q.generic_hom(d, e)))
            hom((1, 0), (1, 0)) = 1
            hom((1, 0), (0, 1)) = 0
            hom((1, 0), (1, 1)) = 0
            hom((0, 1), (1, 0)) = 0
            hom((0, 1), (0, 1)) = 1
            hom((0, 1), (1, 1)) = 1
            hom((1, 1), (1, 0)) = 1
            hom((1, 1), (0, 1)) = 0
            hom((1, 1), (1, 1)) = 0

        """
        d = self._coerce_dimension_vector(d)
        e = self._coerce_dimension_vector(e)

        return self.euler_form(d, e) + self.generic_ext(d, e)

    """
    Harder--Narasimhan types
    """

    @cached_method
    def _all_harder_narasimhan_types(self, d, theta, denom=sum, sorted=False):
        r"""Returns the list of all Harder--Narasimhan types of d.

        INPUT:

        - ``d`` -- dimension vector

        - ``theta` -- stability parameter

        - ``denom`` -- the denominator function (default: sum)

        OUTPUT: list of Harder--Narasimhan types

        EXAMPLES:

        The Harder--Narasimhan types for the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: d = (2, 3)
            sage: theta = (3, -2)
            sage: Q._all_harder_narasimhan_types(d, theta)
            [((1, 0), (1, 1), (0, 2)),
             ((1, 0), (1, 2), (0, 1)),
             ((1, 0), (1, 3)),
             ((1, 1), (1, 2)),
             ((2, 0), (0, 3)),
             ((2, 1), (0, 2)),
             ((2, 2), (0, 1)),
             ((2, 3),)]

        .. NOTE ::

        This is a method of Quiver so that its results can be more efficiently cached.
        See :meth:`QuiverModuli.all_harder_narasimhan_types()` for the better location
        to call it from.

        """
        d = self._coerce_dimension_vector(d)
        theta = self._coerce_vector(theta)

        ds = self.all_subdimension_vectors(d, proper=True, nonzero=True)
        ds = filter(
            lambda e: self.slope(e, theta, denom=denom)
            > self.slope(d, theta, denom=denom),
            ds,
        )
        ds = filter(
            lambda e: self.has_semistable_representation(e, theta, denom=denom),
            ds,
        )
        ds = list(ds)

        if sorted:
            ds.sort(key=(lambda e: self.slope(e, theta, denom=denom)))

        all_types = []
        for e in ds:
            for estar in filter(
                lambda fstar: self.slope(e, theta, denom=denom)
                > self.slope(fstar[0], theta, denom=denom),
                self._all_harder_narasimhan_types(d - e, theta, denom=denom),
            ):
                all_types.append((e,) + estar)

        if self.has_semistable_representation(d, theta, denom=denom):
            all_types.append((d,))

        return all_types

    """
    (Semi-)stability
    """

    def canonical_stability_parameter(self, d):
        r"""
        Returns the canonical stability parameter for ``d``

        INPUT:

        - ``d``: dimension vector

        OUTPUT: canonical stability parameter

        The canonical stability parameter is given by
        :math:`\langle d,-\rangle - \langle -,d\rangle`.

        EXAMPLES:

        Our usual example of the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: Q.canonical_stability_parameter((2, 3))
            (9, -6)

        For the 5-subspace quiver::

            sage: Q = SubspaceQuiver(5)
            sage: Q.canonical_stability_parameter((1, 1, 1, 1, 1, 2))
            (2, 2, 2, 2, 2, -5)

        It takes vertex labels (if present) into account::

            sage: Q = Quiver.from_string("foo---bar", forget_labels=False)
            sage: Q.canonical_stability_parameter((2, 3))
            {'bar': -6, 'foo': 9}

        EXAMPLES:

        Canonical stability parameter for the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q, d = GeneralizedKroneckerQuiver(3), (2, 3)
            sage: Q.canonical_stability_parameter(d)
            (9, -6)

        This method also works with vertex labels::

            sage: from quiver import *
            sage: Q = Quiver.from_string("foo---bar", forget_labels=False)
            sage: d = {"foo": 2, "bar": 3}
            sage: Q.canonical_stability_parameter(d)
            {'bar': -6, 'foo': 9}
        """
        d = self._coerce_dimension_vector(d)
        theta = vector(d) * (-self.euler_matrix().transpose() + self.euler_matrix())

        if self.__has_vertex_labels():
            return dict(zip(self.vertices(), theta))
        return theta

    def has_semistable_representation(self, d, theta=None, denom=sum):
        r"""Checks if there is a ``theta``-semistable of dimension vector ``d``

        INPUT:

        - ``d``: dimension vector

        - ``theta`` (default: canonical stability parameter): stability parameter

        OUTPUT: whether there is a ``theta``-semistable of dimension vector ``d``

        By MR1162487_ a dimension vector `d` admits a :math:`\theta`-semi-stable
        representation if and only if :math:`\mu_{\theta}(e) \leq \mu_{\theta}(d)` for
        all generic subdimension vectors `e` of `d`.

        .. _MR1162487: https://mathscinet.ams.org/mathscinet/relay-station?mr=1162487

        EXAMPLES:

        Semistables for the :math:`\mathrm{A}_2` quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(1)
            sage: Q.has_semistable_representation((1, 1), (1, -1))
            True
            sage: Q.has_semistable_representation((2, 2), (1, -1))
            True
            sage: Q.has_semistable_representation((1, 2), (1, -1))
            False
            sage: Q.has_semistable_representation((0, 0), (1, -1))
            True

        Semistables for the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.has_semistable_representation((2, 3))
            True
            sage: Q.has_semistable_representation((1, 4), (-3, 2))
            False

        """
        if theta is None:
            theta = self.canonical_stability_parameter(d)

        d = self._coerce_dimension_vector(d)
        theta = self._coerce_vector(theta)

        return all(
            self.slope(e, theta, denom=denom) <= self.slope(d, theta, denom=denom)
            for e in self.all_generic_subdimension_vectors(d, nonzero=True)
        )

    def has_stable_representation(self, d, theta=None, denom=sum):
        r"""
        Checks if there is a ``theta``-stable representation of ``d``

        INPUT:

        - ``d``: dimension vector

        - ``theta`` (default: canonical stability parameter): stability parameter

        OUTPUT: whether there is a ``theta``-stable of dimension vector ``d``

        By MR1162487_ `d` admits a theta-stable representation if and only if
        :math:`\mu_{\theta}(e) < \mu_{\theta}(d)` for all proper generic subdimension
        vectors :math:`e` of :math:`d`.

        .. _MR1162487: https://mathscinet.ams.org/mathscinet/relay-station?mr=1162487

        EXAMPLES:

        Stables for the :math:`\mathrm{A}_2` quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(1)
            sage: theta = (1, -1)
            sage: Q.has_stable_representation((1, 1), theta)
            True
            sage: Q.has_stable_representation((2, 2), theta)
            False
            sage: Q.has_stable_representation((0, 0), theta)
            False

        Stables for the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: d = (2, 3)
            sage: theta = Q.canonical_stability_parameter(d)
            sage: Q.has_stable_representation(d, theta)
            True
            sage: Q.has_stable_representation(d)
            True

        """
        if theta is None:
            theta = self.canonical_stability_parameter(d)

        d = self._coerce_dimension_vector(d)
        theta = self._coerce_vector(theta)

        if d == self._coerce_dimension_vector(self.zero_vector()):
            return False

        return all(
            self.slope(e, theta, denom=denom) < self.slope(d, theta, denom=denom)
            for e in self.all_generic_subdimension_vectors(d, proper=True, nonzero=True)
        )

    """
    Canonical decomposition
    """

    def canonical_decomposition(self, d):
        r"""
        Computes the canonical decomposition of a dimension vector.

        INPUT:

        - ``d``: dimension vector

        OUTPUT: canonical decomposition as list of dimension vectors

        The canonical decomposition of a dimension vector `d` is the unique
        decomposition :math:`d = e_1 + e_2 + ... + e_k` such that
        :math:`e_1, e_2, ..., e_k` are such that for all
        :math:`i \neq j, \mathrm{ext}(e_i, e_j) = \mathrm{ext}(e_j, e_i) = 0`.

        The general representation of dimension vector `d` is isomorphic to the direct
        sum of representations of dimension vectors :math:`e_1, e_2, ..., e_k`.

        EXAMPLES:

        Canonical decomposition of the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q.canonical_decomposition((2, 3))
            [(2, 3)]
            sage: for d in Q.all_subdimension_vectors((5, 5)):
            ....:     print(Q.canonical_decomposition(d))
            [(0, 0)]
            [(0, 1)]
            [(0, 1), (0, 1)]
            [(0, 1), (0, 1), (0, 1)]
            [(0, 1), (0, 1), (0, 1), (0, 1)]
            [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1)]
            [(1, 0)]
            [(1, 1)]
            [(1, 2)]
            [(1, 3)]
            [(0, 1), (1, 3)]
            [(0, 1), (0, 1), (1, 3)]
            [(1, 0), (1, 0)]
            [(2, 1)]
            [(2, 2)]
            [(2, 3)]
            [(2, 4)]
            [(2, 5)]
            [(1, 0), (1, 0), (1, 0)]
            [(3, 1)]
            [(3, 2)]
            [(3, 3)]
            [(3, 4)]
            [(3, 5)]
            [(1, 0), (1, 0), (1, 0), (1, 0)]
            [(1, 0), (3, 1)]
            [(4, 2)]
            [(4, 3)]
            [(4, 4)]
            [(4, 5)]
            [(1, 0), (1, 0), (1, 0), (1, 0), (1, 0)]
            [(1, 0), (1, 0), (3, 1)]
            [(5, 2)]
            [(5, 3)]
            [(5, 4)]
            [(5, 5)]
        """
        return self.__canonical_decomposition(d)

    @cached_method(key=lambda self, d: self._coerce_vector(d))
    def __canonical_decomposition(self, d):
        r"""The cacheable implementation of :meth:`Quiver.canonical_decomposition`

        EXAMPLES:

        Canonical decomposition of `(5, 3)` for the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(2)
            sage: Q._Quiver__canonical_decomposition((5, 5))
            [(1, 1), (1, 1), (1, 1), (1, 1), (1, 1)]
        """
        d = self._coerce_dimension_vector(d)

        ds = self.all_generic_subdimension_vectors(d, proper=True, nonzero=True)
        for e in ds:
            if d - e in ds:
                return self.canonical_decomposition(e) + self.canonical_decomposition(
                    d - e
                )
        return [d]

    def dimension_nullcone(self, d):
        r"""
        Returns the dimension of the nullcone

        The nullcone is the set of all nilpotent representations.

        INPUT:

        - ``d`` -- dimension vector

        OUTPUT: dimension of the nullcone

        EXAMPLES:

        The usual example of the 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: Q.dimension_nullcone((2, 3))
            18

        """
        d = self._coerce_dimension_vector(d)

        if self.is_acyclic():
            return d * self.adjacency_matrix() * d
        else:
            raise NotImplementedError()

    def first_hochschild_cohomology(self):
        r"""
        Compute the dimension of the first Hochschild cohomology

        This uses the formula of Happel from Proposition 1.6 in MR1035222_.
        One needs the quiver to be acyclic for this, otherwise it is not necessarily
        finite-dimensional.

        EXAMPLES:

        The first Hochschild cohomology of the `m`-th generalized Kronecker quiver
        is the dimension of :math:`\mathrm{PGL}_{m+1}`::

            sage: from quiver import *
            sage: GeneralizedKroneckerQuiver(3).first_hochschild_cohomology()
            8

        The first Hochschild cohomology vanishes if and only if the quiver is a tree::

            sage: from quiver import *
            sage: SubspaceQuiver(7).first_hochschild_cohomology()
            0

        .. _MR1035222: https://mathscinet.ams.org/mathscinet/relay-station?mr=1035222
        """
        assert self.is_acyclic(), "the quiver needs to be acyclic"

        return (
            1
            - self.number_of_vertices()
            + sum(
                len(self.graph().all_paths(a[0], a[1], use_multiedges=True))
                for a in self.graph().edges()
            )
        )
