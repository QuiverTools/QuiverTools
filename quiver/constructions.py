from sage.combinat.root_system.dynkin_diagram import DynkinDiagram
from sage.combinat.root_system.weyl_characters import WeylCharacterRing
from sage.matrix.special import block_diagonal_matrix, random_matrix, zero_matrix
from sage.misc.prandom import randint
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ

from quiver import Quiver


def disjoint_union(Q1, Q2):
    """Returns the disjoint union of two quivers Q1 and Q2.

    EXAMPLES:

    We construct the disjoint union of 2 generalized Kronecker quivers::

        sage: from quiver import *
        sage: Q1 = GeneralizedKroneckerQuiver(3)
        sage: Q2 = GeneralizedKroneckerQuiver(4)
        sage: Q = disjoint_union(Q1,Q2)
        sage: Q
        disjoint union of 3-Kronecker quiver and 4-Kronecker quiver

    """

    if Q1.get_custom_name() and Q2.get_custom_name():
        name = (
            "disjoint union of " + Q1.get_custom_name() + " and " + Q2.get_custom_name()
        )
    else:
        name = None

    return Quiver(
        block_diagonal_matrix(
            Q1.adjacency_matrix(), Q2.adjacency_matrix(), subdivide=False
        ),
        name=name,
    )


"""Special quivers"""


def GeneralizedKroneckerQuiver(m: int):
    r"""
    Return the generalized Kronecker quiver

    The generalized Kronecker quiver has two vertices and $m$ arrows
    from the first to the second vertex.

    INPUT:

    - ``m`` -- integer; number of arrows in the quiver

    OUTPUT: the generalized Kronecker quiver as Quiver instance

    TESTS::

    The generalized Kronecker quiver is as claimed::

        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: Q.number_of_vertices()
        2
        sage: Q.number_of_arrows()
        3

    """
    return Quiver([[0, m], [0, 0]], name="{}-Kronecker quiver".format(m))


def KroneckerQuiver(m: int = 2):
    r"""
    Return the Kronecker quiver

    The Kronecker quiver has two vertices and 2 arrow from the first
    to the second vertex. If the optional parameter ``m`` is specified
    we construct the generalized Kronecker quiver on ``m`` arrows;

    INPUT:

    - ``m`` -- integer (default: `2`; number of arrows in the quiver

    OUTPUT: the Kronecker quiver as Quiver instance

    TESTS::

    The Kronecker quiver is as claimed::

        sage: from quiver import *
        sage: Q = KroneckerQuiver()
        sage: Q.number_of_vertices()
        2
        sage: Q.number_of_arrows()
        2

    If we specify the number of arrows, we construct a generalized Kronecker quiver::

        sage: KroneckerQuiver(3) == GeneralizedKroneckerQuiver(3)
        True

    """

    return GeneralizedKroneckerQuiver(m)


def ThreeVertexQuiver(m12: int, m13: int, m23: int):
    r"""
    Constructs a 3-vertex quiver, with and `m_{i,j} arrows from `i` to `j`.

    Thus it is always an acyclic quiver.

    INPUT:

    - ``m12`` -- integer; number of arrows from 1 to 2

    - ``m13`` -- integer; number of arrows from 1 to 3

    - ``m23`` -- integer; number of arrows from 2 to 3

    EXAMPLES::

    A 3-vertex quiver with 5 arrows::

        sage: from quiver import *
        sage: Q = ThreeVertexQuiver(2, 2, 1); Q
        an acyclic 3-vertex quiver of type (2, 2, 1)
        sage: Q.number_of_arrows()
        5

    """

    Q = Quiver(
        [[0, m12, m13], [0, 0, m23], [0, 0, 0]],
        name="an acyclic 3-vertex quiver of type {}".format(
            (m12, m13, m23),
        ),
    )
    return Q


def LoopQuiver(m: int):
    r"""
    Return the quiver with 1 vertex and `m` loops.

    This is a synonym for :func:`GeneralizedJordanQuiver`.

    .. SEEALSO::

        :func:`GeneralizedJordanQuiver`

    """
    Q = GeneralizedJordanQuiver(m)
    Q.rename("{}-loop quiver".format(m))

    return Q


def JordanQuiver(m: int = 1):
    r"""
    Return the generalized Jordan quiver with `m` loops, where `m` has default value 1.

    .. SEEALSO::

        :func:`GeneralizedJordanQuiver`

    """

    return GeneralizedJordanQuiver(m)


def GeneralizedJordanQuiver(m: int):
    r"""
    Return the generalized Jordan quiver with `m` loops

    INPUT:

    - ``m`` -- integer; the number of loops in the generalized Jordan quiver

    OUTPUT: the generalized Jordan quiver

    EXAMPLES::

    The generalized Jordan quiver has 1 vertex and `m` arrows::

        sage: from quiver import *
        sage: Q = GeneralizedJordanQuiver(7)
        sage: Q.number_of_vertices()
        1
        sage: Q.number_of_arrows()
        7

    """
    Q = Quiver([[m]], name="generalized Jordan quiver with {} loops".format(m))

    if m == 1:
        Q.rename("Jordan quiver")

    return Q


def SubspaceQuiver(m: int):
    r"""
    Return the subspace quiver with `m` sources

    The sources are labelled `1,\ldots,m` and the sink is `m+1`; there are is one arrow
    from every source to the sink.

    INPUT:

    - ``m`` -- integer; the number of sources in the subspace quiver

    OUTPUT: the subspace quiver with `m` sources

    EXAMPLES:

    The subspace quiver with `m` sources has `m` arrows and `m+1` vertices::

        sage: from quiver import *
        sage: Q = SubspaceQuiver(6)
        sage: Q.number_of_vertices()
        7
        sage: Q.number_of_arrows()
        6

    The subspace quiver with 2 sources is also a 3-vertex quiver::

        sage: SubspaceQuiver(2) == ThreeVertexQuiver(0, 1, 1)
        True

    """
    M = zero_matrix(ZZ, m + 1)
    for i in range(m):
        M[i, m] = 1

    Q = Quiver(M, "{}-subspace quiver".format(m))

    return Q


def ThickenedSubspaceQuiver(m, k):
    r"""
    Return the thickened subspace quiver with `m` sources

    The sources are labelled `1,\ldots,m` and the sink is `m+1`; there are are `k` arrows
    from every source to the sink.

    - ``m`` -- integer; the number of sources in the subspace quiver

    - ``k`` -- integer; the number arrows from a source to the sink

    OUTPUT: the subspace quiver with `m` sources and `k` arrows from each source

    EXAMPLES:

    The `k`-thickened subspace quiver with `m` sources has `km` arrows and `m+1` vertices::

        sage: from quiver import *
        sage: Q = ThickenedSubspaceQuiver(6, 2)
        sage: Q.number_of_vertices()
        7
        sage: Q.number_of_arrows()
        12

    The `k`-thickened subspace quiver with 2 sources is also a 3-vertex quiver::

        sage: ThickenedSubspaceQuiver(2, 6) == ThreeVertexQuiver(0, 6, 6)
        True

    """
    Q = GeneralizedSubspaceQuiver(m, [k] * m)
    Q.rename(
        "thickened subspace quiver with {} sources and multiplicity {}".format(m, k)
    )

    return Q


def GeneralizedSubspaceQuiver(m, K):
    r"""
    Return the generalized subspace quiver with `m` sources and multiplicities `K`

    The sources are labelled `1,\ldots,m` and the sink is `m+1`; there are are `K_i` arrows
    from every source `i=1,\ldots,m` to the sink.

    - ``m`` -- integer; the number of sources in the subspace quiver

    - ``K`` -- list of integers; the number arrows from the `i`th source to the sink

    OUTPUT: the subspace quiver with `m` sources and `K_i` arrows from each source

    EXAMPLES:

    The generalized subspace quiver with `m` sources and multiplicities `K` has `\sum_{i=1}^mK_i`
    arrows and `m+1` vertices::

        sage: from quiver import *
        sage: Q = GeneralizedSubspaceQuiver(6, (1, 2, 3, 4, 5, 6))
        sage: Q.number_of_vertices()
        7
        sage: Q.number_of_arrows()
        21

    The generalized subspace quiver with 2 sources is also a 3-vertex quiver::

        sage: GeneralizedSubspaceQuiver(2, (2, 3)) == ThreeVertexQuiver(0, 2, 3)
        True

    """
    assert len(K) == m

    M = zero_matrix(ZZ, m + 1)
    for i in range(m):
        M[i, m] = K[i]

    Q = Quiver(M, name="a generalized {}-subspace quiver".format(m))

    return Q


def DynkinQuiver(T):
    r"""
    Return the Dynkin quiver of type `T`

    The type `T` is to be taken as in the Sage method `DynkinDiagram`,
    and the quiver is oriented lexigraphically in the vertices of the diagram.

    INPUT:

    - ``T``: a Dynkin type, as documented in the Sage method `DynkinDiagram`

    OUTPUT: the Dynkin quiver with lexicographic ordering on the vertices

    EXAMPLES:

    The `\mathrm{A}_2` quiver is the generalized Kronecker quiver with 1 arrow::

        sage: from quiver import *
        sage: DynkinQuiver("A2") == GeneralizedKroneckerQuiver(1)
        True

    The Dynkin quiver `\mathrm{D}_4` is a different orientation of the 3-subspace quiver::

        sage: DynkinQuiver("D4") == SubspaceQuiver(3)
        False

    We can also consider disconnected Dynkin quivers::

        sage: Q = DynkinQuiver("A3xA4")
        sage: Q.is_connected()
        False

    """
    M = DynkinDiagram(T).adjacency_matrix()

    # orient the quiver lexicographically
    for i in range(M.nrows()):
        for j in range(i):
            M[i, j] = 0

    return Quiver(M, "Dynkin quiver of type {}".format(T))


def ExtendedDynkinQuiver(T):
    r"""
    Return the Dynkin quiver of type `T`

    The type `T` is a string which can be passed onto the Sage method `DynkinDiagram`,
    and the quiver is oriented lexigraphically in the vertices of the diagram.
    The special vertex thus comes first.

    INPUT:

    - ``T`` -- string: a Dynkin type, as documented in the Sage method `DynkinDiagram`

    OUTPUT: the extended Dynkin quiver with lexicographic ordering on the vertices

    EXAMPLES:

    The extended `\mathrm{A}_1` quiver is the Kronecker quiver::

        sage: from quiver import *
        sage: ExtendedDynkinQuiver("A1") == KroneckerQuiver()
        True

    """
    # the adjacency matrix of an extended Dynkin diagram ignores multiple arrows
    # so we modify the Cartan matrix
    D = WeylCharacterRing(T).extended_dynkin_diagram()
    M = -D.cartan_matrix()

    # orient the quiver lexicographically
    for i in range(M.nrows()):
        M[i, i] = 0
        for j in range(i):
            M[i, j] = 0

    return Quiver(M, "Extended Dynkin quiver of type {}".format(T))


def CyclicQuiver(n):
    r"""
    Return the cyclic quiver on `n` vertices

    This is the quiver with `n` vertices and `n` arrows from `i` to `i+1`
    for `i=1,\ldots,n`, with `n+1=1`.

    INPUT:

    - ``n``-- integer; the number of vertices (and arrows)

    OUTPUT: cyclic quiver on `n` vertices

    EXAMPLES:

    The doubled Dynkin quiver of type `\mathrm{A}_2` is the cyclic quiver on 2 vertices::

        sage: from quiver import *
        sage: CyclicQuiver(2) == DynkinQuiver("A2").doubled_quiver()
        True

    """
    M = zero_matrix(n)
    for i in range(n):
        M[i, (i + 1) % n] = 1

    return Quiver(M, "Cyclic quiver on {} vertices".format(n))


def BipartiteQuiver(m, n):
    r"""
    Return the bipartite quiver with `m` sources and `n` sinks

    This is the quiver with `m+n` vertices, having 1 arrow from each of the first `m` vertices
    to the each of the last `n` vertices.

    INPUT:

    - ``m`` -- non-negative integer; number of sources

    - ``n`` -- non-negative integer; number of sinks

    OUTPUT: bipartite quiver with `m` sources and `n` sinks

    EXAMPLES:

    When `m=n=1` we get the `\mathrm{A}_2` quiver::

        sage: from quiver import *
        sage: BipartiteQuiver(1, 1) == DynkinQuiver("A2")
        True

    When `m=2` and `n=1` we get a 3-vertex quiver::

        sage: BipartiteQuiver(2, 1) == ThreeVertexQuiver(0, 1, 1)
        True

    """

    M = zero_matrix(m + n)
    for i in range(m):
        for j in range(n):
            M[i, m + j] = 1

    return Quiver(M, "bipartite quiver on {} sources and {} sinks".format(m, n))


# Sampling and testing methods


def RandomQuiver(vertices, arrow_bound=10, acyclic=False, connected=True):
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
            adjacency = random_matrix(ZZ, vertices, vertices, x=0, y=arrow_bound)

            if acyclic:
                # upper triangular matrix
                for i in range(vertices):
                    for j in range(i, vertices):
                        adjacency[j, i] = 0

            acceptable = Quiver(
                adjacency
            ).is_connected()  # unnecessary overhead in defining Quiver object
    elif not connected:
        adjacency = random_matrix(ZZ, vertices, vertices, x=0, y=arrow_bound)

        if acyclic:
            # upper triangular matrix
            for i in range(vertices):
                for j in range(i, vertices):
                    adjacency[j, i] = 0

    return Quiver(adjacency)


def RandomDimensionVector(quiver, positive=False, upper_bound=10):
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
    return vector(
        [randint(lower_bound, upper_bound) for i in range(quiver.number_of_vertices())]
    )


def RandomStability(quiver, bound=10):
    """Returns a random stability condition for the given quiver.
    Inputs:
        - quiver: a Quiver object;
        - bound: upper and lower bound on the entries. Defaults to 10.
    """
    # what other features should this have?

    return vector(
        [randint(-bound // 2, bound) for i in range(quiver.number_of_vertices())]
    )
