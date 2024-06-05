from quiver import *


def disjoint_union(Q1, Q2):
    """Returns the disjoint union of two quivers Q1 and Q2.

    EXAMPLES:

    We construct the disjoint union of 2 generalized Kronecker quivers::

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
    Q = Quiver(matrix([[0, m], [0, 0]]), name=str(m) + "-Kronecker quiver")
    # TODO do Q.rename here
    return Q


# TODO if optional parameter is given, call GeneralizedKroneckerQuiver
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


# TODO Pieter wants to change this order
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
        sage: Q = ThreeVertexQuiver(2, 2, 1)
        sage: Q.number_of_arrows()
        5

    """

    Q = Quiver(
        matrix([[0, m12, m13], [0, 0, m23], [0, 0, 0]]),
        name="An acyclic 3-vertex quiver",
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
    # TODO do rename
    # Q.rename("{}-loop quiver".format(m))

    return Q


def JordanQuiver(m: int = 1):
    r"""
    Return the generalized Jordan quiver with `m` loops, where `m` has default value 1.

    .. SEEALSO::

        :func:`GeneralizedJordanQuiver`

    """

    Q = GeneralizedJordanQuiver(m)

    # TODO do rename
    # if m == 1:
    #    Q.rename("Jordan quiver")

    return Q


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
    Q = Quiver(matrix([[m]]), name="generalized Jordan quiver with {} loops".format(m))

    return Q


def SubspaceQuiver(m: int):
    r"""
    Return the subspace quiver with `m` sources

    The sources are labelled `1,\ldots,m` and the sink is `m+1`; there are is one arrow
    from every source to the sink.

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
    A = zero_matrix(ZZ, m + 1)
    for i in range(m):
        A[i, m] = 1

    Q = Quiver(A, name=str(m) + "-subspace quiver")
    # TODO do Q.rename here
    return Q


def ThickenedSubspaceQuiver(m, k):
    r"""
    Return the thickened subspace quiver with `m` sources

    The sources are labelled `1,\ldots,m` and the sink is `m+1`; there are are `k` arrows
    from every source to the sink.

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
    Q = GeneralizedSubspaceQuiver(m, [k]*m)
    # TODO do Q.rename here
    return Q


def GeneralizedSubspaceQuiver(m, K):
    r"""
    Return the generalized subspace quiver with `m` sources and multiplicities `K`

    The sources are labelled `1,\ldots,m` and the sink is `m+1`; there are are `K_i` arrows
    from every source `i=1,\ldots,m` to the sink.

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

    A = zero_matrix(ZZ, m + 1)
    # I'm sure you can do this without a for loop
    for i in range(m):
        A[i, m] = K[i]

    Q = Quiver(A, name="A generalized " + str(m) + "-subspace quiver")
    # TODO do Q.rename here
    return Q


def DynkinQuiver(Tn):
    r"""Returns the Dynkin quiver of type Tn. Uses the standard Sagemath implementation of Dynkin diagrams."""
    # use https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/root_system/dynkin_diagram.html
    # TODO: this constructor calls the adjacency_matrix() method many times. Should we call it once and remove lower triangular entries?

    # parse the string Tn
    T = Tn[:-1]
    n = int(Tn[-1])

    return Quiver(
        matrix(
            n,
            n,
            lambda i, j: DynkinDiagram(Tn).adjacency_matrix()[i, j] if i < j else 0,
        ),
        "Dynkin quiver of type " + Tn,
    )


def ExtendedDynkinQuiver(T):
    # TODO implement this
    # TODO orientation: have a default (so for A and D: linear, type E?) but make it possible to change the orientation
    raise NotImplementedError()


def CyclicQuiver(n):
    return ExtendedDynkinQuiver(["A", n])


def BipartiteQuiver(m, n):
    # TODO implement this
    # m is number of sources, n is number of sinks
    raise NotImplementedError()


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
