from sage.matrix.constructor import matrix

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
        r""""Returns the (necessarily symmetric) adjacency matrix of the underlying graph of the quiver.

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
        return (A^n == zero_matrix(ZZ, n))

    def is_connected(self): 
        r""""Returns whether the underlying graph of the quiver is connected or not.

        OUTPUT: Statement truth value as Bool.

        EXAMPLES:

        The 4-Krönecker quiver::
        
            sage: load("quiver.py")
            sage: K = Quiver( matrix(  [[0, 4],
            ....:                       [0, 0]]))
            sage: K.is_connected()
            True

        The doubled 1-Krönecker quiver::
        
            sage: load("quiver.py")
            sage: C1 = Quiver(matrix(  [[0,1],
            ....:                       [1,0]]))
            sage: C1.is_connected()
            True 

        The 3-loop point quiver::

            sage: load("quiver.py")
            sage: L = Quiver(matrix([[3]]))
            sage: L.is_connected()
            True

        The A_10 quiver::

            sage: load("quiver.py")
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

        sage: load("quiver.py")
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
    Basic representation-theoretical properties of the quiver
    """

    def Euler_matrix(self):
        r"""Returns the Euler matrix of the quiver.
        
        OUTPUT: Sage matrix.
        """
        return matrix.identity(self.number_of_vertices()) - self.adjacency_matrix()

    def Euler_form(self, x, y):
        r"""The Euler bilinear form of the quiver.
        
        INPUT: 
        - ``x`` -- vector of integers
        - ``y`` -- vector of integers

        OUTPUT: the multiplication of ``x * self.adjacency_matrix() * y`` as an  Int.
        
        """
        assert (x.length() == self.number_of_vertices() and y.length() == self.number_of_vertices())
        return x * self.Euler_matrix() * y

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

    def canonical_stability_parameter(self,d):
        """The canonical stability parameter is given by <d,_> - <_,d>"""
        E = self.Euler_matrix()
        return d * (-self.Euler_matrix().transpose() + E)

    def has_semistable_representation(self, d, theta, algorithm="reineke"):
        """Checks if there is a theta-semistable representation of dimension vector d."""
        assert algorithm == "reineke"

        # TODO implement this
        if algorithm == "reineke":
            raise NotImplementedError()

    def has_stable_representation(self, d, theta, algorithm="king"):
        """Checks if there is a theta-stable representation of dimension vector d."""
        # TODO implement this
        # https://mathscinet.ams.org/mathscinet-getitem?mr=1315461
        if algorithm == "king":
            raise NotImplementedError()
        # TODO implement this
        # al stands for Adriaenssens--Le Bruyn
        # https://mathscinet.ams.org/mathscinet-getitem?mr=1972892
        if algorithm == "al":
            raise NotImplementedError()

    # TODO dimension vectors should have .is_stable(), .is_amply_stable()?
    def is_amply_stable(self, d, theta):
        """Checks if d is amply stable of theta"""
        # Section 4 of https://arxiv.org/pdf/1410.0466.pdf
        raise NotImplementedError()


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
        # https://arxiv.org/pdf/math/9911014.pdf
        # in Derksen--Weyman's https://mathscinet.ams.org/mathscinet-getitem?mr=1930979 it is claimed that there is a second Schofield algorithm
        # (they do cite the wrong Schofield preprint though...)
        elif algorithm == "schofield-2":
            raise NotImplementedError()

    def harder_narasimhan_stratification(self, d, theta, denominator=sum):
        # TODO what to return?
        # list of the Harder-Narasimhan types?
        # denominator default being sum is total dimension, there are variations possible
        # and the strata will be different!
        # https://mathscinet.ams.org/mathscinet-getitem?mr=1974891
        raise NotImplementedError()

    def in_fundamental_domain(self, d):
        # see e.g. page 3 of https://arxiv.org/pdf/2303.08522.pdf
        raise NotImplementedError()

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



# TODO convention for generator functions is capitalise them?

def GeneralizedKroneckerQuiver(m):
    """
    The generalized Kronecker quiver has two vertices and $m$ arrows from the
    first to the second.

    TESTS::

        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: Q.number_of_vertices()
        2
        sage: Q.number_of_arrows()
        3

    """
    Q = Quiver(matrix([[0, m], [0, 0]]), name = str(m)+"-Kronecker quiver")
    # TODO do Q.rename here
    return Q


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
    # TODO implement this
    return None


def DynkinQuiver( T):
    # TODO implement this
    # TODO orientation: have a default (so for A and D: linear, type E?) but make it possible to change the orientation
    # use https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/root_system/dynkin_diagram.html
    raise NotImplementedError()


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
