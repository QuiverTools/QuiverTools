class Quiver:

    """
    A quiver is represented by its adjacency matrix (a_ij) in M_{n x n}(N) where Q_0 = {1,...,n} and a_{ij} is the number of arrows i --> j.

    Variables:
    adjacencyMatrix
    name = None
    """

    def __init__(self, adjacencyMatrix, name=None):
        assert (adjacencyMatrix.nrows() == adjacencyMatrix.ncols())
        assert all([all([(a >= 0) for a in rows]) for rows in [list(rows) for rows in list(adjacencyMatrix)]])
        # Should we raise an exception/error instead?
        self._adjacencyMatrix = adjacencyMatrix
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
        output += "adjacency matrix:\n"+str(self._adjacencyMatrix)
        return output

    def adjacency_matrix(self):
        return self._adjacencyMatrix

    def number_of_vertices(self):
        return self.adjacency_matrix().nrows()

    def thin_dimension_vector(self):
        return vector([1 for i in range(self.number_of_vertices())])

    def number_of_arrows(self):
        thin = self.thin_dimension_vector()
        return thin*self.adjacency_matrix()*thin

    def Euler_matrix(self):
        return (matrix.identity(self.number_of_vertices()) - self.adjacency_matrix())

    def Euler_form(self, x, y):
        assert ((x.length() == self.number_of_vertices()) & (y.length() == self.number_of_vertices()))
        return x*self.Euler_matrix()*y

    def canonical_stability_parameter(self,d):
        """The canonical stability parameter is given by <d,_> - <_,d>"""
        E = self.Euler_matrix()
        return d*(-E.transpose()+E)

    def is_acyclic(self):
        """A quiver is acyclic iff its adjacency matrix is nilpotent."""
        A = self.adjacency_matrix()
        n = self.number_of_vertices()
        return (A^n == zero_matrix(ZZ,n))

    def opposite_quiver(self):
        """The opposite quiver is given by the transpose of the adjacency matrix of the original quiver."""
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

    def allows_semi_stable_representations(self,d):
        """Checks if there is a semi-stable representation of dimension vector d."""
        """Still needs to be implemented!"""
        return True

    def allows_stable_representations(self,d):
        """Checks if there is a semi-stable representation of dimension vector d."""
        """Still needs to be implemented!"""
        return True


def generalized_Kronecker(m):
    """The generalized Kronecker quiver has two vertices 1,2 and m arrows 1 --> 2.

    TESTS::

        sage: Q = generalized_Kronecker(3)
        sage: Q.number_of_vertices()
        2
        sage: Q.number_of_arrows()
        3

    """
    Q = Quiver(matrix([[0, m], [0, 0]]), name = str(m)+"-Kronecker quiver")
    # TODO do Q.rename here
    return Q

def three_vertex_quiver(m12, m13, m23):
    """An acyclic quiver with 3 vertices and mij many arrows i --> j for 1 <= i < j <= 3."""
    Q = Quiver(matrix([[0,m12,m13],[0,0,m23],[0,0,0]]), name = "An acyclic 3-vertex quiver")
    # TODO do Q.rename here
    return Q

def loop_quiver(cls, m):
    """A quiver with one vertex and m arrows."""
    Q = Quiver(matrix([[m]]), name = str(m)+"-loop quiver")
    # TODO do Q.rename here
    return Q

def Jordan():
    Q = loop_quiver(1)
    # TODO do Q.rename here
    return Q

def subspace_quiver(cls, m):
    """A quiver with m sources 1,...,m and one sink m+1; one arrow from every source to the sink."""
    A = zero_matrix(ZZ, m + 1)
    for i in range(m):
        A[i, m] = 1

    Q = Quiver(A, name = str(m)+"-subspace quiver")
    # TODO do Q.rename here
    return Q
