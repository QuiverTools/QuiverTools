# Don't I first have to import quiver_class.sage? That doesn't work for me somehow. The only way I can make it work is load both quiver_class.sage and this file in sage and go from there. There must be a way to import this but I'm too stupid.

"""Defines how permutations are multiplied."""
Permutations.options(mult='r2l')

class Quiver_moduli:

    """
    A quiver moduli space is implemented by storing a quiver Q, a dimension vector d and a stability parameter theta. We want to distinguish between moduli spaces and moduli stacks (stacky = False/True) and if we're considering the stable or the semi-stable version.
    We do not require theta*d = 0 but work with slope stability instead.

    Variables:
    quiver
    dimensionVector
    stabilityParameter
    stacky = False
    version = 'sst'
    """

    def __init__(self, quiver, dimensionVector, stabilityParameter, stacky=False, version='sst'):
        assert ( (quiver.number_of_vertices() == dimensionVector.length()) & (quiver.number_of_vertices() == stabilityParameter.length()) )
        self._quiver = quiver
        self._dimensionVector = dimensionVector
        self._stabilityParameter = stabilityParameter
        self._stacky = stacky
        assert (version in ['st','sst'])
        self._version = version

    # I would like something along these lines, but I can't seem to import quiver_class.sage
    # @classmethod
    # def canonical_stability(cls, quiver, dimensionVector, stacky=False, version='sst'):
    #     assert (quiver.number_of_vertices() == dimensionVector.length())
    #     theta = quiver_class.canonical_stability_parameter()
    #     return cls(quiver,dimensionVector,theta,stacky=stacky,version=version)

    def __repr__(self):
        output = ""
        if (self._version == 'sst'):
            output += "Semi-stable "
        else:
            output += "Stable "
        output += "quiver moduli "
        if self._stacky:
            output += "stack "
        else:
            output += "space "
        output += "with:\n"+"Q = "+str(self._quiver)+"\n"+"d = "+str(self._dimensionVector)+"\n"+"theta = "+str(self._stabilityParameter)
        return output

    def quiver(self):
        return self._quiver

    def dimension_vector(self):
        return self._dimensionVector

    def stability_parameter(self):
        return self._stabilityParameter

    def slope(self, e):
        """The slope of e is defined as theta*e/(sum_i e_i). We need to ensure that e is non-negative and at least one entry is positive."""
        assert (e.length() == self._dimensionVector.length())
        assert all([(ei >= 0) for ei in e])
        assert any([(ei > 0) for ei in e])
        theta = self.stability_parameter()
        return (theta*e)/(sum(list(e)))

    """
    Notation for explanations:
    G = G_d = prod_{i in Q_0} GL_{d_i}
    T = maximal torus of diagonal matrices
    PG = G/G_m
    PT = T/G_m maximal torus of PT
    W = Weyl group of T in G = Weyl group of PT in PG
      = prod_{i in Q_0} S_{d_i}
    R = bigoplus_{a in Q_1} Hom(k^{d_{s(a)}},k^{d_{t(a)}})
    R^{sst}, R^{st} semi-stable/stable locus
    """

    def Weyl_group(self):
        """Returns the Weyl group W both as a group and as a list of its factors."""
        d = self.dimension_vector()
        asList = [Permutations(di).list() for di in d]
        return {"Group" : cartesian_product(asList),
        "List" : asList}

    def longest_Weyl_group_element(self):
        """Returns the longest Weyl group element as a list of permutations."""
        d = self.dimension_vector()
        return [Permutation(range(di,0,-1)) for di in d]

    """Static methods"""

    """The following methods are general. They do not depend on the specific situation."""

    @staticmethod
    def left_permutation_action_on_list(permut,li):
        """Assumes that the permutation and the list have the same length
        Imitates the following. If we have a permutation p and a set {t_1,...,t_n}
        then we want {t_{p^{-1}(1)},...,t_{p^{-1}(n)}}."""
        n = len(li)
        assert (n == permut.size())
        return list(map(lambda i: li[permut.inverse()[i]-1], range(n)))

    @staticmethod
    def left_permutation_action_on_polynomial(permut,f,alphabet):
        """Computes f(t_{p(1)},...,t_{p(n)})"""
        d = dict(zip(alphabet,left_permutation_action_on_list(permut,alphabet)))
        return f.subs(d)

    @staticmethod
    def divided_difference(i,f,alphabet):
        """Computes (f-(fs_i))/(t_{i+1}-t_i)"""
        n = len(alphabet)
        assert (i < n)
        reflection = Permutations(n).simple_reflection(i)
        return (f-left_permutation_action_on_polynomial(reflection,f,alphabet))/(alphabet[i-1]-alphabet[i])
