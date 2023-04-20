from quiver import *


"""Defines how permutations are multiplied."""
Permutations.options(mult='r2l')


"""
a brainstorm

What about an abstract base class QuiverModuli which deals with all the things common to moduli spaces and moduli stacks?
There'd be abstract methods for things like dimension.

Then we'd have implementations in
- QuiverModuliSpace
- QuiverModuliStack

This avoids the ugly distinction between the two.
The stack would also allow _not_ specifying stable or semistable,
whereas for the quiver moduli space this is a necessity.


Something like enumerating Harder-Narasimhan strata is then a method of QuiverModuliStack?

Something like computing Betti numbers is then only implemented for QuiverModuliSpace.

"""

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
    @classmethod
    def canonical_stability(cls, quiver, dimensionVector, stacky=False, version='sst'):
        assert (quiver.number_of_vertices() == dimensionVector.length())
        theta = quiver.canonical_stability_parameter(d)
        return cls(quiver,dimensionVector,theta,stacky=stacky,version=version)

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

    def dimension(self):
        Q = self.quiver()
        d = self.dimension_vector()
        if self._stacky: # we're looking at the moduli stack
            """dim [R^{(s)st}/G] = dim R^{(s)st} - dim G
            this is -<d,d> if the (semi-)stable locus is non-empty"""
            nonEmpty = ((self._version == 'sst') and Q.allows_semi_stable_representations(d)) or  ((self._version == 'st') and Q.allows_stable_representations(d))
            if nonEmpty:
                return -Q.Euler_form(d,d)
            else:
                return float('NaN')
        else: # we're dealing with the moduli space
            if Q.allows_stable_representations(d):
                """If there are stable representations, then both the stable and the semi-stable moduli space have dimension 1-<d,d>"""
                return 1 - Q.Euler_form(d,d)
            else: # no stables
                if (self.version == 'st'):
                    """The stable moduli space is empty so we return NaN"""
                    return float('NaN')
                else:
                    if Q.allows_semi_stable_representations(d):
                        """This is the case which I don't know how to deal with: there are semi-stables but no stables. Then I think the dimension can be determined using the etale local structure. But I don't quite know how (nor how to implement it)."""
                        return float('NaN')
                    else: # no semi-stables
                        return float('NaN')

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
