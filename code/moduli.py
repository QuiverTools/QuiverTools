from quiver import *
from sage.combinat.permutation import *


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


from abc import ABC, abstractmethod

class QuiverModuli(ABC):
    @abstractmethod
    def __init__(self, Q, d, theta, condition):
        assert (Q.number_of_vertices() == d.length() & Q.number_of_vertices() == theta.length())
        assert (condition in ["semistable", "stable"])
        self._Q = Q
        self._d = d
        self._theta = theta
        self._condition = condition

    def quiver(self):
        return self._Q

    def dimension_vector(self):
        return self._d

    def stability_parameter(self):
        return self._theta

    def is_nonempty(self):
        if self._condition == "stable":
            return self._Q.has_stable_representation(self._d, self._theta)
        elif self._condition == "semistable":
            return self._Q.has_semistable_representation(self._d, self._theta)

    @abstractmethod
    def dimension(self):
        pass

    @abstractmethod
    def is_smooth():
        pass

class QuiverModuliSpace(QuiverModuli):

    def __init__(self, Q, d, theta, condition="stable"):
        QuiverModuli.__init__(self, Q, d, theta, condition)
        self._condition = condition # TODO better name than 'condition' or 'version'?

    def __repr__(self):
        return "A "+self._condition+" quiver moduli space with:\n"+ "Q = "+str(self._Q)+"\n"+ "d = "+str(self._d)+"\n"+ "theta = "+str(self._theta)

    def dimension(self):
        """Computes the dimension of the moduli space M^{theta-(s)st}(Q,d)."""

        """This involves several cases:
        * if there are theta-stable representations then dim M^{theta-sst}(Q,d) = M^{theta-st}(Q,d) = 1 - <d,d>
        * if there are no theta-stable representations then dim M^{theta-st}(Q,d) = -Infinity (by convention) and dim M^{theta-sst} = max_tau dim S_tau, the maximum of the dimension of all Luna strata."""

        """
        EXAMPLES

        The A2-quiver:
        sage: load("quiver.py")
        sage: load("moduli.py")
        sage: Q = GeneralizedKroneckerQuiver(1)
        sage: theta = vector([1,-1])
        sage: d = vector([1,1])
        sage: X = QuiverModuliSpace(Q,d,theta,condition="stable")
        sage: X.dimension()
        0
        sage: X = QuiverModuliSpace(Q,d,theta,condition="semistable")
        sage: X.dimension()
        0
        sage: d = vector([2,2])
        sage: X = QuiverModuliSpace(Q,d,theta,condition="stable")
        sage: X.dimension()
        -Infinity
        sage: X = QuiverModuliSpace(Q,d,theta,condition="semistable")
        sage: X.dimension()
        0

        The Kronecker quiver:
        sage: load("quiver.py")
        sage: load("moduli.py")
        sage: Q = GeneralizedKroneckerQuiver(2)
        sage: theta = vector([1,-1])
        sage: X = QuiverModuliSpace(Q,d,theta,condition="stable")
        sage: X.dimension()
        1
        sage: X = QuiverModuliSpace(Q,d,theta,condition="semistable")
        sage: X.dimension()
        1
        sage: d = vector([2,2])
        sage: X = QuiverModuliSpace(Q,d,theta,condition="stable")
        sage: X.dimension()
        -Infinity
        sage: X = QuiverModuliSpace(Q,d,theta,condition="semistable")
        sage: X.dimension()
        2

        The 3-Kronecker quiver:
        sage: load("quiver.py")
        sage: load("moduli.py")
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: d = vector([2,3])
        sage: X = QuiverModuliSpace(Q,d,theta,condition="semistable")
        sage: X.dimension()
        6
        sage: d = vector([3,3])
        sage: X = QuiverModuliSpace(Q,d,theta,condition="semistable")
        sage: X.dimension()
        10
        sage: d = vector([1,3])
        sage: X = QuiverModuliSpace(Q,d,theta,condition="stable")
        sage: X.dimension()
        0
        sage: d = vector([1,4])
        sage: X = QuiverModuliSpace(Q,d,theta,condition="stable")
        sage: X.dimension()
        -Infinity
        sage: X = QuiverModuliSpace(Q,d,theta,condition="semistable")
        sage: X.dimension()
        -Infinity
        """

        if self._Q.has_stable_representation(self._d, self._theta):
            # if there are stable representations then both the stable and
            # the semi-stable moduli space have dimension `1-<d,d>`
            return 1 - self._Q.euler_form(self._d, self._d)
        else:
            # Stable locus is empty
            if self._condition == "semistable":
                if self._Q.has_semistable_representation(self._d, self._theta):
                    # In this case the dimension is given by the maximum of the dimensions of the Luna strata
                    allLunaTypes = self._Q.all_luna_types(self._d,self._theta)
                    return max([self.dimension_of_luna_stratum(tau) for tau in allLunaTypes])
                else:
                    # I somehow like the convention that the dimension of the empty set is -Infinity
                    return -oo
            else:
                # self._condition == "stable"
                return -oo

    def codimension_of_harder_narasimhan_stratum_in_representation_space(self,dstar,denominator=sum):
        """Computes the codimension of the HN stratum R_{d^*}^HN inside R_d."""
        """The codimension of the HN stratum of d^* = (d^1,...,d^s) is given by - sum_{k < l} <d^k,d^l>"""
        # This check takes a long time. Shall we do it nonetheless?
        assert self._Q.is_harder_narasimhan_type(dstar,self._theta,denominator=denominator)
        s = len(dstar)
        return -sum([self._Q.euler_form(dstar[k],dstar[l]) for k in range(s-1) for l in range(k+1,s)])

    def dimension_of_luna_stratum(self,tau):
        """Computes the dimension of the Luna stratum S_tau."""
        """The dimension of the Luna stratum of tau = [(d^1,p^1),...,(d^s,p^s)] is sum_k l(p^k)(1 - <d^k,d^k>) where for a partition p = (n_1,...,n_l), the length l(p) is l, i.e. the number of rows."""

        """
        EXAMPLES
        The Kronecker quiver:
        sage: load("quiver.py")
        sage: load("moduli.py")
        sage: Q = GeneralizedKroneckerQuiver(2)
        sage: d = vector([2,2])
        sage: theta = vector([1,-1])
        sage: L = Q.all_luna_types(d,theta)
        sage: L
        [[((1, 1), [2])], [((1, 1), [1, 1])]]
        sage: X = QuiverModuliSpace(Q,d,theta,condition="semistable")
        sage: [X.dimension_of_luna_stratum(tau) for tau in L]
        [1, 2]
        """
        # This check takes a long time. Shall we do it nonetheless?
        assert self._Q.is_luna_type(tau,self._theta)
        return sum([len(dn[1])*(1-self._Q.euler_form(dn[0],dn[0])) for dn in tau])

    def betti_numbers(self):
        raise NotImplementedError()

    def is_smooth(self):
        # if theta-coprime then you can shortcut everything
        # if theta != 0 reduce to theta = 0 using https://mathscinet.ams.org/mathscinet-getitem?mr=1972892 (Adriaenssens--Le Bruyn)
        # if theta = 0, then use https://mathscinet.ams.org/mathscinet-getitem?mr=1929191 (Bocklandt)
        raise NotImplementedError()


class QuiverModuliStack(QuiverModuli):

    def __init__(self, Q, d, theta, condition="stable"):
        QuiverModuli.__init__(self, Q, d, theta)
        self._condition = condition # TODO better name than 'condition' or 'version'?

    def __repr__(self):
        return "A "+self._condition+" quiver moduli stack with:\n"+ "Q = "+str(self._Q)+"\n"+ "d = "+str(self._d)+"\n"+ "theta = "+str(self._theta)

    def dimension(self):
        """dim [R^{(s)st}/G] = dim R^{(s)st} - dim G
        this is -<d,d> if the (semi-)stable locus is non-empty"""
        if (self._condition == "stable" and self._Q.has_stable_representation(self._d, self._theta)) or (self._condition == "semistable" and self._Q.has_semistable_representation(self._d, self._theta)):
            return -self._Q.euler_form(self._d,self._d)
        else:
            return -oo

    def is_smooth(self):
        # TODO think about the empty case, should it be smooth?
        return True


class SmoothModel:
    """How about this: instead of a separate class SmoothModel, we could define a method framed_moduli_space(self,n) inside the class QuiverModuliSpace which returns another quiver moduli space. After all, it is itself a quiver moduli space."""
    def __init__(self):
        pass

    def betti_numbers(self):
        raise NotImplementedError()


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
            nonEmpty = ((self._version == 'sst') and Q.has_semi_stable_representations(d)) or  ((self._version == 'st') and Q.has_stable_representations(d))
            if nonEmpty:
                return -Q.Euler_form(d,d)
            else:
                return float('NaN')
        else: # we're dealing with the moduli space
            if Q.has_stable_representations(d):
                """If there are stable representations, then both the stable and the semi-stable moduli space have dimension 1-<d,d>"""
                return 1 - Q.Euler_form(d,d)
            else: # no stables
                if (self.version == 'st'):
                    """The stable moduli space is empty so we return NaN"""
                    return float('NaN')
                else:
                    if Q.has_semi_stable_representations(d):
                        """This is the case which I don't know how to deal with: there are semi-stables but no stables. Then I think the dimension can be determined using the etale local structure. But I don't quite know how (nor how to implement it)."""
                        return float('NaN')
                    else: # no semi-stables
                        return float('NaN')

    # def slope(self, e):
    #     """The slope of e is defined as theta*e/(sum_i e_i). We need to ensure that e is non-negative and at least one entry is positive."""
    #     assert (e.length() == self._dimensionVector.length())
    #     assert all([(ei >= 0) for ei in e])
    #     assert any([(ei > 0) for ei in e])
    #     theta = self.stability_parameter()
    #     return (theta*e)/(sum(list(e)))

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
