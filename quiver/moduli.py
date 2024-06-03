from quiver import *
from sage.combinat.permutation import *
import copy


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
    def __init__(self, Q, d, theta, denominator=sum, condition="semistable"):
        n = Q.number_of_vertices()
        assert (Q.number_of_vertices() == d.length() & Q.number_of_vertices() == theta.length())
        assert (condition in ["semistable", "stable"])
        assert all([denominator(Q.simple_root(i)) > 0 for i in range(1,n+1)])
        self._Q = Q
        self._d = d
        self._theta = theta
        self._denominator = denominator
        self._condition = condition

    def quiver(self):
        return self._Q

    def dimension_vector(self):
        return self._d

    def stability_parameter(self):
        return self._theta
    
    def denominator(self):
        return self._denominator

    def is_nonempty(self):
        if self._condition == "stable":
            return self._Q.has_stable_representation(self._d, self._theta)
        elif self._condition == "semistable":
            return self._Q.has_semistable_representation(self._d, self._theta)
        
    """
    HN business
    """
        
    def all_harder_narasimhan_types(self):
        r"""Returns the list of all HN types.

        OUTPUT: list of list of vectors of Ints
        """

        """A Harder--Narasimhan (HN) type of d with respect to theta is a sequence d^* = (d^1,...,d^s) of dimension vectors such that
        * d^1 + ... + d^s = d
        * mu_theta(d^1) > ... > mu_theta(d^s)
        * Every d^k is theta-semi-stable."""

        """
        EXAMPLES:

        The 3-Kronecker quiver:
        sage: from quiver import *
        sage: Q, d, theta = GeneralizedKroneckerQuiver(3), vector([2,3]), vector([3,-2])
        sage: X = QuiverModuliSpace(Q, d, theta)
        sage: X.all_harder_narasimhan_types()
        [[(1, 0), (1, 1), (0, 2)],
        [(1, 0), (1, 2), (0, 1)],
        [(1, 0), (1, 3)],
        [(1, 1), (1, 2)],
        [(2, 0), (0, 3)],
        [(2, 1), (0, 2)],
        [(2, 2), (0, 1)],
        [(2, 3)]]
        sage: Y = QuiverModuliSpace(Q, d, -theta)
        sage: Y.all_harder_narasimhan_types()
        [[(0, 3), (2, 0)]]

        A 3-vertex quiver:
        sage: from quiver import *
        sage: Q, d = ThreeVertexQuiver(2,3,4), vector([2,3,2])
        sage: theta = Q.canonical_stability_parameter(d)
        sage: Z = QuiverModuliSpace(Q, d, theta)
        sage: Z.all_harder_narasimhan_types()
        [[(0, 1, 0), (2, 0, 1), (0, 2, 1)],
        [(0, 1, 0), (1, 2, 1), (1, 0, 1)],
        [(0, 1, 0), (2, 1, 1), (0, 1, 1)],
        [(0, 1, 0), (2, 2, 1), (0, 0, 1)],
        [(0, 1, 0), (2, 2, 2)],
        [(1, 0, 0), (0, 1, 0), (1, 0, 1), (0, 2, 1)],
        [(1, 0, 0), (0, 1, 0), (1, 1, 1), (0, 1, 1)],
        [(1, 0, 0), (0, 1, 0), (1, 2, 1), (0, 0, 1)],
        [(1, 0, 0), (0, 1, 0), (1, 2, 2)],
        [(1, 0, 0), (0, 2, 0), (1, 0, 1), (0, 1, 1)],
        [(1, 0, 0), (0, 2, 0), (1, 1, 1), (0, 0, 1)],
        [(1, 0, 0), (0, 2, 0), (1, 1, 2)],
        [(1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)],
        [(1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 1, 2)],
        [(1, 0, 0), (1, 1, 0), (0, 2, 0), (0, 0, 2)],
        [(1, 0, 0), (1, 1, 0), (0, 2, 1), (0, 0, 1)],
        [(1, 0, 0), (1, 1, 0), (0, 2, 2)],
        [(1, 0, 0), (0, 3, 0), (1, 0, 1), (0, 0, 1)],
        [(1, 0, 0), (0, 3, 0), (1, 0, 2)],
        [(1, 0, 0), (1, 1, 1), (0, 2, 1)],
        [(1, 0, 0), (1, 2, 0), (0, 1, 0), (0, 0, 2)],
        [(1, 0, 0), (1, 2, 0), (0, 1, 1), (0, 0, 1)],
        [(1, 0, 0), (1, 2, 0), (0, 1, 2)],
        [(1, 0, 0), (0, 3, 1), (1, 0, 1)],
        [(1, 0, 0), (1, 2, 1), (0, 1, 1)],
        [(1, 0, 0), (1, 3, 1), (0, 0, 1)],
        [(1, 0, 0), (1, 3, 2)],
        [(0, 2, 0), (1, 1, 1), (1, 0, 1)],
        [(0, 2, 0), (2, 0, 1), (0, 1, 1)],
        [(0, 2, 0), (2, 1, 1), (0, 0, 1)],
        [(0, 2, 0), (2, 1, 2)],
        [(1, 1, 0), (0, 1, 0), (1, 0, 1), (0, 1, 1)],
        [(1, 1, 0), (0, 1, 0), (1, 1, 1), (0, 0, 1)],
        [(1, 1, 0), (0, 1, 0), (1, 1, 2)],
        [(1, 1, 0), (0, 2, 0), (1, 0, 1), (0, 0, 1)],
        [(1, 1, 0), (0, 2, 0), (1, 0, 2)],
        [(1, 1, 0), (1, 0, 1), (0, 2, 1)],
        [(1, 1, 0), (1, 1, 1), (0, 1, 1)],
        [(1, 1, 0), (1, 2, 0), (0, 0, 2)],
        [(1, 1, 0), (1, 2, 1), (0, 0, 1)],
        [(1, 1, 0), (1, 2, 2)],
        [(2, 0, 0), (0, 1, 0), (0, 2, 1), (0, 0, 1)],
        [(2, 0, 0), (0, 1, 0), (0, 2, 2)],
        [(2, 0, 0), (0, 2, 0), (0, 1, 1), (0, 0, 1)],
        [(2, 0, 0), (0, 2, 0), (0, 1, 2)],
        [(2, 0, 0), (0, 2, 1), (0, 1, 1)],
        [(2, 0, 0), (0, 3, 0), (0, 0, 2)],
        [(2, 0, 0), (0, 3, 1), (0, 0, 1)],
        [(2, 0, 0), (0, 3, 2)],
        [(0, 3, 0), (2, 0, 1), (0, 0, 1)],
        [(0, 3, 0), (2, 0, 2)],
        [(1, 2, 0), (0, 1, 0), (1, 0, 1), (0, 0, 1)],
        [(1, 2, 0), (0, 1, 0), (1, 0, 2)],
        [(1, 2, 0), (1, 0, 1), (0, 1, 1)],
        [(1, 2, 0), (1, 1, 1), (0, 0, 1)],
        [(1, 2, 0), (1, 1, 2)],
        [(2, 0, 1), (0, 3, 1)],
        [(2, 1, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)],
        [(2, 1, 0), (0, 1, 0), (0, 1, 2)],
        [(2, 1, 0), (0, 2, 0), (0, 0, 2)],
        [(2, 1, 0), (0, 2, 1), (0, 0, 1)],
        [(2, 1, 0), (0, 2, 2)],
        [(1, 2, 1), (1, 1, 1)],
        [(2, 1, 1), (0, 2, 1)],
        [(2, 2, 0), (0, 1, 0), (0, 0, 2)],
        [(2, 2, 0), (0, 1, 1), (0, 0, 1)],
        [(2, 2, 0), (0, 1, 2)],
        [(1, 3, 1), (1, 0, 1)],
        [(2, 2, 1), (0, 1, 1)],
        [(2, 3, 0), (0, 0, 2)],
        [(2, 3, 1), (0, 0, 1)],
        [(2, 3, 2)]]

        """
        Q, d, theta, denominator = self._Q, self._d, self._theta, self._denominator
        
        subdimensions = all_subdimension_vectors(d)
        subdimensions.sort(key=(lambda e: deglex_key(e, b=max(d)+1)))
        N = len(subdimensions)

        # sstIndexes is the list of indexes of all non-zero semistable subdimension vectors in subdimensions
        sstIndexes, sstSubdims = Q._Quiver__all_semistable_subdimension_vectors_helper(d, theta)

        # idx_diff(j, i) is the index of the difference subdimensions[j]-subdimensions[i] in the list subdimensions
        idx_diff = (lambda j, i: subdimensions.index(subdimensions[j]-subdimensions[i]))

        hn = [[[]] for j in range(N)]

        for j in range(1,N):
            # sstSub is the list of all indexes in subdimensions of semistable non-zero subdimension vectors of subdimensions[j]
            sstSub = list(filter(lambda i: is_subdimension_vector(subdimensions[i], subdimensions[j]), sstIndexes))
            # The HN types which are not of the form (d) are given by (e,f^1,...,f^s) where e is a proper subdimension vector such that mu_theta(e) > mu_theta(d) and (f^1,...,f^s) is a HN type of f = d-e such that mu_theta(e) > mu_theta(f^1) holds.
            hn[j] = [[i]+fstar for i in sstSub for fstar in list(filter(lambda fstar: fstar == [] or slope(subdimensions[i], theta, denominator=denominator) > slope(subdimensions[fstar[0]], theta, denominator=denominator), hn[idx_diff(j, i)]))]

        hn[0] = [[0]]

        return [[subdimensions[r] for r in fstar] for fstar in hn[N-1]]

    def is_harder_narasimhan_type(self, dstar):
        r"""Checks if dstar is a HN type.
        
        INPUT:
        - ``dstar``: list of vectors of Ints

        OUTPUT: statement truth value as Bool
        """

        """
        EXAMPLES:
        sage: from quiver import *
        sage: Q, d, theta = GeneralizedKroneckerQuiver(3), vector([2,3]), vector([1,0])
        sage: X = QuiverModuliSpace(Q, d, theta)
        sage: hn = X.all_harder_narasimhan_types()
        sage: all([X.is_harder_narasimhan_type(dstar) for dstar in hn])
        True
        sage: dstar = [vector([1,0]), vector([1,0]), vector([0,3])]
        sage: X.is_harder_narasimhan_type(dstar)
        False

        """
        Q, d, theta, denominator = self._Q, self._d, self._theta, self._denominator
        
        assert d == sum(dstar)

        if (d == Q.zero_vector()):
            return (dstar == [Q.zero_vector()])
        else:
            sstIndexes, sstSubdims = Q._Quiver__all_semistable_subdimension_vectors_helper(d, theta)
            slopeDecreasing = all([(slope(dstar[i],theta,denominator=denominator) > slope(dstar[i+1],theta,denominator=denominator)) for i in range(len(dstar)-1)])
            semistable = all([e in sstSubdims for e in dstar])
            return (slopeDecreasing and semistable)
        
    def codimension_of_harder_narasimhan_stratum(self, dstar, secure=True):
        """Computes the codimension of the HN stratum of dstar inside the representation variety.
        
        INPUT:
        - ``dstar``: list of vectors of Ints
        - ``secure``: Bool

        OUTPUT: codimension as Int
        """
        # It checks for dstar to be a HN type iff secure == True. This check is slow.
        # Be sure to be dealing with a HN type if you call it with secure == False. This is fast but yields nonsense, if dstar is not a HN type.

        """The codimension of the HN stratum of d^* = (d^1,...,d^s) is given by - sum_{k < l} <d^k,d^l>"""

        """
        EXAMPLES

        The 3-Kronecker quiver
        sage: from quiver import *
        sage: Q, d, theta = GeneralizedKroneckerQuiver(3), vector([2,3]), vector([1,0])
        sage: X = QuiverModuliSpace(Q, d, theta)
        sage: hn = X.all_harder_narasimhan_types(); hn
        [[(1, 0), (1, 1), (0, 2)],
        [(1, 0), (1, 2), (0, 1)],
        [(1, 0), (1, 3)],
        [(1, 1), (1, 2)],
        [(2, 0), (0, 3)],
        [(2, 1), (0, 2)],
        [(2, 2), (0, 1)],
        [(2, 3)]]
        sage: [X.codimension_of_harder_narasimhan_stratum(dstar) for dstar in hn]
        [12, 9, 8, 3, 18, 10, 4, 0]

        """
        Q = self._Q
        
        n = Q.number_of_vertices()
        assert all([e.length() == n for e in dstar])

        if secure:
            assert self.is_harder_narasimhan_type(dstar)

        s = len(dstar)
        return -sum([Q.euler_form(dstar[k],dstar[l]) for k in range(s-1) for l in range(k+1,s)])
    
    def codimension_unstable_locus(self):
        r"""Computes the codimension of the unstable locus inside the representation variety.

        OUTPUT: codimension as Int
        """

        """"
        EXAMPLES:
        sage: from quiver import *
        sage: Q, d, theta = GeneralizedKroneckerQuiver(3), vector([2,3]), vector([1,0])
        sage: X = QuiverModuliSpace(Q, d, theta)
        sage: X.codimension_unstable_locus()
        3
        sage: Q, d = ThreeVertexQuiver(1,6,1), vector([1,6,6])
        sage: theta = Q.canonical_stability_parameter(d)
        sage: X = QuiverModuliSpace(Q, d, theta)
        sage: X.codimension_unstable_locus()
        1
        sage: Q, d, theta = GeneralizedKroneckerQuiver(1), vector([2,3]), vector([1,0])
        sage: X = QuiverModuliSpace(Q, d, theta)
        sage: X.codimension_unstable_locus()
        0

        """
        d = self._d

        hn = list(filter(lambda dstar: dstar != [d], self.all_harder_narasimhan_types()))
        # Note that while the HN types and strata depend on the denominator, the list of all their codimensions does not.

        return min([self.codimension_of_harder_narasimhan_stratum(dstar, secure=False) for dstar in hn])
    
    """
    Luna
    """

    def all_luna_types(self):
        r"""Returns the unordered list of all Luna types of d for theta.

        OUTPUT: list of tuples containing Int-vector and Int 
        """

        """A Luna type of d for theta is an unordered sequence (i.e. multiset) ((d^1,m_1),...,(d^s,m_s)) of dimension vectors d^k and (positive) natural numbers m_k such that
        * m_1d^1 + ... + m_sd^s = d
        * mu_theta(d^k) = mu_theta(d)
        * All d^k admit a theta-stable representation
        """

        """Example: Suppose that d = 3e and e, 2e, d = 3e are the only stable subdimension vectors. Then the Luna types are:
        ((3e,1))
        ((2e,1),(e,1))
        ((e,3))
        ((e,2),(e,1))
        ((e,1),(e,1),(e,1))
        """

        """Therefore we implement it as follows. A Luna type for us is a list [(d^1,p^1),...,(d^s,p^s)] (actually it should be unordered, but that's difficult because vectors are mutable) of dimension vectors d^k and (non-empty) partitions p^k such that
        * |p^1|d^1 + ... + |p^s|d^s = d
        * same
        * same """

        """So in the above example, the Luna types are
        [(3e,[1])]
        [(2e,[1]),(e,[1])]
        [(e,[3])]
        [(e,[2,1])]
        [(e,[1,1,1])]
        """

        """
        EXAMPLES:

        The Kronecker quiver:
        sage: from quiver import *
        sage: Q, d, theta = KroneckerQuiver(), vector([3,3]), vector([1,-1])
        sage: X = QuiverModuliSpace(Q, d, theta)
        sage: X.all_luna_types()
        [[((1, 1), [3])], [((1, 1), [2, 1])], [((1, 1), [1, 1, 1])]]

        The 3-Kronecker quiver:
        sage: from quiver import *
        sage: Q, d, theta = GeneralizedKroneckerQuiver(3), vector([3,3]), vector([1,-1])
        sage: X = QuiverModuliSpace(Q, d, theta)
        sage: X.all_luna_types()
        [[((1, 1), [3])],
        [((1, 1), [2, 1])],
        [((1, 1), [1, 1, 1])],
        [((1, 1), [1]), ((2, 2), [1])],
        [((3, 3), [1])]]

        The 6-subspace quiver:
        sage: from quiver import *
        sage: Q, d, theta = GeneralizedKroneckerQuiver(3), vector([3,3]), vector([1,-1])
        sage: X = QuiverModuliSpace(Q, d, theta)
        sage: X.all_luna_types()
        [[((1, 1), [3])],
        [((1, 1), [2, 1])],
        [((1, 1), [1, 1, 1])],
        [((1, 1), [1]), ((2, 2), [1])],
        [((3, 3), [1])]]

        """
        Q, d, theta, denominator = self._Q, self._d, self._theta, self._denominator

        if (d == Q.zero_vector()):
            return [tuple([Q.zero_vector(),[1]])]
        else: 
            subdims = all_subdimension_vectors(d)
            subdims.sort(key=(lambda e: deglex_key(e, b=max(d)+1)))
            N = len(subdims)
            # slopeIndexes is the list of indexes j such that the slope of e := subdims[j] equals the slope of d (this requires e != 0)
            slopeIndexes = list(filter(lambda j: slope(subdims[j], theta, denominator=denominator) == slope(d, theta, denominator=denominator), range(1,N)))
            # We consider all subdimension vectors which are not zero, whose slope equals the slope of d, and which admit a stable representation
            # They're in deglex order by the way the helper function works.
            stIndexes, stSubdims = Q._Quiver__all_stable_subdimension_vectors_helper(d, theta, denominator=denominator)
            # idx_diff(j, i) is the index of the difference stSubdims[j]-stSubdims[i] in the list stSubdims
            idx_diff = (lambda j, i: subdims.index(subdims[j]-subdims[i]))

            # partialLunaTypes is going to hold all "partial Luna types" of e for every e in stSubdims; a partial luna type of e is an unordered sequence (i.e. multiset) {(e^1,n_1),...,(e^s,n_s)} such that all e^k are distinct, e^1+...+e^s = e and the slopes of all e^k are the same (and thus equal the slope of e).
            partialLunaTypes = [[] for j in range(N)]
            for j in range(N):
                stSub = list(filter(lambda i: is_subdimension_vector(subdims[i], subdims[j]) and i != j, stIndexes))
                for i in stSub:
                    smaller = partialLunaTypes[idx_diff(j,i)]
                    for tau in smaller:
                        # Check if f := stSubdims[i] occurs as a dimension vector in tau.
                        # If so, say of the form (f,n) then remove this occurrence and add (f,n+1)
                        # If not, then add (f,1)
                        tauNew = copy.deepcopy(tau)
                        occurs = False
                        for dn in tauNew:
                            if (dn[0] == i):
                                # We remove dn from tau and add the tuple (e,dn[1]+1) instead
                                tauNew.remove(dn)
                                tauNew.append(tuple([i,dn[1]+1]))
                                occurs = True
                        if (not occurs):
                            tauNew.append(tuple([i,1]))
                        # Now tauNew is a Luna type of e := subdims[j] the desired form
                        # We sort it, because it's supposed to be unordered
                        tauNew.sort()
                        # If tau isn't already contained, then we add it
                        if tauNew not in partialLunaTypes[j]:
                            partialLunaTypes[j] = partialLunaTypes[j] + [tauNew]
                if (j in stIndexes):
                    # If e = subdims[j] is stable then (e,1) is also a Luna type.
                    partialLunaTypes[j] = partialLunaTypes[j] + [[tuple([j,1])]]
        
            partial = partialLunaTypes[N-1]
            allLunaTypes = []
            for tau in partial:
                listOfPartitions = [Partitions(dn[1]).list() for dn in tau]
                Prod = cartesian_product(listOfPartitions).list()
                allLunaTypes = allLunaTypes + [[tuple([subdims[tau[i][0]],p[i]]) for i in range(len(tau))] for p in Prod]
            return allLunaTypes
        
    def is_luna_type(self, tau):
        r"""Checks if tau is a Luna type for theta.
        
        INPUT:
        - ``tau``: list of tuples

        OUTPUT: statement truth value as Bool
        """

        """
        EXAMPLES:
        
        sage: from quiver import *
        sage: Q, d, theta = KroneckerQuiver(), vector([3,3]), vector([1,-1])
        sage: X = QuiverModuliSpace(Q, d, theta)
        sage: l = X.all_luna_types()
        sage: all([X.is_luna_type(tau) for tau in l])
        True

        """
        Q, d, theta, denominator = self._Q, self._d, self._theta, self._denominator

        n = Q.number_of_vertices()
        assert all([dn[0].length() == n for dn in tau])
        assert d == sum([sum(dn[1])*dn[0] for dn in tau])

        if (d == Q.zero_vector()):
            return (tau == [tuple([Q.zero_vector(),[1]])])
        else:
            dstar = [dn[0] for dn in tau]
            stIndexes, stSubdims = Q._Quiver__all_stable_subdimension_vectors_helper(d, theta, denominator=denominator)
            return all([e in stSubdims for e in dstar]) # Note that in particular the zero vector must not lie in dstar
        
    def dimension_of_luna_stratum(self, tau, secure=True):
        r"""Computes the dimension of the Luna stratum S_tau.
        
        INPUT:
        - ``tau``: list of tuples
        - ``secure``: Bool

        OUTPUT: Dimension as Int
        """

        """The dimension of the Luna stratum of tau = [(d^1,p^1),...,(d^s,p^s)] is sum_k l(p^k)(1 - <d^k,d^k>) where for a partition p = (n_1,...,n_l), the length l(p) is l, i.e. the number of rows."""

        """
        EXAMPLES:

        The Kronecker quiver:
        sage: from quiver import *
        sage: Q, d, theta = KroneckerQuiver(), vector([2,2]), vector([1,-1])
        sage: X = QuiverModuliSpace(Q, d, theta)
        sage: L = X.all_luna_types(); L
        [[((1, 1), [2])], [((1, 1), [1, 1])]]
        sage: [X.dimension_of_luna_stratum(tau) for tau in L]
        [1, 2]

        """

        if secure:
            assert self.is_luna_type(tau)
        return sum([len(dn[1])*(1-self._Q.euler_form(dn[0],dn[0])) for dn in tau])
    
    def local_quiver_setting(self, tau, secure=True):
        r"""Returns the local quiver and dimension vector for the given Luna type.
        
        INPUT:
        - ``tau``: list of tuples
        - ``secure``: Bool

        OUTPUT: tuple cosisting of a Quiver object and a vector
        """

        """
        EXAMPLES:

        sage: from quiver import *
        sage: Q, d, theta = GeneralizedKroneckerQuiver(3), vector([2,2]), vector([1,-1])
        sage: X = QuiverModuliSpace(Q, d, theta)
        sage: L = X.all_luna_types(); L
        [[((1, 1), [2])], [((1, 1), [1, 1])], [((2, 2), [1])]]
        sage: Qloc, dloc = X.local_quiver_setting(L[0]); Qloc.adjacency_matrix() , dloc
        ([1], (2))
        sage: Qloc, dloc = X.local_quiver_setting(L[1]); Qloc.adjacency_matrix() , dloc
        (
        [1 1]        
        [1 1], (1, 1)
        )
        sage: Qloc, dloc = X.local_quiver_setting(L[2]); Qloc.adjacency_matrix() , dloc
        ([4], (1))

        """

        if secure:
            assert self.is_luna_type(tau)

        Q = self._Q

        A = matrix([[Q.generic_ext(dp[0], eq[0]) for eq in tau for n in eq[1]] for dp in tau for m in dp[1]])
        Qloc = Quiver(A)
        dloc = vector([m for dp in tau for m in dp[1]])

        return Qloc, dloc
    
    # TODO: The codimension computation requires the dimension of the nullcone. This is hard, it turns out. It can be done with the Hesselink stratification, but I wasn't willing to go thourgh Lieven's treatment of this.
    def __codimension_inverse_image_luna_stratum(self, tau):
        r"""Computes the codimension of pi^{-1}(S_tau) inside R(Q,d) where pi: R(Q,d)^{theta-sst} --> M^{theta-sst}(Q,d) is the semistable quotient map.
        
        INPUT:
        - ``tau``: list of tuples

        OUTPUT: codimension as Int
        """

        """For tau = [(d^1,p^1),...,(d^s,p^s)] the codimension of pi^{-1}(S_tau) is
        
        -<d,d> + sum_{k=1}^s (<d^k,d^k> - l(p^k) + ||p^k||^2) - dim N(Q_tau, d_tau)

        where for a partition p = (n_1,...,n_l), we define ||p||^2 = sum_v n_v^2 and N(Q_tau, d_tau) is the nullcone of the local quiver setting.
        """

        Q, d = self._Q, self._d
        Qtau, dtau = self.local_quiver_setting(tau, secure=False)
        dimNull = Qtau.dimension_nullcone(dtau)
        return -Q.euler_form(d,d)+sum([Q.euler_form(dk[0], dk[0])-len(dk[1])+sum([nkv**2 for nkv in dk[1]]) for dk in tau])-dimNull
    
    def codimension_properly_semistable_locus(self):
        r"""Computes the codimension of R(Q,d)^{theta-sst} \ R(Q,d)^{theta-st} inside R(Q,d).
        
        OUTPUT: codimension as Int
        """

        """The codimension of the properly semistable locus is the minimal codimension of the inverse image of the non-stable Luna strata."""

        Q, d = self._Q, self._d
        L = self.all_luna_types()
        # This is the stable Luna type; remove it if it occurs
        dstable = [tuple([d, [1]])]
        L = list(filter(lambda tau: tau != dstable, L))
        return min([self.__codimension_inverse_image_luna_stratum(tau) for tau in L])
        
    """
    (Semi-)stability
    """

    def semistable_equals_stable(self):
        # TODO: I moved it here because it uses the luna methods. 

        r"""Checks if every theta-semistable representation of dimension vector d is theta-stable

        OUTPUT: statement truth value as Bool
        """

        """
        EXAMPLES:

        The 3-Kronecker quiver:
        sage: from quiver import *
        sage: Q, d, theta = GeneralizedKroneckerQuiver(3), vector([3,3]), vector([1,-1])
        sage: X = QuiverModuliSpace(Q, d, theta)
        sage: X.semistable_equals_stable()
        False
        sage: e = vector([2,3])
        sage: Y = QuiverModuliSpace(Q, e, theta)
        sage: Y.semistable_equals_stable()
        True

        A double framed example as in our vector fields paper
        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3).framed_quiver(vector([1,0])).coframed_quiver(vector([0,0,1]))
        sage: d = vector([1,2,3,1])
        sage: theta = vector([1,300,-200,-1])
        sage: is_coprime_for_stability_parameter(d, theta)
        False
        sage: X = QuiverModuliSpace(Q, d, theta)
        sage: X.semistable_equals_stable()
        True

        """

        """Every theta-semistable representation is theta-stable if and only if there are no Luna types other than (possibly) (d,[1])."""
        Q, d, theta = self._Q, self._d, self._theta

        # As the computation of all Luna types takes so much time, we should first tests if d is theta-coprime
        if is_coprime_for_stability_parameter(d,theta):
            return True
        else:
            # This is probably the fastest way as checking theta-coprimality is fast whereas checking for existence of a semi-stable representation is a bit slower
            if not Q.has_semistable_representation(d,theta):
                return True
            else:
                allLunaTypes = self.all_luna_types()
                genericType = [tuple([d,[1]])]
                if genericType in allLunaTypes:
                    allLunaTypes.remove(genericType)
                return (not allLunaTypes) # This checks if the list is empty
            
    """
    Ample stability
    """
            
    # TODO dimension vectors should have .is_stable(), .is_amply_stable()?
    # Is this comment now obsolete?
    def is_amply_stable(self):
        r"""Checks if d is amply stable for theta.

        OUTPUT: statement truth value as Bool
        """
        
        """By definition, d is theta-amply stable if the codimension of the theta-stable locus inside R(Q,d) is at least 2."""

        # By Prop. 4.1 of https://arxiv.org/pdf/1410.0466.pdf d is amply stable for theta provided that <e,d-e> <= -2 for every proper subdimension vector.
        # But can we find a necessary and sufficient condition?
        # If every theta-semi-stable representation of dimension vector d is theta-stable then theta-ample stability is equivalent to every proper HN stratum having codimension at least 2.
        # I think I can compute the codimension of the non-stable locus in full generality.
        # TODO: It's more difficult than I thought. I think it's doable though.

        """
        EXAMPLES:

        sage: from quiver import *
        sage: Q, d, theta = GeneralizedKroneckerQuiver(3), vector([2,3]), vector([3,-2])
        sage: X = QuiverModuliSpace(Q, d, theta)
        sage: X.is_amply_stable()
        True
        sage: Y = QuiverModuliSpace(Q, d, -theta)
        sage: Y.is_amply_stable()
        False

        A three vertex example from the rigidity paper:
        sage: from quiver import *
        sage: Q, d = ThreeVertexQuiver(1,6,1), vector([1,6,6])
        sage: theta = Q.canonical_stability_parameter(d)
        sage: X = QuiverModuliSpace(Q, d, theta)
        sage: X.is_amply_stable()
        False

        """

        Q, d, theta = self._Q, self._d, self._theta

        # It's currently only possible with this distinction
        if is_coprime_for_stability_parameter(d, theta):
            return self.codimension_unstable_locus() >= 2
        else:
            return min([self.codimension_unstable_locus(), self.codimension_properly_semistable_locus()]) >= 2

    def is_strongly_amply_stable(self):
        r"""Checks if d is strongly amply stable for theta.

        OUTPUT: statement truth value as Bool
        """

        """We call d strongly amply stable for theta if <e,d-e> <= -2 holds for all subdimension vectors e of d which satisfy slope(e) >= slope(d)."""

        """
        EXAMPLES:

        sage: from quiver import *
        sage: Q, d, theta = GeneralizedKroneckerQuiver(3), vector([2,3]), vector([3,-2])
        sage: X = QuiverModuliSpace(Q, d, theta)
        sage: X.is_strongly_amply_stable()
        True

        sage: from quiver import *
        sage: Q, d = ThreeVertexQuiver(5,1,1), vector([4,1,4])
        sage: theta = Q.canonical_stability_parameter(d)
        sage: X = QuiverModuliSpace(Q, d, theta)
        sage: X.is_amply_stable()
        True
        sage: X.is_strongly_amply_stable()
        False

        """
        Q, d, theta, denominator = self._Q, self._d, self._theta, self._denominator

        # All subdimension vectors of d
        es = all_subdimension_vectors(d)
        # Remove 0 and d
        es.remove(Q.zero_vector())
        es.remove(d)
        # Filter out those of bigger slope
        es = list(filter(lambda e: slope(e, theta, denominator=denominator) >= slope(d, theta, denominator=denominator), es))
        return all([Q.euler_form(e, d-e) <= -2 for e in es])
    
    """
    Tautological relations
    """

    def all_forbidden_subdimension_vectors(self):
        r"""Returns the list of all subdimension vectors d' of d for which mu_theta(d') > mu_theta(d) (in the semistable case) or for which mu_theta(d') >= mu_theta(d) (in the stable case).
        
        OUTPUT: list of vectors
        """

        """
        EXAMPLES:

        sage: from quiver import *
        sage: Q, d, theta, condition = GeneralizedKroneckerQuiver(3), vector([3,3]), vector([1,-1]), "semistable"
        sage: X = QuiverModuliSpace(Q, d, theta, condition=condition)
        sage: X.all_forbidden_subdimension_vectors()
        [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)]
        sage: Q, d, theta, condition = GeneralizedKroneckerQuiver(3), vector([3,3]), vector([1,-1]), "stable"
        sage: Y = QuiverModuliSpace(Q, d, theta, condition=condition)
        sage: Y.all_forbidden_subdimension_vectors()
        [(1, 0), (1, 1), (2, 0), (2, 1), (2, 2), (3, 0), (3, 1), (3, 2)]

        """
        Q, d, theta, condition = self._Q, self._d, self._theta, self._condition

        properSubdimensions = list(filter(lambda e: e != d and e != Q.zero_vector(), all_subdimension_vectors(d)))
        
        if condition == "semistable":
            return list(filter(lambda e: slope(e,theta) > slope(d,theta), properSubdimensions))
        elif condition == "stable":
            return list(filter(lambda e: slope(e,theta) >= slope(d,theta), properSubdimensions))
        
    def all_minimal_forbidden_subdimension_vectors(self):
        r"""Returns the list of all minimal forbidden subdimension vectors of d.
        
        OUTPUT: list of vectors
        """

        """Minimality is with respect to the partial order e << d which means e_i <= d_i for every source i, e_j >= d_j for every sink j, and e_k = d_k for every vertex which is neither a source nor a sink."""

        """
        EXAMPLES:

        sage: from quiver import *
        sage: Q, d, theta, condition = GeneralizedKroneckerQuiver(3), vector([3,3]), vector([1,-1]), "semistable"
        sage: X = QuiverModuliSpace(Q, d, theta, condition=condition)
        sage: X.all_minimal_forbidden_subdimension_vectors()
        [(1, 0), (2, 1), (3, 2)]
        sage: Q, d, theta, condition = GeneralizedKroneckerQuiver(3), vector([3,3]), vector([1,-1]), "stable"
        sage: Y = QuiverModuliSpace(Q, d, theta, condition=condition)
        sage: Y.all_minimal_forbidden_subdimension_vectors()
        [(1, 1), (2, 2)]

        """
        Q = self._Q

        forbidden = self.all_forbidden_subdimension_vectors()
        return list(filter(lambda e: not any([Q.division_order(f,e) for f in list(filter(lambda f: f != e, forbidden))]), forbidden))
    
    def __tautological_presentation(self, inRoots=False, chernClasses=None, chernRoots=None):
        r"""Returns the tautological relations in Chern classes (if inRoots == False) or in Chern roots.
        
        INPUT:
        - ``inRoots``: Bool
        - ``chernClasses``: list of Strings
        - ``chernRoots``: list of Strings
        
        OUTPUT: dict
        """

        # TODO
        """Explanation ..."""

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

        """
        EXAMPLES:
        """
        Q, d, theta = self._Q, self._d, self._theta
        n = Q.number_of_vertices()
        a = Q.adjacency_matrix()

        if chernClasses == None:
            chernClasses = ['x%s_%s'%(i,r) for i in range(1,n+1) for r in range(1,d[i-1]+1)]
        if chernRoots == None:
            chernRoots = ['t%s_%s'%(i,r) for i in range(1,n+1) for r in range(1,d[i-1]+1)]
                
        R = PolynomialRing(QQ, chernRoots)

        def generator(R, i, r):
            r"""Returns generator(R, i, r) = t{i+1}_{r+1}."""
            return R.gen(r+sum([d[j] for j in range(i)]))

        """Generators of the tautological ideal regarded upstairs, i.e. in A*([R/T]).
        For a forbidden subdimension vector e of d, the forbidden polynomial in Chern roots is given by prod_{a: i --> j} prod_{r=1}^{e_i} prod_{s=e_j+1}^{d_j} (tj_s - ti_r) = prod_{i,j} prod_{r=1}^{e_i} prod_{s=e_j+1}^{d_j} (tj_s - ti_r)^{a_{ij}}."""
        forbiddenPolynomials = [prod([prod([(generator(R,j,s) - generator(R,i,r))**a[i,j]  for r in range(e[i]) for s in range(e[j],d[j])]) for i in range(n) for j in range(n)]) for e in self.all_minimal_forbidden_subdimension_vectors()]

        if inRoots:
            return {
                "ParentRing" : R,
                "Generators" : lambda i, r: generator(R, i, r),
                "Relations"  : forbiddenPolynomials
            }
        else:
            """delta is the discriminant"""
            delta = prod([prod([generator(R,i,l) - generator(R,i,k) for k in range(d[i]) for l in range(k+1,d[i])]) for i in range(n)])

            """longest is the longest Weyl group element when regarding W as a subgroup of S_{sum d_i}"""
            longest = []
            r = 0
            for i in range(n):
                longest = longest + list(reversed(range(r+1,r+d[i]+1)))
                r += d[i]
            W = Permutations(bruhat_smaller=longest)

            def antisymmetrization(f):
                """The antisymmetrization of f is the symmetrization divided by the discriminant."""
                # I don't want to define W and delta here but globally because then we need to
                # compute it just once. That's probably a bit faster.
                def permute(f, w):
                    return f.subs({R.gen(i): R.gen(w[i] - 1) for i in range(R.ngens())})

                return sum(w.sign() * permute(f, w) for w in W) // delta

            """Schubert basis of A^*([R/T]) over A^*([R/G])"""
            X = SchubertPolynomialRing(ZZ)
            supp = list(filter(lambda i: d[i] > 0, range(n)))
            B = lambda i: [X(p).expand() for p in Permutations(d[i])] 
            Bprime = [[f.parent().hom([generator(R,i,r) for r in range(f.parent().ngens())], R)(f) for f in B(i)] for i in supp]

            def product_lists(L):
                n = len(L)
                assert n > 0
                if n == 1:
                    return L[0]
                else:
                    P = product_lists([L[i] for i in range(n-1)])
                    return [p*l for p in P for l in L[n-1]]

            schubert = product_lists(Bprime)

            """Define A = A*([R/G])."""
            degrees = []
            for i in range(n):
                degrees = degrees+list(range(1,d[i]+1))
            A = PolynomialRing(QQ, chernClasses, order=TermOrder('wdegrevlex', degrees))

            E = SymmetricFunctions(ZZ).e()
            """The Chern classes of U_i on [R/G] are the elementary symmetric functions in the Chern roots ti_1,...,ti_{d_i}."""
            elementarySymmetric = []
            for i in range(n):
                elementarySymmetric = elementarySymmetric + [E([k]).expand(d[i], alphabet=[generator(R,i,r) for r in range(d[i])]) for k in range(1,d[i]+1)]
            """Map xi_r to the r-th elementary symmetric function in ti_1,...,ti_{d_i}."""
            inclusion = A.hom(elementarySymmetric, R)

            """Tautological relations in Chern classes."""
            tautological = [antisymmetrization(b * f) for b in schubert for f in forbiddenPolynomials]
            tautological = [inclusion.inverse_image(g) for g in tautological]

            return { 
                "ParentRing" : A,
                "Generators" : lambda i, r: generator(A, i, r),
                "Relations"  : tautological
            }

    def tautological_relations(self, inRoots=False, chernClasses=None, chernRoots=None):
        r"""Returns the tautological relations in Chern classes (if inRoots == False) or in Chern roots.
        
        INPUT:
        - ``inRoots``: Bool
        - ``chernClasses``: list of Strings
        - ``chernRoots``: list of Strings
        
        OUTPUT: list
        """

        taut = self.__tautological_presentation(inRoots=inRoots, chernClasses=chernClasses, chernRoots=chernRoots)
        return taut["Relations"]

    @abstractmethod
    def dimension(self):
        pass

    @abstractmethod
    def is_smooth(self):
        pass

    @abstractmethod
    def chow_ring(self):
        pass


class QuiverModuliSpace(QuiverModuli):

    def __init__(self, Q, d, theta, denominator=sum, condition="semistable"):
        QuiverModuli.__init__(self, Q, d, theta, denominator=denominator, condition=condition)

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
        sage: from quiver import *
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
        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(2)
        sage: theta = vector([1,-1])
        sage: d = vector([1,1])
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
        sage: from quiver import *
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
                    allLunaTypes = self.all_luna_types()
                    return max([self.dimension_of_luna_stratum(tau) for tau in allLunaTypes])
                else:
                    # I somehow like the convention that the dimension of the empty set is -Infinity
                    return -oo
            else:
                # self._condition == "stable"
                return -oo
    
    def poincare_polynomial(self):
        r"""Returns the Poincare polynomial of the moduli space.
        
        OUTPUT: polynomial in one variable"""

        r"""The Poincare polynomial is defined as $P_X(q) = \sum_{i \geq 0} (-1)^i dim H^i(X;\mathbb{C}) q^{i/2}$. For a quiver moduli space whose dimension vector is $\theta$-coprime, the odd cohomology vanishes and this is a Polynomial in $q$. We use Cor. 6.9 in Reineke's Harder--Narasimhan paper to compute it."""

        """
        EXAMPLES:

        sage: from quiver import *
        sage: Q, d, theta = KroneckerQuiver(), vector([1,1]), vector([1,-1])
        sage: X = QuiverModuliSpace(Q, d, theta)
        sage: X.poincare_polynomial()
        q + 1
        sage: Q, d, theta = GeneralizedKroneckerQuiver(3), vector([2,3]), vector([3,-2])
        sage: X = QuiverModuliSpace(Q, d, theta)
        sage: X.poincare_polynomial()
        q^6 + q^5 + 3*q^4 + 3*q^3 + 3*q^2 + q + 1
        sage: Q, d = SubspaceQuiver(5), vector([1,1,1,1,1,2])
        sage: theta = Q.canonical_stability_parameter(d)
        sage: X = QuiverModuliSpace(Q, d, theta)
        sage: X.poincare_polynomial()
        q^2 + 5*q + 1

        """

        Q, d, theta = self._Q, self._d, self._theta
        assert (is_coprime_for_stability_parameter(d, theta))

        k = FunctionField(QQ,'L')
        L = k.gen(0)
        K = FunctionField(QQ,'q')
        q = K.gen(0)

        X = QuiverModuliStack(Q, d, theta, condition="semistable")
        mot = X.motive()

        f = k.hom(q, K)

        return (1-q)*f(mot)       


    def betti_numbers(self):
        r"""Returns the Betti numbers of the moduli space.
        
        OUTPUT: List of Ints
        """

        """
        EXAMPLES:

        sage: from quiver import *
        sage: from moduli import *
        sage: Q, d, theta = KroneckerQuiver(), vector([1,1]), vector([1,-1])
        sage: X = QuiverModuliSpace(Q, d, theta, condition="semistable")
        sage: X.poincare_polynomial()
        q + 1
        sage: X.betti_numbers()
        [1, 0, 1]
        sage: Q, d, theta = GeneralizedKroneckerQuiver(3), vector([2,3]), vector([1,-1])
        sage: X = QuiverModuliSpace(Q, d, theta, condition="semistable")
        sage: X.betti_numbers()
        [1, 0, 1, 0, 3, 0, 3, 0, 3, 0, 1, 0, 1]

        """
        
        Q, d, theta = self._Q, self._d, self._theta
        assert (is_coprime_for_stability_parameter(d, theta))
        N = self.dimension()

        K = FunctionField(QQ,'q')
        q = K.gen(0)
        L = FunctionField(QQ,'v')
        v = L.gen(0)
        ext = K.hom(v**2, L)
        # p is the prime place of the DVR associated with v
        p = v.zeros()[0]

        f = ext(self.poincare_polynomial())
        betti = [f.evaluate(p)]
        for i in range(2*N):
            f = (f - f.evaluate(p))/v 
            betti = betti + [f.evaluate(p)]

        return betti


    def is_smooth(self):
        # if theta-coprime then you can shortcut everything
        # if theta != 0 reduce to theta = 0 using https://mathscinet.ams.org/mathscinet-getitem?mr=1972892 (Adriaenssens--Le Bruyn)
        # if theta = 0, then use https://mathscinet.ams.org/mathscinet-getitem?mr=1929191 (Bocklandt)
        if (self._condition == "stable") or (is_coprime_for_stability_parameter(self._d, self._theta)):
            return True
        else:
            raise NotImplementedError()
    
    def picard_rank(self):
        """Computes the Picard rank of the moduli space for known cases."""
        # TODO this should really be a check for theta belonging to the canonical chamber, rather than being equal to the canonical stability.
        # Comment: The Picard rank is also n-1 if theta is not in the canonical chamber (with d theta-coprime and theta-amply stable). This is because Pic(X) is isomorphic to A^1(X) (X is smooth) and the latter is free of rank n-1 by theh tautological presentation.
        if is_coprime_for_stability_parameter(self._d,self._theta) & self._Q.is_amply_stable(self._Q,self._d,self._theta):
            return self._Q.number_of_vertices() - 1
        else:
            raise NotImplementedError()
        
    def index(self):
        """Computes the index of the moduli space for known cases, i.e., the largest integer dividing the canonical divisor in Pic."""
        # TODO this should really be a check for theta belonging to the canonical chamber, rather than being equal to the canonical stability.
        if self._theta == self._Q.canonical_stability_parameter(self._d) & is_coprime_for_stability_parameter(self._d,self._theta) & self._Q.is_amply_stable(self._Q,self._d,self._theta):
            return gcd(self._theta.list())
        else:
            raise NotImplementedError()

    def mukai_inequality_holds(self):
        # TODO ample stability for the canonical stability parameter should be an attribute of the object, so that it is only computed once. Verbatim for many other attributes.
        return 1 - self._Q.euler_form(self._d, self._d) >= self.picard_rank() * (self.index() - 1)
    

    def chow_ring(self, chi=None, chernClasses=None):
        r"""Returns the Chow ring of the moduli space.
        
        INPUT:
        - ``chi``: vector of Ints
        - ``chernClasses``: list of Strings

        OUTPUT: ring
        """

        """
        EXAMPLES:

        The Kronecker quiver:
        sage: from quiver import *
        sage: Q, d, theta = KroneckerQuiver(), vector([1,1]), vector([1,-1])
        sage: X = QuiverModuliSpace(Q, d, theta, condition="semistable")
        sage: chi = vector([1,0])
        sage: A = X.chow_ring(chi=chi)
        sage: I = A.defining_ideal()
        sage: [I.normal_basis(i) for i in range(X.dimension()+1)]
        [[1], [x2_1]]


        The 3-Kronecker quiver:
        sage: from quiver import *
        sage: Q, d, theta = GeneralizedKroneckerQuiver(3), vector([2,3]), vector([3,-2])
        sage: X = QuiverModuliSpace(Q, d, theta, condition="semistable")
        sage: chi = vector([-1,1])
        sage: A = X.chow_ring(chi=chi)
        sage: I = A.defining_ideal()
        sage: [I.normal_basis(i) for i in range(X.dimension()+1)]
        [[1],
        [x2_1],
        [x1_2, x2_1^2, x2_2],
        [x2_1^3, x2_1*x2_2, x2_3],
        [x2_1^2*x2_2, x2_2^2, x2_1*x2_3],
        [x2_2*x2_3],
        [x2_3^2]]

        The 5-subspaces quiver:
        sage: from quiver import *
        sage: Q, d, theta = SubspaceQuiver(5), vector([1,1,1,1,1,2]), vector([2,2,2,2,2,-5])
        sage: X = QuiverModuliSpace(Q, d, theta, condition="semistable")
        sage: chi = vector([-1,-1,-1,-1,-1,3])
        sage: A = X.chow_ring(chi=chi)
        sage: I = A.defining_ideal()
        sage: [I.normal_basis(i) for i in range(X.dimension()+1)]
        [[1], [x2_1, x3_1, x4_1, x5_1, x6_1], [x6_2]]

        """
        Q, d, theta = self._Q, self._d, self._theta
        n = Q.number_of_vertices()

        # This implementation only works if d is theta-coprime which implies that d is indivisible. 
        assert is_coprime_for_stability_parameter(d, theta)

        if chernClasses == None:
            chernClasses = ['x%s_%s'%(i,r) for i in range(1,n+1) for r in range(1,d[i-1]+1)]

        taut = self._QuiverModuli__tautological_presentation(inRoots=False, chernClasses=chernClasses)
        A, generator, rels = taut["ParentRing"], taut["Generators"], taut["Relations"]

        if chi == None:
            [g,m] = extended_gcd(d.list())
            chi = vector(m)

        # TODO assert that chi has integer entries.
        """Make sure that chi has weight one, i.e. provides a retraction for X*(PG) --> X*(G)."""
        assert chi*d == 1

        I = A.ideal(rels) + A.ideal(sum([chi[i]*generator(i,0) for i in range(n)]))

        return QuotientRing(A, I, names=chernClasses)

    def chern_class_line_bundle(self, eta, chernClasses=None):
        """Returns the first Chern class of the line bundle L(eta) = bigotimes_{i in Q_0} det(U_i)^{-eta_i} where eta is a character of PG_d."""

        A = self.chow_ring(chi=None, chernClasses=chernClasses)
        n = self._Q.number_of_vertices()
        d = self._d

        return -sum([eta[i]*A.gen(sum([d[j] for j in range(i)])) for i in range(n)])

    def chern_character_line_bundle(self, eta, chernClasses=None):
        """Computes the Chern character of L(eta).
        The Chern character of a line bundle L with first Chern class x is given by e^x = 1 + x + x^2/2 + x^3/6 + ..."""

        A = self.chow_ring(chi=None, chernClasses=chernClasses)
        N = self.dimension()
        x = self.chern_class_line_bundle(eta, chernClasses=chernClasses)
        return sum([x**i/factorial(i) for i in range(N+1)])


    def total_chern_class_universal(self, i, chi, chernClasses=None):
        """Gives the total Chern class of the universal bundle U_i(chi)."""

        A = self.chow_ring(chi, chernClasses=chernClasses)
        d = self._d

        return 1 + sum([A.gen(r + sum([d[j] for j in range(i-1)])) for r in range(d[i-1])])


    def point_class(self, chi=None, chernClasses=None):
        """Returns the point class as an expression in Chern classes of the U_i(chi)."""

        """The point class is given as the homogeneous component of degree dim X of the expression prod_{a in Q_1} c(U_{t(a)})^{d_{s(a)}} / (prod_{i in Q_0} c(U_i)^{d_i})"""

        """
        EXAMPLES

        P^7 as a quiver moduli space of a generalized Kronecker quiver:
        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(8)
        sage: d = vector([1,1])
        sage: theta = vector([1,-1])
        sage: X = QuiverModuliSpace(Q,d,theta,condition="semistable")
        sage: chi = vector([1,0])
        sage: X.point_class(chi,chernClasses=['o','h'])
        h^7

        Our favorite 6-fold:
        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: d = vector([2,3])
        sage: theta = vector([3,-2])
        sage: X = QuiverModuliSpace(Q,d,theta,condition="semistable")
        sage: chi = vector([-1,1])
        sage: X.point_class(chi,chernClasses=['x1','x2','y1','y2','y3'])
        y3^2

        A moduli space of the 5-subspaces quiver; it agrees with the blow-up of P^2 in 4 points in general position:
        sage: from quiver import *
        sage: Q = SubspaceQuiver(5)
        sage: d = vector([1,1,1,1,1,2])
        sage: theta = vector([2,2,2,2,2,-5])
        sage: X = QuiverModuliSpace(Q,d,theta,condition="semistable")
        sage: chi = vector([-1,-1,-1,-1,-1,3])
        sage: X.point_class(chi,chernClasses=['x1','x2','x3','x4','x5','y','z'])
        1/2*z

        """

        Q, d = self._Q, self._d
        n = Q.number_of_vertices()
        a = Q.adjacency_matrix()
        N = self.dimension()

        A = self.chow_ring(chi=chi, chernClasses=chernClasses)
        pi = A.cover() # The quotient map
        sect = A.lifting_map() # A choice of a section of pi

        if chi == None:
            [g,m] = extended_gcd(d.list())
            chi = vector(m)

        numerator = prod([self.total_chern_class_universal(j+1, chi, chernClasses=chernClasses)**(d*a.column(j)) for j in range(n)])
        denominator = prod([self.total_chern_class_universal(i+1, chi, chernClasses=chernClasses)**d[i] for i in range(n)])

        quotient = numerator/denominator

        return pi(sect(quotient).homogeneous_components()[N])
    
    def degree(self, eta=None, chernClasses=None):
        r"""Computes the degree of the ample line bundle given by eta."""
        # TODO: Need check for ampleness first

        if eta == None:
            eta = self._Q.canonical_stability_parameter(self._d)

        A = self.chow_ring(chi=None, chernClasses=chernClasses)
        N = self.dimension()
        c = self.chern_class_line_bundle(eta, chernClasses=chernClasses)
        p = self.point_class(chernClasses=chernClasses)
        return c**N/p

    def todd_class(self):
        """The Todd class of X is the Todd class of the tangent bundle. For quiver moduli it computes as
        td(X) = (prod_{a:i->j in Q_1} prod_{p=1}^{d_j} prod_{q=1}^{d_i} Q(t_{j,q} - t_{i,p}))/(prod_{i in Q_0} prod_{p,q=1}^{d_i} Q(t_{i,q} - t_{i,p}))
        """

        def todd_generating_series(t,n):
            """We call the series Q(t) = t/(1-e^{-t}) the Todd generating series. The function computes the terms of this series up to degree n."""
            B = [bernoulli(i) for i in range(n+1)]
            return sum([(-1)^i*B[i]/factorial(i)*t^i for i in range(n+1)])

        def truncate(f,n):
            """Takes an element in a graded ring and discards all homogeneous components of degree > n"""
            hom = f.homogeneous_components()
            keyList = [i for i in hom]
            return sum([hom[i] for i in filter(lambda i: i <= n, keyList)])

        raise NotImplementedError()

    # TODO: This is maybe too specific.
    # def diagonal(self, chi=None):
    #     """Computes the class of the diagonal in the Chow ring of X x X where X is the quiver moduli space."""
    #     """It is given by the homogeneous component of degree dim X = 1 - <d,d> of the expression c(F)/C(E), where E = bigoplus_{i in Q_0} U_i^vee boxtimes U_i and F = bigoplus_{a in Q_1} U_{s(a)}^vee boxtimes U_{t(a)} = bigoplus_{i,j in Q_0} (U_i^vee boxtimes U_j)^{a_ij}."""

    #     """
    #     EXAMPLES

    #     P^2 as a quiver moduli space:
    #     sage: from quiver import *
    #     sage: Q = GeneralizedKroneckerQuiver(3)
    #     sage: d = vector([1,1])
    #     sage: theta = vector([1,-1])
    #     sage: X = QuiverModuliSpace(Q,d,theta,condition="semistable")
    #     sage: X.diagonal()
    #     x1_1^2 + x1_1*y1_1 + y1_1^2

    #     """

    #     Q, d, theta = self._Q, self._d, self._theta
    #     n = Q.number_of_vertices()
    #     N = self.dimension()
    #     a = Q.adjacency_matrix()

    #     di = self._QuiverModuli__tautological_presentation()
    #     A = di["Generators"]
    #     I = di["Relations"] + A.ideal(chi)

    #     chernClasses1 = ['x%s_%s'%(i,r) for i in range(1,n+1) for r in range(1,d[i-1]+1)]
    #     chernClasses2 = ['y%s_%s'%(i,r) for i in range(1,n+1) for r in range(1,d[i-1]+1)]
    #     chernClasses = chernClasses1+chernClasses2

    #     AxA = PolynomialRing(QQ,chernClasses)
    #     inclusion1 = A.hom(chernClasses1,AxA)
    #     inclusion2 = A.hom(chernClasses2,AxA)
    #     B = QuotientRing(AxA,inclusion1(I) + inclusion2(I),names=chernClasses)

    #     pi = B.cover() # The quotient map AxA --> B
    #     sect = B.lifting_map() # A choice of a section of pi

    #     chernRoots1 = ['t%s_%s'%(i,r) for i in range(1,n+1) for r in range(1,d[i-1]+1)]
    #     chernRoots2 = ['u%s_%s'%(i,r) for i in range(1,n+1) for r in range(1,d[i-1]+1)]
    #     chernRoots = chernRoots1+chernRoots2
    #     RxR = PolynomialRing(QQ,chernRoots)

    #     def generatorRxR1(i,r):
    #         """Returns generatorRxR1(i,r) = t{i+1}_{r+1}."""
    #         return RxR.gen(r + sum([d[j] for j in range(i)]))

    #     def generatorRxR2(i,r):
    #         """Returns generatorRxR2(i,r) = u{i+1}_{r+1}."""
    #         return RxR.gen(sum([d[j] for j in range(n)]) + r + sum([d[j] for j in range(i)]))

    #     E = SymmetricFunctions(ZZ).e()
    #     elementarySymmetric1 = []
    #     elementarySymmetric2 = []
    #     for i in range(n):
    #         elementarySymmetric1 = elementarySymmetric1 + [E([k]).expand(d[i], alphabet=[generatorRxR1(i,r) for r in range(d[i])]) for k in range(1,d[i]+1)]
    #         elementarySymmetric2 = elementarySymmetric2 + [E([k]).expand(d[i], alphabet=[generatorRxR2(i,r) for r in range(d[i])]) for k in range(1,d[i]+1)]
    #     elementarySymmetric = elementarySymmetric1 + elementarySymmetric2
    #     """Map xi_r to the r-th elementary symmetric function in ti_1,...,ti_{d_i} and yi_r to the same in ui_1,...,ui_{d_i}."""
    #     inclusion = AxA.hom(elementarySymmetric, RxR)

    #     def total_chern_class_boxproduct(i,j):
    #         """Computes the total Chern class of U_i^vee boxtimes U_j"""
    #         c = prod([(1-generatorRxR1(i,r)+generatorRxR2(j,s)) for r in range(d[i]) for s in range(d[j])])
    #         return pi(inclusion.inverse_image(c))

    #     numerator = prod([total_chern_class_boxproduct(i,j)**a[i,j] for i in range(n) for j in range(n)])
    #     denominator = prod([total_chern_class_boxproduct(i,i) for i in range(n)])
    #     quotient = numerator/denominator

    #     return pi(sect(quotient).homogeneous_components()[N])


class QuiverModuliStack(QuiverModuli):

    def __init__(self, Q, d, theta, denominator=sum, condition="semistable"):
        QuiverModuli.__init__(self, Q, d, theta, denominator=denominator, condition=condition)

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
    
    def motive(self):
        r"""Gives an expression for the motive of the semistable moduli stack in an appropriate localization of K_0(Var)"""

        """
        EXAMPLES:

        sage: from quiver import *
        sage: Q, d, theta = LoopQuiver(0), vector([2]), vector([0])
        sage: X = QuiverModuliStack(Q, d, theta, condition="semistable")
        sage: X.motive()
        1/(L^4 - L^3 - L^2 + L)
        sage: Q, d, theta = LoopQuiver(1), vector([2]), vector([0])
        sage: X = QuiverModuliStack(Q, d, theta, condition="semistable")
        sage: X.motive()
        L^3/(L^3 - L^2 - L + 1)
        sage: Q, d, theta = GeneralizedKroneckerQuiver(3), vector([2,3]), vector([3,-2])
        sage: X = QuiverModuliStack(Q, d, theta, condition="semistable")
        sage: X.motive()
        (-L^6 - L^5 - 3*L^4 - 3*L^3 - 3*L^2 - L - 1)/(L - 1)

        """

        # Only for semistable. For stable, we don't know what the motive is. It's not pure in general.
        assert self._condition == "semistable"

        Q, d, theta = self._Q, self._d, self._theta
        n = Q.number_of_vertices()

        if theta == Q.zero_vector():
            K = FunctionField(QQ,'L')
            L = K.gen(0)
            num = L**(-Q.tits_form(d))
            den = prod([prod([(1-L**(-nu)) for nu in range(1,d[i]+1)]) for i in range(n)])
            return num/den
        else:
            I = all_subdimension_vectors(d)
            I = list(filter(lambda e: e != Q.zero_vector() and e != d , I))
            I = list(filter(lambda e: slope(e, theta) > slope(d, theta), I))
            I = I + [Q.zero_vector(), d]
            I.sort(key=(lambda e: deglex_key(e, b=max(d)+1)))

            K = FunctionField(QQ,'L')
            L = K.gen(0)

            # Now define a matrix T of size NxN whose entry at position (i,j) is L^<e-f,e>*mot(f-e) if e = I[i] is a subdimension vector of f = I[j] and 0 otherwise
            mot = lambda e: QuiverModuliStack(Q, e, Q.zero_vector(), condition="semistable").motive()
            N = len(I)
            T = matrix(K, N)
            for i in range(N):
                for j in range(i,N):
                    e, f = I[i], I[j]
                    if is_subdimension_vector(e, f):
                        T[i,j] = L**(Q.euler_form(e-f, e))*mot(f-e)

            # Solve system of linear equations T*x = e_N and extract entry 0 of the solution x.
            y = vector([0 for i in range(N)])
            y[N-1] = 1
            x = T.solve_right(y)

            return x[0]
        
    def chow_ring(self, chernClasses=None):
        r"""Returns the Chow ring of the quotient stack.
        
        INPUT:
        - ``chernClasses``: list of Strings

        OUTPUT: ring
        """

        Q, d = self._Q, self._d
        n = Q.number_of_vertices()

        if chernClasses == None:
            chernClasses = ['x%s_%s'%(i,r) for i in range(1,n+1) for r in range(1,d[i-1]+1)]

        taut = self._QuiverModuli__tautological_presentation(inRoots=False, chernClasses=chernClasses)
        A, generator, rels = taut["ParentRing"], taut["Generators"], taut["Relations"]

        return QuotientRing(A, A.ideal(rels), names=chernClasses)



class SmoothModel:
    """How about this: instead of a separate class SmoothModel, we could define a method framed_moduli_space(self,n) inside the class QuiverModuliSpace which returns another quiver moduli space. After all, it is itself a quiver moduli space."""
    def __init__(self):
        pass

    def betti_numbers(self):
        raise NotImplementedError()

"""Auxiliary methods:"""

def extended_gcd(x):
    """Computes the gcd and the Bezout coefficients of a list of integers."""
    # This exists for two integers but there seems to be no implementation for more than one.
    # That's astonishing.

    n = len(x)
    if n == 1:
        return [x,[1]]
    if n == 2:
        (g,a,b) = xgcd(x[0],x[1])
        return [g,[a,b]]
    if n > 2:
        (g,a,b) = xgcd(x[0],x[1])
        y = [g]+[x[i] for i in range(2,n)]
        [d,c] = extended_gcd(y)
        m = [c[0]*a,c[0]*b]+[c[i] for i in range(1,n-1)]
        return [d,m]

