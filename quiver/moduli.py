from itertools import combinations_with_replacement, product

from sage.arith.misc import bernoulli, factorial, gcd, xgcd
from sage.combinat.partition import Partitions
from sage.combinat.permutation import Permutations
from sage.combinat.schubert_polynomial import SchubertPolynomialRing
from sage.combinat.sf.sf import SymmetricFunctions
from sage.matrix.constructor import matrix
from sage.misc.misc_c import prod
from sage.modules.free_module_element import vector
from sage.rings.function_field.constructor import FunctionField
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.term_order import TermOrder
from sage.rings.quotient_ring import QuotientRing
from sage.rings.rational_field import QQ
from sage.structure.element import Element

from quiver import Quiver

"""Defines how permutations are multiplied."""
Permutations.options(mult="r2l")


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


class QuiverModuli(Element):
    def __init__(self, Q, d, theta=None, denom=sum, condition="semistable"):
        r"""Constructor for an abstract quiver moduli space

        This base class contains everything that is common between
        - quiver moduli spaces, i.e., varieties
        - quiver moduli stacks

        INPUT:

        - ``Q`` -- quiver

        - ``d`` --- dimension vector

        - ``theta`` -- stability parameter (default: canonical stability parameter)

        - ``denom`` -- denominator for slope stability (default: ``sum``), needs to be
          effective on the simple roots

        - ``condition`` -- whether to include all semistables, or only stables
          (default: "semistable")

        See :class:`QuiverModuliSpace` and :class:`QuiverModuliStack` for more details.

        EXAMPLES:

        We can instantiate an abstract quiver moduli space::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: X = QuiverModuli(Q, [2, 3])
            sage: X
            abstract moduli of semistable representations, with
            Q = 3-Kronecker quiver
            d = [2, 3]
            theta = (9, -6)

        It has functionality common to both varieties and stacks, i.e., when it really
        concerns something involving the representation variety::

            sage: X.all_harder_narasimhan_types()
            [((1, 0), (1, 1), (0, 2)),
             ((1, 0), (1, 2), (0, 1)),
             ((1, 0), (1, 3)),
             ((1, 1), (1, 2)),
             ((2, 0), (0, 3)),
             ((2, 1), (0, 2)),
             ((2, 2), (0, 1)),
             ((2, 3),)]

        But things like dimension depend on whether we consider it as a variety or as
        a stack, and thus these are not implemented::

            sage: X.dimension()
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        if theta is None:
            theta = Q.canonical_stability_parameter(d)

        assert Q._is_dimension_vector(d), "`d` is not a dimension vector of `Q`"
        assert Q._is_vector(theta), "`theta` is not a stability parameter for `Q`"
        assert condition in ["semistable", "stable"]
        # TODO this effectivity condition needs to be documented, and maybe be part of Quiver?
        assert all(
            denom(Q._coerce_dimension_vector(Q.simple_root(i))) > 0
            for i in Q.vertices()
        ), "denominator needs to be effective"

        self._Q = Q
        self._d = d
        self._theta = theta
        self._denom = denom
        self._condition = condition

    def _repr_(self):
        r"""
        Give a shorthand string presentation for an abstract quiver moduli space
        """
        if self.get_custom_name():
            return self.get_custom_name()

        output = "abstract moduli of {} representations, with".format(self._condition)
        output += "\nQ = {}\nd = {}\ntheta = {}".format(
            self._Q.repr(), self._d, self._theta
        )

        return output

    def repr(self):
        return self._repr_()

    def quiver(self):
        return self._Q

    def dimension_vector(self):
        return self._d

    def stability_parameter(self):
        return self._theta

    def denominator(self):
        return self._denom

    def is_nonempty(self) -> bool:
        if self._condition == "stable":
            return self._Q.has_stable_representation(self._d, self._theta)
        if self._condition == "semistable":
            return self._Q.has_semistable_representation(self._d, self._theta)

    """
    HN business
    """

    def all_harder_narasimhan_types(self, proper=False, sorted=False):
        r"""
        Returns the list of all HN types.

        A Harder--Narasimhan (HN) type of :math:`d` with respect to :math:`\theta`
        is a sequence :math:`d^* = (d^1,...,d^s)` of dimension vectors such that

        - :math:`d^1 + ... + d^s = d`
        - :math:`\mu_{\theta}(d^1) > ... > \mu_{\theta}(d^s)`
        - Every :math:`d^k` is :math:`\theta`-semi-stable.

        INPUT:

        - ``proper`` -- (default: False) whether to exclude the HN type corresponding to the stable locus

        - ``sorted`` -- (default: False) whether to sort the HN-types according to the given slope

        OUTPUT: list of tuples of dimension vectors encoding HN-types

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, [2, 3])
            sage: X.all_harder_narasimhan_types()
            [((1, 0), (1, 1), (0, 2)),
             ((1, 0), (1, 2), (0, 1)),
             ((1, 0), (1, 3)),
             ((1, 1), (1, 2)),
             ((2, 0), (0, 3)),
             ((2, 1), (0, 2)),
             ((2, 2), (0, 1)),
             ((2, 3),)]
            sage: X.all_harder_narasimhan_types(proper=True)
            [((1, 0), (1, 1), (0, 2)),
             ((1, 0), (1, 2), (0, 1)),
             ((1, 0), (1, 3)),
             ((1, 1), (1, 2)),
             ((2, 0), (0, 3)),
             ((2, 1), (0, 2)),
             ((2, 2), (0, 1))]
            sage: d = [2, 3]
            sage: theta = -Q.canonical_stability_parameter(d)
            sage: Y = QuiverModuliSpace(Q, d, theta)
            sage: Y.all_harder_narasimhan_types()
            [((0, 3), (2, 0))]

        A 3-vertex quiver::

            sage: from quiver import *
            sage: Q = ThreeVertexQuiver(2, 3, 4)
            sage: d = [2, 3, 2]
            sage: Z = QuiverModuliSpace(Q, [2, 3, 2])
            sage: Z.all_harder_narasimhan_types()
            [((0, 1, 0), (1, 2, 1), (1, 0, 1)),
             ((0, 1, 0), (2, 0, 1), (0, 2, 1)),
             ((0, 1, 0), (2, 1, 1), (0, 1, 1)),
             ((0, 1, 0), (2, 2, 1), (0, 0, 1)),
             ((0, 1, 0), (2, 2, 2)),
             ((0, 2, 0), (1, 1, 1), (1, 0, 1)),
             ((0, 2, 0), (2, 0, 1), (0, 1, 1)),
             ((0, 2, 0), (2, 1, 1), (0, 0, 1)),
             ((0, 2, 0), (2, 1, 2)),
             ((0, 3, 0), (2, 0, 1), (0, 0, 1)),
             ((0, 3, 0), (2, 0, 2)),
             ((1, 0, 0), (0, 1, 0), (1, 0, 1), (0, 2, 1)),
             ((1, 0, 0), (0, 1, 0), (1, 1, 1), (0, 1, 1)),
             ((1, 0, 0), (0, 1, 0), (1, 2, 1), (0, 0, 1)),
             ((1, 0, 0), (0, 1, 0), (1, 2, 2)),
             ((1, 0, 0), (0, 2, 0), (1, 0, 1), (0, 1, 1)),
             ((1, 0, 0), (0, 2, 0), (1, 1, 1), (0, 0, 1)),
             ((1, 0, 0), (0, 2, 0), (1, 1, 2)),
             ((1, 0, 0), (0, 3, 0), (1, 0, 1), (0, 0, 1)),
             ((1, 0, 0), (0, 3, 0), (1, 0, 2)),
             ((1, 0, 0), (0, 3, 1), (1, 0, 1)),
             ((1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)),
             ((1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 1, 2)),
             ((1, 0, 0), (1, 1, 0), (0, 2, 0), (0, 0, 2)),
             ((1, 0, 0), (1, 1, 0), (0, 2, 1), (0, 0, 1)),
             ((1, 0, 0), (1, 1, 0), (0, 2, 2)),
             ((1, 0, 0), (1, 1, 1), (0, 2, 1)),
             ((1, 0, 0), (1, 2, 0), (0, 1, 0), (0, 0, 2)),
             ((1, 0, 0), (1, 2, 0), (0, 1, 1), (0, 0, 1)),
             ((1, 0, 0), (1, 2, 0), (0, 1, 2)),
             ((1, 0, 0), (1, 2, 1), (0, 1, 1)),
             ((1, 0, 0), (1, 3, 1), (0, 0, 1)),
             ((1, 0, 0), (1, 3, 2)),
             ((1, 1, 0), (0, 1, 0), (1, 0, 1), (0, 1, 1)),
             ((1, 1, 0), (0, 1, 0), (1, 1, 1), (0, 0, 1)),
             ((1, 1, 0), (0, 1, 0), (1, 1, 2)),
             ((1, 1, 0), (0, 2, 0), (1, 0, 1), (0, 0, 1)),
             ((1, 1, 0), (0, 2, 0), (1, 0, 2)),
             ((1, 1, 0), (1, 0, 1), (0, 2, 1)),
             ((1, 1, 0), (1, 1, 1), (0, 1, 1)),
             ((1, 1, 0), (1, 2, 0), (0, 0, 2)),
             ((1, 1, 0), (1, 2, 1), (0, 0, 1)),
             ((1, 1, 0), (1, 2, 2)),
             ((1, 2, 0), (0, 1, 0), (1, 0, 1), (0, 0, 1)),
             ((1, 2, 0), (0, 1, 0), (1, 0, 2)),
             ((1, 2, 0), (1, 0, 1), (0, 1, 1)),
             ((1, 2, 0), (1, 1, 1), (0, 0, 1)),
             ((1, 2, 0), (1, 1, 2)),
             ((1, 2, 1), (1, 1, 1)),
             ((1, 3, 1), (1, 0, 1)),
             ((2, 0, 0), (0, 1, 0), (0, 2, 1), (0, 0, 1)),
             ((2, 0, 0), (0, 1, 0), (0, 2, 2)),
             ((2, 0, 0), (0, 2, 0), (0, 1, 1), (0, 0, 1)),
             ((2, 0, 0), (0, 2, 0), (0, 1, 2)),
             ((2, 0, 0), (0, 2, 1), (0, 1, 1)),
             ((2, 0, 0), (0, 3, 0), (0, 0, 2)),
             ((2, 0, 0), (0, 3, 1), (0, 0, 1)),
             ((2, 0, 0), (0, 3, 2)),
             ((2, 0, 1), (0, 3, 1)),
             ((2, 1, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)),
             ((2, 1, 0), (0, 1, 0), (0, 1, 2)),
             ((2, 1, 0), (0, 2, 0), (0, 0, 2)),
             ((2, 1, 0), (0, 2, 1), (0, 0, 1)),
             ((2, 1, 0), (0, 2, 2)),
             ((2, 1, 1), (0, 2, 1)),
             ((2, 2, 0), (0, 1, 0), (0, 0, 2)),
             ((2, 2, 0), (0, 1, 1), (0, 0, 1)),
             ((2, 2, 0), (0, 1, 2)),
             ((2, 2, 1), (0, 1, 1)),
             ((2, 3, 0), (0, 0, 2)),
             ((2, 3, 1), (0, 0, 1)),
             ((2, 3, 2),)]
        """

        d = self._Q._coerce_dimension_vector(self._d)
        theta = self._Q._coerce_vector(self._theta)

        all_types = self._Q.all_hn_types(d, theta, denom=self._denom, sorted=sorted)
        if proper and (d,) in all_types:
            all_types.remove((d,))

        return all_types

    def is_harder_narasimhan_type(self, dstar) -> bool:
        r"""
        Checks if ``dstar`` is a HN type.

        A Harder--Narasimhan (HN) type of `d` with respect to :math:`\theta`
        is a sequence :math:`d^* = (d^1,...,d^s)` of dimension vectors such that

        - :math:`d^1 + ... + d^s = d`
        - :math:`\mu_{\theta}(d^1) > ... > \mu_{\theta}(d^s)`
        - Every :math:`d^k` is theta-semi-stable.

        INPUT:

        - ``dstar`` -- list of vectors of Ints

        OUTPUT: statement truth value as Bool

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, [2, 3], [1, 0])
            sage: HNs = X.all_harder_narasimhan_types()
            sage: all(X.is_harder_narasimhan_type(dstar) for dstar in HNs)
            True
            sage: dstar = [[1, 0], [1, 0], [0, 3]]
            sage: X.is_harder_narasimhan_type(dstar)
            False
            sage: X.is_harder_narasimhan_type([Q.zero_vector()])
            False

        """
        # setup shorthand
        Q, d, theta, denom = (
            self._Q,
            self._d,
            self._theta,
            self._denom,
        )

        dstar = list(map(lambda di: Q._coerce_dimension_vector(di), dstar))

        # first condition: sum to dimension vector
        if Q._coerce_dimension_vector(d) != sum(dstar):
            return False

        # second condition: decreasing slopes
        if not all(
            (
                Q.slope(dstar[i], theta, denom=denom)
                > Q.slope(dstar[i + 1], theta, denom=denom)
            )
            for i in range(len(dstar) - 1)
        ):
            return False

        if not all(
            Q.has_semistable_representation(di, theta, denom=denom) for di in dstar
        ):
            return False

        return True

    def codimension_of_harder_narasimhan_stratum(self, dstar, secure=False):
        """
        Computes the codimension of the HN stratum of ``dstar``
        inside the representation variety.

        INPUT:

        - ``dstar`` -- list of vectors of Ints
        - ``secure`` -- (default: False): Bool

        OUTPUT: codimension as Int
        # TODO
        # It checks for dstar to be a HN type iff secure == True. This check is slow.
        # Be sure to be dealing with a HN type if you call it with secure == False. This is fast but yields nonsense, if dstar is not a HN type.

        The codimension of the HN stratum of :math:`d^* = (d^1,...,d^s)` is given by

        .. MATH::

            - sum_{k < l} <d^k,d^l>

        EXAMPLES

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q, d = GeneralizedKroneckerQuiver(3), vector([2,3])
            sage: theta = vector([1,0])
            sage: X = QuiverModuliSpace(Q, d, theta)
            sage: hn = X.all_harder_narasimhan_types(); hn
            [((1, 0), (1, 1), (0, 2)),
             ((1, 0), (1, 2), (0, 1)),
             ((1, 0), (1, 3)),
             ((1, 1), (1, 2)),
             ((2, 0), (0, 3)),
             ((2, 1), (0, 2)),
             ((2, 2), (0, 1)),
             ((2, 3),)]
            sage: [X.codimension_of_harder_narasimhan_stratum(dstar) for dstar in hn]
            [12, 9, 8, 3, 18, 10, 4, 0]

        """
        Q = self._Q

        assert all(Q._is_dimension_vector(di) for di in dstar)

        if secure:
            assert self.is_harder_narasimhan_type(dstar)

        return -sum(
            Q.euler_form(dstar[k], dstar[l])
            for k in range(len(dstar) - 1)
            for l in range(k + 1, len(dstar))
        )

    def codimension_unstable_locus(self):
        r"""
        Computes the codimension of the unstable locus
        inside the representation variety.

        OUTPUT: codimension as Int

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, [2, 3], [1, 0])
            sage: X.codimension_unstable_locus()
            3

        A 3-vertex quiver::

            sage: Q = ThreeVertexQuiver(1,6,1)
            sage: X = QuiverModuliSpace(Q, [1, 6, 6])
            sage: X.codimension_unstable_locus()
            1

        The Kronecker quiver::

            sage: Q = GeneralizedKroneckerQuiver(1)
            sage: X = QuiverModuliSpace(Q, [2, 3], [1, 0])
            sage: X.codimension_unstable_locus()
            0

        """
        HNs = self.all_harder_narasimhan_types(proper=True)

        # note that while the HN types and strata depend on the denominator
        # the list of all their codimensions does not
        return min(
            self.codimension_of_harder_narasimhan_stratum(dstar, secure=False)
            for dstar in HNs
        )

    """
    Luna
    """

    def all_luna_types(self):
        r"""
        Returns the unordered list of all Luna types of d for theta.

        OUTPUT: list of tuples containing Int-vector and Int

        A Luna type of d for theta is an unordered sequence (i.e. multiset)
        :math:`((d^1,m_1),...,(d^s,m_s))` of dimension vectors
        :math:`d^k` and positive integers :math:`m_k` such that
        
        - :math:`m_1d^1 + ... + m_sd^s = d`
        - :math:`\mu_{\theta}(d^k) = \mu_{\theta}(d)`
        - All :math:`d^k` admit a :math:`\theta`-stable representation

        Example: Suppose that `d = 3e` and `e, 2e, d = 3e`
        are the only stable subdimension vectors.
        Then the Luna types are
        
        .. MATH::

            \begin{aligned}
                ((3e,1)) \\
                ((2e,1),(e,1))\\
                ((e,3))\\
                ((e,2),(e,1))\\
                ((e,1),(e,1),(e,1)).
            \end{aligned}

        We implement it as follows.
        
        A Luna type for us is a dictionary
        ``{d^1: p_1^1,..., d^s: p_s^1,..., d_s: p_s^t}``
        of dimension vectors d^k and non-empty partitions p^k such that

        .. MATH::

            |p_1^1|d^1 + ... + |p_s^t|d^s = d

        So in the above example, the Luna types are::

            {3e: [1]}
            {2e: [1], e: [1]}
            {e: [3]}
            {e: [2, 1]}
            {e: [1, 1, 1]}

        EXAMPLES:

        The Kronecker quiver::

            sage: from quiver import *
            sage: Q, d, theta = KroneckerQuiver(), vector([3,3]), vector([1,-1])
            sage: X = QuiverModuliSpace(Q, d, theta)
            sage: X.all_luna_types()
            [{(1, 1): [3]}, {(1, 1): [2, 1]}, {(1, 1): [1, 1, 1]}]

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q, d = GeneralizedKroneckerQuiver(3), vector([3,3])
            sage: theta = vector([1,-1])
            sage: X = QuiverModuliSpace(Q, d, theta)
            sage: X.all_luna_types()
            [{(3, 3): [1]},
             {(1, 1): [1], (2, 2): [1]},
             {(1, 1): [3]},
             {(1, 1): [2, 1]},
             {(1, 1): [1, 1, 1]}]

        The zero vector::

            sage: from quiver import *
            sage: Q, d, theta = KroneckerQuiver(), vector([0,0]), vector([1,-1])
            sage: X = QuiverModuliSpace(Q, d, theta)
            sage: X.all_luna_types()
            [{(0, 0): [1]}]

        """
        # setup shorthand
        Q, d, theta, denom = (
            self._Q,
            self._d,
            self._theta,
            self._denom,
        )

        d = Q._coerce_dimension_vector(d)

        if d == Q.zero_vector():
            # Q.zero_vector() can't be hashed a priori
            z = Q._coerce_vector(Q.zero_vector())
            return [{z: [1]}]

        same_slope = filter(
            lambda e: Q.slope(e, theta, denom=denom) == Q.slope(d, theta, denom=denom),
            Q.all_subdimension_vectors(d, nonzero=True, forget_labels=True),
        )
        same_slope = list(
            filter(
                lambda e: Q.has_stable_representation(e, theta, denom=denom),
                same_slope,
            )
        )

        bound = (sum(d) / min(sum(e) for e in same_slope)).ceil()
        luna_types = []
        for i in range(1, bound + 1):
            for tau in combinations_with_replacement(same_slope, i):
                if sum(tau) == d:
                    # from tau we build all possible Luna types
                    partial = {}
                    for taui in tuple(tau):
                        if taui in partial.keys():
                            partial[taui] += 1
                        else:
                            partial[taui] = 1

                    # partial has the form
                    # {d^1: Partitions(p^1), ..., d^s: Partitions(p^s)}
                    for key in partial.keys():
                        partial[key] = Partitions(partial[key]).list()

                    new_types = [
                        dict(zip(partial.keys(), values))
                        for values in product(*partial.values())
                    ]
                    luna_types += new_types
        return luna_types

    def is_luna_type(self, tau) -> bool:
        r"""
        Checks if tau is a Luna type for theta.

        INPUT:

        - ``tau`` -- dictionary with dimension vectors as keys and lists of ints as values

        OUTPUT: whether ``tau`` is a Luna type.

        EXAMPLES:

        The Kronecker quiver::

            sage: from quiver import *
            sage: Q, d, theta = KroneckerQuiver(), vector([3,3]), vector([1,-1])
            sage: X = QuiverModuliSpace(Q, d, theta)
            sage: l = X.all_luna_types()
            sage: all(X.is_luna_type(tau) for tau in l)
            True

        The 3-Kronecker quiver with zero vector::

            sage: from quiver import *
            sage: Q, d, theta = KroneckerQuiver(), vector([0,0]), vector([1,-1])
            sage: X = QuiverModuliSpace(Q, d, theta)
            sage: d.set_immutable()
            sage: X.is_luna_type({d: [1]})
            True

        """
        Q, d, theta, denom = (
            self._Q,
            self._d,
            self._theta,
            self._denom,
        )

        d = Q._coerce_dimension_vector(d)

        n = Q.number_of_vertices()
        assert all(len(dn) == n for dn in tau.keys())
        assert d == sum(k * dim for k in tau.keys() for dim in tau[k])

        if d == Q.zero_vector():
            # Q.zero_vector() can't be hashed a priori
            z = Q._coerce_vector(Q.zero_vector())
            return tau == {z: [1]}

        return all(
            Q.slope(key, theta, denom=denom) == Q.slope(d, theta, denom=denom)
            and Q.has_semistable_representation(key, theta, denom=denom)
            for key in tau.keys()
        )

    def dimension_of_luna_stratum(self, tau, secure=True):
        r"""
        Computes the dimension of the Luna stratum :math:`S_\tau`.

        INPUT:

        - ``tau`` -- list of tuples
        - ``secure`` -- Bool

        OUTPUT: Dimension as Int

        The dimension of the Luna stratum of
        ``tau = {d^1: p^1,...,d^s: p^s}`` is
        :math:`\sum_k l(p^k)(1 - <\langle d^k,d^k\rangle)`,
        where for a partition :math:`p = (n_1,...,n_l)`,
        the length `l(p)` is `l`, i.e. the number of summands.

        EXAMPLES:

        The Kronecker quiver::

            sage: from quiver import *
            sage: Q, d, theta = KroneckerQuiver(), vector([2,2]), vector([1,-1])
            sage: X = QuiverModuliSpace(Q, d, theta)
            sage: L = X.all_luna_types(); L
            [{(1, 1): [2]}, {(1, 1): [1, 1]}]
            sage: [X.dimension_of_luna_stratum(tau) for tau in L]
            [1, 2]

        """

        if secure:
            assert self.is_luna_type(tau)
        return sum(len(tau[dn]) * (1 - self._Q.euler_form(dn, dn)) for dn in tau.keys())

    def local_quiver_setting(self, tau, secure=True):
        r"""
        Returns the local quiver and dimension vector for the given Luna type.

        INPUT:

        - ``tau`` -- list of tuples
        - ``secure`` -- Bool

        OUTPUT: tuple consisting of a Quiver object and a vector

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q, d = GeneralizedKroneckerQuiver(3), vector([2,2])
            sage: theta = vector([1,-1])
            sage: X = QuiverModuliSpace(Q, d, theta)
            sage: L = X.all_luna_types(); L
            [{(2, 2): [1]}, {(1, 1): [2]}, {(1, 1): [1, 1]}]
            sage: Qloc, dloc = X.local_quiver_setting(L[0]);
            sage: Qloc.adjacency_matrix() , dloc
            ([4], (1))
            sage: Qloc, dloc = X.local_quiver_setting(L[1]);
            sage: Qloc.adjacency_matrix() , dloc
            ([1], (2))
            sage: Qloc, dloc = X.local_quiver_setting(L[2]);
            sage: Qloc.adjacency_matrix() , dloc
            (
            [1 1]
            [1 1], (1, 1)
            )

        """

        if secure:
            assert self.is_luna_type(tau)

        Q = self._Q

        # a word of caution to the future maintainer: Python dictionaries
        # iterate over the keys in the order they were inserted. This ensures
        # that the following is a well-defined adjacency matrix. Python sets
        # DO NOT have this property.
        # TODO this looks fishy; why should this choice of order matter, as long as
        # it is the same when defining dloc?
        A = matrix(
            [
                [Q.generic_ext(dp, eq) for eq in tau.keys() for n in tau[eq]]
                for dp in tau.keys()
                for m in tau[dp]
            ]
        )
        Qloc = Quiver(A)
        dloc = vector(dim for dp in tau.keys() for dim in tau[dp])

        return Qloc, dloc

    # TODO: The codimension computation requires the dimension of the nullcone. This is hard, it turns out. It can be done with the Hesselink stratification, but I wasn't willing to go thourgh Lieven's treatment of this.
    def _codimension_inverse_image_luna_stratum(self, tau):
        r"""
        Computes the codimension of :math:`\pi^{-1}(S_{tau})`
        inside `R(Q,d)` where :math:`\pi: R(Q,d)^{\theta-sst} --> M^{\theta-sst}(Q,d)`
        is the semistable quotient map.

        INPUT:

        - ``tau``: list of tuples

        OUTPUT: codimension as Int

        For ``tau = {d^1: p^1,...,d^s: p^s}``
        the codimension of :math:`\pi^{-1}(S_{tau})` is

        .. MATH::

            -\langle d,d \rangle + \sum_{k=1}^s
            (\langle d^k,d^k\rangle - l(p^k) + ||p^k||^2) -
            \mathrm{dim} N(Q_{tau}, d_{tau}),

        where for a partition :math:`p = (n_1,...,n_l)`, we define
        :math:`||p||^2 = \sum_v n_v^2`
        and :math:`N(Q_{tau}, d_{tau})` is the nullcone of the local quiver setting.
        """
        # setup shorthand
        Q, d = self._Q, self._d

        Qtau, dtau = self.local_quiver_setting(tau, secure=False)
        dimNull = Qtau.dimension_nullcone(dtau)
        return (
            -Q.euler_form(d, d)
            + sum(
                [
                    Q.euler_form(dk[0], dk[0])
                    - len(dk[1])
                    + sum([nkv**2 for nkv in dk[1]])
                    for dk in tau
                ]
            )
            - dimNull
        )

    def codimension_properly_semistable_locus(self):
        r"""
        Computes the codimension of :math:`R^{\theta-sst}(Q,d)
        \setminus R^{\theta-st}(Q,d)` inside :math:`R(Q,d)`.

        OUTPUT: codimension as Int

        The codimension of the properly semistable locus
        is the minimal codimension of the inverse image
        of the non-stable Luna strata."""

        L = self.all_luna_types()
        # This is the stable Luna type; remove it if it occurs
        dstable = [tuple([self._d, [1]])]
        L = list(filter(lambda tau: tau != dstable, L))
        return min([self._codimension_inverse_image_luna_stratum(tau) for tau in L])

    """
    (Semi-)stability
    """

    def semistable_equals_stable(self):
        r"""
        Checks whether every semistable representation is stable
        for the given stability parameter.

        Every :math:`\theta`-semistable representation is
        :math:`\theta`-stable if and only if
        there are no Luna types other than (possibly) ``{d: [1]}``.

        OUTPUT: whether every theta-semistable representation is theta-stable
        for ``self._theta``

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q, d = GeneralizedKroneckerQuiver(3), vector([3,3])
            sage: theta = vector([1,-1])
            sage: X = QuiverModuliSpace(Q, d, theta)
            sage: X.semistable_equals_stable()
            False
            sage: e = vector([2,3])
            sage: Y = QuiverModuliSpace(Q, e, theta)
            sage: Y.semistable_equals_stable()
            True

        A double framed example as in our vector fields paper::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q = Q.framed_quiver([1, 0]).coframed_quiver([0, 0, 1])
            sage: d = [1, 2, 3, 1]
            sage: theta = [1, 300, -200, -1]
            sage: Q.is_theta_coprime(d, theta)
            False
            sage: X = QuiverModuliSpace(Q, d, theta)
            sage: X.semistable_equals_stable()
            True

        """

        # setup shorthand
        Q, d, theta = self._Q, self._d, self._theta
        denom = self._denom
        d = Q._coerce_dimension_vector(d)

        # the computation of all Luna types takes so much time
        # thus we should first tests if d is theta-coprime
        if Q.is_theta_coprime(d, theta):
            return True

        # this is probably the fastest way as checking theta-coprimality is fast
        # whereas checking for existence of a semi-stable representation is a bit slower

        if not Q.has_semistable_representation(d, theta, denom=denom):
            return True
        else:
            allLunaTypes = self.all_luna_types()
            genericType = {d: [1]}
            if genericType in allLunaTypes:
                allLunaTypes.remove(genericType)
            return not allLunaTypes  # This checks if the list is empty

    """
    Ample stability
    """

    # TODO reimplement this with HN strata computation.
    def is_amply_stable(self) -> bool:
        r"""Checks if the dimension vector is amply stable for the stability parameter

        By definition, a dimension vector `d` is :math:`\theta`-amply stable if the
        codimension of the :math:`\theta`-semistable locus
        inside `R(Q,d)` is at least 2.

        OUTPUT: whether the data for the quiver moduli space is amply stable

        EXAMPLES:

        3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: QuiverModuliSpace(Q, [2, 3]).is_amply_stable()
            True
            sage: QuiverModuliSpace(Q, [2, 3], [-3, 2]).is_amply_stable()
            False

        A three-vertex example from the rigidity paper::

            sage: Q = ThreeVertexQuiver(1, 6, 1)
            sage: QuiverModuliSpace(Q, [1, 6, 6]).is_amply_stable()
            False

        """
        HNs = self.all_harder_narasimhan_types(proper=True)
        return (
            min(
                self.codimension_of_harder_narasimhan_stratum(dstar, secure=False)
                for dstar in HNs
            )
            >= 2
        )

    def is_strongly_amply_stable(self) -> bool:
        r"""Checks if the dimension vector is strongly amply stable for the stability
        parameter

        We call `d` strongly amply stable for :math:`\theta` if
        :math:`\langle e,d-e\rangle \leq -2`
        holds for all subdimension vectors :math:`e` of :math:`d` which satisfy
        :math:`\mu_{\theta}(e) >= \mu_{\theta}(d)`.

        OUTPUT: whether the data for the quiver moduli space is strongly amply stable

        EXAMPLES:

        3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: QuiverModuliSpace(Q, [2, 3]).is_strongly_amply_stable()
            True

        A 3-vertex quiver::

            sage: from quiver import *
            sage: Q = ThreeVertexQuiver(5, 1, 1)
            sage: X = QuiverModuliSpace(Q, [4, 1, 4])
            sage: X.is_amply_stable()
            True
            sage: X.is_strongly_amply_stable()
            False

        """
        # setup shorthand
        Q, d, theta, denom = (
            self._Q,
            self._d,
            self._theta,
            self._denom,
        )
        d = Q._coerce_dimension_vector(d)

        # subdimension vectors of smaller slope
        slope = Q.slope(d, theta=theta, denom=denom)
        es = filter(
            lambda e: Q.slope(e, theta=theta, denom=denom) >= slope,
            Q.all_subdimension_vectors(
                d, proper=True, nonzero=True, forget_labels=True
            ),
        )

        return all(Q.euler_form(e, d - e) <= -2 for e in es)

    """
    Methods related to Teleman quantization
    """

    def harder_narasimhan_weight(self, harder_narasimhan_type):
        r"""
        Returns the Teleman weight of a Harder-Narasimhan type
        """
        # setup shorthand
        Q, theta, denom = self._Q, self._theta, self._denom
        HN = harder_narasimhan_type

        return -sum(
            [
                # TODO can we make this cleaner-looking?
                # = unordered tuples without repetition?
                (
                    Q.slope(HN[s], theta, denom=denom)
                    - Q.slope(HN[t], theta, denom=denom)
                )
                * Q.euler_form(HN[s], HN[t])
                for s in range(len(HN) - 1)
                for t in range(s + 1, len(HN))
            ]
        )

    def all_weight_bounds(self, as_dict=False):
        r"""
        Returns the list of all weights appearing in Teleman quantization.

        For each HN type, the 1-PS lambda acts on :math:`\det(N_{S/R}|_Z)`
        with a certain weight. Teleman quantization gives a numerical condition
        involving these weights to compute cohmology on the quotient.

        INPUT:

        - ``as_dict`` -- (default: False) when True it will give a dict whose keys are
          the HN-types and whose values are the weights

        EXAMPLES:

        The 6-dimensional 3-Kronecker example::

            sage: from quiver import *
            sage: X = QuiverModuliSpace(KroneckerQuiver(3), [2, 3])
            sage: X.all_weight_bounds()
            [135, 100, 90, 15/2, 270, 100, 30]
            sage: X.all_weight_bounds(as_dict=True)
            {((1, 0), (1, 1), (0, 2)): 135,
             ((1, 0), (1, 2), (0, 1)): 100,
             ((1, 0), (1, 3)): 90,
             ((1, 1), (1, 2)): 15/2,
             ((2, 0), (0, 3)): 270,
             ((2, 1), (0, 2)): 100,
             ((2, 2), (0, 1)): 30}

        """
        # this is only relevant on the unstable locus
        HNs = self.all_harder_narasimhan_types(proper=True)

        weights = map(lambda HN: self.harder_narasimhan_weight(HN), HNs)

        if as_dict:
            return dict(zip(HNs, weights))

        return list(weights)

    def if_rigidity_inequality_holds(self) -> bool:
        r"""

        OUTPUT: whether the rigidity inequality holds on the given moduli

        If the weights of the 1-PS lambda on :math:`\det(N_{S/R}|_Z)` for each HN type
        are all strictly larger than the weights of the tensors of the universal bundles
        :math:`U_i^\vee \otimes U_j`,
        then the resulting moduli space is infinitesimally rigid.

        EXAMPLES:

            sage: from quiver import *
            sage: X = QuiverModuliSpace(KroneckerQuiver(3), [2, 3])
            sage: X.if_rigidity_inequality_holds()
            True
            sage: X = QuiverModuliSpace(ThreeVertexQuiver(1, 6, 1), [1, 6, 6])
            sage: X.if_rigidity_inequality_holds()
            False

        """
        # setup shorthand
        Q, theta, denom = self._Q, self._theta, self._denom

        weights = self.all_weight_bounds()

        # we compute the maximum weight of the tensors of the universal bundles
        # this is only relevant on the unstable locus
        HNs = self.all_harder_narasimhan_types(proper=True)

        tensor_weights = list(
            map(
                lambda HN: Q.slope(HN[0], theta, denom=denom)
                - Q.slope(HN[-1], theta, denom=denom),
                HNs,
            )
        )

        return all(weights[i] > tensor_weights[i] for i in range(len(HNs)))

    """
    Tautological relations
    """

    def _all_forbidden_subdimension_vectors(self):
        r"""Returns the list of all forbidden subdimension vectors

        These are the dimension vectors `d'` of d for which

        - :math:`\mu_{\theta}(d') > \mu_{\theta}(d)` (in the semistable case)
        - or for which :math:`\mu_{\theta}(d') >= \mu_{\theta}(d)` (in the stable case).

        OUTPUT: list of forbidden subdimension vectors vectors

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, [3, 3], [1, -1], condition="semistable")
            sage: X._all_forbidden_subdimension_vectors()
            [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)]
            sage: X = QuiverModuliSpace(Q, [3, 3], [1, -1], condition="stable")
            sage: X._all_forbidden_subdimension_vectors()
            [(1, 0), (1, 1), (2, 0), (2, 1), (2, 2), (3, 0), (3, 1), (3, 2)]

        """
        # setup shorthand
        Q, d, theta, condition = self._Q, self._d, self._theta, self._condition

        es = Q.all_subdimension_vectors(d, proper=True, nonzero=True)

        # TODO need for denominator?
        slope = Q.slope(d, theta)
        if condition == "semistable":
            return list(filter(lambda e: Q.slope(e, theta) > slope, es))
        elif condition == "stable":
            return list(filter(lambda e: Q.slope(e, theta) >= slope, es))

    # TODO make it private?
    def all_minimal_forbidden_subdimension_vectors(self):
        r"""Returns the list of all `minimal` forbidden subdimension vectors

        Minimality is with respect to the partial order `e << d` which means
        :math:`e_i \leq d_i` for every source `i`, :math:`e_j \geq d_j`
        for every sink `j`, and :math:`e_k = d_k` for every vertex which is neither
        a source nor a sink. See also :meth:`Quiver.division_order`.

        OUTPUT: list of minimal forbidden dimension vectors

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, [3, 3], [1, -1], condition="semistable")
            sage: X.all_minimal_forbidden_subdimension_vectors()
            [(1, 0), (2, 1), (3, 2)]
            sage: Y = QuiverModuliSpace(Q, [3, 3], [1, -1], condition="stable")
            sage: Y.all_minimal_forbidden_subdimension_vectors()
            [(1, 1), (2, 2)]

        """
        # setup shorthand
        Q = self._Q

        forbidden = self._all_forbidden_subdimension_vectors()

        def is_minimal(e):
            return not any(
                Q.division_order(f, e)
                for f in list(filter(lambda f: f != e, forbidden))
            )

        return list(filter(is_minimal, forbidden))

    def __tautological_presentation(
        self, inRoots=False, chernClasses=None, chernRoots=None
    ):
        r"""Returns the tautological relations in Chern classes
        (if ``inRoots == False``) or in Chern roots.

        INPUT:

        - ``inRoots`` -- Bool
        - ``chernClasses`` -- list of Strings
        - ``chernRoots`` -- list of Strings

        OUTPUT: dict

        # TODO
        Explanation ...

        Notation for explanations:
        G = G_d = prod_{i in Q_0} GL_{d_i}
        T = maximal torus of diagonal matrices
        PG = G/G_m
        PT = T/G_m maximal torus of PT
        W = Weyl group of T in G = Weyl group of PT in PG
          = prod_{i in Q_0} S_{d_i}
        R = bigoplus_{a in Q_1} Hom(k^{d_{s(a)}},k^{d_{t(a)}})
        R^{sst}, R^{st} semi-stable/stable locus

        EXAMPLES:

        # TODO
        """
        # setup shorthand
        Q, d = (
            self._Q,
            self._d,
        )

        if chernClasses is None:
            chernClasses = [
                "x%s_%s" % (i, r)
                for i in range(Q.number_of_vertices())
                for r in range(1, d[i] + 1)
            ]
        if chernRoots is None:
            chernRoots = [
                "t%s_%s" % (i, r)
                for i in range(Q.number_of_vertices())
                for r in range(1, d[i] + 1)
            ]

        R = PolynomialRing(QQ, chernRoots)

        def generator(R, i, r):
            r"""Returns generator(R, i, r) = t{i+1}_{r+1}."""
            return R.gen(r + sum([d[j] for j in range(i)]))

        r"""Generators of the tautological ideal regarded upstairs, i.e. in A*([R/T]).
        For a forbidden subdimension vector e of d, the forbidden polynomial in Chern
        roots is given by :math:`\prod_{a: i \to j} \prod_{r=1}^{e_i}
        \prod_{s=e_j+1}^{d_j} (tj_s - ti_r) =
        \prod_{i,j} \prod_{r=1}^{e_i} \prod_{s=e_j+1}^{d_j} (tj_s - ti_r)^{a_{ij}}."""
        forbiddenPolynomials = [
            prod(
                [
                    prod(
                        [
                            (generator(R, j, s) - generator(R, i, r))
                            ** Q.adjacency_matrix()[i, j]
                            for r in range(e[i])
                            for s in range(e[j], d[j])
                        ]
                    )
                    for i in range(Q.number_of_vertices())
                    for j in range(Q.number_of_vertices())
                ]
            )
            for e in self.all_minimal_forbidden_subdimension_vectors()
        ]

        if inRoots:
            return {
                "ParentRing": R,
                "Generators": lambda i, r: generator(R, i, r),
                "Relations": forbiddenPolynomials,
            }
        else:
            """delta is the discriminant"""
            delta = prod(
                [
                    prod(
                        [
                            generator(R, i, l) - generator(R, i, k)
                            for k in range(d[i])
                            for l in range(k + 1, d[i])
                        ]
                    )
                    for i in range(Q.number_of_vertices())
                ]
            )

            """longest is the longest Weyl group element
            when regarding W as a subgroup of S_{sum d_i}"""
            longest = []
            r = 0
            for i in range(Q.number_of_vertices()):
                longest = longest + list(reversed(range(r + 1, r + d[i] + 1)))
                r += d[i]
            W = Permutations(bruhat_smaller=longest)

            def antisymmetrization(f):
                """The antisymmetrization of f is the symmetrization
                divided by the discriminant."""

                # I don't want to define W and delta here but globally because then we need to
                # compute it just once. That's probably a bit faster.
                def permute(f, w):
                    return f.subs({R.gen(i): R.gen(w[i] - 1) for i in range(R.ngens())})

                return sum(w.sign() * permute(f, w) for w in W) // delta

            """Schubert basis of A^*([R/T]) over A^*([R/G])"""
            X = SchubertPolynomialRing(ZZ)
            supp = list(filter(lambda i: d[i] > 0, range(Q.number_of_vertices())))

            def B(i):
                return [X(p).expand() for p in Permutations(d[i])]

            Bprime = [
                [
                    f.parent().hom(
                        [generator(R, i, r) for r in range(f.parent().ngens())], R
                    )(f)
                    for f in B(i)
                ]
                for i in supp
            ]

            # TODO is this not something already implemented? if not, explain what it does!
            def product_lists(L):
                n = len(L)
                assert n > 0
                if n == 1:
                    return L[0]
                else:
                    P = product_lists([L[i] for i in range(n - 1)])
                    return [p * l for p in P for l in L[n - 1]]

            schubert = product_lists(Bprime)

            """Define A = A*([R/G])."""
            degrees = []
            for i in range(Q.number_of_vertices()):
                degrees = degrees + list(range(1, d[i] + 1))
            A = PolynomialRing(QQ, chernClasses, order=TermOrder("wdegrevlex", degrees))

            E = SymmetricFunctions(ZZ).e()
            """The Chern classes of U_i on [R/G] are the elementary symmetric functions
            in the Chern roots ti_1,...,ti_{d_i}."""
            elementarySymmetric = []
            for i in range(Q.number_of_vertices()):
                elementarySymmetric = elementarySymmetric + [
                    E([k]).expand(
                        d[i],
                        alphabet=[generator(R, i, r) for r in range(d[i])],
                    )
                    for k in range(1, d[i] + 1)
                ]
            """Map xi_r to the r-th elementary symmetric function
            in ti_1,...,ti_{d_i}."""
            inclusion = A.hom(elementarySymmetric, R)

            """Tautological relations in Chern classes."""
            tautological = [
                antisymmetrization(b * f)
                for b in schubert
                for f in forbiddenPolynomials
            ]
            tautological = [inclusion.inverse_image(g) for g in tautological]

            return {
                "ParentRing": A,
                "Generators": lambda i, r: generator(A, i, r),
                "Relations": tautological,
            }

    def tautological_relations(self, inRoots=False, chernClasses=None, chernRoots=None):
        r"""
        Returns the tautological relations in
        Chern classes (if inRoots == False) or in Chern roots.

        INPUT:

        - ``inRoots`` -- Bool
        - ``chernClasses`` -- list of Strings
        - ``chernRoots`` -- list of Strings

        OUTPUT: list
        """

        taut = self.__tautological_presentation(
            inRoots=inRoots, chernClasses=chernClasses, chernRoots=chernRoots
        )
        return taut["Relations"]

    def dimension(self) -> int:
        raise NotImplementedError()

    def is_smooth(self) -> bool:
        raise NotImplementedError()

    def chow_ring(self):
        raise NotImplementedError()


class QuiverModuliSpace(QuiverModuli):
    def __init__(self, Q, d, theta=None, denom=sum, condition="semistable"):
        r"""Constructor for a quiver moduli space

        This is the quiver moduli space as a variety.

        INPUT:

        - ``Q`` -- quiver

        - ``d`` --- dimension vector

        - ``theta`` -- stability parameter (default: canonical stability parameter)

        - ``denom`` -- denominator for slope stability (default: ``sum``), needs to be
          effective on the simple roots

        - ``condition`` -- whether to include all semistables, or only stables
          (default: "semistable")

        EXAMPLES:

        An example::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: QuiverModuliSpace(Q, [2, 3])
            moduli space of semistable representations, with
            Q = 3-Kronecker quiver
            d = [2, 3]
            theta = (9, -6)

        """
        QuiverModuli.__init__(
            self,
            Q,
            d,
            theta=theta,
            denom=denom,
            condition=condition,
        )

    def _repr_(self):
        r"""
        Give a shorthand string presentation for the quiver moduli space
        """
        if self.get_custom_name():
            return self.get_custom_name()

        output = "moduli space of {} representations, with".format(self._condition)
        output += "\nQ = {}\nd = {}\ntheta = {}".format(
            self._Q.repr(), self._d, self._theta
        )

        return output

    def repr(self):
        return self._repr_()

    def dimension(self):
        r"""
        Computes the dimension of the moduli space :math:`M^{\theta-(s)st}(Q,d)`.

        This involves several cases:

        - If there are :math:`\theta`-stable representations then
          :math:`\mathrm{dim} M^{\theta-sst}(Q,d) =
          M^{\theta-st}(Q,d) = 1 - \langle d,d\rangle`;
        - if there are no :math:`\theta`-stable representations then
          :math:`\mathrm{dim} M^{\theta-st}(Q,d) = -\infty` by convention,
          and we define :math:`\mathrm{dim} M^{\theta-sst} =
          \mathrm{max}_{\tau} \{\mathrm{dim} S_{\tau}\}`,
          the maximum of the dimension of all Luna strata.

        EXAMPLES

        The A2-quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(1)
            sage: X = QuiverModuliSpace(Q, [1, 1], condition="stable")
            sage: X.dimension()
            0
            sage: X = QuiverModuliSpace(Q, [1, 1], condition="semistable")
            sage: X.dimension()
            0
            sage: X = QuiverModuliSpace(Q, [2, 2], condition="stable")
            sage: X.dimension()
            -Infinity
            sage: X = QuiverModuliSpace(Q, [2, 2], condition="semistable")
            sage: X.dimension()
            0

        The Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(2)
            sage: X = QuiverModuliSpace(Q, [1, 1], [1, -1], condition="stable")
            sage: X.dimension()
            1
            sage: X = QuiverModuliSpace(Q, [1, 1], [1, -1], condition="semistable")
            sage: X.dimension()
            1
            sage: X = QuiverModuliSpace(Q, [2, 2], [1, -1], condition="stable")
            sage: X.dimension()
            -Infinity
            sage: X = QuiverModuliSpace(Q, [2, 2], [1, -1], condition="semistable")
            sage: X.dimension()
            2

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, [2, 3], condition="semistable")
            sage: X.dimension()
            6
            sage: X = QuiverModuliSpace(Q, [3, 3],condition="semistable")
            sage: X.dimension()
            10
            sage: X = QuiverModuliSpace(Q, [1, 3],condition="stable")
            sage: X.dimension()
            0
            sage: X = QuiverModuliSpace(Q, [1, 4],condition="stable")
            sage: X.dimension()
            -Infinity
            sage: X = QuiverModuliSpace(Q, [1, 4],condition="semistable")
            sage: X.dimension()
            -Infinity

        """
        # setup shorthand
        Q, d, theta = (
            self._Q,
            self._d,
            self._theta,
        )

        # if there are stable representations then both the stable and
        # the semi-stable moduli space have dimension `1-<d,d>`
        if Q.has_stable_representation(d, theta):
            return 1 - Q.euler_form(d, d)

        # stable locus is empty
        if self._condition == "stable":
            return -Infinity

        # we care about the semistable locus
        if Q.has_semistable_representation(d, theta):
            # in this case the dimension is given by
            # the maximum of the dimensions of the Luna strata
            return max(
                self.dimension_of_luna_stratum(tau) for tau in self.all_luna_types()
            )

        # semistable locus is also empty
        return -Infinity

    def poincare_polynomial(self):
        r"""
        Returns the Poincare polynomial of the moduli space.

        OUTPUT: polynomial in one variable
        # TODO allow a user-supplied ring?

        The Poincare polynomial is defined as

        .. MATH::
            P_X(q) = \sum_{i \geq 0} (-1)^i \mathrm{dim} H^i(X;\mathbb{C}) q^{i/2}.

        For a quiver moduli space whose dimension vector is
        :math:`\theta`-coprime, the odd cohomology vanishes
        and this is a Polynomial in :math:`q`.
        We use Cor. 6.9 in Reineke's Harder--Narasimhan paper to compute it.

        EXAMPLES:

        Some Kronecker quivers::

            sage: from quiver import *
            sage: Q = KroneckerQuiver()
            sage: X = QuiverModuliSpace(Q, [1, 1])
            sage: X.poincare_polynomial()
            q + 1
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, [2, 3])
            sage: X.poincare_polynomial()
            q^6 + q^5 + 3*q^4 + 3*q^3 + 3*q^2 + q + 1
            sage: Q = SubspaceQuiver(5)
            sage: X = QuiverModuliSpace(Q, [1, 1, 1, 1, 1, 2])
            sage: X.poincare_polynomial()
            q^2 + 5*q + 1

        """
        # setup shorthand
        Q, d, theta = self._Q, self._d, self._theta

        assert Q.is_theta_coprime(d, theta)

        k = FunctionField(QQ, "L")
        K = FunctionField(QQ, "q")
        q = K.gen(0)
        f = k.hom(q, K)

        X = QuiverModuliStack(Q, d, theta, condition="semistable")

        return (1 - q) * f(X.motive())

    def betti_numbers(self):
        r"""
        Returns the Betti numbers of the moduli space.

        OUTPUT: List of Ints

        EXAMPLES:

        Some Kronecker quivers::

            sage: from quiver import *
            sage: Q, d, theta = KroneckerQuiver(), vector([1,1]), vector([1,-1])
            sage: X = QuiverModuliSpace(Q, d, theta, condition="semistable")
            sage: X.poincare_polynomial()
            q + 1
            sage: X.betti_numbers()
            [1, 0, 1]
            sage: Q, d = GeneralizedKroneckerQuiver(3), vector([2,3])
            sage: theta = vector([1,-1])
            sage: X = QuiverModuliSpace(Q, d, theta, condition="semistable")
            sage: X.betti_numbers()
            [1, 0, 1, 0, 3, 0, 3, 0, 3, 0, 1, 0, 1]

        """

        assert self._Q.is_theta_coprime(self._d, self._theta)
        N = self.dimension()

        K = FunctionField(QQ, "q")
        L = FunctionField(QQ, "v")
        v = L.gen(0)
        ext = K.hom(v**2, L)
        # p is the prime place of the DVR associated with v
        p = v.zeros()[0]

        f = ext(self.poincare_polynomial())
        betti = [f.evaluate(p)]
        for i in range(2 * N):
            f = (f - f.evaluate(p)) / v
            betti = betti + [f.evaluate(p)]

        return betti

    def is_smooth(self) -> bool:
        # if theta-coprime then you can shortcut everything
        # if theta != 0 reduce to theta = 0 using https://mathscinet.ams.org/mathscinet-getitem?mr=1972892 (Adriaenssens--Le Bruyn)
        # if theta = 0, then use https://mathscinet.ams.org/mathscinet-getitem?mr=1929191 (Bocklandt)
        if (self._condition == "stable") or (
            self._Q.is_theta_coprime(self._d, self._theta)
        ):
            return True
        else:
            raise NotImplementedError()

    def is_projective(self) -> bool:
        raise NotImplementedError()

    def picard_rank(self):
        """Computes the Picard rank of the moduli space for known cases."""
        # setup shorthand
        Q, d, theta = self._Q, self._d, self._theta
        # TODO requires smooth and projective?

        if Q.is_theta_coprime(d, theta) and Q.is_amply_stable(Q, d, theta):
            return Q.number_of_vertices() - 1
        else:
            # TODO if smooth: compute the Betti numbers and return b_2
            raise NotImplementedError()

    def index(self):
        """Computes the index of the moduli space for known cases,
        i.e., the largest integer dividing the canonical divisor in Pic."""
        # TODO this should really be a check for theta belonging to the canonical chamber, rather than being equal to the canonical stability.
        # setup shorthand
        Q, d, theta = self._Q, self._d, self._theta
        if (
            # TODO at the very least check for multiple of canonical
            theta == Q.canonical_stability_parameter(d)
            and Q.is_theta_coprime(d, theta)
            and Q.is_amply_stable(Q, d, theta)
        ):
            # TODO what if theta is rescaled?
            return gcd(Q._to_vector(theta))
        else:
            raise NotImplementedError()

    def mukai_inequality_holds(self):
        # TODO ample stability for the canonical stability parameter should be an attribute of the object, so that it is only computed once. Verbatim for many other attributes.
        # setup shorthand
        Q, d = self._Q, self._d

        return 1 - Q.tits_form(d) >= self.picard_rank() * (self.index() - 1)

    def chow_ring(self, chi=None, chernClasses=None):
        r"""
        Returns the Chow ring of the moduli space.

        INPUT:

        - ``chi`` -- vector of Ints
        - ``chernClasses`` -- list of Strings

        OUTPUT: ring

        EXAMPLES:

        The Kronecker quiver::

            sage: from quiver import *
            sage: Q, d, theta = KroneckerQuiver(), vector([1,1]), vector([1,-1])
            sage: X = QuiverModuliSpace(Q, d, theta, condition="semistable")
            sage: chi = vector([1,0])
            sage: A = X.chow_ring(chi=chi)
            sage: I = A.defining_ideal()
            sage: [I.normal_basis(i) for i in range(X.dimension()+1)]
            [[1], [x1_1]]


        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q, d = GeneralizedKroneckerQuiver(3), vector([2,3])
            sage: theta = vector([3,-2])
            sage: X = QuiverModuliSpace(Q, d, theta, condition="semistable")
            sage: chi = vector([-1,1])
            sage: A = X.chow_ring(chi=chi)
            sage: I = A.defining_ideal()
            sage: [I.normal_basis(i) for i in range(X.dimension()+1)]
            [[1],
            [x1_1],
            [x0_2, x1_1^2, x1_2],
            [x1_1^3, x1_1*x1_2, x1_3],
            [x1_1^2*x1_2, x1_2^2, x1_1*x1_3],
            [x1_2*x1_3],
            [x1_3^2]]

        The 5-subspaces quiver::

            sage: from quiver import *
            sage: Q, d = SubspaceQuiver(5), vector([1,1,1,1,1,2])
            sage: theta = vector([2,2,2,2,2,-5])
            sage: X = QuiverModuliSpace(Q, d, theta, condition="semistable")
            sage: chi = vector([-1,-1,-1,-1,-1,3])
            sage: A = X.chow_ring(chi=chi)
            sage: I = A.defining_ideal()
            sage: [I.normal_basis(i) for i in range(X.dimension()+1)]
            [[1], [x1_1, x2_1, x3_1, x4_1, x5_1], [x5_2]]

        """
        Q, d, theta = self._Q, self._d, self._theta
        n = Q.number_of_vertices()

        # This implementation only works if d is theta-coprime
        # which implies that d is indivisible.
        assert Q.is_theta_coprime(d, theta)

        if chernClasses is None:
            chernClasses = [
                "x%s_%s" % (i, r) for i in range(n) for r in range(1, d[i] + 1)
            ]

        taut = self._QuiverModuli__tautological_presentation(
            inRoots=False, chernClasses=chernClasses
        )
        A, generator, rels = taut["ParentRing"], taut["Generators"], taut["Relations"]

        if chi is None:
            [g, m] = extended_gcd(d.list())
            chi = vector(m)

        # TODO assert that chi has integer entries.
        """Make sure that chi has weight one, i.e.,
        provides a retraction for X*(PG) --> X*(G)."""
        assert chi * d == 1

        I = A.ideal(rels) + A.ideal(sum([chi[i] * generator(i, 0) for i in range(n)]))

        return QuotientRing(A, I, names=chernClasses)

    def chern_class_line_bundle(self, eta, chernClasses=None):
        r"""
        Returns the first Chern class of the line bundle
        :math:`L(\eta) = \bigotimes_{i \in Q_0} \det(U_i)^{-\eta_i}`,
        where :math:`\eta` is a character of :math:`PG_d`."""

        A = self.chow_ring(chi=None, chernClasses=chernClasses)
        n = self._Q.number_of_vertices()
        d = self._d

        return -sum([eta[i] * A.gen(sum([d[j] for j in range(i)])) for i in range(n)])

    def chern_character_line_bundle(self, eta, chernClasses=None):
        r"""
        Computes the Chern character of L(eta).

        The Chern character of a line bundle `L` with first Chern class `x`
        is given by :math:`e^x = 1 + x + \frac{x^2}{2} + \frac{x^3}{6} + \dots`
        """

        N = self.dimension()
        x = self.chern_class_line_bundle(eta, chernClasses=chernClasses)
        return sum([x**i / factorial(i) for i in range(N + 1)])

    def total_chern_class_universal(self, i, chi, chernClasses=None):
        """Gives the total Chern class of the universal bundle U_i(chi)."""

        A = self.chow_ring(chi, chernClasses=chernClasses)
        d = self._d

        return 1 + sum(
            [A.gen(r + sum([d[j] for j in range(i - 1)])) for r in range(d[i - 1])]
        )

    def point_class(self, chi=None, chernClasses=None):
        r"""
        Returns the point class as an expression in Chern classes of the
        :math:`U_i` (``chi``).

        The point class is given as the homogeneous component of degree
        :math:`\mathrm{dim} X` of the expression

        .. MATH::

            \prod_{a \in Q_1} c(U_{t(a)})^{d_{s(a)}} / (\prod_{i \in Q_0} c(U_i)^{d_i})

        EXAMPLES

        :math:`\mathbb{P}^7` as a quiver moduli space
        of a generalized Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(8)
            sage: d = vector([1,1])
            sage: theta = vector([1,-1])
            sage: X = QuiverModuliSpace(Q,d,theta,condition="semistable")
            sage: chi = vector([1,0])
            sage: X.point_class(chi,chernClasses=['o','h'])
            h^7

        Our favorite 6-fold::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: d = vector([2,3])
            sage: theta = vector([3,-2])
            sage: X = QuiverModuliSpace(Q,d,theta,condition="semistable")
            sage: chi = vector([-1,1])
            sage: X.point_class(chi,chernClasses=['x1','x2','y1','y2','y3'])
            y3^2

        A moduli space of the 5-subspaces quiver;
        it agrees with the blow-up of :math:`\mathbb{P}^2` in 4 points
        in general position::

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
        pi = A.cover()  # The quotient map
        sect = A.lifting_map()  # A choice of a section of pi

        if chi is None:
            [g, m] = extended_gcd(d.list())
            chi = vector(m)

        my_numerator = prod(
            [
                self.total_chern_class_universal(j + 1, chi, chernClasses=chernClasses)
                ** (d * a.column(j))
                for j in range(n)
            ]
        )
        my_denom = prod(
            [
                self.total_chern_class_universal(i + 1, chi, chernClasses=chernClasses)
                ** d[i]
                for i in range(n)
            ]
        )

        quotient = my_numerator / my_denom

        return pi(sect(quotient).homogeneous_components()[N])

    def degree(self, eta=None, chernClasses=None):
        r"""Computes the degree of the ample line bundle given by eta."""
        # TODO: Need check for ampleness first

        if eta is None:
            eta = self._Q.canonical_stability_parameter(self._d)

        N = self.dimension()
        c = self.chern_class_line_bundle(eta, chernClasses=chernClasses)
        p = self.point_class(chernClasses=chernClasses)

        return c**N / p

    def todd_class(self):
        r"""
        The Todd class of `X` is the Todd class of the tangent bundle.

        For quiver moduli it computes as

        .. MATH::

            td(X) =
            (\prod_{a:i \to j \in Q_1} \prod_{p=1}^{d_j} \prod_{q=1}^{d_i} Q(t_{j,q} -
            t_{i,p}))/(prod_{i \in Q_0} \prod_{p,q=1}^{d_i} Q(t_{i,q} - t_{i,p}))
        """

        def todd_generating_series(t, n):
            r"""
            We call the series :math:`Q(t) = t/(1-e^{-t})` the Todd generating series.

            The function computes the terms of this series up to degree n."""
            B = [bernoulli(i) for i in range(n + 1)]
            return sum([(-1) ^ i * B[i] / factorial(i) * t ^ i for i in range(n + 1)])

        def truncate(f, n):
            r"""
            Takes an element in a graded ring and discards
            all homogeneous components of degree > n"""
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
    def __init__(self, Q, d, theta=None, denom=sum, condition="semistable"):
        r"""Constructor for a quiver moduli stack

        This is the quiver moduli space as a stack.

        INPUT:

        - ``Q`` -- quiver

        - ``d`` --- dimension vector

        - ``theta`` -- stability parameter (default: canonical stability parameter)

        - ``denom`` -- denominator for slope stability (default: ``sum``), needs to be
          effective on the simple roots

        - ``condition`` -- whether to include all semistables, or only stables
          (default: "semistable")

        EXAMPLES:

        An example::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: X = QuiverModuliStack(Q, [2, 3])

        """
        QuiverModuli.__init__(self, Q, d, theta=theta, denom=denom, condition=condition)

    def __repr__(self):
        return (
            "A "
            + self._condition
            + " quiver moduli stack with:\n"
            + "Q = "
            + str(self._Q)
            + "\n"
            + "d = "
            + str(self._d)
            + "\n"
            + "theta = "
            + str(self._theta)
        )

    def dimension(self):
        r"""
        Computes the dimension of the moduli stack :math:`[R^{(s)st}/G]`.

        .. MATH::

            dim [R^{(s)st}/G] = dim R^{(s)st} - dim G

        The dimension turns out to be :math:`-\langle d,d\rangle`
        if the (semi-)stable locus is non-empty"""
        # setup shorthand
        Q, d, theta = self._Q, self._d, self._theta

        if self._condition == "stable" and Q.has_stable_representation(d, theta):
            return -Q.euler_form(d, d)
        # TODO is this one correct? we need to check for existence of a stable I think?
        if self._condition == "semistable" and Q.has_semistable_representation(
            d, theta
        ):
            return -Q.euler_form(d, d)
        else:
            return -Infinity

    def is_smooth(self) -> bool:
        # TODO think about the empty case, should it be smooth?
        return True

    def motive(self):
        r"""Gives an expression for the motive of the semistable moduli stack
        in an appropriate localization of K_0(Var)

        # TODO more explanation

        EXAMPLES:

        Loop quivers::

            sage: from quiver import *
            sage: Q, d, theta = LoopQuiver(0), vector([2]), vector([0])
            sage: X = QuiverModuliStack(Q, d, theta, condition="semistable")
            sage: X.motive()
            1/(L^4 - L^3 - L^2 + L)
            sage: Q, d, theta = LoopQuiver(1), vector([2]), vector([0])
            sage: X = QuiverModuliStack(Q, d, theta, condition="semistable")
            sage: X.motive()
            L^3/(L^3 - L^2 - L + 1)

        The 3-Kronecker quiver::

            sage: Q, d = GeneralizedKroneckerQuiver(3), vector([2,3])
            sage: theta = vector([3,-2])
            sage: X = QuiverModuliStack(Q, d, theta, condition="semistable")
            sage: X.motive()
            (-L^6 - L^5 - 3*L^4 - 3*L^3 - 3*L^2 - L - 1)/(L - 1)

        """

        # Only for semistable.
        # For stable, we don't know what the motive is. It's not pure in general.
        assert self._condition == "semistable"
        # TODO well: if we have stable == semistable then we can also compute it!

        # setup shorthand
        Q, d, theta = self._Q, self._d, self._theta

        # TODO allow some other ring?
        K = FunctionField(QQ, "L")
        L = K.gen(0)

        # TODO coercion needs to be checked here
        if theta == Q.zero_vector():
            num = L ** (-Q.tits_form(d))
            den = prod(
                [
                    prod([(1 - L ** (-nu)) for nu in range(1, d[i] + 1)])
                    for i in range(Q.number_of_vertices())
                ]
            )
            return num / den
        else:
            # TODO use proper=True, nonzero=True, or maybe not?
            # in any case, the next 6 lines are an atrocity
            I = Q.all_subdimension_vectors(d)
            I = list(filter(lambda e: e != Q.zero_vector() and e != d, I))
            I = list(filter(lambda e: Q.slope(e, theta) > Q.slope(d, theta), I))
            I = I + [Q.zero_vector(), d]
            I = [Q._coerce_dimension_vector(e) for e in I]
            # TODO I believe max(d) on a dict should give the wrong result
            I.sort(key=(lambda e: Q._deglex_key(e, b=max(d) + 1)))

            # Now define a matrix T of size NxN whose entry at position (i,j) is
            # L^<e-f,e>*mot(f-e) if e = I[i] is a subdimension vector of f = I[j]
            # and 0 otherwise
            # TODO it's bad to have a function motive inside a motive method
            def motive(e):
                return QuiverModuliStack(
                    Q, e, Q.zero_vector(), condition="semistable"
                ).motive()

            N = len(I)
            T = matrix(K, N)
            for i in range(N):
                for j in range(i, N):
                    e, f = I[i], I[j]
                    if Q.is_subdimension_vector(e, f):
                        T[i, j] = L ** (Q.euler_form(e - f, e)) * motive(f - e)

            # Solve system of linear equations T*x = e_N
            # and extract entry 0 of the solution x.
            y = vector([0 for i in range(N)])
            y[N - 1] = 1
            x = T.solve_right(y)

            return x[0]

    def chow_ring(self, chernClasses=None):
        r"""Returns the Chow ring of the quotient stack.

        INPUT:

        - ``chernClasses``: list of Strings

        OUTPUT: ring
        """
        # setup shorthand
        Q, d = self._Q, self._d

        # TODO there is very similar code earlier
        if chernClasses is None:
            chernClasses = [
                "x%s_%s" % (i, r)
                for i in range(Q.number_of_vertices())
                for r in range(1, d[i] + 1)
            ]

        taut = self._QuiverModuli__tautological_presentation(
            inRoots=False, chernClasses=chernClasses
        )
        A, _, rels = taut["ParentRing"], taut["Generators"], taut["Relations"]

        return QuotientRing(A, A.ideal(rels), names=chernClasses)


class SmoothModel:
    """How about this: instead of a separate class SmoothModel,
    we could define a method framed_moduli_space(self,n)
    inside the class QuiverModuliSpace which returns another quiver moduli space.
    After all, it is itself a quiver moduli space."""

    def __init__(self):
        pass

    def betti_numbers(self):
        raise NotImplementedError()


"""Auxiliary methods:"""


def extended_gcd(x):
    """Computes the gcd and the Bezout coefficients of a list of integers."""
    # This exists for two integers but there seems to be
    # no implementation for more than one.
    # That's astonishing.

    n = len(x)
    if n == 1:
        return [x, [1]]
    if n == 2:
        (g, a, b) = xgcd(x[0], x[1])
        return [g, [a, b]]
    if n > 2:
        (g, a, b) = xgcd(x[0], x[1])
        y = [g] + [x[i] for i in range(2, n)]
        [d, c] = extended_gcd(y)
        m = [c[0] * a, c[0] * b] + [c[i] for i in range(1, n - 1)]
        return [d, m]
