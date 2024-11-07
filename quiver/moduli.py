from itertools import combinations_with_replacement, product

from sage.arith.misc import bernoulli, factorial, gcd, xgcd
from sage.combinat.partition import Partitions
from sage.combinat.permutation import Permutations
from sage.combinat.schubert_polynomial import SchubertPolynomialRing
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.tuple import UnorderedTuples
from sage.matrix.constructor import matrix
from sage.misc.misc_c import prod
from sage.modules.free_module_element import vector, zero_vector
from sage.rings.function_field.constructor import FunctionField
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.term_order import TermOrder
from sage.rings.quotient_ring import QuotientRing
from sage.rings.rational_field import QQ
from sage.structure.element import Element

from . import Quiver

"""Defines how permutations are multiplied."""
Permutations.options(mult="r2l")


class QuiverModuli(Element):
    def __init__(self, Q, d, theta=None, denom=sum, condition="semistable"):
        r"""
        Constructor for an abstract quiver moduli space.

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
            sage: X = QuiverModuli(Q, (2, 3))
            sage: X
            abstract moduli of semistable representations, with
            - Q = 3-Kronecker quiver
            - d = (2, 3)
            - θ = (9, -6)

        It has functionality common to both varieties and stacks, i.e., when it really
        concerns something involving the representation variety::

            sage: X.all_harder_narasimhan_types()
            (((1, 0), (1, 1), (0, 2)),
             ((1, 0), (1, 2), (0, 1)),
             ((1, 0), (1, 3)),
             ((1, 1), (1, 2)),
             ((2, 0), (0, 3)),
             ((2, 1), (0, 2)),
             ((2, 2), (0, 1)),
             ((2, 3),))

        But things like dimension depend on whether we consider it as a variety or as
        a stack, and thus these are not implemented::

            sage: X.dimension()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if theta is None:
            theta = Q.canonical_stability_parameter(d)

        assert Q._is_dimension_vector(d), "``d`` needs to be a dimension vector"
        assert Q._is_vector(theta), "`theta` needs to be a stability parameter"
        assert condition in [
            "semistable",
            "stable",
        ], "condition needs to be (semi)stable"
        assert all(
            denom(Q._coerce_dimension_vector(Q.simple_root(i))) > 0
            for i in Q.vertices()
        ), "denominator needs to be effective"
        assert (
            Q._coerce_dimension_vector(d) * Q._coerce_vector(theta) == 0
        ), "for the moment we require that `theta(d) == 0`"

        self._Q = Q
        self._d = d
        self._theta = theta
        self._denom = denom
        self._condition = condition

    def __repr_helper(self, description):
        r"""
        Standard format for shorthand string presentation.

        EXAMPLES:

        A Kronecker moduli space with non-standard description::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: X = QuiverModuli(Q, (2, 3))
            sage: print(X._QuiverModuli__repr_helper("Kronecker moduli space"))
            Kronecker moduli space of semistable representations, with
            - Q = 3-Kronecker quiver
            - d = (2, 3)
            - θ = (9, -6)
        """
        output = "{} of {} representations, with".format(description, self._condition)
        output += "\n- Q = {}\n- d = {}\n- θ = {}".format(
            self._Q.repr(), self._d, self._theta
        )

        return output

    def _repr_(self):
        r"""
        Give a shorthand string presentation for an abstract quiver moduli space.

        EXAMPLES:

        A Kronecker moduli space::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: QuiverModuli(Q, (2, 3))
            abstract moduli of semistable representations, with
            - Q = 3-Kronecker quiver
            - d = (2, 3)
            - θ = (9, -6)
        """
        if self.get_custom_name():
            return self.get_custom_name()

        return self.__repr_helper("abstract moduli")

    def repr(self):
        r"""
        Give a shorthand string presentation for an abstract quiver moduli space.

        EXAMPLES:

        A Kronecker moduli space::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: QuiverModuli(Q, (2, 3))
            abstract moduli of semistable representations, with
            - Q = 3-Kronecker quiver
            - d = (2, 3)
            - θ = (9, -6)
        """
        return self._repr_()

    def to_space(self):
        r"""
        Make the abstract quiver moduli a variety.

        This is an explicit way of casting an abstract :class:`QuiverModuli`
        to a :class:`QuiverModuliSpace`.

        EXAMPLES:

        From an abstract quiver moduli to a space::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: X = QuiverModuli(Q, (2, 3))
            sage: X.to_space()
            moduli space of semistable representations, with
            - Q = 3-Kronecker quiver
            - d = (2, 3)
            - θ = (9, -6)

        From a stack to a space::

            sage: X = QuiverModuliStack(Q, (2, 3))
            sage: X.to_space()
            moduli space of semistable representations, with
            - Q = 3-Kronecker quiver
            - d = (2, 3)
            - θ = (9, -6)
        """
        return QuiverModuliSpace(
            self._Q, self._d, self._theta, self._denom, self._condition
        )

    def to_stack(self):
        r"""
        Make the abstract quiver moduli a stack.

        This is an explicit way of casting an abstract :class:`QuiverModuli`
        to a :class:`QuiverModuliStack`.

        EXAMPLES:

        From an abstract quiver moduli to a stack::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: X = QuiverModuli(Q, (2, 3))
            sage: X.to_stack()
            moduli stack of semistable representations, with
            - Q = 3-Kronecker quiver
            - d = (2, 3)
            - θ = (9, -6)

        From a space to a stack::

            sage: X = QuiverModuliSpace(Q, (2, 3))
            sage: X.to_stack()
            moduli stack of semistable representations, with
            - Q = 3-Kronecker quiver
            - d = (2, 3)
            - θ = (9, -6)
        """
        return QuiverModuliStack(
            self._Q, self._d, self._theta, self._denom, self._condition
        )

    def quiver(self):
        r"""
        Returns the quiver of the moduli space.

        OUTPUT: the underlying quiver as an instance of the :class:`Quiver` class

        EXAMPLES:

        The quiver of a Kronecker moduli space::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuli(Q, (2, 3))
            sage: Q == X.quiver()
            True
        """
        return self._Q

    def dimension_vector(self):
        r"""
        Returns the dimension vector of the moduli space.

        OUTPUT: the dimension vector

        EXAMPLES:

        The dimension vector of a Kronecker moduli space::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuli(Q, (2, 3))
            sage: X.dimension_vector()
            (2, 3)

        The dimension vector is stored in the same format as it was given::

            sage: Q = Quiver.from_string("foo---bar", forget_labels=False)
            sage: X = QuiverModuli(Q, {"foo": 2, "bar": 3})
            sage: X.dimension_vector()
            {'bar': 3, 'foo': 2}
        """
        return self._d

    def stability_parameter(self):
        r"""
        Returns the stability parameter of the moduli space.

        OUTPUT: the stability parameter

        EXAMPLES:

        The stability parameter of a Kronecker moduli space::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuli(Q, (2, 3), (3, -2))
            sage: X.stability_parameter()
            (3, -2)

        The stability parameter is stored in the same format as it was given::

            sage: Q = Quiver.from_string("foo---bar", forget_labels=False)
            sage: d, theta = {"foo": 2, "bar": 3}, {"foo": 3, "bar": -2}
            sage: X = QuiverModuliSpace(Q, d, theta);
            sage: X.stability_parameter()
            {'bar': -2, 'foo': 3}
        """
        return self._theta

    def denominator(self):
        r"""
        Returns the denominator of the slope function :math:`\mu_{\theta}`.

        OUTPUT: the denominator as a function

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (2, 3))
            sage: X.denominator()
            <built-in function sum>
        """
        return self._denom

    def is_nonempty(self) -> bool:
        r"""
        Checks if the moduli space is nonempty.

        OUTPUT: whether there exist stable/semistable representations, according
        to the condition

        EXAMPLES:

        The 3-Kronecker quiver for `d = (2, 3)` has stable representations::

            sage: from quiver import *
            sage: Q, d = GeneralizedKroneckerQuiver(3), (2, 3)
            sage: X = QuiverModuliSpace(Q, d, condition="stable"); X.is_nonempty()
            True

        The Jordan quiver does not have stable representations, but it has semistable
        ones::

            sage: Q = JordanQuiver()
            sage: X = QuiverModuliSpace(Q, [3], condition="stable")
            sage: X.is_nonempty()
            False
            sage: X = QuiverModuliSpace(Q, [3], condition="semistable")
            sage: X.is_nonempty()
            True
        """
        if self._condition == "stable":
            return self._Q.has_stable_representation(self._d, self._theta)
        if self._condition == "semistable":
            return self._Q.has_semistable_representation(self._d, self._theta)

    def is_theta_coprime(self) -> bool:
        r"""
        Checks whether the combination of `d` and `theta` is coprime.

        This just calls :meth:`Quiver.is_theta_coprime` for the data defining the
        moduli space.

        EXAMPLES:

        A coprime example::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: QuiverModuliSpace(Q, (2, 3)).is_theta_coprime()
            True

        And a non-example::

            sage: QuiverModuliSpace(Q, (3, 3)).is_theta_coprime()
            False
        """
        return self._Q.is_theta_coprime(self._d, self._theta)

    """
    Harder--Narasimhan stratification
    """

    def all_harder_narasimhan_types(self, proper=False, sorted=False):
        r"""
        Returns the list of all Harder--Narasimhan types.

        A Harder--Narasimhan (HN) type of `d` with respect to :math:`\theta`
        is a sequence :math:`{\bf d}^* = ({\bf d}^1,...,{\bf d}^s)` of dimension vectors
        such that

        - :math:`{\bf d}^1 + ... + {\bf d}^s = {\bf d}`
        - :math:`\mu_{\theta}({\bf d}^1) > ... > \mu_{\theta}({\bf d}^s)`
        - Every :math:`{\bf d}^k` is :math:`\theta`-semistable.

        INPUT:

        - ``proper`` -- (default: False) whether to exclude the HN-type corresponding
          to the stable locus

        - ``sorted`` -- (default: False) whether to sort the HN-types according to the
          given slope

        OUTPUT: list of tuples of dimension vectors encoding HN-types

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (2, 3))
            sage: X.all_harder_narasimhan_types()
            (((1, 0), (1, 1), (0, 2)),
             ((1, 0), (1, 2), (0, 1)),
             ((1, 0), (1, 3)),
             ((1, 1), (1, 2)),
             ((2, 0), (0, 3)),
             ((2, 1), (0, 2)),
             ((2, 2), (0, 1)),
             ((2, 3),))
            sage: X.all_harder_narasimhan_types(proper=True)
            (((1, 0), (1, 1), (0, 2)),
             ((1, 0), (1, 2), (0, 1)),
             ((1, 0), (1, 3)),
             ((1, 1), (1, 2)),
             ((2, 0), (0, 3)),
             ((2, 1), (0, 2)),
             ((2, 2), (0, 1)))
            sage: d = (2, 3)
            sage: theta = -Q.canonical_stability_parameter(d)
            sage: Y = QuiverModuliSpace(Q, d, theta)
            sage: Y.all_harder_narasimhan_types()
            (((0, 3), (2, 0)),)

        A 3-vertex quiver::

            sage: from quiver import *
            sage: Q = ThreeVertexQuiver(2, 3, 4)
            sage: Z = QuiverModuliSpace(Q, (2, 3, 2))
            sage: Z.all_harder_narasimhan_types()
            (((0, 1, 0), (1, 2, 1), (1, 0, 1)),
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
             ((2, 3, 2),))
        """
        d = self._Q._coerce_dimension_vector(self._d)
        theta = self._Q._coerce_vector(self._theta)

        all_types = self._Q._all_harder_narasimhan_types(
            d, theta, denom=self._denom, sorted=sorted
        )

        if proper:
            all_types = tuple(ds for ds in all_types if ds != (d,))

        return all_types

    def is_harder_narasimhan_type(self, dstar) -> bool:
        r"""
        Checks if ``dstar`` is a Harder--Narasimhan type.

        A Harder--Narasimhan (HN) type of :math`{\bf d}` with respect to :math:`\theta`
        is a sequence :math:`{\bf d}^* = ({\bf d}^1,...,{\bf d}^s)` of dimension vectors
        such that

        - :math:`{\bf d}^1 + ... + {\bf d}^s = {\bf d}`
        - :math:`\mu_{\theta}({\bf d}^1) > ... > \mu_{\theta}({\bf d}^s)`
        - Every :math:`{\bf d}^k` is :math:`\theta`-semistable.

        INPUT:

        - ``dstar`` -- list of dimension vectors

        OUTPUT: whether ``dstar`` is a valid HN type for the moduli space

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (2, 3))
            sage: HNs = X.all_harder_narasimhan_types()
            sage: all(X.is_harder_narasimhan_type(dstar) for dstar in HNs)
            True
            sage: dstar = [(1, 0), (1, 0), (0, 3)]
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

        assert all(
            Q._is_dimension_vector(di) for di in dstar
        ), "elements of ``dstar`` need to be dimension vectors"

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

        # third condition
        if not all(
            Q.has_semistable_representation(di, theta, denom=denom) for di in dstar
        ):
            return False

        return True

    def codimension_of_harder_narasimhan_stratum(self, dstar, secure=False):
        r"""
        Computes the codimension of the HN stratum of ``dstar``
        inside the representation variety :math:`R(Q,{\bf d})`.

        INPUT:

        - ``dstar`` -- the HN type as a list of dimension vectors

        - ``secure`` -- whether to first check it is an HN-type (default: False)

        OUTPUT: codimension as an integer

        By default, the method does not check if ``dstar`` is a valid HN type.
        This can be enabled by passing ``secure=True``.

        The codimension of the HN stratum of
        :math:`{\bf d}^* = ({\bf d}^1,...,{\bf d}^s)` is given by

        .. MATH::

            - \sum_{k < l} \langle {\bf d}^k,{\bf d}^l\rangle

        INPUT:

        - ``dstar`` -- list of dimension vectors

        - ``secure`` -- whether to check ``dstar`` is an HN-type (default: False)

        OUTPUT: codimension of the HN-stratum

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (2, 3))
            sage: HNs = X.all_harder_narasimhan_types()
            sage: [X.codimension_of_harder_narasimhan_stratum(dstar) for dstar in HNs]
            [12, 9, 8, 3, 18, 10, 4, 0]

        """
        Q = self._Q

        assert all(
            Q._is_dimension_vector(di) for di in dstar
        ), "elements of ``dstar`` need to be dimension vectors"

        if secure:
            assert self.is_harder_narasimhan_type(dstar), "``dstar`` must be HN-type"

        return -sum(
            Q.euler_form(dstar[k], dstar[l])
            for k in range(len(dstar) - 1)
            for l in range(k + 1, len(dstar))
        )

    def codimension_unstable_locus(self):
        r"""
        Computes codimension of the unstable locus inside the representation variety.

        This is the minimum of the codimensions of the proper Harder--Narasimhan strata
        of the representation variety.

        OUTPUT: codimension of the unstable locus

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (2, 3))
            sage: X.codimension_unstable_locus()
            3

        A 3-vertex quiver::

            sage: Q = ThreeVertexQuiver(1, 6, 1)
            sage: X = QuiverModuliSpace(Q, (1, 6, 6))
            sage: X.codimension_unstable_locus()
            1

        The :math:`\mathrm{A}_2` quiver is of finite type::

            sage: Q = GeneralizedKroneckerQuiver(1)
            sage: X = QuiverModuliSpace(Q, (2, 3))
            sage: X.codimension_unstable_locus()
            0

        """
        HNs = self.all_harder_narasimhan_types(proper=True)

        # note that while the HN types and strata depend on the denominator
        # the maximum of their codimensions does not
        return min(
            self.codimension_of_harder_narasimhan_stratum(dstar, secure=False)
            for dstar in HNs
        )

    """
    Luna
    """

    def all_luna_types(self, exclude_stable=False):
        r"""
        Returns the unordered list of all Luna types of ``d`` for ``theta``.

        INPUT:

        - ``exclude_stable`` -- whether to exclude the stable Luna type ``{d: [1]}``
          (default: False)

        OUTPUT: the list of all the Luna types as dictionaries.

        The Luna stratification of the representation variety concerns the étale-local
        structure of the moduli space of semistable quiver representations. It is
        studied in MR1972892_, and for more details one is referred there.

        .. _MR1972892: https://mathscinet.ams.org/mathscinet/relay-station?mr=1972892

        A Luna type of :math:`{\bf d}` for :math:`\theta` is an unordered sequence
        :math:`(({\bf d}^1,m_1),...,({\bf d}^s,m_s))` of pairs of dimension vectors
        :math:`{\bf d}^k` and positive integers :math:`m_k` such that

        - :math:`m_1{\bf d}^1 + ... + m_s{\bf d}^s = {\bf d}`,
        - :math:`\mu_{\theta}({\bf d}^k) = \mu_{\theta}({\bf d})`, and
        - all the :math:`{\bf d}^k` admit a :math:`\theta`-stable representation.

        Note that a pair :math:`({\bf d}^i, m_i)`
        can appear multiple times in a Luna type, and the same dimension vector
        :math:`{\bf d}^i` can appear coupled with different integers.

        IMPLEMENTATION:

        Here a Luna type is a dictionary
        ``{d^1: p^1, ... d^s: p^s}``
        whose keys are dimension vectors :math:`{\bf d}^k` and values are non-empty
        lists of positive integers
        ``p^k = [p_{k, 1}, ..., p_{k, t_k}]``.

        The corresponding Luna type is then the unordered sequence of tuples

        .. MATH::

            ({\bf d}^1, p_{1, 1}), \dots, ({\bf d}^1, p_{1, t_1}), \dots
            ({\bf d}^s, p_{s, 1}), \dots, ({\bf d}^s, p_{s, t_s}),

        such that

        .. MATH::

            (p_{1, 1} + \dots + p_{1, t_1}) \cdot {\bf d}^1 + \dots +
            (p_{s, 1} + \dots + p_{s, t_s}) \cdot {\bf d}^s = {\bf d}.

        ALGORITHM:

        The way we compute the Luna types of a quiver moduli space is taken from
        Section 4 in MR2511752_.

        .. _MR2511752: https://mathscinet.ams.org/mathscinet/relay-station?mr=2511752

        EXAMPLES:

        The Kronecker quiver::

            sage: from quiver import *
            sage: Q = KroneckerQuiver()
            sage: X = QuiverModuliSpace(Q, (3, 3), (1, -1))
            sage: X.all_luna_types()
            [{(1, 1): [3]}, {(1, 1): [2, 1]}, {(1, 1): [1, 1, 1]}]

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (3, 3), (1, -1))
            sage: X.all_luna_types()
            [{(3, 3): [1]},
             {(1, 1): [1], (2, 2): [1]},
             {(1, 1): [3]},
             {(1, 1): [2, 1]},
             {(1, 1): [1, 1, 1]}]

        The zero vector::

            sage: from quiver import *
            sage: Q = KroneckerQuiver()
            sage: X = QuiverModuliSpace(Q, (0, 0), (1, -1))
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

        # we will build all possible Luna types from the bottom up
        Ls = []

        # start with all subdimension vectors
        ds = Q.all_subdimension_vectors(d, nonzero=True, forget_labels=True)
        # look for subdimension vectors with the same slope as ``d``
        # and which admit a stable representation:
        # this encodes the second and third condition in the definition
        same_slope = filter(
            lambda e: Q.slope(e, theta, denom=denom) == Q.slope(d, theta, denom=denom)
            and Q.has_stable_representation(e, theta, denom=denom),
            ds,
        )
        same_slope = list(same_slope)

        # bounds how long a Luna type can be
        bound = (sum(d) / min(sum(e) for e in same_slope)).ceil()

        for i in range(1, bound + 1):
            for tau in combinations_with_replacement(same_slope, i):
                # first condition is not satisfied
                if not sum(tau) == d:
                    continue

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

                # we add all possible Luna types we can build to our list
                Ls += [
                    dict(zip(partial.keys(), values))
                    for values in product(*partial.values())
                ]

        stable = {d: [1]}
        if exclude_stable and stable in Ls:
            Ls.remove(stable)

        return Ls

    def is_luna_type(self, tau) -> bool:
        r"""
        Checks if ``tau`` is a Luna type.

        INPUT:

        - ``tau`` -- Luna type encoded by a dictionary of multiplicities indexed by
          dimension vectors

        OUTPUT: whether ``tau`` is a Luna type

        For a description of Luna types, see :meth:`all_luna_types`.

        EXAMPLES:

        The Kronecker quiver::

            sage: from quiver import *
            sage: Q = KroneckerQuiver()
            sage: X = QuiverModuliSpace(Q, (3, 3), (1, -1))
            sage: Ls = X.all_luna_types()
            sage: all(X.is_luna_type(tau) for tau in Ls)
            True

        The 3-Kronecker quiver with zero vector::

            sage: from quiver import *
            sage: Q = KroneckerQuiver()
            sage: X = QuiverModuliSpace(Q, (0, 0), (1, -1))
            sage: X.is_luna_type({Q.zero_vector(): [1]})
            True
        """
        Q, d, theta, denom = (
            self._Q,
            self._d,
            self._theta,
            self._denom,
        )

        d = Q._coerce_dimension_vector(d)

        assert all(
            Q._is_dimension_vector(dk) for dk in tau.keys()
        ), "elements of ``tau`` need to be dimension vectors"

        if d == Q.zero_vector():
            # Q.zero_vector() can't be hashed a priori
            z = Q._coerce_vector(Q.zero_vector())
            return tau == {z: [1]}

        # we check the 3 conditions in that order
        return d == sum(sum(m) * dk for (dk, m) in tau.items()) and all(
            Q.slope(dk, theta, denom=denom) == Q.slope(d, theta, denom=denom)
            and Q.has_semistable_representation(dk, theta, denom=denom)
            for dk in tau.keys()
        )

    def dimension_of_luna_stratum(self, tau, secure=True):
        r"""
        Computes the dimension of the Luna stratum :math:`S_\tau`.

        INPUT:

        - ``tau`` -- Luna type encoded by a dictionary of multiplicities indexed by
          dimension vectors

        - ``secure`` -- whether to first check it is a Luna type (default: False)

        OUTPUT: dimension of the corresponding Luna stratum

        The dimension of the Luna stratum of ``tau = {d^1: p^1,...,d^s: p^s}`` is

        .. MATH::

            \sum_k l(p^k)(1 - \langle {\bf d}^k,{\bf d}^k\rangle),

        where for a partition :math:`p = (n_1,...,n_l)`,
        the length `l(p)` is `l`, i.e., the number of summands.

        EXAMPLES:

        The Kronecker quiver::

            sage: from quiver import *
            sage: Q = KroneckerQuiver()
            sage: X = QuiverModuliSpace(Q, (2, 2), (1, -1))
            sage: Ls = X.all_luna_types(); Ls
            [{(1, 1): [2]}, {(1, 1): [1, 1]}]
            sage: [X.dimension_of_luna_stratum(tau) for tau in Ls]
            [1, 2]
        """
        if secure:
            assert self.is_luna_type(tau), "``tau`` needs to be a Luna type"

        return sum(len(tau[di]) * (1 - self._Q.euler_form(di, di)) for di in tau.keys())

    def local_quiver_setting(self, tau, secure=True):
        r"""
        Returns the local quiver and dimension vector for the given Luna type.

        The local quiver describes the singularities of a moduli space,
        and is introduced and studied in studied in MR1972892_.

        .. _MR1972892: https://mathscinet.ams.org/mathscinet/relay-station?mr=1972892

        INPUT:

        - ``tau`` -- Luna type encoded by a dictionary of multiplicities indexed by
          dimension vectors

        - ``secure`` -- whether to first check it is a Luna type (default: False)

        OUTPUT: tuple consisting of a Quiver object and a dimension vector

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (2, 2), (1, -1))
            sage: Ls = X.all_luna_types(); Ls
            [{(2, 2): [1]}, {(1, 1): [2]}, {(1, 1): [1, 1]}]
            sage: Qloc, dloc = X.local_quiver_setting(Ls[0]);
            sage: Qloc.adjacency_matrix(), dloc
            ([4], (1))
            sage: Qloc, dloc = X.local_quiver_setting(Ls[1]);
            sage: Qloc.adjacency_matrix(), dloc
            ([1], (2))
            sage: Qloc, dloc = X.local_quiver_setting(Ls[2]);
            sage: Qloc.adjacency_matrix(), dloc
            (
            [1 1]
            [1 1], (1, 1)
            )
        """
        if secure:
            assert self.is_luna_type(tau), "``tau`` needs to be a Luna type"

        Q = self._Q

        # we use the order of vertices provided by ``tau.keys()`` for Qloc and dloc
        A = matrix(
            [
                [Q.generic_ext(dp, eq) for eq in tau.keys() for n in tau[eq]]
                for dp in tau.keys()
                for m in tau[dp]
            ]
        )
        Qloc = Quiver(A)
        dloc = vector(m for dp in tau.keys() for m in tau[dp])

        return Qloc, dloc

    def _codimension_inverse_image_luna_stratum(self, tau):
        r"""
        Computes the codimension of the preimage of the Luna stratum

        This is the codimension of :math:`\pi^{-1}(S_{tau})`
        inside:math:`R(Q,{\bf d})` where

        .. MATH::

            \pi\colon R(Q,{\bf d})^{\theta{\rm-sst}}\to M^{\theta{\rm-sst}}(Q,{\bf d})

        is the semistable quotient map.

        INPUT:

        - ``tau`` -- Luna type encoded by a dictionary of multiplicities indexed by
          dimension vectors

        OUTPUT: the codimension of the inverse image of the Luna stratum

        For ``tau = {d^1: p^1,...,d^s: p^s}``
        the codimension of :math:`\pi^{-1}(S_{tau})` is

        .. MATH::

            -\langle {\bf d},{\bf d} \rangle + \sum_{k=1}^s
            (\langle {\bf d}^k,{\bf d}^k\rangle - l(p^k) + ||p^k||^2) -
            \dim N(Q_{tau}, {\mathbf{d}_{tau}),

        where for a partition :math:`p = (n_1,...,n_l)`, we define
        :math:`||p||^2 = \sum_v n_v^2`
        and :math:`N(Q_{\tau}, d_{\tau})` is the nullcone of the local quiver setting.

        This is currently not working properly because we cannot compute the dimension
        of the nullcone::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (3, 3))
            sage: Ls = X.all_luna_types()
            sage: X._codimension_inverse_image_luna_stratum(Ls[0])
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        # setup shorthand
        Q, d = self._Q, self._d

        Qtau, dtau = self.local_quiver_setting(tau, secure=False)
        return (
            -Q.euler_form(d, d)
            + sum(
                [
                    Q.euler_form(dk, dk) - sum(m) + sum([nkv**2 for nkv in m])
                    for (dk, m) in tau.items()
                ]
            )
            - Qtau.dimension_nullcone(dtau)
        )

    def codimension_properly_semistable_locus(self):
        r"""
        Computes the codimension of :math:`R^{\theta\rm-sst}(Q,{\bf d})
        \setminus R^{\theta\rm-st}(Q,{\bf d})` inside :math:`R(Q,{\bf d})`.

        OUTPUT: codimension of the properly semistable locus

        The codimension of the properly semistable locus
        is the minimal codimension of the inverse image
        of the non-stable Luna strata.

        EXAMPLES:

        If the semistable locus is the stable locus the codimension is -Infinity::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: QuiverModuliSpace(Q, (2, 3)).codimension_properly_semistable_locus()
            -Infinity

        This is currently not working properly because we cannot compute the dimension
        of the nullcone::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: QuiverModuliSpace(Q, (3, 3)).codimension_properly_semistable_locus()
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        Ls = self.all_luna_types(exclude_stable=True)
        codimensions = [self._codimension_inverse_image_luna_stratum(tau) for tau in Ls]

        return min(codimensions, default=-Infinity)

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

        OUTPUT: whether every theta-semistable representation is :math:`\theta`-stable

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (3, 3))
            sage: X.semistable_equals_stable()
            False
            sage: Y = QuiverModuliSpace(Q, (2, 3))
            sage: Y.semistable_equals_stable()
            True

        A double framed example as in arXiv.2311.17004_::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: Q = Q.framed_quiver((1, 0)).coframed_quiver((0, 0, 1))
            sage: d = (1, 2, 3, 1)
            sage: theta = (1, 300, -200, -1)
            sage: X = QuiverModuliSpace(Q, d, theta)
            sage: X.is_theta_coprime()
            False
            sage: X.semistable_equals_stable()
            True

        .. _arXiv.2311.17004: https://doi.org/10.48550/arXiv.2311.17004
        """
        # setup shorthand
        Q, d, theta, denom = self._Q, self._d, self._theta, self._denom
        d = Q._coerce_dimension_vector(d)

        # the computation of all Luna types takes so much time
        # thus we should first tests if ``d`` is ``theta``-coprime
        if self.is_theta_coprime():
            return True

        # this is probably the fastest way as checking theta-coprimality is fast
        # whereas checking for existence of a semi-stable representation
        # is a bit slower
        if not Q.has_semistable_representation(d, theta, denom=denom):
            return True
        else:
            Ls = self.all_luna_types(exclude_stable=True)
            return not Ls  # this checks if the list is empty

    """
    Ample stability
    """

    def is_amply_stable(self) -> bool:
        r"""Checks if the dimension vector is amply stable for the stability parameter

        By definition, a dimension vector :math`{\bf d}` is :math:`\theta`-amply stable
        if the codimension of the :math:`\theta`-semistable locus
        inside:math:`R(Q,{\bf d})` is at least 2.

        OUTPUT: whether the data for the quiver moduli space is amply stable

        EXAMPLES:

        3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: QuiverModuliSpace(Q, (2, 3)).is_amply_stable()
            True
            sage: QuiverModuliSpace(Q, (2, 3), [-3, 2]).is_amply_stable()
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

        We call :math:`{\bf d}` strongly amply stable for :math:`\theta` if
        :math:`\langle{\bf e},{\bf d}-{\bf e}\rangle \leq -2`
        holds for all subdimension vectors :math:`{\bf e}` of :math:`{\bf d}` for which
        :math:`\mu_{\theta}({\bf e})\geq\mu_{\theta}({\bf d})`.

        OUTPUT: whether the data for the quiver moduli space is strongly amply stable

        EXAMPLES:

        3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: QuiverModuliSpace(Q, (2, 3)).is_strongly_amply_stable()
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

        INPUT:

        - ``harder_narasimhan_type`` -- list of vectors of Ints

        OUTPUT: weight as a fraction

        The weight of a Harder-Narasimhan type :math:`{\bf d}^*`
        is the weight of the associated 1-PS :math:`\lambda` acting on
        :math:`\det(N_{S/R})^{\vee}|_Z`, where `S` is the
        corresponding Harder--Narasimhan stratum.

        .. SEEALSO:: :meth:`all_weight_bounds`, :meth:`if_rigidity_inequality_holds`

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (2, 3))
            sage: HN = X.all_harder_narasimhan_types(proper=True)
            sage: {dstar: X.harder_narasimhan_weight(dstar) for dstar in HN}
            {((1, 0), (1, 1), (0, 2)): 135,
             ((1, 0), (1, 2), (0, 1)): 100,
             ((1, 0), (1, 3)): 90,
             ((1, 1), (1, 2)): 15/2,
             ((2, 0), (0, 3)): 270,
             ((2, 1), (0, 2)): 100,
             ((2, 2), (0, 1)): 30}
        """
        # setup shorthand
        Q, theta, denom = self._Q, self._theta, self._denom
        HN = harder_narasimhan_type

        return -sum(
            [
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

        For each HN type, the 1-PS lambda acts on :math:`\det(N_{S/R}^{\vee}|_Z)`
        with a certain weight. Teleman quantization gives a numerical condition
        involving these weights to compute cohomology on the quotient.

        INPUT:

        - ``as_dict`` -- (default: False) when True it will give a dict whose keys are
          the HN-types and whose values are the weights

        EXAMPLES:

        The 6-dimensional 3-Kronecker example::

            sage: from quiver import *
            sage: X = QuiverModuliSpace(KroneckerQuiver(3), (2, 3))
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

        weights = map(lambda dstar: self.harder_narasimhan_weight(dstar), HNs)

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

        Kronecker moduli satisfy the rigidity inequality::

            sage: from quiver import *
            sage: X = QuiverModuliSpace(KroneckerQuiver(3), (2, 3))
            sage: X.if_rigidity_inequality_holds()
            True

        The following 3-vertex example does not (however, it is rigid by other means)::

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
                lambda dstar: Q.slope(dstar[0], theta, denom=denom)
                - Q.slope(dstar[-1], theta, denom=denom),
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
        Q, d, theta, denom, condition = (
            self._Q,
            self._d,
            self._theta,
            self._denom,
            self._condition,
        )

        es = Q.all_subdimension_vectors(d, proper=True, nonzero=True)

        slope = Q.slope(d, theta, denom=denom)

        if condition == "semistable":
            return list(filter(lambda e: Q.slope(e, theta, denom=denom) > slope, es))
        elif condition == "stable":
            return list(filter(lambda e: Q.slope(e, theta, denom=denom) >= slope, es))

    def _all_minimal_forbidden_subdimension_vectors(self):
        r"""Returns the list of all `minimal` forbidden subdimension vectors

        Minimality is with respect to the partial order :math`e\ll d` which means
        :math:`e_i \leq d_i` for every source `i`, :math:`e_j \geq d_j`
        for every sink `j`, and :math:`e_k = d_k` for every vertex which is neither
        a source nor a sink. See also :meth:`Quiver.division_order`.

        OUTPUT: list of minimal forbidden dimension vectors

        EXAMPLES:

        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (3, 3))
            sage: X._all_minimal_forbidden_subdimension_vectors()
            [(1, 0), (2, 1), (3, 2)]
            sage: Y = QuiverModuliSpace(Q, (3, 3), condition="stable")
            sage: Y._all_minimal_forbidden_subdimension_vectors()
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

    def __generator(self, R, i, r):
        r"""
        Returns the appropriate generator of R.

        This is a repeatedly used helper function in dealing with Chow rings.

        EXAMPLES:

        We index some generators of the polynomial ring defining the Chow ring of the
        Kronecker 6-fold::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (2, 3))
            sage: chi = (-1, 1)
            sage: A = X.chow_ring(chi=chi, classes=["x1", "x2", "y1", "y2", "y3"])
            sage: R = A.defining_ideal().ring()
            sage: R
            Multivariate Polynomial Ring in x1, x2, y1, y2, y3 over Rational Field
            sage: X._QuiverModuli__generator(R, 0, 0)
            x1
            sage: X._QuiverModuli__generator(R, 1, 2)
            y3
        """
        # setup shorthand
        Q, d = self._Q, self._d
        d = Q._coerce_dimension_vector(d)

        return R.gen(r + sum(d[j] for j in range(i)))

    def tautological_ideal(self, use_roots=False, classes=None, roots=None):
        r"""
        Returns the tautological presentation of the Chow ring of the moduli space.

        INPUT:

        - ``use_roots`` -- (default: False) whether to return the relations in Chern
          roots

        - ``classes`` -- (default: None) optional list of strings to name the Chern
          classes

        - ``roots`` -- (default: None) optional list of strings to name the Chern roots

        OUTPUT: ideal of a polynomial ring

        EXAMPLES:

        The tautological ideal for our favourite 6-fold has 9 non-zero generators::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (2, 3))
            sage: len(X.tautological_ideal().gens())
            9

        """

        return self.__tautological_ideal_helper(
            use_roots=use_roots, classes=classes, roots=roots
        )["ideal"]

    def __tautological_ideal_helper(self, use_roots=False, classes=None, roots=None):
        r"""
        Helper function for the tautological ideal.

        INPUT:

        - ``use_roots`` -- (default: False) whether to return the relations in Chern
          roots
        - ``classes`` -- (default: None) optional list of strings to name the Chern
          classes
        - ``roots`` -- (default: None) optional list of strings to name the Chern roots

        OUTPUT: dictionary with keys "ideal", "inclusion" and "ambient_ring"

        EXAMPLES:

        The tautological ideal for our favourite 6-fold has 9 non-zero generators::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (2, 3))
            sage: X._QuiverModuli__tautological_ideal_helper()["ambient_ring"]
            Multivariate Polynomial Ring in t0_1, t0_2, t1_1, t1_2, t1_3
            over Rational Field

        .. SEEALSO:: :meth:`tautological_ideal`
        """
        # setup shorthand
        Q, d = self._Q, self._d
        d = Q._coerce_dimension_vector(d)

        if classes is None:
            classes = [
                "x{}_{}".format(i, r)
                for i in range(Q.number_of_vertices())
                for r in range(1, d[i] + 1)
            ]

        if roots is None:
            roots = [
                "t{}_{}".format(i, r)
                for i in range(Q.number_of_vertices())
                for r in range(1, d[i] + 1)
            ]

        assert len(classes) == sum(d), "number of classes must be number of generators"
        assert len(roots) == sum(d), "number of roots must be number of generators"

        R = PolynomialRing(QQ, roots)

        r"""Generators of the tautological ideal regarded upstairs, i.e. in A*([R/T]).
        For a forbidden subdimension vector e of d, the forbidden polynomial in Chern
        roots is given by :math:`\prod_{a: i \to j} \prod_{r=1}^{e_i}
        \prod_{s=e_j+1}^{d_j} (tj_s - ti_r) =
        \prod_{i,j} \prod_{r=1}^{e_i} \prod_{s=e_j+1}^{d_j} (tj_s - ti_r)^{a_{ij}}."""
        forbidden_polynomials = [
            prod(
                prod(
                    (self.__generator(R, j, s) - self.__generator(R, i, r))
                    ** Q.adjacency_matrix()[i, j]
                    for r in range(e[i])
                    for s in range(e[j], d[j])
                )
                for i in range(Q.number_of_vertices())
                for j in range(Q.number_of_vertices())
            )
            for e in self._all_minimal_forbidden_subdimension_vectors()
        ]

        # the user wants to have the ideal in `R`
        if use_roots:
            return {"ideal": R.ideal(forbidden_polynomials), "ambient_ring": R}

        # delta is the discriminant: precomputed for antisymmetrization(f)
        delta = prod(
            prod(
                self.__generator(R, i, l) - self.__generator(R, i, k)
                for k in range(d[i])
                for l in range(k + 1, d[i])
            )
            for i in range(Q.number_of_vertices())
        )

        # longest is the longest Weyl group element
        # regarding W as a subgroup of S_{sum d_i}
        longest = []
        r = 0
        for i in range(Q.number_of_vertices()):
            longest = longest + list(reversed(range(r + 1, r + d[i] + 1)))
            r += d[i]

        # Weyl group: precomputed for antisymmetrization(f)
        W = Permutations(bruhat_smaller=longest)

        def antisymmetrization(f):
            r"""The antisymmetrization of a polynomial `f` is the symmetrization
            divided by the discriminant."""

            def permute(f, w):
                return f.subs({R.gen(i): R.gen(w[i] - 1) for i in range(R.ngens())})

            return sum(w.sign() * permute(f, w) for w in W) // delta

        # we construct the Schubert basis of CH^*([R/T]) over CH^*([R/G])
        X = SchubertPolynomialRing(ZZ)

        def B(i):
            return [X(p).expand() for p in Permutations(d[i])]

        Bprime = [
            [
                f.parent().hom(
                    [self.__generator(R, i, r) for r in range(f.parent().ngens())],
                    R,
                )(f)
                for f in B(i)
            ]
            for i in Q.support(d)
        ]

        # take a list of lists of elements of a ring and multiply them recursively
        # multiplying each element of the first list with the products of the remaining
        # there might be a more Pythonic way of doing this, but it'll do for now
        def product_lists(L):
            n = len(L)
            assert n > 0
            if n == 1:
                return L[0]
            else:
                P = product_lists([L[i] for i in range(n - 1)])
                return [p * l for p in P for l in L[n - 1]]

        schubert = product_lists(Bprime)

        # define A = CH^*([R/G])
        degrees = []
        for i in range(Q.number_of_vertices()):
            degrees = degrees + list(range(1, d[i] + 1))
        A = PolynomialRing(QQ, classes, order=TermOrder("wdegrevlex", degrees))

        E = SymmetricFunctions(ZZ).e()
        """The Chern classes of U_i on [R/G] are the elementary symmetric functions
        in the Chern roots ti_1,...,ti_{d_i}."""
        elementarySymmetric = []
        for i in range(Q.number_of_vertices()):
            elementarySymmetric = elementarySymmetric + [
                E([k]).expand(
                    d[i],
                    alphabet=[self.__generator(R, i, r) for r in range(d[i])],
                )
                for k in range(1, d[i] + 1)
            ]
        """Map xi_r to the r-th elementary symmetric function
        in ti_1,...,ti_{d_i}."""
        inclusion = A.hom(elementarySymmetric, R)

        """Tautological relations in Chern classes."""
        tautological = [
            antisymmetrization(b * f) for b in schubert for f in forbidden_polynomials
        ]
        tautological = [inclusion.inverse_image(g) for g in tautological]

        # get rid of zeroes
        tautological = [f for f in tautological if f]

        return {
            "ideal": A.ideal(tautological),
            "ambient_ring": R,
            "inclusion": inclusion,
        }

    def dimension(self) -> int:
        r"""
        Returns the dimension of the moduli space.

        Abstract method, see the concrete implementations for details.

        .. SEEALSO::

            - :meth:`QuiverModuliSpace.dimension`
            - :meth:`QuiverModuliStack.dimension`

        EXAMPLES:

        This is not implemented as it is ambiguous: it depends on whether we consider
        it as a variety or as a stack::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: QuiverModuli(Q, (2, 3)).dimension()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()

    def is_smooth(self) -> bool:
        r"""
        Checks if the moduli space is smooth.

        Abstract method, see the concrete implementations for details.

        .. SEEALSO::

            - :meth:`QuiverModuliSpace.is_smooth`
            - :meth:`QuiverModuliStack.is_smooth`

        This is not implemented as it is ambiguous: it depends on whether we consider
        it as a variety or as a stack::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: QuiverModuli(Q, (2, 3)).is_smooth()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()

    def chow_ring(self):
        r"""
        Returns the Chow ring of the moduli space.

        Abstract method, see the concrete implementations for details.

        .. SEEALSO::

            - :meth:`QuiverModuliSpace.chow_ring`
            - :meth:`QuiverModuliStack.chow_ring`

        EXAMPLES:

        This is not implemented as it is ambiguous: it depends on whether we consider
        it as a variety or as a stack::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: QuiverModuli(Q, (2, 3)).is_smooth()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
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
            sage: QuiverModuliSpace(Q, (2, 3))
            moduli space of semistable representations, with
            - Q = 3-Kronecker quiver
            - d = (2, 3)
            - θ = (9, -6)

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

        EXAMPLES:

        A Kronecker moduli space::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: QuiverModuliSpace(Q, (2, 3))
            moduli space of semistable representations, with
            - Q = 3-Kronecker quiver
            - d = (2, 3)
            - θ = (9, -6)

        """
        if self.get_custom_name():
            return self.get_custom_name()

        return super()._QuiverModuli__repr_helper("moduli space")

    def repr(self):
        r"""
        Give a shorthand string presentation for the quiver moduli space

        EXAMPLES:

        A Kronecker moduli space::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: QuiverModuliSpace(Q, (2, 3))
            moduli space of semistable representations, with
            - Q = 3-Kronecker quiver
            - d = (2, 3)
            - θ = (9, -6)

        """
        return self._repr_()

    def dimension(self):
        r"""
        Computes the dimension of the moduli space :math:`M^{\theta-(s)st}(Q,{\bf d})`.

        This involves several cases:

        - If there are :math:`\theta`-stable representations then
          :math:`\dim M^{\theta\rm-sst}(Q,{\bf d}) =
          M^{\theta-st}(Q,{\bf d}) = 1 - \langle {\bf d},{\bf d}\rangle`;
        - if there are no :math:`\theta`-stable representations then
          :math:`\dim M^{\theta-st}(Q,{\bf d}) = -\infty` by convention,
          and we define :math:`\dim M^{\theta\rm\rm-sst} =
          \mathrm{max}_{\tau} \{\dim S_{\tau}\}`,
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
            sage: X = QuiverModuliSpace(Q, (2, 3), condition="semistable")
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

        The Jordan quiver::

            sage: QuiverModuliSpace(JordanQuiver(1), (0,)).dimension()
            0
            sage: X = QuiverModuliSpace(JordanQuiver(1), (0,), condition="stable")
            sage: X.dimension()
            -Infinity
            sage: QuiverModuliSpace(JordanQuiver(1), (1,)).dimension()
            1
            sage: QuiverModuliSpace(JordanQuiver(1), (2,)).dimension()
            2
            sage: QuiverModuliSpace(JordanQuiver(1), (3,)).dimension()
            3
            sage: QuiverModuliSpace(JordanQuiver(1), (4,)).dimension()
            4

        Some generalized Jordan quivers::

            sage: QuiverModuliSpace(JordanQuiver(2), (0,)).dimension()
            0
            sage: QuiverModuliSpace(JordanQuiver(2), (1,)).dimension()
            2
            sage: QuiverModuliSpace(JordanQuiver(2), (2,)).dimension()
            5
            sage: QuiverModuliSpace(JordanQuiver(2), (3,)).dimension()
            10
            sage: QuiverModuliSpace(JordanQuiver(2), (4,)).dimension()
            17

        More generalized Jordan quivers::

            sage: QuiverModuliSpace(JordanQuiver(3), (0,)).dimension()
            0
            sage: QuiverModuliSpace(JordanQuiver(3), (1,)).dimension()
            3
            sage: QuiverModuliSpace(JordanQuiver(3), (2,)).dimension()
            9
            sage: QuiverModuliSpace(JordanQuiver(3), (3,)).dimension()
            19
            sage: QuiverModuliSpace(JordanQuiver(3), (4,)).dimension()
            33

        """
        # setup shorthand
        Q, d, theta = (
            self._Q,
            self._d,
            self._theta,
        )

        # the zero dimension vector only has the zero representation which is semistable
        # but not stable
        if Q._coerce_dimension_vector(d) == Q.zero_vector():
            if self._condition == "semistable":
                return 0
            return -Infinity

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

        OUTPUT: Poincaré polynomial in the variable ``q``

        The Poincare polynomial is defined as

        .. MATH::
            P_X(q) = \sum_{i \geq 0} (-1)^i \dim{\rm H}^i(X;\mathbb{C}) q^{i/2}

        For a quiver moduli space whose dimension vector is
        :math:`\theta`-coprime, the odd cohomology vanishes
        and this is a polynomial in :math:`q`.

        ALGORITHM:

        Corollary 6.9 in MR1974891_.

        .. _MR1974891: https://mathscinet.ams.org/mathscinet/relay-station?mr=1974891

        EXAMPLES:

        Some Kronecker quivers::

            sage: from quiver import *
            sage: Q = KroneckerQuiver()
            sage: X = QuiverModuliSpace(Q, (1, 1))
            sage: X.poincare_polynomial()
            q + 1
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (2, 3))
            sage: X.poincare_polynomial()
            q^6 + q^5 + 3*q^4 + 3*q^3 + 3*q^2 + q + 1
            sage: Q = SubspaceQuiver(5)
            sage: X = QuiverModuliSpace(Q, (1, 1, 1, 1, 1, 2))
            sage: X.poincare_polynomial()
            q^2 + 5*q + 1
        """
        # setup shorthand
        Q, d, theta = self._Q, self._d, self._theta
        d = Q._coerce_dimension_vector(d)
        theta = Q._coerce_vector(theta)

        assert self.is_theta_coprime(), "need coprime"

        k = FunctionField(QQ, "L")
        K = FunctionField(QQ, "q")
        q = K.gen(0)
        f = k.hom(q, K)

        X = QuiverModuliStack(Q, d, theta, condition="semistable")

        P = (1 - q) * f(X.motive())

        assert P.denominator() == 1, "must live in the polynomial ring"

        return P.numerator()

    def betti_numbers(self):
        r"""
        Returns the Betti numbers of the moduli space.

        OUTPUT: Betti numbers of the moduli space

        ALGORITHM:

        Corollary 6.9 in MR1974891_.

        .. _MR1974891: https://mathscinet.ams.org/mathscinet/relay-station?mr=1974891

        EXAMPLES:

        Some Kronecker quivers::

            sage: from quiver import *
            sage: Q = KroneckerQuiver()
            sage: X = QuiverModuliSpace(Q, (1, 1), condition="semistable")
            sage: X.poincare_polynomial()
            q + 1
            sage: X.betti_numbers()
            [1, 0, 1]
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (2, 3), condition="semistable")
            sage: X.betti_numbers()
            [1, 0, 1, 0, 3, 0, 3, 0, 3, 0, 1, 0, 1]

        """
        # setup shorthand
        Q, d, theta = self._Q, self._d, self._theta
        d = Q._coerce_dimension_vector(d)
        theta = Q._coerce_vector(theta)

        assert self.is_theta_coprime(), "need coprime"

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
        r"""
        Returns whether the moduli space is smooth.

        This is easy if the condition is ``"stable"``, because this moduli space
        is always smooth. In the ``"semistable"`` case there is an algorithm,
        by combining the work of Adriaenssens--Le Bruyn and Bocklandt,
        which is currently not implemented.

        EXAMPLES:

        Some 3-Kronecker example::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: QuiverModuliSpace(Q, (2, 3)).is_smooth()
            True
            sage: QuiverModuliSpace(Q, (2, 3), condition="stable").is_smooth()
            True
            sage: QuiverModuliSpace(Q, (3, 3), condition="stable").is_smooth()
            True
            sage: QuiverModuliSpace(Q, (3, 3)).is_smooth()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        # stable locus is always smooth
        if self._condition == "stable":
            return True

        # if we have semistables, it is more subtle
        # this guarantees smoothness without an expensive calculation
        if self._Q.is_theta_coprime(self._d, self._theta):
            return True
        # also guarantees smoothness
        if self.semistable_equals_stable():
            return True

        # need to combine the local quivers from Adriaenssens--Le Bruyn
        # with Bocklandt's criterion for smoothness
        # see https://github.com/QuiverTools/QuiverTools/issues/24
        raise NotImplementedError()

    def semisimple_moduli_space(self):
        r"""
        Return the moduli space with ``theta`` replaced by zero.

        This is the moduli space of semisimple representations for the same quiver
        and the same dimension vector.

        EXAMPLES:

        For an acyclic quiver this moduli space is a point::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (2, 3))
            sage: X.semisimple_moduli_space().dimension()
            0

        For a quiver with oriented cycles we get an affine variety::

            sage: Q = JordanQuiver(2)
            sage: X = QuiverModuliSpace(Q, (3,))
            sage: X.dimension()
            10
        """
        # setup shorthand
        Q, d = (
            self._Q,
            self._d,
        )

        return QuiverModuliSpace(Q, d, theta=Q.zero_vector())

    def is_projective(self) -> bool:
        r"""
        Check whether the moduli space is projective

        EXAMPLES:

        For acyclic quivers the semistable moduli space is always projective::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: QuiverModuliSpace(Q, (2, 3)).is_projective()
            True

        If we have strictly semistable representations, then the stable moduli space
        is only quasiprojective but not projective::

            sage: QuiverModuliSpace(Q, (3, 3), condition="stable").is_projective()
            False

        In pathological cases we can have that the affine moduli space of semisimples
        is reduced to a point, and the projective-over-affine becomes projective::

            sage: Q = CyclicQuiver(3)
            sage: QuiverModuliSpace(Q, (2, 0, 2)).is_projective()
            True

        For the zero dimension vector we get either a point or an empty space, which is
        always projective::

            sage: Q = KroneckerQuiver(3)
            sage: QuiverModuliSpace(Q, (0, 0)).is_projective()
            True
            sage: QuiverModuliSpace(Q, (0, 0), condition="stable").is_projective()
            True

        """
        # setup shorthand
        Q, condition = self._Q, self._condition

        # in the acyclic case the semistable moduli space is always projective
        # the stable moduli space is projective if semistability is stability
        if Q.is_acyclic():
            if condition == "semistable":
                return True
            if condition == "stable":
                return self.semistable_equals_stable()
        # so now Q has oriented cycles: the moduli space is projective-over-affine
        # if we have semistable, or quasiprojective-over-affine is we have stable
        # it suffices that the affine is just a point then
        if condition == "semistable":
            return self.semisimple_moduli_space().dimension() <= 0
        if condition == "stable":
            return (
                self.semisimple_moduli_space().dimension() <= 0
                and self.semistable_equals_stable()
            )

    def picard_rank(self):
        r"""
        Computes the Picard rank of the moduli space.

        We compute this as the Betti number :math:`\mathrm{b}_2`.

        EXAMPLES:

        Kronecker moduli are rank 1::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: QuiverModuliSpace(Q, (2, 3)).picard_rank()
            1

        """
        assert self.is_smooth and self.is_projective(), "must be smooth and projective"

        return self.betti_numbers()[2]

    def index(self):
        r"""
        Computes the index of the moduli space

        The index is the largest integer dividing the canonical divisor in Pic.
        For now this is only implemented for the canonical stability condition.

        EXAMPLES:

        The usual 3-Kronecker example::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: QuiverModuliSpace(Q, (2, 3)).index()
            3

        Subspace quiver moduli have index 1::

            sage: Q = SubspaceQuiver(7)
            sage: QuiverModuliSpace(Q, (1, 1, 1, 1, 1, 1, 1, 2)).index()
            1
        """
        # setup shorthand
        Q, d, theta = self._Q, self._d, self._theta

        if (
            theta == Q.canonical_stability_parameter(d)
            and self.is_theta_coprime()
            and self.is_amply_stable()
        ):
            return gcd(Q._coerce_vector(theta))

        raise NotImplementedError()

    def chow_ring(self, chi=None, classes=None):
        r"""
        Returns the Chow ring of the moduli space.

        For a given datum :math:`(Q, {\bf d}, \theta)` such that
        :math:`Q` is acyclic and :math:`{\bf d}` is :math:`\theta`-coprime,
        the Chow ring of the moduli space of quiver representations
        is described in MR3318266_ and arXiv.2307.01711_.

        .. _MR3318266: https://mathscinet.ams.org/mathscinet-getitem?mr=3318266
        .. _arXiv.2307.01711: https://doi.org/10.48550/arXiv.2307.01711

        Let

        .. MATH::

            R = \bigotimes{i \in Q_0} \mathbb{Q}[x_{i, 1}, \dots, x_{i,d_i}]

        Let :math:`e_{i, j}` be the elementary symmetric function of degree :math:`j`
        in :math:`d_i` variables, and let :math:`\xi_{i, j}` be
        :math:`e_{i, j}(x_{i, 1},\dots,x_{i, d_i})`.
        We denote by :math:`A` the ring of invariants

        .. MATH::

            A := R^{S_{\bf d}} = \mathbb{Q}[\xi_{i, j}],

        where :math:`S_{\bf d} = \prod_{i \in Q_0} S_{{\bf d}_i}` acts by permuting
        the variables.

        The ring :math:`\operatorname{CH}(M^{\theta-st}(Q,{\bf d}))` is a quotient
        of `A` by two types of relations:
        a single linear relation, given by the choice of linearization upon which
        the universal bundles are constructed, and the so-called
        tautological relations, which we define below.

        The *linear relation* given by the linearization `a` is the identity
        :math:`\sum_{i \in Q_0} a_i c_1(U_i) = 0` in :math:`A`.

        A subdimension vector :math:`{\bf e}` of :math:`{\bf d}` is said to be
        "forbidden" if :math:`\mu_{\theta}({\bf e}) > \mu_{\theta}({\bf d})`.
        One actually only needs to consider forbidden dimension vectors that are minimal
        with respect to a certain partial order, see :meth:`Quiver.division_order`.

        We define the *tautological ideal* :math:`I_{\rm taut}` of `R` as the ideal
        generated by the polynomials

        .. MATH::

            \prod_{a\in Q_1}\prod_{k=1}^{e_{s(a)}}
            \prod_{\ell=d_{t(a)}+1}^{d_{t(a)}}
            \left( x_{t(a),\ell}-x_{s(a),k} \right),

        for every forbidden subdimension vector `e` of `d`.

        The tautological relations in `A` are then given by the image of
        :math:`I_{\rm taut}` under the `antisymmetrization` map

        .. MATH::

            \rho : R \to A: \frac{1}{\delta}
            \sum_{\sigma \in S_{\bf d}} sign(\sigma) \sigma \cdot f,

        where :math:`\delta` is the discriminant
        :math:`\prod_{i\in Q_0}\prod_{1\leq k<\ell\leq d_i}(x_{i,\ell}-x_{i,k})`.

        The Chow ring :math:`\operatorname{CH}(M^{\theta\rm-st}(Q,{\bf d}))` is then
        the quotient of `A` by :math:`(\sum_{i\in Q_0} a_i c_1(U_i)) + \rho(I_{taut})`.

        INPUT:

        - ``chi`` -- choice of linearization, we need that :math:`\chi({\bf d})=1`

        - ``classes`` -- list of generators for the polynomial ring (default: None)

        OUTPUT: ring

        EXAMPLES:

        The Kronecker quiver::

            sage: from quiver import *
            sage: Q= KroneckerQuiver()
            sage: X = QuiverModuliSpace(Q, (1, 1))
            sage: chi = (1, 0)
            sage: A = X.chow_ring(chi=chi)
            sage: I = A.defining_ideal()
            sage: [I.normal_basis(i) for i in range(X.dimension()+1)]
            [[1], [x1_1]]


        The 3-Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (2, 3))
            sage: chi = (-1, 1)
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

        The 5-subspace quiver::

            sage: from quiver import *
            sage: Q, d = SubspaceQuiver(5), (1, 1, 1, 1, 1, 2)
            sage: theta = (2, 2, 2, 2, 2, -5)
            sage: X = QuiverModuliSpace(Q, d, theta, condition="semistable")
            sage: chi = (-1, -1, -1, -1, -1, 3)
            sage: A = X.chow_ring(chi=chi)
            sage: I = A.defining_ideal()
            sage: [I.normal_basis(i) for i in range(X.dimension()+1)]
            [[1], [x1_1, x2_1, x3_1, x4_1, x5_1], [x5_2]]

        The ideal Chow ring for our favourite 6-fold has 10 generators, 9 from the
        tautological ideal, and 1 linear relation::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (2, 3))
            sage: chi = (-1, 1)
            sage: R = X.chow_ring(chi=chi);
            sage: R.ambient()
            Multivariate Polynomial Ring in x0_1, x0_2, x1_1, x1_2, x1_3
            over Rational Field
            sage: len(R.defining_ideal().gens())
            10

        """
        Q, d, theta = self._Q, self._d, self._theta
        n = Q.number_of_vertices()

        d = Q._coerce_dimension_vector(d)
        theta = Q._coerce_vector(theta)

        # this implementation only works if d is theta-coprime
        # which implies that d is indivisible.
        assert Q.is_theta_coprime(d, theta), "need coprime"

        def extended_gcd(x):
            r"""
            Computes the gcd and the Bezout coefficients of a list of integers.

            This exists for two integers but seemingly not for more than two.
            """
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

        # if a linearization is not given we compute one here
        if chi is None:
            [g, m] = extended_gcd(d.list())
            chi = vector(m)

        chi = Q._coerce_vector(chi)

        # make sure that chi has weight one, i.e., provides a retraction for
        # X*(PG) --> X*(G).
        assert chi * d == 1

        tautological = self.tautological_ideal(use_roots=False, classes=classes)

        A = tautological.ring()
        linear = A.ideal(
            sum(chi[i] * self._QuiverModuli__generator(A, i, 0) for i in range(n))
        )

        return QuotientRing(A, tautological + linear, names=classes)

    def chern_class_line_bundle(self, eta, classes=None):
        r"""
        Returns the first Chern class of the line bundle

        .. MATH::

            L(\eta) = \bigotimes_{i \in Q_0} \det(U_i)^{-\eta_i},

        where :math:`\eta` is a character of :math:`PG_d`.

        INPUT:

        - ``eta`` -- character of :math:`PG_d` as vector in :math:`\mathbb{Z}Q_0`

        EXAMPLES:

        On the Kronecker 6-fold we can take the canonical line bundle, which we can
        see to have index 3::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (2, 3))
            sage: eta = Q.canonical_stability_parameter((2, 3))
            sage: X.chern_class_line_bundle(eta)
            -3*x1_1bar
        """
        # setup shorthand
        Q, d = self._Q, self._d
        d = Q._coerce_dimension_vector(d)

        A = self.chow_ring(chi=None, classes=classes)

        return -sum(
            eta[i] * A.gen(sum(d[j] for j in range(i)))
            for i in range(Q.number_of_vertices())
        )

    def chern_character_line_bundle(self, eta, classes=None):
        r"""
        Computes the Chern character of L(eta).

        The Chern character of a line bundle `L` with first Chern class `x`
        is given by :math:`e^x = 1 + x + \frac{x^2}{2} + \frac{x^3}{6} + \dots`

        EXAMPLES:

        On the Kronecker 6-fold the canonical line bundle has the following Chern
        character::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (2, 3))
            sage: eta = Q.canonical_stability_parameter((2, 3))
            sage: X.chern_character_line_bundle(eta)
            4617/80*x1_3bar^2 - 1539/40*x1_2bar*x1_3bar + 81/8*x1_1bar^2*x1_2bar
            - 27/8*x1_2bar^2 - 27/4*x1_1bar*x1_3bar - 9/2*x1_1bar^3 + 9/2*x1_1bar^2
            - 3*x1_1bar + 1
        """
        x = self.chern_class_line_bundle(eta, classes=classes)

        return sum(x**i / factorial(i) for i in range(self.dimension() + 1))

    def total_chern_class_universal(self, i, chi, classes=None):
        r"""
        Gives the total Chern class of the universal bundle :math:`U_i(chi)`.

        EXAMPLES:

        The two summands for the Kronecker 6-fold::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (2, 3))
            sage: chi = (-1, 1)
            sage: X.total_chern_class_universal(0, chi)
            x0_2bar + 2*x1_1bar + 1
            sage: X.total_chern_class_universal(1, chi)
            x0_2bar + x1_1bar + 1
        """
        # setup shorthand
        Q, d = self._Q, self._d
        d = Q._coerce_dimension_vector(self._d)

        A = self.chow_ring(chi, classes=classes)

        return 1 + sum(
            A.gen(r + sum(d[j] for j in range(i - 1))) for r in range(d[i - 1])
        )

    def point_class(self, chi=None, classes=None):
        r"""
        Returns the point class as an expression in Chern classes of the
        :math:`U_i` (``chi``).

        INPUT:

        - ``chi`` -- linearization of the universal bundles (default: None)

        The point class is given as the homogeneous component of degree
        :math:`\dim X` of the expression

        .. MATH::

            \prod_{a \in Q_1} c(U_{t(a)})^{d_{s(a)}} / (\prod_{i \in Q_0} c(U_i)^{d_i})

        EXAMPLES

        :math:`\mathbb{P}^7` as a quiver moduli space
        of a generalized Kronecker quiver::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(8)
            sage: X = QuiverModuliSpace(Q, (1, 1))
            sage: chi = (1, 0)
            sage: X.point_class(chi, classes=["o", "h"])
            h^7

        Our favorite 6-fold::

            sage: from quiver import *
            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (2, 3))
            sage: chi = (-1, 1)
            sage: X.point_class(chi, classes=["x1", "x2", "y1", "y2", "y3"])
            y3^2

        A moduli space of the 5-subspace quiver;
        it agrees with the blow-up of :math:`\mathbb{P}^2` in 4 points
        in general position::

            sage: from quiver import *
            sage: Q = SubspaceQuiver(5)
            sage: theta = (2, 2, 2, 2, 2, -5)
            sage: X = QuiverModuliSpace(Q, (1, 1, 1, 1, 1, 2))
            sage: chi = (-1, -1, -1, -1, -1, 3)
            sage: X.point_class(chi, classes=['x1', 'x2', 'x3', 'x4', 'x5', 'y', 'z'])
            1/2*z

        If we don't specify ``chi`` a default that will still work is used, but the
        results do depend on it::

            sage: X.point_class(classes=['x1', 'x2', 'x3', 'x4', 'x5', 'y', 'z'])
            -1/3*y^2
        """
        # setup shorthand
        Q, d = self._Q, self._d
        d = Q._coerce_dimension_vector(d)

        A = self.chow_ring(chi=chi, classes=classes)
        section = A.lifting_map()  # a choice of a section of pi

        p = prod(
            self.total_chern_class_universal(j + 1, chi, classes=classes)
            ** (d * Q.adjacency_matrix().column(j))
            for j in range(Q.number_of_vertices())
        )
        q = prod(
            self.total_chern_class_universal(i + 1, chi, classes=classes) ** d[i]
            for i in range(Q.number_of_vertices())
        )

        quotient = p / q

        pi = A.cover()  # the quotient map

        return pi(section(quotient).homogeneous_components()[self.dimension()])

    def degree(self, eta=None, classes=None):
        r"""
        Computes the degree of the line bundle given by eta.

        INPUT:

        - ``eta`` -- class of line bundle (default: anticanonical line bundle)

        - ``classes`` -- variables to be used (default: None)

        EXAMPLES:

        Rederive a calculation from arXiv.2307.01711_::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: d = (2, 3)
            sage: X = QuiverModuliSpace(Q, d)
            sage: eta = Q.canonical_stability_parameter(d)
            sage: eta = eta / 3
            sage: X.degree(eta)
            57

        .. _arXiv.2307.01711: https://doi.org/10.48550/arXiv.2307.01711

        """
        if eta is None:
            eta = self._Q.canonical_stability_parameter(self._d)

        c = self.chern_class_line_bundle(eta, classes=classes)
        p = self.point_class(classes=classes)

        return c ** self.dimension() / p

    def todd_class(self, chi=None, classes=None):
        r"""
        The Todd class of `X` is the Todd class of the tangent bundle.

        INPUT:

        - ``chi`` -- linearization of the universal bundles (default: None)
        - ``classes`` -- variables to be used (default: None)

        OUTPUT: the Todd class as an element of the Chow ring


        The Todd class is computed in arXiv.2307.01711_. It is given by the formula

        .. MATH::

            td(X) =
            (\prod_{a:i \to j \in Q_1} \prod_{p=1}^{d_j} \prod_{q=1}^{d_i} Q(t_{j,q} -
            t_{i,p}))/(\prod_{i \in Q_0} \prod_{p,q=1}^{d_i} Q(t_{i,q} - t_{i,p}))

        EXAMPLES:

        An example from arXiv.2307.01711_::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (2, 3))
            sage: X.todd_class()
            x1_3bar^2 - 77/60*x1_2bar*x1_3bar + 823/1440*x1_1bar^2*x1_2bar -
            823/4320*x1_2bar^2 - 257/4320*x1_1bar*x1_3bar - 17/32*x1_1bar^3 -
            15/32*x1_3bar + 5/12*x0_2bar + x1_1bar^2 - 3/2*x1_1bar + 1

        .. _arXiv.2307.01711: https://doi.org/10.48550/arXiv.2307.01711
        """

        def todd_Q(t, n):
            r"""
            We call the series :math:`Q(t) = t/(1-e^{-t})` the Todd generating series.
            The function computes the terms of this series up to degree n.
            We use this instead of the more conventional notation `Q` to avoid a
            clash with the notation for the quiver.
            """
            return sum(
                (-1) ** i * (bernoulli(i) * t**i) / factorial(i) for i in range(n + 1)
            )

        def truncate(f, n):
            r"""
            Takes an element in a graded ring and discards all homogeneous components
            of degree > n
            """
            components = f.homogeneous_components()
            keys = [i for i in components]

            return sum(components[i] for i in filter(lambda i: i <= n, keys))

        # setup shorthand
        Q, d = self._Q, self._d
        A = self.chow_ring(chi=chi, classes=classes)
        taut = self._QuiverModuli__tautological_ideal_helper(
            use_roots=False, classes=classes, roots=None
        )
        R, inclusion = taut["ambient_ring"], taut["inclusion"]

        n = self.dimension()

        def short_t(i, p):
            r"""
            Shorthand for the generators of the ambient ring
            from which the Chow ring is constructed
            """
            return self._QuiverModuli__generator(R, i, p)

        num = 1
        den = 1

        # truncating after each step
        # massively cuts runtime
        for a in Q.arrows():
            i, j = a
            for p in range(d[i]):
                for q in range(d[j]):
                    num *= todd_Q(short_t(j, q) - short_t(i, p), n)
                    num = truncate(num, n)

        for i in range(Q.number_of_vertices()):
            for p in range(d[i]):
                for q in range(d[i]):
                    den *= todd_Q(short_t(i, q) - short_t(i, p), n)
                    den = truncate(den, n)

        num = inclusion.inverse_image(num)
        den = inclusion.inverse_image(den)

        # return an element in the Chow ring
        return A(num) / A(den)

    def integral(self, L, chi=None, classes=None):
        r"""
        Integrates the Todd class against an element of the Chow ring.

        INPUT:

        - ``L`` -- element of the Chow ring
        - ``chi`` -- linearization of the universal bundles (default: None)
        - ``classes`` -- variables to be used (default: None)

        OUTPUT: the integral of :math:`td(X) \cdot L` over the moduli space

        EXAMPLES:

        The integral of :math:`\mathcal{O}(i)` on the projective line for some `i`::

            sage: from quiver import *
            sage: Q = KroneckerQuiver()
            sage: X = QuiverModuliSpace(Q, (1, 1))
            sage: L = X.chern_character_line_bundle((1, -1))
            sage: [X.integral(L ** i) for i in range(-5, 5)]
            [-4, -3, -2, -1, 0, 1, 2, 3, 4, 5]

        Hilbert series for the 3-Kronecker quiver as in arXiv.2307.01711_::

            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliSpace(Q, (2, 3))
            sage: L = Q.canonical_stability_parameter((2, 3)) / 3
            sage: O = X.chern_character_line_bundle(L)
            sage: [X.integral(O ** i) for i in range(5)]
            [1, 20, 148, 664, 2206]

        .. _arXiv.2307.01711: https://doi.org/10.48550/arXiv.2307.01711
        """

        integrand = (
            (self.todd_class(chi=chi, classes=classes) * L)
            .lift()
            .homogeneous_components()
        )

        if self.dimension() not in integrand.keys():
            return 0

        return integrand[self.dimension()] / self.point_class(chi=chi, classes=classes)


class QuiverModuliStack(QuiverModuli):
    def __init__(self, Q, d, theta=None, denom=sum, condition="semistable"):
        r"""
        Constructor for a quiver moduli stack.

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
            sage: X = QuiverModuliStack(Q, (2, 3))

        """
        QuiverModuli.__init__(self, Q, d, theta=theta, denom=denom, condition=condition)

    def _repr_(self):
        r""".
        Give a shorthand string presentation for the quiver moduli stack

        EXAMPLES:

        A Kronecker moduli stack::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: QuiverModuliStack(Q, (2, 3))
            moduli stack of semistable representations, with
            - Q = 3-Kronecker quiver
            - d = (2, 3)
            - θ = (9, -6)

        """
        if self.get_custom_name():
            return self.get_custom_name()

        return super()._QuiverModuli__repr_helper("moduli stack")

    def repr(self):
        r"""
        Give a shorthand string presentation for a quiver moduli stack.

        EXAMPLES:

        A Kronecker moduli spac::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: QuiverModuliStack(Q, (2, 3))
            moduli stack of semistable representations, with
            - Q = 3-Kronecker quiver
            - d = (2, 3)
            - θ = (9, -6)

        """
        return self._repr_()

    def dimension(self):
        r"""
        Computes the dimension of the moduli stack :math:`[R^{(s)st}/G]`.

        This is the dimension of a quotient stack, thus we use

        .. MATH::

            dim [R^{{\rm (s)st}}/G] = dim R^{{\rm (s)st}} - dim G

        The dimension turns out to be :math:`-\langle d,d\rangle`
        if the (semi-)stable locus is non-empty.

        EXAMPLES:

        The dimension of a moduli space of stable is off by one from the moduli stack
        because of the generic stabilizer being 1-dimensional::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: X = QuiverModuli(Q, (2, 3))
            sage: X.to_stack().dimension()
            5
            sage: X.to_space().dimension()
            6
        """
        # setup shorthand
        Q, d = self._Q, self._d

        if self.is_nonempty():
            return -Q.euler_form(d, d)
        else:
            return -Infinity

    def is_smooth(self) -> bool:
        r"""
        Return whether the stack is smooth.

        The stack is a quotient of a smooth variety, thus it is always smooth.

        EXAMPLES:

        Nothing interesting to see here::

            sage: from quiver import *
            sage: QuiverModuliSpace(KroneckerQuiver(3), (2, 3)).is_smooth()
            True
        """
        return True

    def motive(self):
        r"""Gives an expression for the motive of the semistable moduli stack

        This really lives inside an appropriate localization of :math:`K_0(Var)`,
        but it only involves the Lefschetz class.

        EXAMPLES:

        Loop quivers::

            sage: from quiver import *
            sage: Q = LoopQuiver(0)
            sage: X = QuiverModuliStack(Q, (2,), (0,))
            sage: X.motive()
            1/(L^4 - L^3 - L^2 + L)
            sage: Q = LoopQuiver(1)
            sage: X = QuiverModuliStack(Q, (2,), (0,))
            sage: X.motive()
            L^3/(L^3 - L^2 - L + 1)

        The 3-Kronecker quiver::

            sage: Q = GeneralizedKroneckerQuiver(3)
            sage: X = QuiverModuliStack(Q, (2, 3))
            sage: X.motive()
            (-L^6 - L^5 - 3*L^4 - 3*L^3 - 3*L^2 - L - 1)/(L - 1)
        """
        # only for semistable.
        # for stable, we don't know what the motive is: it's not pure in general.
        assert self._condition == "semistable"

        # setup shorthand
        Q, d, theta = self._Q, self._d, self._theta
        d = Q._coerce_dimension_vector(d)

        d = Q._coerce_dimension_vector(d)
        theta = Q._coerce_vector(theta)

        K = FunctionField(QQ, "L")
        L = K.gen(0)

        if theta == Q.zero_vector():
            return L ** (-Q.tits_form(d)) / prod(
                prod(1 - L ** (-nu) for nu in range(1, d[i] + 1))
                for i in range(Q.number_of_vertices())
            )

        # start with all subdimension vectors
        ds = Q.all_subdimension_vectors(d, proper=True, nonzero=True)
        # only consider those of greater slope
        ds = list(filter(lambda e: Q.slope(e, theta) > Q.slope(d, theta), ds))
        # put zero and ``d`` back in and sort them conveniently
        ds = ds + [Q.zero_vector(), d]
        ds.sort(key=(lambda e: Q._deglex_key(e, b=max(d) + 1)))

        # Now define a matrix T of size NxN whose entry at position (i,j) is
        # L^<e-f,e>*mot(f-e) if e = I[i] is a subdimension vector of f = I[j]
        # and 0 otherwise
        T = matrix(K, len(ds))
        for i, j in UnorderedTuples(range(len(ds)), 2):
            e, f = ds[i], ds[j]
            if not Q.is_subdimension_vector(e, f):
                continue

            T[i, j] = (
                L ** (Q.euler_form(e - f, e))
                * QuiverModuliStack(
                    Q, f - e, Q.zero_vector(), condition="semistable"
                ).motive()
            )

        # solve system of linear equations T*x = e_N
        # and extract entry 0 of the solution x.
        y = zero_vector(len(ds))
        y[len(ds) - 1] = 1
        x = T.solve_right(y)

        return x[0]

    def chow_ring(self, classes=None):
        r"""Returns the Chow ring of the quotient stack.

        INPUT:

        - ``classes`` -- variables to be used (default: None)

        EXAMPLES:

        The Chow ring of the stack defining the Kronecker 6-fold has as its defining
        ideal the tautological ideal::

            sage: from quiver import *
            sage: Q = KroneckerQuiver(3)
            sage: X = QuiverModuliStack(Q, (2, 3))
            sage: X.tautological_ideal() == X.chow_ring().defining_ideal()
            True
        """

        tautological = self.tautological_ideal(use_roots=False, classes=classes)

        return QuotientRing(tautological.ring(), tautological, names=classes)
