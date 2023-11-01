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
                    allLunaTypes = self._Q.all_luna_types(self._d,self._theta)
                    return max([self.dimension_of_luna_stratum(tau) for tau in allLunaTypes])
                else:
                    # I somehow like the convention that the dimension of the empty set is -Infinity
                    return -oo
            else:
                # self._condition == "stable"
                return -oo

    def dimension_of_luna_stratum(self,tau):
        """Computes the dimension of the Luna stratum S_tau."""
        """The dimension of the Luna stratum of tau = [(d^1,p^1),...,(d^s,p^s)] is sum_k l(p^k)(1 - <d^k,d^k>) where for a partition p = (n_1,...,n_l), the length l(p) is l, i.e. the number of rows."""

        """
        EXAMPLES
        The Kronecker quiver:
        sage: from quiver import *
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
        # The Luna stratification is a stratification of the semistable moduli space
        assert self._condition == "semistable"
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

    def tautological_presentation(self, chi=None, chernClasses=None):
        """Returns the Chow ring of the moduli space in terms of generators and relations."""

        """chi is a character of weight one, i.e. it provides a retraction of the inclusion X*(PG) --> X*(G). For this chi*d = 1 must hold. This requires d to be indivisible. If chi=None is given the chi is computed by the extended Euclidean algorithm."""
        """chernClasses is a list of alphanumeric strings that denote the Chern classes of the universal bundles U_i(chi) in the Chow ring. If None is given then they are of the form xi_r."""

        """
        OUTPUT

        A dict
        { "ChowRing" : The Chow ring of the moduli space
        "Generators" : The Chow ring of [R/G],
        "Relations"  : The tautological ideal}
        """

        """
        EXAMPLES

        The Kronecker quiver:
        sage: from quiver import *
        sage: Q = KroneckerQuiver()
        sage: d = vector([1,1])
        sage: theta = vector([1,-1])
        sage: chi = vector([1,0])
        sage: X = QuiverModuliSpace(Q,d,theta,condition="semistable")
        sage: di = X.tautological_presentation(chi)
        sage: I = di["Relations"]
        sage: N = X.dimension()
        sage: [I.normal_basis(i) for i in range(N+1)]
        [[1], [c2_1]]

        The 3-Kronecker quiver:
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: d = vector([2,3])
        sage: theta = vector([3,-2])
        sage: chi = vector([-1,1])
        sage: X = QuiverModuliSpace(Q,d,theta,condition="semistable")
        sage: di = X.tautological_presentation(chi)
        sage: I = di["Relations"]
        sage: N = X.dimension()
        sage: [I.normal_basis(i) for i in range(N+1)]
        [[1],
         [c2_1],
         [c1_2, c2_1^2, c2_2],
         [c2_1^3, c2_1*c2_2, c2_3],
         [c2_1^2*c2_2, c2_2^2, c2_1*c2_3],
         [c2_2*c2_3],
         [c2_3^2]]

        The 5-subspaces quiver:
        sage: Q = SubspaceQuiver(5)
        sage: d = vector([1,1,1,1,1,2])
        sage: theta = vector([2,2,2,2,2,-5])
        sage: chi = vector([-1,-1,-1,-1,-1,3])
        sage: X = QuiverModuliSpace(Q,d,theta,condition="semistable")
        sage: di = X.tautological_presentation(chi)
        sage: I = di["Relations"]
        sage: N = X.dimension()
        sage: [I.normal_basis(i) for i in range(N+1)]
        [[1], [c2_1, c3_1, c4_1, c5_1, c6_1], [c6_2]]

        """

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

        Q, d, theta = self._Q, self._d, self._theta
        # The Chow group has a ring structure only if the space is smooth.
        # Our algorithm works only in the case that sst=st.
        assert is_coprime_for_stability_parameter(d,theta)
        n = Q.number_of_vertices()
        a = Q.adjacency_matrix()

        if chi == None:
            [g,m] = extended_gcd(d.list())
            chi = vector(m)

        # TODO assert that chi has integer entries.
        """Make sure that chi has weight one, i.e. provides a retraction for X*(PG) --> X*(G)."""
        assert chi*d == 1

        if chernClasses == None:
            chernClasses = ['x%s_%s'%(i,r) for i in range(1,n+1) for r in range(1,d[i-1]+1)]

        """This is the Chow ring of the quotient stack [R/T]. The generators ti_r denote the Chern roots of the universal bundles U_i."""
        R = PolynomialRing(QQ,['t%s_%s'%(i,r) for i in range(1,n+1) for r in range(1,d[i-1]+1)])

        def generator(R,i,r):
            """Returns generator(R,i,r) = t{i+1}_{r+1}."""
            return R.gen(r+sum([d[j] for j in range(i)]))

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
        B = [[X(p).expand() for p in Permutations(d[i])] for i in range(n)]
        Bprime = cartesian_product([[f.parent().hom([generator(R,i,r) for r in range(f.parent().ngens())], R)(f) for f in B[i]] for i in range(n)])
        schubert = [prod([bi for bi in b]) for b in Bprime]

        """Generators of the tautological ideal regarded upstairs, i.e. in A*([R/T])."""
        minimalForbiddenSubdimensionVectors = Q.all_minimal_forbidden_subdimension_vectors(d,theta)
        """For a forbidden subdimension vector e of d, the forbidden polynomial in Chern roots is given by prod_{a: i --> j} prod_{r=1}^{e_i} prod_{s=e_j+1}^{d_j} (tj_s - ti_r) = prod_{i,j} prod_{r=1}^{e_i} prod_{s=e_j+1}^{d_j} (tj_s - ti_r)^{a_{ij}}."""
        forbiddenPolynomials = [prod([prod([(generator(R,j,s) - generator(R,i,r))**a[i,j]  for r in range(e[i]) for s in range(e[j],d[j])]) for i in range(n) for j in range(n)]) for e in minimalForbiddenSubdimensionVectors]

        """Define A = A*([R/G])."""
        degrees = []
        for i in range(n):
            degrees = degrees+list(range(1,d[i]+1))
        A = PolynomialRing(QQ, ['c%s_%s'%(i,r) for i in range(1,n+1) for r in range(1,d[i-1]+1)], order=TermOrder('wdegrevlex', degrees))

        E = SymmetricFunctions(ZZ).e()
        """The Chern classes of U_i on [R/G] are the elementary symmetric functions in the Chern roots ti_1,...,ti_{d_i}."""
        elementarySymmetric = []
        for i in range(n):
            elementarySymmetric = elementarySymmetric + [E([k]).expand(d[i], alphabet=[generator(R,i,r) for r in range(d[i])]) for k in range(1,d[i]+1)]
        """Map xi_r to the r-th elementary symmetric function in ti_1,...,ti_{d_i}."""
        inclusion = A.hom(elementarySymmetric, R)

        """Definition of the tautological ideal."""
        tautological = [antisymmetrization(b * f) for b in schubert for f in forbiddenPolynomials]
        tautological = A.ideal([inclusion.inverse_image(g) for g in tautological])

        linear = A.ideal(sum([chi[i]*generator(A,i,0) for i in range(n)]))
        I = linear + tautological

        return { "ChowRing" : QuotientRing(A,I,names=chernClasses),
        "Generators" : A,
        "Relations" : I
        }

    def chow_ring(self, chi=None, chernClasses=None):
        """Returns the Chow ring of the moduli space."""

        """
        EXAMPLES

        The Kronecker quiver:
        sage: from quiver import *
        sage: Q = KroneckerQuiver()
        sage: d = vector([1,1])
        sage: theta = vector([1,-1])
        sage: chi = vector([1,0])
        sage: X = QuiverModuliSpace(Q,d,theta,condition="semistable")
        sage: A = X.chow_ring(chi=chi,chernClasses=['o','h'])
        sage: A.inject_variables()
        Defining o, h
        sage: o
        0
        sage: h
        h
        """

        di = self.tautological_presentation(chi=chi, chernClasses=chernClasses)
        return di["ChowRing"]

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

        raise NotImplementedError

    def diagonal(self, chi=None):
        """Computes the class of the diagonal in the Chow ring of X x X where X is the quiver moduli space."""
        """It is given by the homogeneous component of degree dim X = 1 - <d,d> of the expression c(F)/C(E), where E = bigoplus_{i in Q_0} U_i^vee boxtimes U_i and F = bigoplus_{a in Q_1} U_{s(a)}^vee boxtimes U_{t(a)} = bigoplus_{i,j in Q_0} (U_i^vee boxtimes U_j)^{a_ij}."""

        """
        EXAMPLES

        P^2 as a quiver moduli space:
        sage: from quiver import *
        sage: Q = GeneralizedKroneckerQuiver(3)
        sage: d = vector([1,1])
        sage: theta = vector([1,-1])
        sage: X = QuiverModuliSpace(Q,d,theta,condition="semistable")
        sage: X.diagonal()
        x1_1^2 + x1_1*y1_1 + y1_1^2

        """

        Q, d, theta = self._Q, self._d, self._theta
        n = Q.number_of_vertices()
        N = self.dimension()
        a = Q.adjacency_matrix()

        di = self.tautological_presentation(chi=chi)
        A = di["Generators"]
        I = di["Relations"]

        chernClasses1 = ['x%s_%s'%(i,r) for i in range(1,n+1) for r in range(1,d[i-1]+1)]
        chernClasses2 = ['y%s_%s'%(i,r) for i in range(1,n+1) for r in range(1,d[i-1]+1)]
        chernClasses = chernClasses1+chernClasses2

        AxA = PolynomialRing(QQ,chernClasses)
        inclusion1 = A.hom(chernClasses1,AxA)
        inclusion2 = A.hom(chernClasses2,AxA)
        B = QuotientRing(AxA,inclusion1(I) + inclusion2(I),names=chernClasses)

        pi = B.cover() # The quotient map AxA --> B
        sect = B.lifting_map() # A choice of a section of pi

        chernRoots1 = ['t%s_%s'%(i,r) for i in range(1,n+1) for r in range(1,d[i-1]+1)]
        chernRoots2 = ['u%s_%s'%(i,r) for i in range(1,n+1) for r in range(1,d[i-1]+1)]
        chernRoots = chernRoots1+chernRoots2
        RxR = PolynomialRing(QQ,chernRoots)

        def generatorRxR1(i,r):
            """Returns generatorRxR1(i,r) = t{i+1}_{r+1}."""
            return RxR.gen(r + sum([d[j] for j in range(i)]))

        def generatorRxR2(i,r):
            """Returns generatorRxR2(i,r) = u{i+1}_{r+1}."""
            return RxR.gen(sum([d[j] for j in range(n)]) + r + sum([d[j] for j in range(i)]))

        E = SymmetricFunctions(ZZ).e()
        elementarySymmetric1 = []
        elementarySymmetric2 = []
        for i in range(n):
            elementarySymmetric1 = elementarySymmetric1 + [E([k]).expand(d[i], alphabet=[generatorRxR1(i,r) for r in range(d[i])]) for k in range(1,d[i]+1)]
            elementarySymmetric2 = elementarySymmetric2 + [E([k]).expand(d[i], alphabet=[generatorRxR2(i,r) for r in range(d[i])]) for k in range(1,d[i]+1)]
        elementarySymmetric = elementarySymmetric1 + elementarySymmetric2
        """Map xi_r to the r-th elementary symmetric function in ti_1,...,ti_{d_i} and yi_r to the same in ui_1,...,ui_{d_i}."""
        inclusion = AxA.hom(elementarySymmetric, RxR)

        def total_chern_class_boxproduct(i,j):
            """Computes the total Chern class of U_i^vee boxtimes U_j"""
            c = prod([(1-generatorRxR1(i,r)+generatorRxR2(j,s)) for r in range(d[i]) for s in range(d[j])])
            return pi(inclusion.inverse_image(c))

        numerator = prod([total_chern_class_boxproduct(i,j)**a[i,j] for i in range(n) for j in range(n)])
        denominator = prod([total_chern_class_boxproduct(i,i) for i in range(n)])
        quotient = numerator/denominator

        return pi(sect(quotient).homogeneous_components()[N])


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


# class Quiver_moduli:
#
#     """
#     A quiver moduli space is implemented by storing a quiver Q, a dimension vector d and a stability parameter theta. We want to distinguish between moduli spaces and moduli stacks (stacky = False/True) and if we're considering the stable or the semi-stable version.
#     We do not require theta*d = 0 but work with slope stability instead.
#
#     Variables:
#     quiver
#     dimensionVector
#     stabilityParameter
#     stacky = False
#     version = 'sst'
#     """
#
#     def __init__(self, quiver, dimensionVector, stabilityParameter, stacky=False, version='sst'):
#         assert ( (quiver.number_of_vertices() == dimensionVector.length()) & (quiver.number_of_vertices() == stabilityParameter.length()) )
#         self._quiver = quiver
#         self._dimensionVector = dimensionVector
#         self._stabilityParameter = stabilityParameter
#         self._stacky = stacky
#         assert (version in ['st','sst'])
#         self._version = version
#
#     # I would like something along these lines, but I can't seem to import quiver_class.sage
#     @classmethod
#     def canonical_stability(cls, quiver, dimensionVector, stacky=False, version='sst'):
#         assert (quiver.number_of_vertices() == dimensionVector.length())
#         theta = quiver.canonical_stability_parameter(d)
#         return cls(quiver,dimensionVector,theta,stacky=stacky,version=version)
#
#     def __repr__(self):
#         output = ""
#         if (self._version == 'sst'):
#             output += "Semi-stable "
#         else:
#             output += "Stable "
#         output += "quiver moduli "
#         if self._stacky:
#             output += "stack "
#         else:
#             output += "space "
#         output += "with:\n"+"Q = "+str(self._quiver)+"\n"+"d = "+str(self._dimensionVector)+"\n"+"theta = "+str(self._stabilityParameter)
#         return output
#
#     def quiver(self):
#         return self._quiver
#
#     def dimension_vector(self):
#         return self._dimensionVector
#
#     def stability_parameter(self):
#         return self._stabilityParameter
#
#     def dimension(self):
#         Q = self.quiver()
#         d = self.dimension_vector()
#         if self._stacky: # we're looking at the moduli stack
#             """dim [R^{(s)st}/G] = dim R^{(s)st} - dim G
#             this is -<d,d> if the (semi-)stable locus is non-empty"""
#             nonEmpty = ((self._version == 'sst') and Q.has_semi_stable_representations(d)) or  ((self._version == 'st') and Q.has_stable_representations(d))
#             if nonEmpty:
#                 return -Q.Euler_form(d,d)
#             else:
#                 return float('NaN')
#         else: # we're dealing with the moduli space
#             if Q.has_stable_representations(d):
#                 """If there are stable representations, then both the stable and the semi-stable moduli space have dimension 1-<d,d>"""
#                 return 1 - Q.Euler_form(d,d)
#             else: # no stables
#                 if (self.version == 'st'):
#                     """The stable moduli space is empty so we return NaN"""
#                     return float('NaN')
#                 else:
#                     if Q.has_semi_stable_representations(d):
#                         """This is the case which I don't know how to deal with: there are semi-stables but no stables. Then I think the dimension can be determined using the etale local structure. But I don't quite know how (nor how to implement it)."""
#                         return float('NaN')
#                     else: # no semi-stables
#                         return float('NaN')
#
#     # def slope(self, e):
#     #     """The slope of e is defined as theta*e/(sum_i e_i). We need to ensure that e is non-negative and at least one entry is positive."""
#     #     assert (e.length() == self._dimensionVector.length())
#     #     assert all([(ei >= 0) for ei in e])
#     #     assert any([(ei > 0) for ei in e])
#     #     theta = self.stability_parameter()
#     #     return (theta*e)/(sum(list(e)))
#
#     """
#     Notation for explanations:
#     G = G_d = prod_{i in Q_0} GL_{d_i}
#     T = maximal torus of diagonal matrices
#     PG = G/G_m
#     PT = T/G_m maximal torus of PT
#     W = Weyl group of T in G = Weyl group of PT in PG
#       = prod_{i in Q_0} S_{d_i}
#     R = bigoplus_{a in Q_1} Hom(k^{d_{s(a)}},k^{d_{t(a)}})
#     R^{sst}, R^{st} semi-stable/stable locus
#     """
#
#     def Weyl_group(self):
#         """Returns the Weyl group W both as a group and as a list of its factors."""
#         d = self.dimension_vector()
#         asList = [Permutations(di).list() for di in d]
#         return {"Group" : cartesian_product(asList),
#         "List" : asList}
#
#     def longest_Weyl_group_element(self):
#         """Returns the longest Weyl group element as a list of permutations."""
#         d = self.dimension_vector()
#         return [Permutation(range(di,0,-1)) for di in d]
#
#     """Static methods"""
#
#     """The following methods are general. They do not depend on the specific situation."""
#
#     @staticmethod
#     def left_permutation_action_on_list(permut,li):
#         """Assumes that the permutation and the list have the same length
#         Imitates the following. If we have a permutation p and a set {t_1,...,t_n}
#         then we want {t_{p^{-1}(1)},...,t_{p^{-1}(n)}}."""
#         n = len(li)
#         assert (n == permut.size())
#         return list(map(lambda i: li[permut.inverse()[i]-1], range(n)))
#
#     @staticmethod
#     def left_permutation_action_on_polynomial(permut,f,alphabet):
#         """Computes f(t_{p(1)},...,t_{p(n)})"""
#         d = dict(zip(alphabet,left_permutation_action_on_list(permut,alphabet)))
#         return f.subs(d)
#
#     @staticmethod
#     def divided_difference(i,f,alphabet):
#         """Computes (f-(fs_i))/(t_{i+1}-t_i)"""
#         n = len(alphabet)
#         assert (i < n)
#         reflection = Permutations(n).simple_reflection(i)
#         return (f-left_permutation_action_on_polynomial(reflection,f,alphabet))/(alphabet[i-1]-alphabet[i])
