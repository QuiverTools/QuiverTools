"""input: Nothing, everything hard-coded for 3-Kronecker, d = (2,3) and theta = (3,-2)
output: Chow ring of the semi-stable moduli space in terms of generators and relations."""

d = [2,3]
theta = [3,-2]

def canonical_stability():
    return 3*[3,-2]

Permutations.options(mult='r2l')
#RR.<tt1,tt2,ss1,ss2,ss3,xx1,xx2,yy1,yy2,yy3> = PolynomialRing(QQ, order=TermOrder('wdegrevlex', (1,1,1,1,1,1,2,1,2,3)))
#II = RR.ideal([xx1-tt1-tt2, xx2-tt1*tt2, yy1-ss1-ss2-ss3, yy2-ss1*ss2-ss1*ss3-ss2*ss3, yy3-ss1*ss2*ss3])
#R.<t1,t2,s1,s2,s3,x1,x2,y1,y2,y3> = QuotientRing(RR,II)
# The above doesn't work because subs() is not defined in a quotient ring.
# We could find a workaround though, if necessary, using the quotient map.

"""t1,t2 are the Chern roots of the bundle U_1 and s1,s2,s3 are the Chern roots
of the bundle U_2. They generate the Chow ring of the quotient stack [R/T]."""
R.<t1,t2,s1,s2,s3> = PolynomialRing(QQ)
Delta = (t2-t1)*(s2-s1)*(s3-s1)*(s3-s2)
t = [t1,t2]
s = [s1,s2,s3]

"""We give names to the elementary symmetric functions in the t's and the s's
but they live in the ring A^*([R/T]). We need to change the ring to obtain the
subring which is generated by them."""
# There seems to be no direct way to define a subring in sage. You can do a
# workaround though using homomorphisms.
e1, e2 = t1+t2, t1*t2
d1, d2, d3 = s1+s2+s3, s1*s2+s1*s3+s2*s3, s1*s2*s3

"""The Chow ring of [R/G] is the subring of A^*([R/T]) generated by the elementary symmetric functions."""
A.<xx1,xx2,yy1,yy2,yy3> = PolynomialRing(QQ, order=TermOrder('wdegrevlex', (1,2,1,2,3)))

"""Now we define the homomorphism A --> R which sends xx_r to e_r and yy_s to d_s."""
inclusion = A.hom([e1,e2,d1,d2,d3],R)

"""A^*([R/T]) is free as an A^*([R/G])-module. A basis is given by all Demazure operators applied to the following element X."""
X = t1*s1^2*s2

"""W is the Weyl group of T inside G."""
S2 = Permutations(2).list()
S3 = Permutations(3).list()
W = cartesian_product([S2,S3])
w0 = (Permutation([2,1]), Permutation([3,2,1]))

"""Now to the 'forbidden' polynomials. They live in A^*([R/T]).
They should be defined as a function of the forbidden subdimension vector
in the general implementation. The general formula is
prod_{a in Q_1} prod_{r=1}^{d'_{s(a)}} prod_{s=d'_{t(a)}+1}^{d_{t(a)}} (t_{t(a),s} - t_{s(a),r}).
The forbidden subdimension vectors are (1,1) and (2,2)."""
f11 = (s2-t1)^3*(s3-t1)^3
f22 = (s3-t1)^3*(s3-t2)^3

"""The following methods are general. They do not depend on the specific situation."""

def left_permutation_action_on_list(permut,li):
    """Assumes that the permutation and the list have the same length
    Imitates the following. If we have a permutation p and a set {t_1,...,t_n}
    then we want {t_{p^{-1}(1)},...,t_{p^{-1}(n)}}."""
    n = len(li)
    #permutedList = [None]*n
    #for i in range(n):
    #    permutedList[i] = li[permut[i]-1]
    #result = {li[i] : li[permut[i]-1] for i in range(n)}
    #return permutedList
    return list(map(lambda i: li[permut.inverse()[i]-1], range(n)))

def left_permutation_action_on_polynomial(permut,f,alphabet):
    """Computes f(t_{p(1)},...,t_{p(n)})"""
    d = dict(zip(alphabet,left_permutation_action_on_list(permut,alphabet)))
    return f.subs(d)

def divided_difference(i,f,alphabet):
    """Computes (f-(fs_i))/(t_{i+1}-t_i)"""

    n = len(alphabet)
    # Assumes that i<n
    reflection = Permutations(n).simple_reflection(i)
    return (f-left_permutation_action_on_polynomial(reflection,f,alphabet))/(alphabet[i-1]-alphabet[i])

"""The methods from here on use the specific setup."""

def Demazure_operator(w,f):
    """Iterated application of the divided difference operators"""

    reducedt = w[0].reduced_word()
    reducedt.reverse()
    reduceds = w[1].reduced_word()
    reduceds.reverse()
    #g = reduce(lambda i,j: divided_difference(j,divided_difference(i,f,[t1,t2]),[t1,t2]), reducedt)
    g = f
    for i in reducedt:
        g = divided_difference(i,g,[t1,t2])
    for j in reduceds:
        g = divided_difference(j,g,[s1,s2,s3])
    return g

def Schubert(w):
    """The Schubert polynomial S_w is in this context defined as S_w = d_{w_0w^{-1}}X, where d is the Demazure operator and w_0 is the longest Weyl group element."""

    u = (w0[0]*w[0].inverse(), w0[1]*w[1].inverse())
    return Demazure_operator(u,X)


def left_Weyl_group_action_on_polynomial(w,f):
    """Computes the left action of Weyl group element w in W
    on f(t_1,t_2,s_1,s_2,s_3)"""

    # First let w act on {t_1,t_2} and {s_1,s_2,s_3} from the left
    wAppliedTot = left_permutation_action_on_list(w[0],t)
    wAppliedTos = left_permutation_action_on_list(w[1],s)
    # Make dicts {t_i:t_w(i)} and similarly for s
    #dictOft = dict(zip(t,wAppliedTot))
    #dictOfs = dict(zip(s,wAppliedTos))
    #d = dict(dictOft,**dictOfs)
    d = dict(zip(t,wAppliedTot))
    d.update(dict(zip(s,wAppliedTos)))

    # substitute t_i = t_w(i) in f
    return f.subs(d)

def sign(w):
    """Sign of Weyl group element w"""
    return Permutation(w[0]).sign()*Permutation(w[1]).sign()


def antisymmetrization(f):
    """Computes sum_{w in W} sign(w)*(fw)"""
    #rho = 0
    #for w in W:
    #    p1 = w[0]
    #    p2 = w[1]
    #    d = dict(right_action_on_list(p1,[t1,t2]),**right_action_on_list(p2,[s1,s2,s3]))
    #    rho = rho + Permutation(p1).sign()*Permutation(p2).sign()*f.subs(d)
    #return rho
    return sum(list(map(lambda w: sign(w)*left_Weyl_group_action_on_polynomial(w,f), W)))

def symmetrization(f):
    return antisymmetrization(f)/Delta

def schubert_basis():
    """Consists of all Schubert polynomials."""
    return list(map(Schubert,W))

def tautological_upstairs():
    """Auxiliary method. These are the tautological relations already, but they
    live in the wrong ring."""
    indexSet = cartesian_product([schubert_basis(),[f11,f22]])
    return list(map(lambda  pairOfPolys: symmetrization(pairOfPolys[0]*pairOfPolys[1]), indexSet))

def tautological():
    return list(map(inclusion.inverse_image,tautological_upstairs()))

#def ChowRing():
#    """Computation of the Chow ring.
#    It is the quotient by the ideal generated by the tautological relations and the
#    (non-canonical) linear relation xx1-yy1."""
#    li = tautological()
#    li.append(xx1-yy1)
#    I = A.ideal(li)
#    return QuotientRing(A,I)

"""Note to self: If we define R.<x,y,z> = PoylnomialRing(QQ), I = R.ideal(...) and S = QuotientRing(R,I) then the images of the generators x,y,z are by default labeled S.0, S.1, S.2. Can we use this to give names to the Chern classes downstairs?"""

"""Computation of the Chow ring.
It is the quotient by the ideal generated by the tautological relations and the
(non-canonical) linear relation xx1-yy1."""
li = tautological()
li.append(xx1-yy1)
I = A.ideal(li)
ChowRing.<x1,x2,y1,y2,y3> = QuotientRing(A,I)
pi = ChowRing.cover()
sect = ChowRing.lifting_map() # A choice of a section of pi

def point_class():
    """Point class. The computation of the point class comes from my notes. We have to take the homogeneous component of degree 6 = 1-<d,d> of the quotient quot of the two total Chern classes of the bundles F and E.
    When evaluating quot we see that the highest term (degree 6) is y3^2.
    This means that the point class is y3^2."""
    ce = (1-x1+x2)^2*(1-y1+y2-y3)^3
    cf = (1-x1+x2)^9
    quot = cf/ce
    return pi(sect(quot).homogeneous_components()[6])

def ample_Generator():
    """The generator of the ample cone is L(theta) = det(U_1)^{-3} otimes det(U_2)^2, so its first Chern class is -3x1 + 2y1 = -y1."""
    ampleGenerator = -y1
    return ampleGenerator

def Chern_class_of_anticanonical_bundle():
    """The anti-canonical bundle is L(3*theta), so its first Chern class is -3y1."""
    anticanonicalBundle = -3*y1
    return anticanonicalBundle

def degree(eta):
    """Computes the degree of the line bundle L(eta)."""
    chernClass = -eta[0]*x1-eta[1]*y1
    return chernClass^6/point_class()

def degree_of_ample_Generator():
    return ample_Generator()^6/point_class()

def degree_of_anticanonical_bundle():
    return Chern_class_of_anticanonical_bundle()^6/point_class()

def Todd_generating_series(t,n):
    """We call the series
    Q(t) = t/(1-e^{-t}) the Todd generating series. The function computes the terms of this series up to degree n."""
    B = [bernoulli(i) for i in range(n+1)]
    return sum([(-1)^i*B[i]/factorial(i) for i in range(n+1)])

def Todd_class():
    n = 6 # Dimension of moduli space
    """The values numerator and denominator are computed in R, i.e. in the Chow ring of [R/T] and then the inverse image is taken under the inclusion A --> R, where A is the Chow ring of [R/G]"""
    numerator = inclusion.inverse_image(prod([[Todd_generating_series(q-p,n) for q in [s1,s2,s3]] for p in [t1,t2]]))
    denominator = inclusion.inverse_image(prod([[Todd_generating_series(q-p,n) for q in [t1,t2]] for p in [t1,t2]])*prod([[Todd_generating_series(q-p,n) for q in [s1,s2,s3]] for p in [s1,s2,s3]]))
    """The todd class of the tangent bundle is this:"""
    return numerator/denominator
    



"""Other nice methods:"""

def topological_Euler_characteristic():
    """This is the alternating sum of the Betti numbers. As all even cohomology groups vanish in this case and the cohomology is algebraic, it agrees with the vector space dimension of the Chow ring."""
    return I.vector_space_dimension()

def monomial_basis():
    """The function returns a basis in monomials in Chern classes of universal bundles of the Chow ring of the moduli space."""
    basisUpstairs = [I.normal_basis(d) for d in range(7)]
    return [list(map(pi,basisUpstairs[d])) for d in range(7)]

def even_betti_numbers():
    """Returns the even Betti numbers of the moduli space (note that the odd cohomology vanishes). As the cohomology is algebraic, h^{2i} agrees with dim A^i."""
    return [len(I.normal_basis(d)) for d in range(7)]

def Poincare_polynomial():
    """The Poincare polynomial is defined as sum_i (-1)^i h^i q^i. As the cohomology in odd degrees vanishes, it is a Polynomial in q^2 with non-negative coefficients."""
    QT.<q> = PolynomialRing(QQ)
    return sum([even_betti_numbers()[d]*q^(2*d) for d in range(7)])

# ----- old
#I = R.ideal([x1-y1, 3*x2^2-3*x2*y2+y2^2-y1*y3, (3*x2-2*y2)*y3, x2^3-y1*y2*y3+y3^2, -4*x2*y1+y1^3+3*y3, 3*x2^2-x2*y1^2, 3*x2^2+x2*y2-y1^2*y2, x2*y1*y2-3*y2*y3, 3*y1^5-5*y2*y3, x2^3-x2*y1*y3])
