"""input: Nothing, everything hard-coded for 3-Kronecker, d = (2,3) and theta = (3,-2)
output: Chow ring of the semi-stable moduli space in terms of generators and relations."""

Permutations.options(mult='r2l')
#RR.<tt1,tt2,ss1,ss2,ss3,xx1,xx2,yy1,yy2,yy3> = PolynomialRing(QQ, order=TermOrder('wdegrevlex', (1,1,1,1,1,1,2,1,2,3)))
#II = RR.ideal([xx1-tt1-tt2, xx2-tt1*tt2, yy1-ss1-ss2-ss3, yy2-ss1*ss2-ss1*ss3-ss2*ss3, yy3-ss1*ss2*ss3])
#R.<t1,t2,s1,s2,s3,x1,x2,y1,y2,y3> = QuotientRing(RR,II)
# The above doesn't work because subs() is not defined in a quotient ring.
# We could find a workaround though, if necessary, usind the quotient map.

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

"""The Chow ring of [R/G] is the subring of A^*([R/T]) generated by the elementary
symmetric functions."""
A.<xx1,xx2,yy1,yy2,yy3> = PolynomialRing(QQ, order=TermOrder('wdegrevlex', (1,2,1,2,3)))

"""Now we define the homomorphism A --> R which sends xx_r to e_r and yy_s to d_s."""
inclusion = A.hom([e1,e2,d1,d2,d3],R)

"""A^*([R/T]) is free as an A^*([R/G])-module. A basis is given
by all Demazure operators applied to the following element X."""
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

"""Computation of the Chow ring.
It is the quotient by the ideal generated by the tautological relations and the
(non-canonical) linear relation xx1-yy1."""
li = tautological()
li.append(xx1-yy1)
I = A.ideal(li)
ChowRing.<x1,x2,y1,y2,y3> = QuotientRing(A,I)

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
    """The Schubert polynomial S_w is in this context defined as S_w = d_{w_0w^{-1}}X,
    where d is the Demazure operator and w_0 is the longest Weyl group element."""

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

# ----- old
#I = R.ideal([x1-y1, 3*x2^2-3*x2*y2+y2^2-y1*y3, (3*x2-2*y2)*y3, x2^3-y1*y2*y3+y3^2, -4*x2*y1+y1^3+3*y3, 3*x2^2-x2*y1^2, 3*x2^2+x2*y2-y1^2*y2, x2*y1*y2-3*y2*y3, 3*y1^5-5*y2*y3, x2^3-x2*y1*y3])
