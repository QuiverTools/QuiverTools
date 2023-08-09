"""We compute the class of the diagonal inside the Chow ring A^*(X x X) where X is the moduli space of the 3-Kronecker quiver with d = [2,3] and theta = [3,-2]."""

"""t1,t2 are the Chern roots of the bundle p_1^*U_1 and s1,s2,s3 are the Chern roots
of the bundle p_1^*U_2 while tt1,tt2 and ss1,ss2,ss3 are the Chern roots of the bundles p_2^*U_1 and p_2^*U_2. They generate the Chow ring of the quotient stack [R/T] x [R/T]."""
R.<t1,t2,s1,s2,s3,tt1,tt2,ss1,ss2,ss3> = PolynomialRing(QQ)
Delta = (t2-t1)*(s2-s1)*(s3-s1)*(s3-s2)*(tt2-tt1)*(ss2-ss1)*(ss3-ss1)*(ss3-ss2)
t = [t1,t2]
s = [s1,s2,s3]
tt = [tt1,tt2]
ss = [ss1,ss2,ss3]
