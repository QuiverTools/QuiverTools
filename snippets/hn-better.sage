"""
Corollary 5.5 of Reineke's paper: Euler characteristic appears
- if non-zero there exists a semistable?
- no issues with stable != semistable?

=> avoid recursion!

below is the code taken from the Hodge diamond cutter, which implements Corollary 6.9 (which requires d coprime)
"""

def non_empty(Q, d):
    # solve Ax=b for A upper triangular via back substitution
    def solve(A, b):
        assert A.is_square() and A.nrows() == len(b)

        n = len(b) - 1
        x = [0] * (n + 1)

        # start
        x[n] = b[n] / A[n, n]

        # induct
        for i in range(n - 1, -1, -1):
            x[i] = (b[i] - sum([A[i, j] * x[j] for j in range(i + 1, n + 1)])) / A[i, i]

        return x

    @cached_function
    def GL(n):
        r"""Cardinality of general linear group $\\mathrm{GL}_n(\\mathbb{F}_v)$"""
        return prod([v**n - v**k for k in range(n)])

    # not caching this one seems faster
    # @cached_function
    def Rd(Q, d):
        """Cardinality of ``R_d`` from Definition 3.1"""
        return v**sum([d[i] * d[j] * Q[i, j] for i, j in Q.dict()])

    @cached_function
    def Gd(d):
        """Cardinality of ``G_d`` from Definition 3.1"""
        return prod([GL(di) for di in d])

    def Id(Q, d, mu):
        """Returns the indexing set from Corollary 5.5

        These are the dimension vectors smaller than ``d`` whose slope is bigger
        than ``d``, together with the zero dimension vector and ``d`` itself.
        """
        # all possible dimension vectors e <= d
        E = cartesian_product([range(di + 1) for di in d])

        # predicate from Corollary 5.5, E[0] is the zero dimension vector
        return [E[0]] + list(filter(lambda e: mu(e) > mu(d), E[1:])) + [d]

    def Td(Q, d, mu):
        """Returns the upper triangular transfer matrix from Corollary 5.5"""
        # Euler form
        chi = matrix.identity(len(d)) - Q
        # indexing set for the transfer matrix
        I = Id(Q, d, mu)
        # make them vectors now so that we only do it once
        I = list(map(vector, I))

        def entry(Q, e, f):
            """Entry of the transfer matrix, as per Corollary 6.9"""
            fe = f - e

            if all(fei >= 0 for fei in fe):
                return v**(-fe * chi * e) * Rd(Q, fe) / Gd(fe)
            return 0

        T = matrix(K, len(I), len(I))

        for i, Ii in enumerate(I):
            for j in range(i, len(I)):  # upper triangular
                T[i, j] = entry(Q, Ii, I[j])

        return T


