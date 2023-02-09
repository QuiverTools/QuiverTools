"""
input: Q quiver (in reality only the Euler form is needed?)
       dimension vector alpha
output: canonical decomposition
        ordered list of pairs of multiplicity r_i and dimension vector alpha_i

P1: r_i\geq 0
P2: alpha_i a Schur root
P3: alpha_i left orthogonal to alpha_j
P4: r_i=1 when alpha_i imaginary and not isotropic
P5: for all i<j: < alpha_j,alpha_i > >= 0

# how to check whether something is a Schur root
- use Proposition 11.2.7 (due to Schofield) to check for existence of a stable of that dimension vector for the canonical stability condition
    this check is something due to King's original paper
- by Lemma 11.2.8 (2): we can assume indivisibility to check Schur root?


# how to check left orthogonality
    observe that you can _start_ with the exceptional collection given by simples, and thus you have P3
    then the operations somehow preserve the left orthogonality?


Comments:
    - original Schofield paper assumes Q acyclic in Section 4, but not in the other sections?
    - Derksen--Weyman in Compositio assume Q acyclic
"""
