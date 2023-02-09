"""
How to compute generic ext for two dimension vectors a, b?

Option 1: Theorem 5.4 of Schofield's General representations of quivers
    ext(a, b) = max { -<a'',b> } where a'' ranges of generic subdimension vectors of b

Option 2: it is enough to check that a is a generic subdimension vector of a+b for _ext-vanishing_
    (not for general computation)
    = Theorem 3.3 of Schofield's General representations of quivers

    (implemented in code/schur.sage for now)
"""
