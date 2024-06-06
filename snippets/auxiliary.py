"""The following methods are general. They do not depend on the specific situation."""


def left_permutation_action_on_list(permut, li):
    """Assumes that the permutation and the list have the same length
    Imitates the following. If we have a permutation p and a set {t_1,...,t_n}
    then we want {t_{p^{-1}(1)},...,t_{p^{-1}(n)}}."""
    n = len(li)
    return list(map(lambda i: li[permut.inverse()[i] - 1], range(n)))


def left_permutation_action_on_polynomial(permut, f, alphabet):
    """Computes f(t_{p(1)},...,t_{p(n)})"""
    d = dict(zip(alphabet, left_permutation_action_on_list(permut, alphabet)))
    return f.subs(d)


def divided_difference(i, f, alphabet):
    """Computes (f-(fs_i))/(t_{i+1}-t_i)"""

    n = len(alphabet)
    # Assumes that i<n
    reflection = Permutations(n).simple_reflection(i)
    return (f - left_permutation_action_on_polynomial(reflection, f, alphabet)) / (
        alphabet[i - 1] - alphabet[i]
    )
