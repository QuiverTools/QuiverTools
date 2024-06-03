# quiver via Euler matrix
Q = [[1, -3], [0, 1]] # 3-Kronecker quiver

d = (2,3)
Theta = (9,-6)
#Theta = (d[1], -d[0])

def mu(Theta, d):
    return vector(Theta) * vector(d) / sum(d)

def euler(d, e):
    return vector(d) * matrix(Q) * vector(e)

# check whether a dimension vector is semistable
def is_semistable(d):
    D = HN_types(d, strict=True)
    # TODO this is not right: needs to be list of HN-types, now we treat it as a single HN-type?!

    return all([euler(D[k], D[l]) != 0 for k in range(len(D)) for l in range(k + 1, len(D))])

# subdimension vectors
# always exclude zero
# exclude d if asked for
def subvectors(d, strict=False):
    result = cartesian_product([range(di + 1) for di in d])
    return [e for e in result if e != (0,)*len(d) and (not strict or (e != d))]

def weight(D):
    return tuple(map(sum, zip(*D)))

# is there a nice Pythonic way for this?
def pairs(L):
    return [(L[i], L[j]) for i in range(len(L)) for j in range(i + 1, len(L))]

# enumerate HN-types of a dimension vector
def HN_types(d, strict=False):
    types = []

    # the set of candidate d^1's
    # TODO only semistable weights should be used
    # question: should we check whether d is a Schur root at some point to get the induction started?
    # TODO there is some infinite recursion going on...
    candidates = [[d1] for d1 in subvectors(d, strict=strict)] # if is_semistable(d1)]

    while candidates:
        # intermediate HN-type
        D = candidates.pop()
        e = weight(D)

        # if this happens we have found a HN-type
        if e == d:
            types.append(D)
            continue

        difference = tuple(di - ei for (di, ei) in zip(d, e))
        # considering candidate d^{k+1}'s
        for dkp in subvectors(difference):
            # the two necessary conditions, checked in order of expected computational cost
            if mu(Theta, D[-1]) <= mu(Theta, dkp): continue
            # TODO fix missing semistability condition!
            #if not is_semistable(dkp): continue

            candidates.append(D + [dkp])

    return types

# codimension of HN-type, per Proposition 3.4
def codim(D):
    return - sum([euler(dr, ds) for (dr, ds) in pairs(D)])


# Hans' notation k_s is the product of
# - mu(d^i)
# - the minimal integer for which mu(d^i) becomes an integer (and not merely rational) for all d^i in the HN-type
L = HN_types(d)
print(L)


# we don't rescale by m for now
# so not really the width!
def width(D):
    #print([(mu(Theta, ds), mu(Theta, dr), euler(dr, ds)) for (dr, ds) in pairs(D)])
    return sum((mu(Theta, ds) - mu(Theta, dr)) * euler(dr, ds) for (dr, ds) in pairs(D))


for D in L:
    print("HN-type is {}".format(D))
    print("  codimension is {}".format(codim(D)))
    print("  width is {}".format(width(D)))
    w = width(D)
    print("  checking the condition from Conjecture 6.20")
    for (dr, ds) in pairs(D):
        print("    {}".format(mu(Theta, dr) - mu(Theta, ds)))
        print("    strict inequality {}".format(mu(Theta, dr) - mu(Theta, ds) < w))
    print("")
