# constructor for a quiver from string à la msinvar

def construct_from_strings(input):
    arrows = [[int(v) for v in chain.split('-')] for chain in input] # list of lists of strings
    n = max(max(chain) for chain in arrows)

    out = Matrix([[0]*n for i in range(n)])
    for chain in arrows:
        for i in range(len(chain)-1):
            out[chain[i]-1, chain[i+1]-1] += 1
    return out