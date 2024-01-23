"""
the main functionality for v1:

Quiver:
- generic subdimension vectors
- canonical decomposition

QuiverModuli:
- Harder-Narasimhan types
- Luna types

QuiverModuliSpace:
- Betti numbers
"""

class Quiver:
    _adjacency
    _name
    """
    TODO figure out how Sage really deals with these things, starting from
    - https://doc.sagemath.org/html/en/reference/structure/sage/structure/sage_object.html
    - https://doc.sagemath.org/html/en/reference/structure/sage/structure/element.html
    - the rename method
    """

    def __init__(self, M, name=None)
    # M should be also allowed to be a directed multigraph
    # in fact, I think there should be a getter .graph() which returns the underlying directed multigraph

    def adjacency_matrix(self)
    def underlying_graph(self) # currently returns an undirected graph, make sure to distinguish between directed and undirected!

    # TODO things which are purely graph-theoretic: relegate them to .graph()?
    # i.e., is_sink, is_source, outdegree, indegree, is_connected, is_acyclic, number_of_arrows, number_of_vertices@

    def oppositive_quiver(self)
    def double_quiver(self)

    def euler_matrix(self)
    def euler_form(self, d, e)

    # TODO what about thin_dimension_vector and simple_root? these shouldn't be methods of Quiver
    # maybe have them as helper function?

    def support(self, d)

    def in_fundamental_domain(self, d)
    # maybe call it fundamental_domain_contains? then Q.fundamental_domain_contains(d) reads nice

    def canonical_stability_parameter(self, d)

    def is_generic_subdimension_vector(self, e, d)
    def generic_subdimension_vectors(self, d)

    # don't include generic_ext_vanishing?

    def is_schur_root(self, d)

    def canonical_decomposition(self, d, algorithm="")
    # TODO figure out the different algorithms


class QuiverModuli:

    def __init__(self, d, theta, denominator=sum)
    """
    TODO it should really be possible to initialize this without a theta, defaulting to canonical stability condition?
    also, zero stability should be easy to use
    """
    # Do we want to have denominator as a variable of the class? It doesn't change the moduli space/stack. 
    # Also we might want to compute the HN strata for different denominators without having to change the moduli space.

    # TODO should this be a public method?
    def all_slope_decreasing_sequences(self)

    def has_stable_representation(self, algorithm="generic")
    # other algorithm should be called `local`, generic: Schofield + King, local: Adriaenssens--Le Bruyn
    def has_semistable_representation(self, algorithm="generic")
    # is there another algorithm? is there some variation of Adriaenssens--Le Bruyn?

    def is_amply_stable(self)
    def is_strongly_amply_stable(self)


    # TODO can / should we abbreviate to HN? or hn?
    def harder_narasimhan_types(self)
    # TODO should this be a public method?
    def is_HN_type(self, dstar)
    def HN_codimension(self, dstar)

    def luna_types(self)
    # TODO should this be a public method?
    def is_luna_type(self, tau)

    def semistable_equals_stable(self, algorithm="")

    # TODO this should have a more descriptive name!
    def partial_order(self, d, e)

class QuiverModuliSpace(QuiverModuli):
    # TODO default is semistable, False means stable?
    def __init__(self, d, theta, denominator=sum, semistable=True)

    # uses has_stable, resp. has_semistable accordingly
    def is_empty(self, algorithm="")

    def dimension(self)

    # TODO work out the interface, assert that it is smooth and projective (and not so by accident)
    def chow_ring(self)


# TODO these should be external I think? they feel too specific to be part of the class
def all_weight_bounds(M)
def rigidity_inequality_satisfied(M)
