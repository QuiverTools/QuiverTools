[![tests](https://github.com/quiver-tools/quiver.tools/actions/workflows/tests.yml/badge.svg)](https://github.com/quiver-tools/quiver.tools/actions)


# Instructions

You can install it using

``sage --pip install -e .``

and then you can just do `from quiver import *` to get started.


# Summary

The goal is to implement algorithms related to
- quiver representations
- moduli spaces of quiver representations

For now we steer clear of quivers with potentials, see https://smzg.github.io/msinvar/index.html for that.

# Reading list

## Urgent

- https://mathscinet.ams.org/mathscinet-getitem?mr=2556134 (!!)
- https://mathscinet.ams.org/mathscinet-getitem?mr=1752774
- https://arxiv.org/abs/2109.08905
- examples from https://nbviewer.org/github/smzg/msinvar/tree/main/notebooks/
  e.g. https://nbviewer.org/github/smzg/msinvar/blob/main/notebooks/Stable%20representations.ipynb is useful (!)

## Less urgent

- https://mathscinet.ams.org/mathscinet-getitem?mr=2483828

# List of algorithms to implement

These are now mentioned in the classes, usually as placeholders.

* King's criterion for the existence of a stable representation of given dimension vector

* Adriaenssens--Le Bruyn's method to determine dimension vectors admitting stable representations

* Reineke's criterion for the existence of a semistable representation of given dimension vector

* Schofield's and Derksen--Weyman's algorithm for Kac's canonical decomposition of a dimension vector
  - Lemma 17 of the Derksen--Weyman Compositio paper has a variation for non-acyclic quivers

  Derksen--Weyman book, Algorithm 11.9.11
  Algorithm 8 in their Compositio paper

* Reineke's computation of Betti numbers for arbitrary~$Q$, $\mathbf{d}$ and~$\theta$
  - https://mathscinet.ams.org/mathscinet-getitem?mr=4000572
      - Hans says: generic means the Euler form is symmetric on the ray (or subspace) of dimension vectors with fixed slope
      - intersection homology appears, and is related to DT invariants
      - see Higgs branch formula in https://www.lpthe.jussieu.fr/~pioline/Computing/CoulombHiggs.pdf

* Engel--Reineke's computation of Betti numbers for smooth models
  - see https://mathscinet.ams.org/mathscinet-getitem?mr=2511752
  it is Corollary 5.3 of https://arxiv.org/pdf/0706.4306.pdf ?

* Reineke's enumeration of the Harder--Narasimhan strata
  (and thus also Kempf--Ness strata)
  including a generalisation for arbitrary denominators in the slope-stability function

* Bocklandt's reduction method to check smoothness of semisimple quiver moduli
  moduli spaces of semisimple quiver representations

* Bocklandt--Van de Weyer's condition for cofreeness of quiver moduli

* Schofield's description of perpendicular categories in terms of a smaller quiver

* window methods for tautological vector bundles using the Kempf--Ness stratification

* point counts 

* computations in Chow ring, à la 
  - https://www.math.sciences.univ-nantes.fr/~sorger/en/chow/
    
    the source code can be found at https://github.com/sagemath/sagetrac-mirror/tree/u/gh-sorger-c/chow/src/sage/schemes/chow

    don't forget to set the tangent bundle, à la https://www.math.sciences.univ-nantes.fr/~sorger/en/chow/html/chow/sage/schemes/chow/scheme.html

    there is also another abandoned Sage implementation at https://github.com/hiepdang/Sage, called Schubert3
  - http://www2.macaulay2.com/Macaulay2/doc/Macaulay2-1.18/share/doc/Macaulay2/Schubert2/html/index.html
  using https://mathscinet.ams.org/mathscinet-getitem?mr=3318266
  
  - also: https://github.com/8d1h/IntersectionTheory Julia implementation, could be a good learning opportunity for everyone involved

* check for fundamental set à la Domokos, page 3, https://arxiv.org/abs/2303.08522

* criterion for amply stable in \S4 of https://arxiv.org/pdf/1410.0466.pdf

# Possible features

- interface with GAP/QPA to compute e.g. dimensions of Hom-spaces for _explicit_ representations (cool idea by the horse Clever Hans)


# Quivers that need a constructor

- generalised Kronecker quiver
- 3 vertex quivers
- Jordan quiver
- subspace quivers
- generalised subspace quivers
- bipartite quivers (easier way to describe it: use I for the sources, J for the sinks, then an IxJ matrix suffices)
- Dynkin quivers
- cyclic quiver


# Window methods

Our rigidity paper gives the formulas for the windows.

What needs to be done is determine the weights for equivariant bundles.
These are all determined using highest weights of GL(d),
and thus we need to compute the Kempf--Ness thingie of Dan (or Constantin)
for a highest weight of GL(d).


# More things to think about?

- Implement the reduction steps of https://arxiv.org/pdf/2303.08522.pdf


# Understanding what's in CoulombHiggs
- QuiverPoincarePolynomial
- QuiverPoincarePolynomialRat


# Other things to consider

* counting stables over finite fields: https://arxiv.org/abs/0708.1259
  implemented in https://smzg.github.io/msinvar/wall_crossing.html?highlight=simple#msinvar.wall_crossing.WallCrossingStructure.simple
  see https://mathscinet.ams.org/mathscinet-getitem?mr=2250021 too (Reineke) (section 7 has examples?)

* Kac polynomial: relates to the previous thing
  - Hua's paper: https://mathscinet.ams.org/mathscinet-getitem?mr=1752774
  - Hua's thesis has more examples: https://mathscinet.ams.org/mathscinet-getitem?mr=2716788
  - what's up with https://arxiv.org/abs/2207.09839 ?
  - see https://arxiv.org/abs/math/0608321 for more algorithmic approach?
    implemented in msinvar (in a opaque way)

* torus localisation (see email of Hans from February 8 2023)
  - see also torus-equivariant Chow rings

* there's some fun results in https://link.springer.com/article/10.1007/s10468-005-8762-y
  these are generalised in https://mathscinet.ams.org/mathscinet-getitem?mr=2511752: noncommutative Hilbert schemes are a particular case of smooth models
  they give a description of the cell decomposition
  Hans and Sergey have a simpler description of those

* tautological presentation of cohomology
  - https://mathscinet.ams.org/mathscinet-getitem?mr=3318266
  - https://mathscinet.ams.org/mathscinet-getitem?mr=1324213
  - Hans has Singular code that works for generalised Kronecker quivers?

* torus-equivariant Chow rings of quiver moduli
  - is https://mathscinet.ams.org/mathscinet-getitem?mr=4155242 algorithmic?

* wall-and-chamber decompositions: how?!
  - Proposition 8.6 of https://scholar.harvard.edu/files/elanakalashnikov/files/short_version.pdf? Hans thinks not
    https://faculty.math.illinois.edu/Macaulay2/doc/Macaulay2-1.19.1/share/doc/Macaulay2/ThinSincereQuivers/html/_cone__System.html thinks otherwise
    or is 8.6 only for quiver
    https://arxiv.org/abs/1508.00539
  - example in Remark 2.2 of https://alco.centre-mersenne.org/item/10.5802/alco.152.pdf
  - is there some abelian-non-abelian correspondence that can be done?
  - Elena's abelianised quiver is called a covering quiver
  
* Python version in https://github.com/marybarker/ThinSincereQuivers/tree/master/python_version

* various genera (in the sense of Hirzebruch), see https://arxiv.org/abs/2306.06660v1 for inspiration


# Other things to not (?) consider

* https://arxiv.org/abs/1104.5698 and https://smzg.github.io/msinvar/higgs_bundles.html
* https://smzg.github.io/msinvar/curve_DT.html and https://arxiv.org/abs/1310.4991
* what's in https://github.com/marybarker/quiver_mutations
* https://www.math.uni-bielefeld.de/~jgeuenich/string-applet/
  dimension of kQ
* https://arxiv.org/abs/2306.06660v1 is an excellent example of how things should be done!
  also, we should be able to implement these calculations in complete generality

# Implementation ideas

- msinvar has a short way of describing quivers it seems? https://smzg.github.io/msinvar/quivers.html
  e.g. `Q=Quiver('1-2-3,1-3', prec=[2,2,2]); Q`

- msinvar has some ideas on stability functions? see https://smzg.github.io/msinvar/stability.html


# Papers with explicit calculations that ought to be done in QuiverTools

* https://arxiv.org/abs/2303.00503
  their terminology of anti-attractor stability refers to the anticanonical stability condition
* https://arxiv.org/pdf/2012.14358.pdf ?
