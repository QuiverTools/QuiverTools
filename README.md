[![docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://sage.quiver.tools)
[![tests](https://github.com/quiver-tools/quiver.tools/actions/workflows/tests.yml/badge.svg)](https://github.com/quiver-tools/quiver.tools/actions)
[![code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

[![DOI](https://zenodo.org/badge/599021304.svg)](https://zenodo.org/doi/10.5281/zenodo.12680224)
[![mybinder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/QuiverTools/mybinder-sage/master)
# QuiverTools

A SageMath package to deal with quivers and moduli of quiver representations.

You can read its documentation as

* [a webpage](https://sage.quiver.tools)
* [a pdf](https://sage.quiver.tools/documentation.pdf)

A more detailed user guide is in the works.

# Instructions

You can install it by running

``sage --pip install git+https://github.com/QuiverTools/QuiverTools.git``

and then you can just do `from quiver import *` to get started.

Alternatively, you can (if you are patient enough to let it get started) run it from your browser in a notebook: 
[![mybinder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/QuiverTools/mybinder-sage/master)

# How to cite QuiverTools

If you have used this code in any way, please consider citing it as explained on [Zenodo](https://zenodo.org/doi/10.5281/zenodo.12680224). You can choose to cite a specific version, or always the latest version. For the latter you can use `doi:10.5281/zenodo.1268022`.

The following BibTeX entry is a good starting point:

```bibtex
@software{quivertools,
  author = {Belmans, Pieter and Franzen, Hans and Petrella, Gianni},
  title = {QuiverTools},
  url = {https://quiver.tools},
  doi = {10.5281/zenodo.1268022},
}
```

which leads to something like

> Pieter Belmans, Hans Franzen and Gianni Petrella. _QuiverTools_. doi:10.5281/zenodo.1268022. url: https://quiver.tools

# Other toolsets

Other algorithmic aspects of quivers, mostly disjoint of QuiverTools, can be found in

* https://folk.ntnu.no/oyvinso/QPA/, which focuses on calculations with actual representations
* https://smzg.github.io/msinvar/index.html, which focuses on DT-invariants

# Authors

* Pieter Belmans (University of Luxembourg)
* Hans Franzen (University of Paderborn)
* Gianni Petrella (University of Luxembourg)

# Funding

We acknowledge the generous support of:

* the Luxembourg National Research Fund (FNR–17113194 and FNR–17953441)
* the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) SFB-TRR 358/1 2023 "Integral Structures in Geometry and Representation Theory" (491392403)
