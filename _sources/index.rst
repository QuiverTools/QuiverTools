***********
QuiverTools
***********

QuiverTools is a SageMath package to deal with quivers and moduli of quiver representations.
Below you can find its documentation.
A more detailed user guide is in the works.

To install it, run

.. code-block::

   sage --pip install git+https://github.com/QuiverTools/QuiverTools.git

and then you can simply run

.. code-block:: python

   from quiver import *

to get started.

For more information, see `https://quiver.tools <https://quiver.tools>`_.

**Authors**

* Pieter Belmans (University of Luxembourg)
* Hans Franzen (University of Paderborn)
* Gianni Petrella (University of Luxembourg)

**Funding**

We acknowledge the generous support of:

* the Luxembourg National Research Fund (FNR–17113194 and FNR–17953441)
* the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) SFB-TRR 358/1 2023 "Integral Structures in Geometry and Representation Theory" (491392403)

Quivers
=======

.. autoclass:: quiver.Quiver
    :members:
    :special-members:
    :member-order: bysource

Moduli spaces
=============

.. autoclass:: quiver.QuiverModuli
    :members:
    :special-members:
    :member-order: bysource

.. autoclass:: quiver.QuiverModuliSpace
    :members:
    :special-members:
    :member-order: bysource

.. autoclass:: quiver.QuiverModuliStack
    :members:
    :special-members:
    :member-order: bysource

Constructing quivers
====================

.. autofunction:: quiver.disjoint_union
.. autofunction:: quiver.GeneralizedKroneckerQuiver
.. autofunction:: quiver.KroneckerQuiver
.. autofunction:: quiver.ThreeVertexQuiver
.. autofunction:: quiver.LoopQuiver
.. autofunction:: quiver.JordanQuiver
.. autofunction:: quiver.GeneralizedJordanQuiver
.. autofunction:: quiver.SubspaceQuiver
.. autofunction:: quiver.ThickenedSubspaceQuiver
.. autofunction:: quiver.GeneralizedSubspaceQuiver
.. autofunction:: quiver.DynkinQuiver
.. autofunction:: quiver.ExtendedDynkinQuiver
.. autofunction:: quiver.CyclicQuiver
.. autofunction:: quiver.BipartiteQuiver
