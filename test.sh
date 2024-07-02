#!/bin/sh

# doctests
sage -t quiver/quiver.py
sage -t quiver/constructions.py
sage -t quiver/moduli.py

# coverage
sage -coverage quiver/quiver.py
sage -coverage quiver/constructions.py
sage -coverage quiver/moduli.py

# sage -coverage.py quiver/constructions.py
# sage -coverage.py quiver/quiver.py
# sage -coverage.py quiver/moduli.py
