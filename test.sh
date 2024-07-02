#!bin/sh

# doctests
sage -t quiver/constructions.py
sage -t quiver/quiver.py
sage -t quiver/moduli.py

# coverage
sage -coverage quiver/constructions.py
sage -coverage quiver/quiver.py
sage -coverage quiver/moduli.py

#linting
echo "black"
black quiver/*.py

echo "isort"
isort --profile black quiver/*.py

echo "flake8"
flake8 --max-line-length 88 quiver/*.py

echo "ruff"
ruff check --ignore=E501,E741 quiver/*.py

# sage -coverage.py quiver/constructions.py
# sage -coverage.py quiver/quiver.py
# sage -coverage.py quiver/moduli.py