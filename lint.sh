#!/bin/sh

echo "black"
black quiver/*.py

echo "isort"
isort --profile black quiver/*.py

echo "flake8"
flake8 --ignore=E741 --max-line-length 88 quiver/*.py

echo "ruff"
ruff check --ignore=E741 quiver/*.py
