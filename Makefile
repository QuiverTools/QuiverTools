install:
	sage -pip install --upgrade -v -e .
	rm -r build
	rm -r quiver.egg-info

test:
	sage -t quiver

coverage:
	sage --coverage quiver

docs: docs/Makefile
	cd docs && make html
	cd docs && make latexpdf

docs-clean:
	cd docs && make clean

lint:
	black quiver
	isort --profile black quiver
	flake8 --extend-ignore=E741 --max-line-length 88 quiver
	ruff check --ignore=E741 quiver

.PHONY: install test coverage docs-clean docs lint
