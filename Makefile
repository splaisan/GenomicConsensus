SHELL = /bin/bash -e

bdist:
	python setup.py build --executable="/usr/bin/env python"
	python setup.py bdist --formats=egg

install:
	python setup.py install

develop:
	python setup.py develop

check:
	# Unit tests
	nosetests tests/unit
	# End-to-end tests
	PATH=`pwd`:$(PATH) cram tests/integration/*.t

test: check
tests: check

doc:
	cd doc; make html

docs: doc


clean:
	-rm -f *.gz *.gff *.fq *.fa *.csv
	-rm -rf dist/ build/ *.egg-info
	-rm -rf doc/_build
	-find . -name "*.pyc" | xargs rm -f


.PHONY: check test tests doc docs clean
