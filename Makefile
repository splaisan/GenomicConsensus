SHELL = /bin/bash -e

bdist:
	python setup.py build --executable="/usr/bin/env python"
	python setup.py bdist --formats=egg

install:
	python setup.py install

develop:
	python setup.py develop

tests:
	# Unit tests
	nosetests tests/unit
	# End-to-end tests
	PATH=`pwd`:$(PATH) cram tests/integration/*.t

extra-tests:
	# Tests that need to be run by Jenkins but are slowing
	# down the development cycle, so aren't run by "tests"
	# target.
	PATH=`pwd`:$(PATH) cram tests/integration/extra/*.t


INTERNAL_UTILS_PATH = /mnt/secondary/Share/Quiver/Tools

internal-tests:
	# Tests that depend on files located on PacBio internal NFS
	# servers, including some utilities (exonerate suite, MuMMer)
	PATH=`pwd`:$(INTERNAL_UTILS_PATH):$(PATH) cram tests/integration/internal/*.t

doc:
	cd doc; make html

docs: doc


clean:
	-rm -f *.gz *.gff *.fq *.fa *.csv
	-rm -rf dist/ build/ *.egg-info
	-rm -rf doc/_build
	-find . -name "*.pyc" | xargs rm -f


.PHONY: check test tests doc docs clean
