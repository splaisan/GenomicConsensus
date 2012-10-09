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
	PATH=`pwd`:$(PATH) cram `find tests/integration ! -name "internal*" -name "*.t"`

internal-tests:
	# Tests that depend on files located on PacBio internal NFS
	# servers

	# FIXME: since the `jenkins` user does not have a shared NFS
	# home directory there is no easy way for us to set up ssh
	# keys and have it work nicely on all SGE exec hosts.  Thus
	# for the moment I am specifying one specific host.  This is
	# obviously problematic.

	PATH=`pwd`:$(PATH) qrsh -l hostname=mp-f082 -V -q secondary -cwd cram tests/integration/internal*.t

doc:
	cd doc; make html

docs: doc


clean:
	-rm -f *.gz *.gff *.fq *.fa *.csv
	-rm -rf dist/ build/ *.egg-info
	-rm -rf doc/_build
	-find . -name "*.pyc" | xargs rm -f


.PHONY: check test tests doc docs clean
