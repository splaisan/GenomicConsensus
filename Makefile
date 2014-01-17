SHELL = /bin/bash -e
INTERNAL_UTILS_PATH = /mnt/secondary/Share/Quiver/Tools

bdist:
	python setup.py build --executable="/usr/bin/env python"
	python setup.py bdist --formats=egg

install:
	python setup.py install

develop:
	python setup.py develop

tests:
	# Unit tests
	nosetests --with-xunit tests/unit
	# End-to-end tests
	PATH=`pwd`:$(PATH) cram tests/cram/*.t

extra-tests:
	# Tests that need to be run by Jenkins but are slowing
	# down the development cycle, so aren't run by "tests"
	# target.
	PATH=`pwd`:$(PATH) cram tests/cram/extra/*.t

internal-tests:
	# Tests that depend on files located on PacBio internal NFS
	# servers, including some utilities (exonerate suite, MuMMer)
	PATH=`pwd`:$(INTERNAL_UTILS_PATH):$(PATH) cram tests/cram/internal/*.t

doc:
	cd doc; make html

clean:
	-rm -rf dist/ build/ *.egg-info
	-rm -rf doc/_build
	-rm -f nosetests.xml
	-find . -name "*.pyc" | xargs rm -f

tags:
	find GenomicConsensus -name "*.py" | xargs etags

pip-install:
	@which pip > /dev/null
	@pip freeze|grep 'GenomicConsensus=='>/dev/null \
      && pip uninstall -y GenomicConsensus \
      || echo -n ''
	@pip install --no-index \
           --install-option="--install-scripts=$(PREFIX)/bin" \
           ./

# Aliases
docs: doc
check: tests
test: tests

.PHONY: check test tests doc docs clean tags
