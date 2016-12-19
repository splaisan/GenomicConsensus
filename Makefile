SHELL = /bin/bash -e
INTERNAL_UTILS_PATH = /mnt/secondary/Share/Quiver/Tools

develop:
	python setup.py develop

tests:
	# Unit tests
	nosetests --with-xunit tests/unit
	# End-to-end tests
	PATH=`pwd`:$(PATH) cram --xunit-file=gc-cram.xml tests/cram/*.t

extra-tests:
	# Tests that need to be run by Jenkins but are slowing
	# down the development cycle, so aren't run by "tests"
	# target.
	PATH=`pwd`:$(PATH) cram --xunit-file=gc-extra-cram.xml tests/cram/extra/*.t

internal-tests:
	# Long running tests that depend on files located on PacBio internal NFS
	# servers, including some utilities (exonerate suite, MuMMer)
	(. /mnt/software/Modules/current/init/bash && \
	 module add mummer/3.23         && \
	 module add exonerate/2.0.0     && \
	 module add blasr/2.3.0         && \
	 module add gfftools/dalexander && \
	 cram --xunit-file=gc-internal-cram.xml tests/cram/internal/*.t)

doc:
	cd doc; make html

clean:
	-rm -rf dist/ build/ *.egg-info
	-rm -rf doc/_build
	-rm -f nosetests.xml
	-find . -name "*.pyc" | xargs rm -f

tags:
	find GenomicConsensus -name "*.py" | xargs etags

# Aliases
docs: doc
check: tests
test: tests

.PHONY: check test tests doc docs clean tags
