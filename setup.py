from setuptools import setup, find_packages
from os.path import join, dirname

# Load __VERSION__ from the GenomicConsensus package that is under
# this directory---do NOT import GenomicConsensus, as importing
# GenomicConsensus may fail if it has not actually been installed yet.
globals = {}
execfile("GenomicConsensus/__init__.py", globals)
__VERSION__ = globals["__VERSION__"]


setup(
    name = 'GenomicConsensus',
    version=__VERSION__,
    author='Pacific Biosciences',
    author_email='devnet@pacificbiosciences.com',
    license=open('LICENSES').read(),
    scripts = ['bin/variantCaller.py',
               'bin/summarizeConsensus.py',
               'bin/plurality',
               'bin/quiver'],
    packages = find_packages(),
    package_data={'GenomicConsensus.quiver': ['resources/*/GenomicConsensus/*.ini']},
    include_package_data=True,
    zip_safe = False,
    install_requires=[
        'pbcore >= 0.6.4',
        'numpy >= 1.6.0',
        'h5py >= 2.0.1',
        'ConsensusCore >= 0.6.1'
        ]
    )
