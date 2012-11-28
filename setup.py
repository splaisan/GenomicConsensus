from setuptools import setup, find_packages

setup(
    name = 'GenomicConsensus',
    version='0.4.0',
    author='Pacific Biosciences',
    author_email='pbiDevNet@pacificbiosciences.com',
    license=open('LICENSES').read(),
    scripts = ['variantCaller.py',
               'summarizeConsensus.py',
               'dumbview.py',
               'plurality',
               'quiver'],
    packages = find_packages(),
    package_data = { "GenomicConsensus.resources" : ["*.json"] },
    zip_safe = False,
    install_requires=[
        'pbcore >= 0.2',
        'numpy >= 1.6.0',
        'h5py >= 1.3.0'
        ]
    )
