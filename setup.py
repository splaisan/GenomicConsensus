from setuptools import setup, find_packages

setup(
    name = 'GenomicConsensus',
    version='0.4.0',
    author='Pacific Biosciences',
    author_email='devnet@pacificbiosciences.com',
    license=open('LICENSES').read(),
    scripts = ['variantCaller.py',
               'summarizeConsensus.py',
               'plurality',
               'quiver'],
    packages = find_packages(),
    zip_safe = False,
    install_requires=[
        'pbcore >= 0.2',
        'numpy >= 1.6.0',
        'h5py >= 1.3.0'
        ]
    )
