from setuptools import setup, find_packages

setup(
    name = 'GenomicConsensus',
    version='0.5.0',
    author='Pacific Biosciences',
    author_email='devnet@pacificbiosciences.com',
    license=open('LICENSES').read(),
    scripts = ['variantCaller.py',
               'summarizeConsensus.py',
               'plurality',
               'quiver'],
    packages = find_packages(),
    package_data={'GenomicConsensus.quiver': ['resources/*.ini']},
    include_package_data=True,
    zip_safe = False,
    install_requires=[
        'pbcore >= 0.5.0',
        'numpy >= 1.6.0',
        'h5py >= 2.0.1'
        ]
    )
