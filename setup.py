from distutils.core import setup


setup(
    name='GMM',
    version='0.3.1',
    author='Drew Conway',
    author_email='drew.conway@nyu.edu',
    packages=['gmm', 'gmm.test'],
    scripts=[],
    url='http://pypi.python.org/pypi/GMM/',
    license='LICENSE.txt',
    description='This package provides a basic framework and supporting functionality for generating network structure network structure using graph motifs.',
    long_description=open('README.txt').read(),
)