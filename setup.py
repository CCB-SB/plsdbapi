from setuptools import setup

setup(
    name='plsdbapi',
    version='0.1.0',
    author='Anna Hartung',
    packages=['plsdbapi'],
    install_requires=['requests', 'pandas'],
    licence='LICENSE.txt',
    description='API for plasmid database PLSDB',
)
