#!/usr/bin/env python

from distutils.core import setup

LONG_DESCRIPTION = \
'''The program reads one or more input FASTA files.
For each file it computes a variety of statistics, and then
prints a summary of the statistics as output.

The goal is to provide a solid foundation for new bioinformatics command line tools,
and is an ideal starting place for new projects.'''


setup(
    name='MC_Star',
    version='0.1.0.0',
    author='TianyiLu',
    author_email='lutianyi9494@icloud.com',
    packages=['MC_Star'],
    package_dir={'MC_Star': 'MC_Star'},
    entry_points={
        'console_scripts': ['MC_Star = MC_Star.MC_Star:main']
    },
    url='https://github.com/TianyiLu9494/MC_Star',
    license='LICENSE',
    description=('A prototypical bioinformatics command line tool'),
    long_description=(LONG_DESCRIPTION),
    install_requires=["biopython"],
)
