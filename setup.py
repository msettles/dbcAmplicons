#!/usr/bin/env python


try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension


editdist = Extension('editdist', sources=['lib/editdist.c'])
trim = Extension('trim', sources=['lib/trim.c'])

config = \
    {
        'description': 'Processing of Illumina double barcoded amplicon projects',
        'author': 'Matt Settles, Alida Gerritsen',
        'url': 'https://github.com/msettles/dbcAmplicons',
        'download_url': 'https://github.com/msettles/dbcAmplicons',
        'author_email': 'settles@ucdavis.edu',
        'version': 'v0.7.0-02082016',
        'install_requires': ['biom-format>=2.1.3'],
        'packages': ['dbcAmplicons'],
        'scripts': ['bin/dbcAmplicons', 'scripts/python/convert2ReadTo4Read.py', 'scripts/python/splitReadsBySample.py', 'scripts/R/reduce_amplicons.R'],
        'name': 'dbcAmplicons',
        "ext_package": 'dbcAmplicons',
        'ext_modules': [editdist, trim]
    }

setup(**config)
