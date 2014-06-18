#!/usr/bin/env python

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

editdist = Extension('editdist', sources = ['lib/editdist.c'])
trim = Extension('trim', sources = ['lib/trim.c'])
config = {
    'description': 'Processing of Illumina double barcoded amplicon projects',
    'author': 'Matt Settles',
    'url': 'https://github.com/msettles/DBC_amplicons/py-dbcAmplicons',
    'download_url': 'https://github.com/msettles/DBC_amplicons/py-dbcAmplicons',
    'author_email': 'msettles@uidaho.edu',
    'version': 'v0.5.5-6182014',
    'install_requires': [],
    'packages': ['dbcAmplicons'],
    'scripts': ['bin/dbcAmplicons', 'scripts/python/convert2ReadTo4Read.py'],
    'name': 'dbcAmplicons',
    "ext_package":'dbcAmplicons',
    'ext_modules': [editdist,trim]
}

setup(**config)

