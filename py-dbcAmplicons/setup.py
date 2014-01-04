#!/usr/bin/env python

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

editdist = Extension('editdist', sources = ['lib/editdist.c'])
config = {
    'description': 'Processing of Illumina double barcoded amplicon projects',
    'author': 'Matt Settles',
    'url': 'https://github.com/msettles/DBC_amplicons/py-dbcAmplicons',
    'download_url': 'https://github.com/msettles/DBC_amplicons/py-dbcAmplicons',
    'author_email': 'msettles@uidaho.edu',
    'version': '0.2',
    'install_requires': [],
    'packages': ['dbcAmplicons'],
    'scripts': ['bin/dbcPreprocess'],
    'name': 'dbcAmplicons',
    "ext_package":'dbcAmplicons',
    'ext_modules': [editdist]
}

setup(**config)

