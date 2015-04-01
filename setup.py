#!/usr/bin/env python


try:
    from setuptools import setup, Extension
    from setuptools.command.build_ext import build_ext as _build_ext
    from setuptools.spawn import spawn
except ImportError:
    from distutils.core import setup, Extension
    from distutils.command.build_ext import build_ext as _build_ext
    from distutils.spawn import spawn

FLASH2DIR = 'included-apps/FLASH'

editdist = Extension('editdist', sources=['lib/editdist.c'])
trim = Extension('trim', sources=['lib/trim.c'])

config = \
    {
        'description': 'Processing of Illumina double barcoded amplicon projects',
        'author': 'Matt Settles, Alida Gerritsen',
        'url': 'https://github.com/msettles/DBC_amplicons/py-dbcAmplicons',
        'download_url': 'https://github.com/msettles/DBC_amplicons/py-dbcAmplicons',
        'author_email': 'msettles@uidaho.edu',
        'version': 'v0.6.0-04012015',
        'install_requires': ['biom-format'],
        'packages': ['dbcAmplicons'],
        'scripts': ['bin/dbcAmplicons', 'included-apps/FLASH/flash2-amp', 'scripts/python/convert2ReadTo4Read.py', 'scripts/python/splitReadsBySample.py', 'scripts/R/reduce_amplicons.R'],
        'name': 'dbcAmplicons',
        "ext_package": 'dbcAmplicons',
        'ext_modules': [editdist, trim]
    }


class dbcBuildExt(_build_ext):
    def run(self):
        flash2cmd = ['bash', '-c', 'cd ' + FLASH2DIR + ' && make']
        spawn(cmd=flash2cmd, dry_run=self.dry_run)
        _build_ext.run(self)

CMDCLASS = {"build_ext": dbcBuildExt}
setup(cmdclass=CMDCLASS, **config)
