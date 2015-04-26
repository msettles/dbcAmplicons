#!/usr/bin/env python


try:
    from setuptools import setup, Extension
#    from setuptools.command.build_ext import build_ext as _build_ext
#    from setuptools import spawn
except ImportError:
    from distutils.core import setup, Extension
#    from distutils.command.build_ext import build_ext as _build_ext
#    from distutils import spawn


editdist = Extension('editdist', sources=['lib/editdist.c'])
trim = Extension('trim', sources=['lib/trim.c'])

config = \
    {
        'description': 'Processing of Illumina double barcoded amplicon projects',
        'author': 'Matt Settles, Alida Gerritsen',
        'url': 'https://github.com/msettles/DBC_amplicons/py-dbcAmplicons',
        'download_url': 'https://github.com/msettles/DBC_amplicons/py-dbcAmplicons',
        'author_email': 'msettles@uidaho.edu',
        'version': 'v0.6.1-04222015',
        'install_requires': ['biom-format>=2.1.3'],
        'packages': ['dbcAmplicons'],
        'scripts': ['bin/dbcAmplicons', 'scripts/python/convert2ReadTo4Read.py', 'scripts/python/splitReadsBySample.py', 'scripts/R/reduce_amplicons.R'],
        'name': 'dbcAmplicons',
        "ext_package": 'dbcAmplicons',
        'ext_modules': [editdist, trim]
    }


#FLASH2DIR = 'included-apps'
#class dbcBuildExt(_build_ext):
#    def run(self):
#        find_flash2 = spawn.find_executable("flash2")
#        if find_flash2 is None:
#            flash2dld = ['bash', '-c', 'wget', '-o',FLASH2DIR + '/master.zip', 'https://github.com/dstreett/FLASH2/archive/master.zip']
#            spawn(cmd=flash2dld, dry_run=self.dry_run)
#            flash2cmd = ['bash', '-c', 'cd ' + FLASH2DIR + ' && make']
#            spawn(cmd=flash2cmd, dry_run=self.dry_run)
#            config['scripts'].extend( 'included-apps/FLASH2/flash2' )
#        _build_ext.run(self)

#CMDCLASS = {"build_ext": dbcBuildExt}
#setup(cmdclass=CMDCLASS, **config)
#
setup(**config)
