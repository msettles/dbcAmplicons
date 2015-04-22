#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sqlite3 as lite
import os
import misc
import urllib
import logging as log
import zipfile as zipf
import itertools

log.basicConfig(format='%(levelname)s:%(message)s', level=log.INFO)

ncbi_data_url = 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip'


def read_from_zip(archive, file):
    """
    Read a file from an NCBI taxonomy zip archive
    Return an iterator already separated
    """
    zfile = zipf.ZipFile(archive, 'r')
    for line in zfile.read(file).splitlines():
        yield line.rstrip('\t|\n').split('\t|\t')


def fetch_ncbi_data(dest_dir='.', overwrite=False, url=ncbi_data_url):
    """
    Download data from NCBI to generate local taxonomy database.

    * dest_dir - directory in which to save output files
                 (created if necessary).
    * overwrite - don't download if False and target of url exists in dest_dir
    * url - url to archive; default is ncbi_data_url

    Returns (fname, downloaded), where fname is the name of the
    downloaded zip archive, and downloaded is True if a new files was
    downloaded, false otherwise.

    see ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_readme.txt
    portion of function borrowed from taxtastic
    https://github.com/fhcrc/taxtastic
    """

    dest_dir = misc.make_sure_path_exists(dest_dir)
    fout = os.path.join(dest_dir, os.path.split(ncbi_data_url)[-1])

    if os.access(fout, os.F_OK) and not overwrite:
        log.info('%s exists; not downloading' % fout)
        downloaded = False
    else:
        # check for writeablility of directory
        # add in progress bar?
        # ability to stream directly?
        log.info('downloading %s to %s' % (ncbi_data_url, fout))
        urllib.urlretrieve(ncbi_data_url, fout)
        downloaded = True

    return (fout, downloaded)


def ncbi_db_load(engine, archive, root_name='root', maxrows=None):
    """
    Load data from zip archive into database identified by con. Data
    is not loaded if target tables already contain data.

    portion of function borrowed from taxtastic
    https://github.com/fhcrc/taxtastic
    """

    try:
        # nodes
        log.info("Inserting nodes")
        rows = read_nodes(
            rows=read_archive(archive, 'nodes.dmp'),
            root_name=root_name,
            ncbi_source_id=1)
        # Add is_valid
        do_insert(engine, 'nodes', rows, maxrows, add=False)

        # names
        log.info("Inserting names")
        rows = read_names(
            rows=read_archive(archive, 'names.dmp'),
            unclassified_regex=UNCLASSIFIED_REGEX)
        do_insert(engine, 'names', rows, maxrows, add=False)

        # merged
        log.info("Inserting merged")
        rows = read_archive(archive, 'merged.dmp')
        rows = (dict(zip(['old_tax_id', 'new_tax_id'], row)) for row in rows)
        do_insert(engine, 'merged', rows, maxrows, add=False)

        fix_missing_primary(engine)

        # Mark names as valid/invalid
        mark_is_valid(engine)
        update_subtree_validity(engine)

    except lite.IntegrityError, err:
        raise IntegrityError(err)


def read_nodes(rows, root_name, ncbi_source_id):
    """
    fields:
      tax_id                           -- node id in GenBank taxonomy database
      parent tax_id                    -- parent node id in GenBank taxonomy database
      rank                             -- rank of this node (superkingdom, kingdom, ...)
      embl code                        -- locus-name prefix; not unique
      division id                      -- see division.dmp file
      inherited div flag  (1 or 0)     -- 1 if node inherits division from parent
      genetic code id                  -- see gencode.dmp file
      inherited GC  flag  (1 or 0)     -- 1 if node inherits genetic code from parent
      mitochondrial genetic code id    -- see gencode.dmp file
      inherited MGC flag  (1 or 0)     -- 1 if node inherits mitochondrial gencode from parent
      GenBank hidden flag (1 or 0)     -- 1 if name is suppressed in GenBank entry lineage
      hidden subtree root flag (1 or 0)-- 1 if this subtree has no sequence data yet
      comments                         -- free-text comments and citations    

    Return an iterator of rows ready to insert into table "nodes".

    * rows - iterator of lists (eg, output from read_archive or read_dmp)
    * root_name - string identifying the root node (replaces NCBI's default).
    """

    keys = 'tax_id parent_id rank embl_code division_id'.split()
    idx = dict((k, i) for i, k in enumerate(keys))
    rank = idx['rank']

    # assume the first row is the root, reassign rank of root to 'root'
    row = rows.next()
    row[rank] = 'root'
    rows = itertools.chain([row], rows)

    colnames = keys + ['source_id']
    for r in rows:
        row = dict(zip(colnames, r))
        assert len(row) == len(colnames)

        # replace whitespace in "rank" with underscore
        row['rank'] = '_'.join(row['rank'].split())
        yield row

"""
con = None


con = lite.connect('test.db')

with con:

    cur = con.cursor()

    cur.execute("DROP TABLE IF EXISTS Cars")
    cur.execute("CREATE TABLE Cars(Id INT, Name TEXT, Price INT)")
    cur.executemany("INSERT INTO Cars VALUES(?, ?, ?)", cars)


cur.execute("DROP TABLE IF EXISTS Cars")
cur.execute("CREATE TABLE Cars(Id INT, Name TEXT, Price INT)")

"""
