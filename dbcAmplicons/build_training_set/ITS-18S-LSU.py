#!/usr/bin/env python

## download eukariotic classification sequences 18S, ITS1, 5.2S, ITS2 28S
#
## NCBI search term
## ((ITS2[Title] OR ITS1[Title] OR 18S[Title] OR 28S[Title] OR 5.8S[Title]) AND 100:2000[Sequence Length]) NOT (Uncultured* OR Unknown* OR Unidentified*)
##
#
## Written by Matt Settles
## Nov 10, 2014
#
#
#
from Bio import SeqIO
from Bio import Entrez
import urllib
import sys
import re
import cStringIO
import time

Entrez.email = "msettles@uidaho.edu"

query = '((ITS2[Title] OR ITS1[Title] OR 18S[Title] OR 28S[Title] OR 5.8S[Title]) AND 100:2000[Sequence Length]) NOT (Uncultured* OR Unknown* OR Unidentified*)'
handle = Entrez.esearch(db="nuccore", term=query, retmax=1000000)
record = Entrez.read(handle)
gi_list = record["IdList"]
len(gi_list)

gi_str = ",".join(gi_list)

gi = gi_list[0]

#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nucleotide&db=taxonomy&id=341926284


handle = Entrez.elink(dbfrom="nucleotide",db="taxonomy",id=gi, retmode="text" )
record = Entrez.read(handle)

taxid = record[0]['LinkSetDb'][0]['Link'][0]['Id']

handle = Entrez.efetch(db="nuccore", id=gi, rettype="gb", retmode="text")
gb = SeqIO.parse(handle, "gb")

'gi|' + gi + '|dbj|' + seq.id + '|taxid|' + taxid + '|org|' + seq.annotations['organism'] + '|'







if len(sys.argv) != 3:
    print 'Usage: download_from_ncbi.py "query" outfile'
    print "where query is quoted, and contains a string such as:"
    print '\t"\'Homo sapiens\'[porgn:__txid9606] AND 5000:50000000[Sequence Length]"'
    sys.exit()

query = sys.argv[1]
outfile = sys.argv[2]

# query  = "'Homo sapiens'[porgn:__txid9606] AND 5000:50000000[Sequence Length])"
# outfile = "tmp.fasta"

startt = time.time()

utils = 'http://www.ncbi.nlm.nih.gov/entrez/eutils'
db = 'nucleotide'
report = "gb"

#Set up the initial search and get statistics about the search results:
esearch = "%s/esearch.fcgi?db=%s&retmax=1&usehistory=y&term=%s" % (utils, db, query)

print "Executing first query:\n\t%s" % esearch

esearch_results_handle = urllib.urlopen(esearch)
tmp = esearch_results_handle.read()
Count = re.search('<Count>.*</Count>', tmp).group(0)[7:-8]
QueryKey = re.search('<QueryKey>.*</QueryKey>', tmp).group(0)[10:-11]
WebEnv = re.search('<WebEnv>.*</WebEnv>', tmp).group(0)[8:-9]
print "Count = ", Count, "QueryKey = ", QueryKey, "WebEnv = ", WebEnv


#Set up for downloading actual data:
outf = open(outfile, 'w')
retmax = 1000
restart = 0
total_recs = 0
lens = {}
while restart < int(Count):
    rquery = utils + "/efetch.fcgi?rettype=" + report + "&retmode=text&retstart=%s" % restart
    rquery += "&retmax=%s" % retmax + '&db=' + db + '&query_key=' + QueryKey
    rquery += "&WebEnv=" + WebEnv

    print "Executing query:\n\t%s" % rquery
    outSIO = cStringIO.StringIO()
    #stream the data into a cStringIO object
    try:
        outSIO.write(urllib.urlopen(rquery).read())
        outSIO.seek(0)
    except Exception as exc:
        retmax = int(retmax/2)
        print "Error on restart=%s" % restart
        print "Producing exception: \n\t" + str(exc)
        print "Retrying with restart=%s, retmax=%s" % (restart, retmax)
        continue

    print "Query successful, parsing records."
    #And parse the cStringIO object with SeqIO
    recs = 0
    retmax = 1000
    for record in SeqIO.parse(outSIO, 'gb'):
        
        SeqIO.write(record, outf, 'gb')
        recs += 1
        lens[record.id] = len(record)
    total_recs += recs
    restart = restart + recs
        del outSIO
        print "Processed %s records." % recs
        print "Processed %s/%s total records at a rate of %s/S" % (total_recs, Count, (total_recs)/(time.time() - startt))
        if restart < int(Count):
            time.sleep(30)

print "Finished processing all records."
print "Expected record count: %s, actual record count: %s" % (Count, total_recs)
print "Elapsed time=%s s" % (time.time() - startt)
outf.close()

#Write out stats
outf = open("download_stats.csv", 'w')
outf.write("seqid,length\n")
for k in lens:
    outf.write("%s,%s\n" % (k, lens[k]))
outf.close()




