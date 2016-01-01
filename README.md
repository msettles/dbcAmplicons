#dbcAmplicons
=============

Analysis of Double Barcoded Illumina Amplicon Data

## Installation:
Minimum requirements for the dbcAmplicons package: Python 2.X is required to run any and all of the
elements in the package. Some applications within the package rely on FLASH or the RDP classifier to
complete, so the user must determine of those external tools are necessary if the entire pipeline will not be run.
###Programs:
External programs used in the pipeline:
FLASH,
RDP classifier
###Windows Systems
dbcAmplicons has not been tested on Windows systems, although Windows Powershell (active in XP systems and later) may be sufficient to run the analysis. It is highly recommended, however, that the user
install a Linux-based virtual machine to run the pipeline.
###Mac OSX
Mac users should first download and install Xcode and then download the app via a terminal by:  
`git clone https://github.com/msettles/dbcAmplicons.git`  
`cd dbcAmplicons`  
`python setup.py install`  
`dbcAmplicons`  
If installation was successful, this will bring up the usage arguments.
###Linux
For Linux users, some python dependencies may be required before the installation is succesful. To install,
open a terminal and type  
`sudo apt-get install python-dev`  
And type the sudo password at the prompt. Then, download the app via a terminal by:  
`git clone https://github.com/msettles/dbcAmplicons.git`  
`cd dbcAmplicons`  
`sudo python setup.py install`  
`dbcAmplicons`  
If installation was successful, this will bring up the usage arguments.

## VERSION
=======

### version 0.3.0
initial release add support for preprocess and splitreads
### version 0.4.0
add support for multiple input files
### version 0.5.0
reconfigure repository and interface
### version 0.5.1
added support for rdp classification
### version 0.5.2
rdp classification and abundance table generation available
### version 0.5.3
added support for project with no primers
### version 0.5.4
added convert2ReadTo4Read.py script in scripts folder, minor edits to main dbcAmplicons application
### version 0.5.5
added support for gzip in a separate process and ability to infer read file names
### version 0.5.6
added support for Fluidgm down stream processing of samples, in scripts/R
### version 0.5.7
increased speed in gzip read, and bug fixes
### version 0.5.8
fixed bug in barcodes mismatches, added two new flags to preprocess keep primer and test
### version 0.5.9
added script to split by sample, appending output files through iterations
### version 0.6.0
added support for output of biom format in abundance tables
### version 0.6.1
added validate app, some pep modification and various tweeks
### version 0.6.2
bugfixes and cleanup, merge customRDP, scripts, mapping with develop, close issue #1, varius bug fixes
### version 0.6.3
addressed issue #11, validate pairs in primer file
### version 0.6.4
modified abundance output, flash output and reporting, fix issue #12
### version 0.6.5
Show error message and download for RDP when not found
### version 0.6.6
Added support for variable length barcodes in convert2ReadTo4Read.py
### version 0.6.7
Added new subroutine screen, using bowtie2 to screen reads into mapped and unmapped sets
