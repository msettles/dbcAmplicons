DBC_amplicons
=============

Analysis of Double Barcoded Illumina Amplicon Data

## Installation:
Minimum requirements for the dbcAmplicons package: Python 2.X is required to run any and all of the
elements in the package. Some applications within the package rely on FLASH or the RDP classifier to
complete, so the user must determine of those external tools are necessary if the entire pipeline will not be run.
###Programs:
External programs used in the pipeline:
FLASH
RDP classifier
###Windows Systems
dbcAmplicons has not been tested on Windows systems, although Windows Powershell (active in XP sys-
tems and later) may be sufficient to run the analysis. It is highly recommended, however, that the user
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
