# MoNvIso
The Modeling eNvironment for Isoforms

# How to use MoNvIso

To install and use MoNvIso you will need:

- Python 3 or beyond
- MODELLER (https://salilab.org/modeller/download_installation.html)

Download the folder *MoNvIso.tar.gz*, unpack it and place the files *master_script.py* *library.py* and *installer.sh* directory where you wish to install MoNvIso.

Run the command *./installer.sh*. This will download the necessary libraries (biopython and psutil), HMMER, COBALT, and the databases of isoforms (canonical and variants) from Uniprot; and it will create an environment called Monviso as a subfolder of the directory where you are installing. It will also create an example of the parameters file and the mutations file.

To run the program, update the mutations file with your genes and mutations, the parameters file (if needed) and then run *python master_script.py* and wait.
