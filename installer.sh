#!/bin/bash

echo ""
echo "Thank you for installing MoNvIso!"
echo ""

check_python(){

python_version=$(python -V)
python_version="${python_version//[!0-9]/}" 
python_version="${python_version:0:1}"

if [ $python_version -lt 3 ]
then
	echo "Please update your Python! You will need Python 3.0 or beyond to run MoNvIso."
	exit 1 

else
	echo "Creating Monviso environment"
	echo ""
	python -m venv $1/Monviso
	cd Monviso
	source bin/activate
fi
}


install_hmmer(){

echo "installing HMMER"
echo ""

wget http://eddylab.org/software/hmmer/hmmer.tar.gz

tar xvzf hmmer.tar.gz

rm hmmer.tar.gz
hmmer_folder=$(ls | grep hmmer)
cd $hmmer_folder
mkdir $1/Monviso/hmmer
./configure --prefix=$1/Monviso/hmmer/
make
make install
cd $1/Monviso
rm $hmmer_folder -r
}

download_cobalt(){

wget -r ftp://ftp.ncbi.nlm.nih.gov/pub/cobalt/executables/LATEST/*x64-linux.tar.gz
mv ftp.ncbi.nlm.nih.gov/pub/cobalt/executables/LATEST/*x64-linux.tar.gz .
rm ftp.ncbi.nlm.nih.gov -r

tar -xzvf ncbi-cobalt-*-linux.tar.gz
rm *.tar.gz
}

download_databases(){

while true; do
    read -p "Do you wish to download Uniprot/Swissprot databases (MoNvIso will look for them)?   " yn
    case $yn in
        [Yy]* ) wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz; wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz; gzip -d uniprot_sprot.fasta.gz; gzip -d uniprot_sprot_varsplic.fasta.gz; rm uniprot*.gz; break;;
        [Nn]* ) break;;
        * ) echo "Please answer yes or no.";;
    esac
done
}

create_parameters(){

cobalt=$(ls | grep *cobalt*)

touch parameters.dat

if [ -f "uniprot_sprot.fasta" ]; then

echo "DB_LOCATION=$1/Monviso/

COBALT_HOME=$1/Monviso/$cobalt/bin/

HMMER_HOME=$1/Monviso/hmmer/bin/
  
RESOLUTION=4.50

SEQID=25

HMM_TO_IMPORT=100

MODEL_CUTOFF=5

PDB_TO_USE=10

NUM_OF_MOD_WT=1

NUM_OF_MOD_MUT=1" > parameters.dat

else
echo "DB_LOCATION=/path/to/db/

COBALT_HOME=$1/Monviso/$cobalt/bin/

HMMER_HOME=$1/Monviso/hmmer/bin/
  
RESOLUTION=4.50

SEQID=25

HMM_TO_IMPORT=100

MODEL_CUTOFF=5

PDB_TO_USE=10 

NUM_OF_MOD_WT=1

NUM_OF_MOD_MUT=1" > parameters.dat
fi

}

create_example_mutations(){

touch mutations.txt

echo 'GRIN1
R       844     C
Ala     349     Thr
Pro     578     Arg
Ser     688     Tyr
Tyr     647     Ser

GRIN2B
E413G
C436R
M1342R
L1424F
PRO1439ALA

AP2M1
R170W
V106L' > mutations.txt

}

download_libraries(){

pip install biopython
pip install psutil
}

cwd=$(pwd)
check_python $cwd
install_hmmer $cwd
download_cobalt
download_databases
create_parameters $cwd
create_example_mutations
download_libraries
echo 'All done! Remember you need to also download MODELLER (https://salilab.org/modeller/download_installation.html)'
