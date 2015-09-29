#!/bin/bash

# since the overall pipeline has been tested with python 2.7 and imports many python modules available in Anaconda (https://www.continuum.io/downloads), I would suggest to install it first

#MAC OS X
wget https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda-2.3.0-MacOSX-x86_64.sh
 
bash Anaconda-2.3.0-MacOSX-x86_64.sh

#LINUX 
wget https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda-2.3.0-Linux-x86_64.sh

bash Anaconda-2.3.0-Linux-x86_64.sh

#export PATH="/path/to/anaconda/bin:$PATH"

#check if python 2.7 is installed

#
echo ""
echo "Check python version... (2.7 required)"
min=$(python -c "import sys; print (sys.version_info[:])[1]")
maj=$(python -c "import sys; print (sys.version_info[:])[0]")
if [[ $maj != 2 ]] || [[ $min != 7 ]]
then
	echo "You need Python2.7 to run this pipeline."
	exit 1
else
	echo "OK."
	echo ""
fi

##################################################################################################################
# limix (https://github.com/PMBio/limix)
#
# limix dependancies:
#
# python libraries (should anaconda installation went well, no need to install them):
# scipy, numpy, pandas, cython
#
# Swig:
# swig 2.0 or higher (only required if you need to recompile C++ interfaces)
#
# GCC:
# gcc version > 4.2.1
#
# LIMIX reccomanded installation
pip install limix

#################################################################################################################
#
# Peer (https://github.com/PMBio/peer/wiki/Installation-instructions)
#
# peer Dependencies:
#
# For PYTHON package (used in the pipeline) of peer:
#
# Python 2.5 or later
# scipy
# numpy
#
# For FULL source package installation:
#
# CMake 2.8 or later
# SWIG 2.0 or later
#

#Install PEER from source (reccomanded)
 
git clone git@github.com:PMBio/peer.git
#create a target directory for the build:
cd peer 
mkdir build

#Then go to the build directory and call CMake

cd build
cmake ./..

#Finally, compile and install PEER
make
make install

#WARNING: By default, peer is not built as an universal binary on OS X. 
#This can be easily changed by adding the following configuration parameter to CMake: -DBUILD_UNIVERSAL=1.
# This will cause peer to be built as an Intel 32/64 bit binary. 
#For more control over the included architectures, edit the file CMakeLists.txt in the root directory of the project.


limix_path=
#
pythonfiles2fix=scripts/geno_preprocessing.py
for f in $pythonfiles2fix; do
    sed -i "s/^limix_path=.*/limix_path=$limix_path/" $f
done

#eqtl


