#!/bin/bash

EPIPELINE_VERSION=
# where is the pipeline stuff
EPIPELINE_SRC_DIR=.
# where the pipeline is going to be installed
EPIPELINE_INS_DIR=



if [ "`uname 2>/dev/null`-" = "Linux-" ]; then
    OS=linux
else
    # there are only linux and mac os :)
    # lets try to also install in mac...
    OS=mac
fi


############################################################
function pprint_msg {
    echo "* $*"
}


function usage {
    echo "Usage: install.sh  -s cloned_directory  -a dir   ";
    echo " -s dir : toplevel source directory";
    echo " -a dir : install  to directory 'dir'. If dir is not given the value of EPIPELINE_DIR will be used (if available).";
    echo " -x software: install/update software.";
}


function python_2_7_installed {
    #check if python 2.7 is installed
    #
    pprint_msg "Check python version... (2.7 required)"
    min=$(python -c "import sys; print (sys.version_info[:])[1]")
    maj=$(python -c "import sys; print (sys.version_info[:])[0]")
    if [[ $maj != 2 ]] || [[ $min != 7 ]] ; then
	pprint_msg "You need Python2.7 to run this pipeline."
	exit 1
    else
	pprint_msg "OK."
    fi
}

# include whatever is necessary to install and compile the programs
function check_dependencies {
    DEVEL_LIBRARIES_REQUIRED="python-devel python cmake python-pip"
    MISSING=0
    pprint_msg "If installation fails then please check if the following libraries are installed:"
    pprint_msg "$DEVEL_LIBRARIES_REQUIRED"
    # Binaries that should be available
    # make is required to...compile make
    BINARIES="wget  git which make bzip2 unzip cmake pip"
    pprint_msg "Checking dependencies..."
    for bin in $BINARIES; do
	PATH2BIN=`which $bin 2> /dev/null`
	if [ "$PATH2BIN-" == "-" ]; then
	    pprint_msg " $bin not found!"
 	    MISSING=1
	else
	    pprint_msg " $bin found: $PATH2BIN"
	fi
    done
    pprint_msg "Checking dependencies...done."
    if [ $MISSING == 1 ]; then
	pprint_msg "ERROR: Unable to proceed"
	exit 1
    fi

}

function download {
    URL=$1
    OFILE=$2
#    if [ $OS == "linux" ]; then	    
     wget  --no-check-certificate -c -nc -L "$URL" -O $OFILE
#    else
	# 
#	curl $URL
#    fi

}

function get_fullpath {

    F=$1
    if [ $OS == "linux" ]; then	    
	readlink -f $F
    else
	# 
	greadlink -f $F
    fi    
}
################################################################
# Software to install 

MAKE_VERSION=4.1
MAKE_FILE=make-${MAKE_VERSION}.tar.gz
MAKE_URL=http://ftp.gnu.org/gnu/make/$MAKE_FILE

function make_install {
    pprint_msg "Installing make..."
    download $MAKE_URL $MAKE_FILE
    tar xzvf $MAKE_FILE
    pushd `echo $MAKE_FILE|sed "s/.tar.gz//"`
    ./configure 
    make -j $J
    cp make $EPIPELINE_DIR/bin
    popd
    pprint_msg "Installing make...done."
}

# since the overall pipeline has been tested with python 2.7 and
# imports many python modules available in Anaconda
# (https://www.continuum.io/downloads), I would suggest to install it
# first
ANACONDA_VERSION=2.3.0
# TODO: add powered by Anaconda to the github repo and initial output of the pipeline
function anaconda_install {
    
    pprint_msg "Installing anaconda..."
    file=Anaconda-$ANACONDA_VERSION-Linux-x86_64.sh
    if [ "$OS" == "mac" ]; then
	file=Anaconda-$ANACONDA_VERSION-MacOSX-x86_64.sh
    fi
    #MAC OS X
    download https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/$file $file
    chmod +x $file
    mkdir -p $EPIPELINE_DIR/anaconda
    ./$file -p $EPIPELINE_DIR/anaconda -b -f
    pprint_msg "Installing anaconda...done."
}
#export PATH="/path/to/anaconda/bin:$PATH"


#####################################################
# limix (https://github.com/PMBio/limix)
function limix_install {
# limix dependancies:
#
# python libraries (should have been already installed by anaconda):
# scipy, numpy, pandas, cython
# Swig:
# swig 2.0 or higher (only required if you need to recompile C++ interfaces)
#
# GCC:
# gcc version > 4.2.1
# LIMIX reccomended installation
    pprint_msg "Installing Limix..."
    pip install limix --user
    pprint_msg "Installing Limix...done."
}

#####################################################################
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
#Install PEER from source (recommended)
function peer_install {
    pprint_msg "Installing PEER..."
    # Install dependencies
    SWIG_VERSION=3.0.7
    SWIG_FILE=swig-$SWIG_VERSION.tar.gz
    SWIG_URL=http://prdownloads.sourceforge.net/swig/$SWIG_FILE
    download $SWIG_URL $SWIG_FILE
    tar xzvf $SWIG_FILE
    pushd swig-$SWIG_VERSION
    ./configure --prefix $EPIPELINE_DIR
    make
    make install
    popd
    #
    git clone git@github.com:PMBio/peer.git
    #create a target directory for the build:
    mkdir -p peer/build
    pushd peer/build 
    if [ "$OS" == "mac" ]; then
	#By default, peer is not built as an universal binary on OS X. 
	#This can be easily changed by adding the following configuration parameter to CMake: -DBUILD_UNIVERSAL=1.
	cmake ./.. -DBUILD_UNIVERSAL=1 -DCMAKE_INSTALL_PREFIX:PATH=$EPIPELINE_DIR
    else
	cmake ./.. -DCMAKE_INSTALL_PREFIX:PATH=$EPIPELINE_DIR
    fi

    #Finally, compile and install PEER
    make
    make install
    popd
    pprint_msg "Installing PEER...done."
}

function fix_paths {
    # TODO: FIX THIS PATH
    limix_path=$EPIPELINE_DIR/
    peer_path=$EPIPELINE_DIR/
    # limix
    pythonfiles2fix="geno_preprocessing.py"
    for f in $pythonfiles2fix; do
	pprint_msg "Fixing limix_path in $f..."
	sed -i "s|^limix_path=.*|limix_path='$limix_path'|" $EPIPELINE_DIR/scripts/$f
	pprint_msg "Fixing limix_path in $f...done."	
    done
}


function epipeline_install {

    pprint_msg "Installing epipeline..."
    mkdir -p  $EPIPELINE_DIR/scripts
    mkdir -p  $EPIPELINE_DIR/bin
    mkdir -p  $EPIPELINE_DIR/include
    mkdir -p  $EPIPELINE_DIR/lib
    cp -a $EPIPELINE_SRC_DIR/scripts/* $EPIPELINE_DIR/scripts
    #cp -a $EPIPELINE_SRC_DIR/bin/* $EPIPELINE_DIR/bin
    #fix_paths
    pprint_msg "Installing epipeline...done."
}

#eqtl

####################################################################
# 
install=everything

while getopts "s:a:x:h"  Option
do
    case $Option in
# update/reinstall
        a ) install=everything;EPIPELINE_INS_DIR=$OPTARG;;# send all output to a log file
        s ) EPIPELINE_SRC_DIR=$OPTARG;;# send all output to a log file
	x ) install=software_install;EPIPELINE_INS_DIR=$EPIPELINE_DIR;SOFTWARE=$OPTARG;;
        h ) usage; exit;;
    esac
done


if [ "$EPIPELINE_INS_DIR-" != "-" ]; then
    export EPIPELINE_DIR=$EPIPELINE_INS_DIR
fi


if [ "$EPIPELINE_DIR-" == "-" ]; then
    echo ERROR: EPIPELINE_DIR directory not defined. > /dev/stderr
    usage
    exit 1    
fi

if [ "$EPIPELINE_SRC_DIR-" = "-" ]; then
    usage
    exit 1    
fi

# Check if env is available
DEF_ENV="/usr/bin/env"
ENV_FP=$DEF_ENV
if [ -x $ENV_FP ]; then
    pprint_msg "env found in $ENV_FP"
else
    ENV_FP="/bin/env"
    if [ -x $ENV_FP ]; then
	pprint_msg "env found in $ENV_FP"
    else
	echo "ERROR: env command not found - please ensure that it is in the PATH"
	exit 1
    fi
fi
# print system info
uname -a
# Check dependencies
check_dependencies
# Full path
EPIPELINE_DIR=`get_fullpath "$EPIPELINE_DIR"`
EPIPELINE_SRC_DIR=`get_fullpath "$EPIPELINE_SRC_DIR"`
pprint_msg "EPIPELINE_DIR=$EPIPELINE_DIR"
pprint_msg "EPIPELINE_SRC_DIR=$EPIPELINE_SRC_DIR"
pprint_msg "PATH=$PATH"


SETUP_FILE=$EPIPELINE_DIR/eqtl_setup.sh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$EPIPELINE_DIR/lib

# Python - where to install the packages
export PYTHONUSERBASE=$EPIPELINE_DIR/anaconda
mkdir -p $PYTHONUSERBASE/lib/python2.7/site-packages

# R
mkdir -p $IRAP_DIR/Rlibs

# 
TMP_DIR=$EPIPELINE_DIR/tmp
mkdir -p $TMP_DIR
pushd $TMP_DIR
pprint_msg "Cleaning up $TMP_DIR..."
rm -rf *
pprint_msg "Cleaning up $TMP_DIR...done."

# useful for debugging
if [ "$install" == "software_install" ]; then
    ${SOFTWARE}_install
    exit 0
fi

# install everything
check_dependencies
python_2_7_installed

# generate the file with the environment variables
cat <<EOF > $SETUP_FILE
export EPIPELINE_DIR=$EPIPELINE_DIR
export PATH=\$EPIPELINE_DIR/bin:\$EPIPELINE_DIR/scripts:\$EPIPELINE_DIR/anaconda/bin:\$PATH
export LD_LIBRARY_PATH=\$EPIPELINE_DIR/lib:\$LD_LIBRARY_PATH:/usr/local/lib
export CFLAGS="-I\$EPIPELINE_DIR/include \$CFLAGS"
export R_LIBS_USER=$EPIPELINE_DIR/Rlibs
export CXXFLAGS="-I\$EPIPELINE_DIR/include -L\$IRAP_DIR/lib \$CXXFLAGS"
export PYTHONUSERBASE=\$EPIPELINE_DIR/anaconda
EOF
source $SETUP_FILE
set -e 
epipeline_install

make_install
anaconda_install
limix_install
peer_install
popd
exit 0

