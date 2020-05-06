#!/bin/bash
# https://www.aithercfd.com/2016/12/03/using-travisci.html
# Updated to use openMPI 4.0
mkdir -p ~/openmpi
if [ "$TRAVIS_OS_NAME" == "osx" ]; then
	cd ~/openmpi
	echo "ls /usr/local"
	ls /usr/local
	echo "ls /usr/local/Cellar"
	ls /usr/local/Cellar
	echo "ls /usr/local/Cellar/4.0.3"
	ls /usr/local/Cellar/4.0.3
	if [ -d "/usr/local/Cellar" ]; then
		echo "Using cached OpenMPI on " $TRAVIS_OS_NAME
		ln -s /usr/local/bin /usr/local/Cellar/4.0.3/bin
		ln -s /usr/local/lib /usr/local/Cellar/4.0.3/lib
		ln -s /usr/local/include /usr/local/Cellar/4.0.3/include
    else
		echo "Installing OpenMPI with homebrew on " $TRAVIS_OS_NAME
		rm -rf ~/openmpi/*
		NUM_CORES=$(sysctl -n hw.ncpu)
		HOMEBREW_MAKE_JOBS=$NUM_CORES brew install open-mpi
	fi
else
	cd ~/openmpi
	if [ -f "$HOME/openmpi/bin/mpirun" ] && [ -f "$HOME/openmpi-4.0.3/config.log" ]; then
		echo "Using cached OpenMPI on " $TRAVIS_OS_NAME
		echo "Configuring OpenMPI"
		cd ~/openmpi-4.0.3
		./configure --prefix=~/openmpi CC=$C_COMPILER CXX=$CXX_COMPILER &> openmpi.configure
	else
		echo "Downloading OpenMPI Source on " $TRAVIS_OS_NAME
		rm -rf ~/openmpi
		rm -rf ~/openmpi-4.0.3
		cd ~
		wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.3.tar.gz
		tar zxf openmpi-4.0.3.tar.gz
		echo "Configuring and building OpenMPI"
		mkdir -p ~/openmpi
		cd ~/openmpi-4.0.3
		./configure --prefix=$HOME/openmpi CC=$C_COMPILER CXX=$CXX_COMPILER &> openmpi.configure
		make -j &> openmpi.make
		make install &> openmpi.install
	fi
	test -n $CC && unset CC
	test -n $CXX && unset CXX
fi
cd $TRAVIS_BUILD_DIR
