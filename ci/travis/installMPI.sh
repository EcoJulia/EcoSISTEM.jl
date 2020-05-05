#!/bin/bash
# https://www.aithercfd.com/2016/12/03/using-travisci.html
# Updated to use openMPI 4.0
mkdir -p ~/openmpi
cd ~/openmpi
if [ "$TRAVIS_OS_NAME" == "osx" ]; then
	echo $PWD
	ls *
	if [ -f "./bin/mpirun" ]; then
		echo "Using cached OpenMPI on " $TRAVIS_OS_NAME
    else
		echo "Installing OpenMPI with homebrew on " $TRAVIS_OS_NAME
		NUM_CORES=$(sysctl -n hw.ncpu)
		HOMEBREW_MAKE_JOBS=$NUM_CORES brew install open-mpi
		cp -r /usr/local/Cellar/open-mpi/4.0.3/* .
    fi
	ln -s /usr/local/bin bin
	ln -s /usr/local/lib lib
	ln -s /usr/local/include include
else
	echo "ls ~/openmpi"
	ls ~/openmpi
	echo "ls ~/openmpi-4.0.3"
	ls ~/openmpi-4.0.3
	if [ -f "./bin/mpirun" ] && [ -f "$HOME/openmpi-4.0.3/config.log" ]; then
		echo "Using cached OpenMPI on " $TRAVIS_OS_NAME
		echo "Configuring OpenMPI"
		cd ../openmpi-4.0.3
		./configure --prefix=~/openmpi CC=$C_COMPILER CXX=$CXX_COMPILER &> openmpi.configure
	else
		echo "Downloading OpenMPI Source on " $TRAVIS_OS_NAME
		rm -rf ~/openmpi
		rm -rf ~/openmpi-4.0.3
		wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.3.tar.gz
		tar zxf openmpi-4.0.3.tar.gz
		echo "Configuring and building OpenMPI"
		mv openmpi-4.0.3 ~
		cd ~/openmpi-4.0.3
		./configure --prefix=$HOME/openmpi CC=$C_COMPILER CXX=$CXX_COMPILER &> openmpi.configure
		make -j &> openmpi.make
		make install &> openmpi.install
	fi
	test -n $CC && unset CC
	test -n $CXX && unset CXX
fi
cd $TRAVIS_BUILD_DIR
