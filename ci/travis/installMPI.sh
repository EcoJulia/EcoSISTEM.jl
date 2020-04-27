#!/bin/bash
# https://www.aithercfd.com/2016/12/03/using-travisci.html
# Updated to use openMPI 4.0

if [ "$TRAVIS_OS_NAME" == "osx" ]; then
    cd openmpi
    if [ -f "bin/mpirun" ]; then
	echo "Using cached OpenMPI"
    else
	echo "Installing OpenMPI with homebrew"
			brew install open-mpi
			ln -s /usr/local/bin bin
			ln -s /usr/local/lib lib
			ln -s /usr/local/include include
    fi
else
    if [ -f "openmpi/bin/mpirun" ] && [ -f "openmpi-4.0.3/config.log" ]; then
			echo "Using cached OpenMPI"
			echo "Configuring OpenMPI"
			cd openmpi-4.0.3
			./configure --prefix=$TRAVIS_BUILD_DIR/openmpi CC=$C_COMPILER CXX=$CXX_COMPILER &> openmpi.configure
		else
			echo "Downloading OpenMPI Source"
			wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.3.tar.gz
			tar zxf openmpi-4.0.3.tar.gz
			echo "Configuring and building OpenMPI"
			cd openmpi-4.0.3
			./configure --prefix=$TRAVIS_BUILD_DIR/openmpi CC=$C_COMPILER CXX=$CXX_COMPILER &> openmpi.configure
			make -j4 &> openmpi.make
			make install &> openmpi.install
			cd ..
    fi
    test -n $CC && unset CC
    test -n $CXX && unset CXX
fi
