cd lammps
rm main.tar.gz
rm -r mlo_lammps-main/
wget https://github.com/moritzobenauer/mlo_lammps/archive/main.tar.gz
tar -zvxf main.tar.gz
cp mlo_lammps-main/patches/* src/
rm -r build/
mkdir build && cd build

# Following:
# https://github.com/PrincetonUniversity/install_lammps/blob/master/01_installing/ins/della/della9_amd_double_prec_aocc_aocl.sh

module purge
module load gcc-toolset/14
module load aocc/5.0.0
module load aocl/aocc/5.0.0
module load openmpi/aocc-5.0.0/4.1.6
FFTW3DIR=/opt/AMD/aocl/aocl-linux-aocc-5.0.0/aocc


cmake3 \
    -D CMAKE_INSTALL_PREFIX=$HOME/.local \
    -D LAMMPS_MACHINE=d9_double_aocc \
    -D ENABLE_TESTING=no \
    -D BUILD_MPI=yes \
    -D BUILD_OMP=yes \
    -D CMAKE_CXX_COMPILER=clang++ \
    -D CMAKE_BUILD_TYPE=Release \
    -D CMAKE_CXX_FLAGS_RELEASE="-Ofast -march=native -DNDEBUG" \
    -D PKG_KSPACE=yes \
    -D FFT=FFTW3 \
    -D FFT_SINGLE=yes \
    -D FFTW3F_INCLUDE_DIR=${FFTW3DIR}/include_LP64 \
    -D FFTW3F_LIBRARY=${FFTW3DIR}/lib_LP64/libfftw3f.so \
    -D PKG_MOLECULE=yes \
    -D PKG_RIGID=yes \
    -D PKG_MC=yes \
    ../cmake

make -j 8
make install

