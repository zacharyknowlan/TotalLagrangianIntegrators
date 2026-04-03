dir="../MFEM/"
mkdir -p "$dir"/

export CFLAGS="$CFLAGS -Wno-misleading-indentation" # Suppress compiler warnings for METIS
wget https://github.com/mfem/tpls/raw/refs/heads/gh-pages/metis-5.1.0.tar.gz -P "$dir"/
tar -xzvf "$dir"/metis-5.1.0.tar.gz -C "$dir"/
cmake -S "$dir"/metis-5.1.0/ -B "$dir"/metis-5.1.0/build/ \
    -DCMAKE_INSTALL_PREFIX="$dir"/metis-5.1.0/build/install/ \
    -DCMAKE_BUILD_TYPE=Release \
    -DGKLIB_PATH="$dir"/metis-5.1.0/GKlib/ \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_C_COMPILER=`which gcc`
cmake --build "$dir"/metis-5.1.0/build/ -j 8 --target install

wget https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v7.11.0.tar.gz -P "$dir"/
tar -xzvf "$dir"/v7.11.0.tar.gz -C "$dir"/
cmake -S "$dir"/SuiteSparse-7.11.0/ -B "$dir"/SuiteSparse-7.11.0/build/ \
    -DCMAKE_INSTALL_PREFIX="$dir"/SuiteSparse-7.11.0/build/install/ \
    -DCMAKE_C_COMPILER=`which gcc` \
    -DCMAKE_CXX_COMPILER=`which g++` \
    -DCMAKE_Fortran_COMPILER=`which gfortran` \
    -DSUITESPARSE_USE_CUDA=OFF \
    -DBUILD_SHARED_LIBS=ON \
    -DBUILD_STATIC_LIBS=OFF \
    -DSUITESPARSE_ENABLE_PROJECTS="suitesparse_config;amd;btf;camd;ccolamd;colamd;cholmod;cxsparse;klu;umfpack;"
cmake --build "$dir"/SuiteSparse-7.11.0/build/ -j 8 --target install

wget https://github.com/mfem/mfem/archive/refs/tags/v4.8.tar.gz -P "$dir"/
tar -xzvf "$dir"/v4.8.tar.gz -C "$dir"/
cmake -S "$dir"/mfem-4.8/ -B "$dir"/mfem-4.8/build/ \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$dir"/mfem-4.8/build/install/ \
    -DCMAKE_CXX_COMPILER=`which g++` \
    -DMFEM_ENABLE_TESTING=OFF \
    -DMFEM_USE_MPI=OFF \
    -DMFEM_USE_CUDA=OFF \
    -DMFEM_USE_METIS_5=YES \
    -DMETIS_DIR="$dir"/metis-5.1.0/build/install/ \
    -DMFEM_USE_SUITESPARSE=ON \
    -DSuiteSparse_UMFPACK_LIBRARY="$dir"/SuiteSparse-7.11.0/build/install/lib/libumfpack.so \
    -DSuiteSparse_UMFPACK_INCLUDE_DIR="$dir"/SuiteSparse-7.11.0/build/install/include/suitesparse/ \
    -DSuiteSparse_KLU_LIBRARY="$dir"/SuiteSparse-7.11.0/build/install/lib/libklu.so \
    -DSuiteSparse_KLU_INCLUDE_DIR="$dir"/SuiteSparse-7.11.0/build/install/include/suitesparse/ \
    -DSuiteSparse_AMD_LIBRARY="$dir"/SuiteSparse-7.11.0/build/install/lib/libamd.so \
    -DSuiteSparse_AMD_INCLUDE_DIR="$dir"/SuiteSparse-7.11.0/build/install/include/suitesparse/ \
    -DSuiteSparse_BTF_LIBRARY="$dir"/SuiteSparse-7.11.0/build/install/lib/libbtf.so \
    -DSuiteSparse_BTF_INCLUDE_DIR="$dir"/SuiteSparse-7.11.0/build/install/include/suitesparse/ \
    -DSuiteSparse_CHOLMOD_LIBRARY="$dir"/SuiteSparse-7.11.0/build/install/lib/libcholmod.so \
    -DSuiteSparse_CHOLMOD_INCLUDE_DIR="$dir"/SuiteSparse-7.11.0/build/install/include/suitesparse/ \
    -DSuiteSparse_COLAMD_LIBRARY="$dir"/SuiteSparse-7.11.0/build/install/lib/libcolamd.so \
    -DSuiteSparse_COLAMD_INCLUDE_DIR="$dir"/SuiteSparse-7.11.0/build/install/include/suitesparse/ \
    -DSuiteSparse_CAMD_LIBRARY="$dir"/SuiteSparse-7.11.0/build/install/lib/libcamd.so \
    -DSuiteSparse_CAMD_INCLUDE_DIR="$dir"/SuiteSparse-7.11.0/build/install/include/suitesparse/ \
    -DSuiteSparse_CCOLAMD_LIBRARY="$dir"/SuiteSparse-7.11.0/build/install/lib/libccolamd.so \
    -DSuiteSparse_CCOLAMD_INCLUDE_DIR="$dir"/SuiteSparse-7.11.0/build/install/include/suitesparse/ \
    -DSuiteSparse_config_LIBRARY="$dir"/SuiteSparse-7.11.0/build/install/lib/libsuitesparseconfig.so \
    -DSuiteSparse_config_INCLUDE_DIR="$dir"/SuiteSparse-7.11.0/build/install/include/suitesparse/
cmake --build "$dir"/mfem-4.8/build/ -j 8 --target install
