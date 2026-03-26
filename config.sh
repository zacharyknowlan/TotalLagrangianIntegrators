cmake -B build -S . \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCMAKE_CXX_COMPILER=`which g++` \
    -DMFEM_ROOT="../../MFEM/mfem-4.8/build/install/"
cmake --build build/ -j 8
