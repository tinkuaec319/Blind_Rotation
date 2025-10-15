Implementation of Blind rotation with reduced key-switch computation and Key merging 
====================================================================================

## Blind rotation with reduced key-switch computation and Key merging for Bootstrapping FHEs
The library contains the implementation of the fully homorphic encryption schemes presented in the paper [Blind rotation with reduced key-switch computation and Key merging ]
by using [OpenFHE_v1.1.1](https://github.com/openfheorg/openfhe-development/releases/tag/v1.1.1) updated by Xiang et al. in [Fast Blind Rotation for Bootstrapping FHEs](https://eprint.iacr.org/2023/1564) 

For optimal performance, we also employ an approximate gadget decomposition and provide improved parameter sets as in LMKCDEY (see `binfhecontext.cpp`).
### Requirements
A C++ compiler, the NTL libraries.

## Run the code
1. Configure, build and compile the project.
```
mkdir build
cd build
cmake -DWITH_NTL=ON .. 
cmake -DCMAKE_POSITION_INDEPENDENT_CODE=ON ..
make 
```
2. In this work we have added two new implementations of blind rotation viz., blind rotation technique proposed by Li et al. in [Faster ntru-based bootstrapping in less than 4 ms.] 
and the proposed blind rotation technique presented in this paper.
 
To run the desired blind rotation technique please uncomment one of the line number 1012 and/or 1015 and/or 1018 based on the desired blind rotation computation needed
accordingly in file [binfhe-base-scheme.cpp] present at CHIFHE_RGSW_KeyMerge/src/binfhe/lib/


2. Run the `boolean-SS` program in `build/bin/examples/binfhe`.
   
Experiment results carried out on a commodity laptop with a $12^{th}$ Gen Intel(R) Core $i5$-12450H  $\times 12$ running at 400MHz-$4400$GHz on a Linux 22.04.5 LTS machine.
 ```
//NTRU based 
1000 times of Xiang et al. blind rotation using parameters P128G for gate bootstrapping:   77 ms

//RGSW based
1000 times of AP blind rotation using parameters STD128 for gate bootstrapping:            105 ms  
1000 times of CGGI blind rotation using parameters STD128 for gate bootstrapping:          80  ms  
1000 times of Lee et al. blind rotation using parameters STD128 for gate bootstrapping:    86  ms  
1000 times of Li et al. blind rotation using parameters  STD128 for gate bootstrapping:    85  ms  
1000 times of Proposed blind rotation using parameters   STD128 for gate bootstrapping:    60 ms   
```

3. We recommend using the following CMake command-line configuration for best performance.
```
cmake -DWITH_NTL=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DNATIVE_SIZE=32 -DWITH_NATIVEOPT=ON -DCMAKE_C_COMPILER=clang-12 -DCMAKE_CXX_COMPILER=clang++-12 -DWITH_OPENMP=OFF -DCMAKE_C_FLAGS="-pthread" -DCMAKE_CXX_FLAGS="-pthread" .. 

//NTRU based 
1000 times of Xiang et al. blind rotation using parameters P128G for gate bootstrapping:   35 ms

//RGSW based
1000 times of AP blind rotation using parameters STD128 for gate bootstrapping:            74 ms
1000 times of CGGI blind rotation using parameters STD128 for gate bootstrapping:          46 ms
1000 times of Lee et al. blind rotation using parameters STD128 for gate bootstrapping:    53 ms
1000 times of Li et al. blind rotation using parameters  STD128 for gate bootstrapping:    52 ms
1000 times of Proposed blind rotation using parameters   STD128 for gate bootstrapping:    38 ms 