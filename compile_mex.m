cd utilities
mex CXXOPTIMFLAGS="-march=corei7-avx -mtune=corei7-avx -O3 -DNDEBUG" ep2d_full.cpp
mex CXXOPTIMFLAGS="-march=corei7-avx -mtune=corei7-avx -O3 -DNDEBUG" ep2d_simple.cpp
mex CXXOPTIMFLAGS="-march=corei7-avx -mtune=corei7-avx -O3 -DNDEBUG" ep2d_guided.cpp
mex CXXOPTIMFLAGS="-march=corei7-avx -mtune=corei7-avx -O3 -DNDEBUG" ep2d_simple_9.cpp
mex CXXOPTIMFLAGS="-march=corei7-avx -mtune=corei7-avx -O3 -DNDEBUG" ep2d_simple_15.cpp
mex CXXOPTIMFLAGS="-march=corei7-avx -mtune=corei7-avx -O3 -DNDEBUG" ep2d_simple_31.cpp
cd ..