mpicc -std=c++14 integration_parallel.cpp
mpirun -n 3 ./a.out
convert out.ppm out.png
