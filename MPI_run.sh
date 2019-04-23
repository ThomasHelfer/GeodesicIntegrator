mpicc -std=c++14 integration_parallel.cpp
mpirun -n 32 ./a.out
for f in Image*.ppm; do
  convert ./"$f" ./"${f%.ppm}.png"
done
