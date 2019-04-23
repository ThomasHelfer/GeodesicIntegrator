#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>  
#include <iostream>
#include <fstream>
#define SIZE 4
#define W 100

int main() {
int numtasks, rank, sendcount, recvcount, source;
int size_per_task; 

 float* sendbuf = (float *)malloc(sizeof(float) * SIZE); 

 MPI_Init(NULL,NULL);
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

 float rank1[SIZE*SIZE];

  for(int i = 0; i < SIZE; i++ ){sendbuf[i] = rank;}

	
  MPI_Gather(sendbuf, 4, MPI_FLOAT, rank1 , 4, MPI_FLOAT, 0, MPI_COMM_WORLD);
  
  
  if(rank == 0)
   for(int i = 0; i < SIZE*SIZE; i++){
    printf("%f\n",rank1[i]);
   }
  
  free(sendbuf);

  MPI_Finalize();
}

