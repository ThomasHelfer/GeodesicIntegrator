/*
 * compiled with gcc 5.4:
 * g++-mp-5 -std=c++14 -o rk4 rk4.cc
 *
 */
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <algorithm> 
#include <ctime>
#include <mpi.h>

#include "tensor.hpp"
#include "schwarzschild.hpp"
#include "rk4.hpp"
#include "render.hpp"

using namespace std;

#define FOR1(IDX) for (int IDX = 0; IDX < 4; ++IDX)
#define FOR2(IDX1,IDX2) FOR1(IDX1) FOR1(IDX2)
#define FOR3(IDX1,IDX2,IDX3) FOR2(IDX1,IDX2) FOR1(IDX3)
#define FOR4(IDX1,IDX2,IDX3,IDX4) FOR2(IDX1,IDX2) FOR2(IDX3,IDX4)


int main(void)
{
	int numtasks, rank, sendcount, recvcount, source;
	int size_per_task;


 	MPI_Init(NULL,NULL);
 	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

 	char name_render[] = "out.ppm";	
	render_black_hole<Black_Hole> rend;
	const double size_x = 20;
	const double size_y = 20;
	
	size_per_task = (int)ceil(H*H/(double)numtasks);
	cout << size_per_task << endl;

	int* pic_local = (int *)malloc(sizeof(int) * size_per_task);

	int* pic = (int *)malloc(sizeof(int) * size_per_task * numtasks);

	const Vec3 center(100.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0);

	rend.picture(pic_local,center, size_x,size_y, rank*size_per_task , size_per_task*(rank+1));	

	MPI_Gather(pic_local, size_per_task, MPI_INT, pic , size_per_task, MPI_INT, 0, MPI_COMM_WORLD);

	rend.render(pic,name_render);
/*
	free(pic_local);
	free(pic);
*/
	 return 0;
}
