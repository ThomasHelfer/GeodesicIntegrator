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


# define H 100

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

	// ============= MPI INIT ================
 	MPI_Init(NULL,NULL);
 	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	size_per_task = (int)ceil(H*H/(double)numtasks);
	// ============= DEFINITIONS ==============

 	char name_render[] = "out.ppm";	
	render_black_hole<Black_Hole> rend;
	const double size_x = 20;
	const double size_y = 20;
	double alpha = 0;
	double t1,t2;	

	int* red_local   = (int *)malloc(sizeof(int) * size_per_task);
	int* green_local = (int *)malloc(sizeof(int) * size_per_task);
	int* blue_local  = (int *)malloc(sizeof(int) * size_per_task);
	int* red   = (int *)malloc(sizeof(int) * size_per_task * numtasks);
	int* green = (int *)malloc(sizeof(int) * size_per_task * numtasks);
	int* blue  = (int *)malloc(sizeof(int) * size_per_task * numtasks);

	// ========== Main CODE =================
	int i = 0;	
	for(alpha = 0; alpha < M_PI/2.; alpha += M_PI/(2.*50.)){

		// ------ Create Filename -----------------
		t1 = MPI_Wtime();
		string imgname="Image_";
		char cbuff[20];
		sprintf (cbuff, "%03d", i);
		imgname.append(cbuff);
		imgname.append(".ppm");	


		// ------ Making Picture (all the hard work is here) --- 
		const Vec3 center(100.0*cos(alpha), 0.0, 100*sin(alpha), 0.0, -1.0*cos(alpha), 0.0, -1.0*sin(alpha), 1.0); // Set Up center and viewing angle
		rend.picture(red_local, green_local, blue_local, center, size_x,size_y,alpha, rank*size_per_task , size_per_task*(rank+1));	
		if(rank==0){ cout << imgname <<  endl;}

		// Draw A circle 
		rend.render_circle(red_local, green_local, blue_local, size_x, size_y, rank*size_per_task , size_per_task*(rank+1));

		// --------- MPI communication ------------
		MPI_Gather(red_local, size_per_task, MPI_INT, red , size_per_task, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Gather(green_local, size_per_task, MPI_INT, green , size_per_task, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Gather(blue_local, size_per_task, MPI_INT, blue , size_per_task, MPI_INT, 0, MPI_COMM_WORLD);

		// ---------- Write out data --------------
		if(rank==0){rend.render(red, green, blue, imgname);}


		t2 = MPI_Wtime();
                if (rank ==0) {cout << "duration per frame " << t2-t1 << "sec" << endl; }
			
		i++;
	}
	// =========== Clean up =================

	free(red_local);
	free(green_local);
	free(blue_local);
 	free(red);
 	free(green);
 	free(blue);

	MPI_Finalize();

	return 0;
}
