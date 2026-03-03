// ---------------
//   Dongzhuo Li
//   March, 2015
// ---------------



#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <math.h>
#include "par_model_source_receiver.h"
#include "eikonal.h"
#include "raytracing.h"
#include <omp.h>
#include <time.h>

using namespace std;

int main() {

	printf ("---------- Ray Tracing with a linearized eikonal equation--------\n");
	printf ("\n");
	printf ("---------- reading input parameters:\n");

	int dims[3] = {0, 0, 0};
	// parameters in Par_file
	int dimension, nx, ny, nz, invnx, invny, invnz, lengthG, \
        ixsrcn, iysrcn, izsrcn, ixrecn,iyrecn,izrecn,iglob,i1,j1,k1,ii,jj,kk,ijk;
	double dx, dy, dz, dzxy, xmin, ymin, zmin, invdx, invdy, invdz, invxmin, invymin, invzmin;
	string velModel_file, attModel_file, source_file, receiver_file;
	double *velModel, *attModel, *T, *TS;
	bool usesecond = false;
	bool usecross = false;
	ofstream timeMap, rayPointsOutPut, time_tstar;
	vector <double *> SourcePoints;
	vector <double *> AllTravelTimeMaps;
	map < int, vector< vector <double> > > ReceiverPoints;


	// ------------------------ Read Parameters & Load Model ----------------------
	read_parafile(dimension, nx, ny, nz, dx, dy, dz, xmin, ymin, zmin, invnx, invny, invnz, \
	              invdx, invdy, invdz, invxmin, invymin, invzmin, velModel_file, attModel_file, source_file, receiver_file);

	dzxy = dz;

	if (dimension == 2) {
		dims[0] = nz; dims[1] = nx;
		lengthG = (invnx - 1) * (invnz - 1);
		velModel = (double *)malloc(dims[0] * dims[1] * sizeof(double));
		T        = (double *)malloc(dims[0] * dims[1] * sizeof(double));
		attModel = (double *)malloc(dims[0] * dims[1] * sizeof(double));  // -wdd
		TS       = (double *)malloc(dims[0] * dims[1] * sizeof(double));  // -wdd
	}
	if (dimension == 3) {
		dims[0] = nz; dims[1] = nx; dims[2] = ny;
		lengthG = (invnx - 1) * (invny - 1) * (invnz - 1);
		velModel = (double *)malloc(dims[0] * dims[1] * dims[2] * sizeof(double));
		T        = (double *)malloc(dims[0] * dims[1] * dims[2] * sizeof(double));
		attModel = (double *)malloc(dims[0] * dims[1] * dims[2] * sizeof(double));  // -wdd
		TS       = (double *)malloc(dims[0] * dims[1] * dims[2] * sizeof(double));  // -wdd		
	}

	velModelLoad(dimension, velModel_file, velModel);
	attModelLoad(dimension, attModel_file, attModel);  // -wdd

	sourceLoad(dimension, source_file, SourcePoints);

	receiverLoad(dimension, receiver_file, ReceiverPoints);
	printf ("---------- finished reading input parameters!\n");


	// --------------------------- Time Map Calculations ----------------------------
	// output time map
	//timeMap.open("./Out/timeMap.txt");
	time_tstar.open("./Out/time_tstar_table");
	printf ("\n");
	printf ("----- Traveltime computation using MSFM3D: Multistencil Fast Marching Method -----\n");
	printf ("\n");
//	printf (":::::::::::::::::::::::::::::: parallel section ::::::::::::::::::::::::::::::\n");
//	printf ("\n");


	//***************** OpenMP parallel region **************************************//
	rayPointsOutPut.open("./Out/raypoints.txt");
	double *temptt;
	double *tempts;  // -wdd
	//	double *rayPoints;

	vector< vector<double> > rayPoints;

// omp_lock_t lck;
// omp_init_lock(&lck);

	// start clock
	double start = omp_get_wtime();

  	#pragma omp parallel 
	{
	    #pragma omp for ordered private(izsrcn, ixsrcn, iysrcn, iglob, temptt, rayPoints) schedule(static,1)
		for (int i = 0; i < SourcePoints.size(); i++) 
		{
			temptt  = (double *)malloc(dims[0] * dims[1] * dims[2] * sizeof(double) );
			tempts  = (double *)malloc(dims[0] * dims[1] * dims[2] * sizeof(double) );   // -wdd
			int ID = omp_get_thread_num();
			if (dimension == 3)
			{
				izsrcn = round((SourcePoints.at(i)[1] - zmin) / dz);
				ixsrcn = round((SourcePoints.at(i)[2] - xmin) / dx);
				iysrcn = round((SourcePoints.at(i)[3] - ymin) / dy);
				#pragma omp critical 
				{
					cout << "parallel, thread id #" << ID << ". time map computing source #" << i + 1 << endl;			
					cout << "izsrcn = " << izsrcn << ". ixsrcn = " << ixsrcn << ". iysrcn= " << iysrcn << endl;
				}
				//msfm3dCpp(temptt, velModel, izsrcn, ixsrcn, iysrcn, dims, usesecond, usecross);
				//msfm3dCpp(temptt, tempts, velModel, attModel, izsrcn, ixsrcn, iysrcn, dims, usesecond, usecross);
				msfm3dCpp(temptt, tempts, velModel, attModel, izsrcn, ixsrcn, iysrcn, dzxy, dims, usesecond, usecross);
				#pragma omp ordered 
			// ------------------------------- Save the time and tstar Table ----------------------------------
                                time_tstar << int(nx) << " " << int(ny) << " " << int(nz) << "\n";
				for (int jj = 0; jj < int(ny); ++jj) {
				    for (int ii = 0; ii < int(nx); ++ii) {
           				 for (int kk = 0; kk < int(nz); ++kk) {
				              //ijk = mindex3(k, i, j, nz, nx);
					      ijk = jj * (nz * nx) + ii * nz + kk; 	      	 
					      time_tstar << temptt[ijk] << " " << tempts[ijk] << endl;
           				 }
       				    }
   				 }
			// ------------------------------- Ray Tracing ----------------------------------
				for (int j = 0; j < ReceiverPoints[int(SourcePoints.at(i)[0])].size(); j++) 
				{
	            				int ID2 = omp_get_thread_num();
                    izrecn = round((ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(1) - zmin) / dz);
                    ixrecn = round((ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(2) - xmin) / dx);
                    iyrecn = round((ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(3) - ymin) / dy);
                    iglob  = dims[0]*dims[1]*iyrecn+dims[0]*ixrecn+izrecn;			
				//	raytracing(dimension, nx, ny, nz, dx, dy, dz, xmin, ymin, zmin, \
				//          invnx, invny, invnz, invdx, invdy, invdz, invxmin, invymin, invzmin, \
				//          SourcePoints.at(i)[2], SourcePoints.at(i)[3], SourcePoints.at(i)[1], \
				//          ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(2), \
				//       ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(3), \
				//         ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(1), \
				//          temptt, rayPoints);
					//cout << "Thread id #" << ID2 << ", raytracing: " << "Source " << i + 1 << ", " << "Receiver " << j + 1 << endl; 
			                rayPointsOutPut << int(SourcePoints.at(i)[0]) << " " << ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(0) \
				              << " " << rayPoints.size() <<" "<< dx*temptt[iglob] << endl;
					for (int k = 0; k < rayPoints.size(); k++)
					{
						rayPointsOutPut << rayPoints.at(k).at(0) << " " << rayPoints.at(k).at(1) << " " << rayPoints.at(k).at(2) << endl;
					}
					rayPoints.clear();
				}
			}
                        free(temptt); free(tempts);
		}
	}
	//***************** end of OpenMP parallel region **************************************//

	// output time map
	// timeMap.close();

	// -------------------------------- Clean Up ------------------------------------
	rayPointsOutPut.close();
	time_tstar.close();
	free(velModel); free(attModel);
	for (int i = 0; i < SourcePoints.size(); i++) {
		delete SourcePoints.at(i);
	}

// show time:
	printf (":::::::::::::::::::::::::::::: done! ::::::::::::::::::::::::::::::\n");
	double end = omp_get_wtime();
	printf("Time taken:: %.2fs\n", end - start);
	return 0;
}
