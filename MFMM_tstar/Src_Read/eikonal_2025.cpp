// ---------------
//   Dongzhuo Li
//   March, 2015
//
//   Modified from the fast marching eikonal solver
//   written by D.Kroon University of Twente (Oct 2010)
// ---------------

#include "math.h"
#include "common.c"

// function prototype
struct CalculateDis2DResult {
    double Tt;
    double TSt;
};
struct CalculateDis3DResult {
    double Tt;
    double TSt;
};

//double CalculateDistance2D(double *T, double *TS, double Fij, double Qij, int *dims, int i, int j, bool usesecond, bool usecross, bool *Frozen);
//CalculateDis2DResult CalculateDistance2D(double *T, double *TS, double Fij, double Qij, int *dims, int i, int j, bool usesecond, bool usecross, bool *Frozen);
CalculateDis2DResult CalculateDistance2D(double *T, double *TS, double Fij, double Qij, int *dims, double dzxy, int i, int j, bool usesecond, bool usecross, bool *Frozen);
//double CalculateDistance3D(double *T, double Fijk, int *dims, int i, int j, int k, bool usesecond, bool usecross, bool *Frozen);
//CalculateDis3DResult CalculateDistance3D(double *T, double *TS, double Fijk, double Qijk, int *dims, int i, int j, int k, bool usesecond, bool usecross, bool *Frozen);
CalculateDis3DResult CalculateDistance3D(double *T, double *TS, double Fijk, double Qijk, int *dims, double dzxy, int i, int j, int k, bool usesecond, bool usecross, bool *Frozen);
double second_derivative(double Txm1, double Txm2, double Txp1, double Txp2);
double second_derivative_tt_ts(double Txm1, double Txm2, double Txp1, double Txp2, double TSxm1, double TSxm2, double TSxp1, double TSxp2);
double averageVelocityAroundSource2D(int bz, int bx, int izsrcn, int ixsrcn, int *dims, double *F);
double averageAttenuationAroundSource2D(int bz, int bx, int izsrcn, int ixsrcn, int *dims, double *Q);  // -wdd
double CalculateTConstant2D(double aveVel, int izsrcn, int ixsrcn, int zloc, int xloc);
double CalculateTSConstant2D(double aveVel, double aveAtt, int izsrcn, int ixsrcn, int zloc, int xloc);  // -wdd
double averageVelocityAroundSource3D(int bz, int bx, int by, int izsrcn, int ixsrcn, int iysrcn, int *dims, double *F);
double averageAttenuationAroundSource3D(int bz, int bx, int by, int izsrcn, int ixsrcn, int iysrcn, int *dims, double *Q); // -wdd
//double CalculateTConstant3D(double aveVel, int izsrcn, int ixsrcn, int iysrcn, int zloc, int xloc, int yloc);
double CalculateTConstant3D(double aveVel, int izsrcn, int ixsrcn, int iysrcn, double dzxy, int zloc, int xloc, int yloc);
//double CalculateTSConstant3D(double aveVel, double aveAtt, int izsrcn, int ixsrcn, int iysrcn, int zloc, int xloc, int yloc);
double CalculateTSConstant3D(double aveVel, double aveAtt, int izsrcn, int ixsrcn, int iysrcn, double dzxy, int zloc, int xloc, int yloc);

//struct CalculateDis2DResult {
//    double Tt;
//    double TSt;
//};

//void msfm2dCpp(double *T, double *F, int izsrcn, int ixsrcn, int *dims, bool usesecond, bool usecross) { // bool-This keyword is a built-in data type in C++ that can 
													 //       hold one of two values: true or false.  -wdd
//void msfm2dCpp(double *T, double *TS, double *F, double *Q, int izsrcn, int ixsrcn, int *dims, bool usesecond, bool usecross) { // bool-This keyword is a built-in data type in C++ that can
void msfm2dCpp(double *T, double *TS, double *F, double *Q, int izsrcn, int ixsrcn, double dzxy, int *dims, bool usesecond, bool usecross) {

	// box around the source - Dongzhuo
	int bz = 2;
	int bx = 2;

	/* Euclidian distance image */
	double *Y; // NOT USED CURRENTLY!!!

	/* Current distance values */
	//double Tt, Ty;
	double Tt, Ty, TSt, TSy;

	/* Matrix containing the Frozen Pixels" */
	bool *Frozen;

	/* Augmented Fast Marching (For skeletonize) */
	bool Ed = 0;
	//  in C++ is used to declare a boolean variable named Ed and initialize it to the value 0. -wdd
	//  This is a data type in C++ used to represent boolean values. It can hold two possible values: true (which is usually 
	//  represented by 1) and false (which is usually represented by 0).  -wdd

	// /* Size of input image */
	// const mwSize *dims_c;
	// mwSize dims[2];

	// /* Size of  SourcePoints array */
	// const mwSize *dims_sp_c;
	// mwSize dims_sp[2];

	/* Number of pixels in image */
	int npixels;  // declares an integer variable named npixels.  -wdd 
	int dims_sp[2] = {2, 1}; // declares and initializes an array of integers named dims_sp. -wdd
				 // dims_sp[0] is set to 2.; dims_sp[1] is set to 1. -wdd

	/* Neighbour list */
	int neg_free;
	int neg_pos;
	double *neg_listv; // declares a pointer to a double variable (or an array of double variables) named neg_listv. -wdd
	double *neg_listq; // -wdd for tstar-'t*' value. -wdd

	double *neg_listx;
	double *neg_listy;
	double *neg_listo;

	int *listproptt;
        int *listpropts;	
	double **listval; // declares a pointer to a pointer of type double, commonly referred to as a "double pointer." - wdd
	double **listtstar;

	/* Neighbours 4x2 */
	int ne[8] = {-1, 1, 0, 0, 0, 0, -1, 1};

	/* Loop variables  */
	int z, k, itt, q;

	/* Current location */
	int i, j, x, y;

	/* Index */
	int IJ_index, XY_index, index;

	npixels = dims[0] * dims[1];  // what is dims, and dims[0], dims[1] ??? -wdd
        // dims[0] = nz; dims[1] = nx;
	// compute average velocity in the small box around the soruce - Dongzhuo
	double aveVel = averageVelocityAroundSource2D(bz, bx, izsrcn, ixsrcn, dims, F);
	double aveAtt = averageAttenuationAroundSource2D(bz, bx, izsrcn, ixsrcn, dims, Q);
	// is a function call in C++ that calculates the average velocity around a specified source in a two-dimensional context. -wdd
	//tp cout << "aveVel = " << aveVel << endl;

	/* Pixels which are processed and have a final distance are frozen */
	Frozen = (bool *)malloc( npixels * sizeof(bool) );
	// sizeof(bool) - returns the size (in bytes) of a single bool variable. The size of a bool is typically 1 byte. -wdd
	// malloc(...) - (memory allocation) is a standard library function in C and C++ that allocates a specified number of bytes in 
	//               memory and returns a pointer to the beginning of the allocated memory block.                    -wdd
	// (bool *) - This is a type cast that converts the void pointer returned by malloc to a pointer of type bool*.  -wdd               
	for (q = 0; q < npixels; q++) {  // q++ - q = q + 1;
		Frozen[q] = 0; // The purpose of Frozen[q] = 0; is to initialize or reset the q-th element of the Frozen array to a defined value, which is 0 in this case.
		T[q] = -1;  // Same with 'Frozen[q] = 0'.
	        TS[q] = -1;
	}
	if (Ed) {
		for (q = 0; q < npixels; q++) {
			Y[q] = -1;
		}
	} //If Ed is true (or equivalently, 1), the code block inside the if statement will be executed. If Ed is false (or 0), the code block will be skipped.  -wdd
        // Although we do not know what are 'T' and 'Y', but we know that they are specified as '-1'.
	/*Free memory to store neighbours of the (segmented) region */
	neg_free = 100000;
	neg_pos = 0;
	neg_listx = (double *)malloc( neg_free * sizeof(double) );
	neg_listy = (double *)malloc( neg_free * sizeof(double) );
	if (Ed) {
		neg_listo = (double *)malloc( neg_free * sizeof(double) );
		for (q = 0; q < neg_free; q++) {
			neg_listo[q] = 0;
		}
	}
	// If Ed is false (or 0), the code block will be skipped.  -wdd

	/* List parameters array */
	//listprop = (int *)malloc(3 * sizeof(int)); // using malloc to allocate memory for an integer array of three elements. -wdd
	listproptt = (int *)malloc(3 * sizeof(int));
	listpropts = (int *)malloc(3 * sizeof(int));
	/* Make jagged list to store a maximum of 2^64 values */
	listval = (double **)malloc( 64 * sizeof(double *) );
	listtstar = (double **)malloc( 64 * sizeof(double *) ); // -wdd
	// malloc(64 * sizeof(double *)) allocates memory for 64 elements, each of which is a pointer to a double. -wdd
	// (double **) casts the returned void * pointer from malloc to double **, which is a pointer to a pointer of type double. -wdd
	/* Initialize parameter list */
	initialize_list(listval, listproptt); // initialize_list, to initialize the 2D array listval and the 1D array listprop. -wdd 
	initialize_list(listtstar, listpropts);  // -wdd

	neg_listv = listval[listproptt[1] - 1]; //assigns neg_listv the pointer to the row in listval specified by the index listprop[1] - 1. -wdd
	neg_listq = listtstar[listpropts[1] - 1];
        
	/*(There are 3 pixel classes:
	 *  - frozen (processed)
	 *  - narrow band (boundary) (in list to check for the next pixel with smallest distance)
	 *  - far (not yet used)
	 */

	/* set all starting points to distance zero and frozen  */
	/* and add all neighbours of the starting points to narrow list  */
	for (z = 0; z < dims_sp[1]; z++) {
		/*starting point  */
		// dims_sp[2] = {2, 1}, thus dims_sp[1] = 1;
		x = izsrcn;
		y = ixsrcn;
		XY_index = x + y * dims[0]; // dims[0] is 'nz'. Thus it denotes the total rows. -wdd
					    // Here we can guess that 'dims[0]' is the total number of gird in z direction. -wdd

		/*Set starting point to frozen and distance to zero  */
		Frozen[XY_index] = 1;
		T[XY_index] = 0;  // what is 'T' ?? -wdd So, 'T' is the traveltime ?
                TS[XY_index] = 0;  // -wdd
		if (Ed) {
			Y[XY_index] = 0; // what is 'Y' ?? -wdd So, 'Y' is the distance ?
		}
	}

	for (z = 0; z < dims_sp[1]; z++) { // dims_sp[2] = {2, 1}, thus dims_sp[1] = 1;
		/*starting point  */
		x = izsrcn;
		y = ixsrcn;
		XY_index = x + y * dims[0];

		/* Add neigbours of starting points  */
		for (k = 0; k < 4; k++) {
			// k = 0, 1, 2, 3. -wdd
			/*Location of neighbour  */
			i = x + ne[k]; j = y + ne[k + 4];
			// ne[8] = { -1, 1, 0, 0, 0, 0, -1, 1} -wdd
			//  thus, i = x + (-1, 1, 0, 0); j = y + (0, 0, -1, 1);
			//  Finally, (-1, 0), (1, 0), (0, -1) and (0, 1).
			//  It is up, down, left, right
			IJ_index = i + j * dims[0]; // thus dims[0] is 'nz'. -wdd
						    // thus j - column number; i - row number;

			/*Check if current neighbour is not yet frozen and inside the
			 *picture  */
			// isntfrozen2d - A function that checks if the element at (i, j) in the 2D logical array is not frozen. -wdd
			if (isntfrozen2d(i, j, dims, Frozen)) {
				Tt = (1 / (max(F[IJ_index], eps)));  // Let's guess what is 'F' ??? -wdd
								     // "F" is the velocity vector for whole region ! -wdd
								     // Thus, 'Tt' is the slowness. -wdd
				TSt = (1 / (max(F[IJ_index] * Q[IJ_index], eps))); // -wdd
				Ty = 1;
				TSy = 1;   // -wdd
				/*Update distance in neigbour list or add to neigbour list */
				if (T[IJ_index] > 0) {
					if (neg_listv[(int)T[IJ_index]] > Tt) {
						listupdate(listval, listproptt, (int)T[IJ_index], Tt);
						listupdate(listtstar, listpropts, (int)TS[IJ_index], TSt);  // -wdd
													  // tstar depends on time. -wdd
					}
					if (Ed) {
						neg_listo[(int)T[IJ_index]] = min(neg_listo[(int)T[IJ_index]], Ty);
					}
				} else {
					/*If running out of memory at a new block  */
					if (neg_pos >= neg_free) {
						// At begnning, neg_free = 100000; neg_pos = 0;  -wdd
						neg_free += 100000;
						neg_listx = (double *)realloc(neg_listx, neg_free * sizeof(double) );
						neg_listy = (double *)realloc(neg_listy, neg_free * sizeof(double) );
						if (Ed) {
							neg_listo = (double *)realloc(neg_listo, neg_free * sizeof(double) );
						}
					}
					list_add(listval, listproptt, Tt);
					list_add(listtstar, listpropts, TSt);  // -wdd
					// listprop - 3: A 1D array (likely int*) that might contain row indices or some values related to the elements of listval. -wdd
					// listval - 64: A 2D array (probably double**) where elements will be updated. -wdd
					// Tt: A double or int value to be added to elements in listval. -wdd
					// If 'Tt' is slowness or not ? -wdd
					neg_listv = listval[listproptt[1] - 1];
					neg_listq = listtstar[listpropts[1] - 1];  // -wdd
					// listprop[1] accesses the second element of the listprop array (since indexing starts at 0). -wdd
					// listprop[1] - 1 calculates the index of the row in listval that neg_listv will point to. -wdd
					// listval[listprop[1] - 1] gives the pointer to the beginning of that row, which is then assigned to neg_listv. -wdd
					neg_listx[neg_pos] = i; // At beginning, neg_pos = 0; -wdd
					neg_listy[neg_pos] = j;
					// 'neg_listx' and 'neg_listy' are used to save the number of row and column for narrow-band. -wdd 
					if (Ed) {
						neg_listo[neg_pos] = Ty;  // At beginning, for i, j, neg_listo = Ty = 1.  -wdd
					}
					T[IJ_index] = neg_pos;  // At beginning, T is '-1'. Then, convert it to 'neg_pos = 0, 1, 2, 3, 4,...' -wdd
					TS[IJ_index] = neg_pos; // -wdd
					neg_pos++;
				}
			}
		}
	}
	/*Loop through all pixels of the image  */
	for (itt = 0; itt < npixels; itt++) {
		/*Get the pixel from narrow list (boundary list) with smallest
		 *distance value and set it to current pixel location  */
		index = list_minimum(listval, listproptt);
		// listval: A 2D array (double** or similar) containing values. -wdd
		// listprop: A 1D array (likely int*) that may specify the rows or specific elements of listval to be searched for the minimum. -wdd
		// Now, we find the index of the minimun value in the narrow band. -wdd
		neg_listv = listval[listproptt[1] - 1];
		neg_listq = listtstar[listpropts[1] - 1]; // -wdd

		/* Stop if pixel distance is infinite (all pixels are processed)  */
		if (IsInf(neg_listv[index])) {
			break;
		}
		if (IsInf(neg_listq[index])) {
			break;
		}	// -wdd
		x = (int)neg_listx[index]; y = (int)neg_listy[index];
		// Based on the 'index' value, we then can find the row (x) and column (y) of this pixel. -wdd

		XY_index = x + y * dims[0];
		Frozen[XY_index] = 1;  // Then, frozen this minimum points.  -wdd
		T[XY_index] = neg_listv[index];  // what is 'neg_listv' ???
		TS[XY_index] = neg_listq[index];  // -wdd
		if (Ed) {
			Y[XY_index] = neg_listo[index];
		} // Here we should identify what are 'neg_listv' and 'neg_listo'. ??? -wdd

		/*Remove min value by replacing it with the last value in the array  */
		list_remove_replace(listval, listproptt, index) ;
		list_remove_replace(listtstar, listpropts, index) ;  // -wdd
		// listval: A 2D array (double**) where elements will be removed or replaced. -wdd
		// listprop: A 1D array (likely int*) containing row indices that we may need to modify. -wdd
		// index: An integer specifying the specific row or element in listval (via listprop) to remove or replace. -wdd
		neg_listv = listval[listproptt[1] - 1];
		neg_listq = listtstar[listpropts[1] - 1];   // -wdd
		if (index < (neg_pos - 1)) {
			neg_listx[index] = neg_listx[neg_pos - 1];
			neg_listy[index] = neg_listy[neg_pos - 1];
			if (Ed) {
				neg_listo[index] = neg_listo[neg_pos - 1];
			}
			T[(int)(neg_listx[index] + neg_listy[index]*dims[0])] = index;
			TS[(int)(neg_listx[index] + neg_listy[index]*dims[0])] = index;  // -wdd
		}
		neg_pos = neg_pos - 1;


		/*Loop through all 4 neighbours of current pixel  */
		for (k = 0; k < 4; k++) {
                        // 'k' - 0, 1, 2, 3. -wdd
			/*Location of neighbour  */
			i = x + ne[k]; j = y + ne[k + 4];
			IJ_index = i + j * dims[0];

			/*Check if current neighbour is not yet frozen and inside the  */
			/*picture  */
			if (isntfrozen2d(i, j, dims, Frozen)) {
				// analyticaly compute the time table near the source - Dongzhuo
				// originally
				// Tt = CalculateDistance2D(T, F[IJ_index], dims, i, j, usesecond, usecross, Frozen);
				if (abs(izsrcn - i) <= bz && abs(ixsrcn - j) <= bx) {
					Tt = CalculateTConstant2D(aveVel, izsrcn, ixsrcn, i, j);
					TSt = CalculateTSConstant2D(aveVel, aveAtt, izsrcn, ixsrcn, i, j); // -wdd
				} else {
					// Tt = CalculateDistance2D(T, F[IJ_index], dims, i, j, usesecond, usecross, Frozen);
					// should multiply the dimension of the grid - Dongzhuo
					//Tt = CalculateDistance2D(T, TS, F[IJ_index], Q[IJ_index], dims, i, j, usesecond, usecross, Frozen);
					CalculateDis2DResult result = CalculateDistance2D(T, TS, F[IJ_index], Q[IJ_index], dims, dzxy, i, j, usesecond, usecross, Frozen);
                                        Tt = result.Tt; TSt = result.TSt;
					//CalculateDis2DResult CalculateDistance2D
				}

				if (Ed) {
					// should multiply the dimension of the grid - Dongzhuo
					//Ty = CalculateDistance2D(Y, 1, dims, i, j, usesecond, usecross, Frozen);
					//Ty = CalculateDistance2D(Y, Y, 1, 1, dims, i, j, usesecond, usecross, Frozen);
					CalculateDis2DResult result = CalculateDistance2D(Y, Y, 1, 1, dims, dzxy, i, j, usesecond, usecross, Frozen);
					Tt = result.Tt; TSt = result.TSt;
				}

				/*Update distance in neigbour list or add to neigbour list */
				IJ_index = i + j * dims[0];
				if ((T[IJ_index] > -1) && T[IJ_index] <= listproptt[0]) {
					if (neg_listv[(int)T[IJ_index]] > Tt) {
						listupdate(listval, listproptt,    (int)T[IJ_index], Tt);
						listupdate(listtstar,  listpropts,    (int)TS[IJ_index], TSt);  // -wdd
					}
					if (Ed) {
						neg_listo[neg_pos] = min(neg_listo[neg_pos], Ty);
					}
				} else {
					/*If running out of memory at a new block */
					if (neg_pos >= neg_free) {
						neg_free += 100000;
						neg_listx = (double *)realloc(neg_listx, neg_free * sizeof(double) );
						neg_listy = (double *)realloc(neg_listy, neg_free * sizeof(double) );
						if (Ed) {
							neg_listo = (double *)realloc(neg_listo, neg_free * sizeof(double) );
						}
					}
					list_add(listval, listproptt, Tt);
					list_add(listtstar, listpropts, TSt);  // -wdd

					neg_listv = listval[listproptt[1] - 1];
					neg_listq = listtstar[listpropts[1] - 1];  // -wdd

					neg_listx[neg_pos] = i; neg_listy[neg_pos] = j;
					if (Ed) {
						neg_listo[neg_pos] = Ty;
					}
					T[IJ_index] = neg_pos;
					TS[IJ_index] = neg_pos; // -wdd
					neg_pos++;
				}
			}
		}

	}
	/* Free memory */
	/* Destroy parameter list */
	destroy_list(listval, listproptt);
	destroy_list(listtstar, listpropts);

	free(neg_listx);
	free(neg_listy);
	if (Ed) {
		free(neg_listo);
	}
	free(Frozen);
}

//void msfm3dCpp(double *T, double *TS, double *F, double *Q, int izsrcn, int ixsrcn, int iysrcn, int *dims, bool usesecond, bool usecross) {
void msfm3dCpp(double *T, double *TS, double *F, double *Q, int izsrcn, int ixsrcn, int iysrcn, double dzxy, int *dims, bool usesecond, bool usecross) {
	/* The input variables */

	// box around the source - Dongzhuo
	int bz = 2;
	int bx = 2;
	int by = 2;

	/* Euclidian distance image */
	double *Y; // NOT USED CURRENTLY!!!

	/* Current distance values */
	double Tt, Ty, TSt, TSy;

	/* Matrix containing the Frozen Pixels" */
	bool *Frozen;

	/* Augmented Fast Marching (For skeletonize) */
	bool Ed = 0;

	/* Size of input image */
	// const mwSize *dims_c;
	// mwSize dims[3];

	/* Size of  SourcePoints array */
	// const mwSize *dims_sp_c;
	// mwSize dims_sp[3];

	/* Number of pixels in image */
	int npixels;
	int dims_sp[2] = {3, 1}; // dims_sp[0] = 3; dims_sp[1] = 1; -wdd

	/* Neighbour list */
	int neg_free;
	int neg_pos;
	double *neg_listv;
	double *neg_listq; // -wdd for tstar-'t*' value. -wdd

	double *neg_listx;
	double *neg_listy;
	double *neg_listz;
	double *neg_listo;

	int *listproptt;
	int *listpropts;
	double **listval;
        double **listtstar;

	/* Neighbours 6x3 */
	int ne[18] = {-1, 0, 0, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, -1, 0, 0, 1};

	/* Loop variables */
	int s, w, itt, q;

	/* Current location */
	int i, j, k, x, y, z;

	/* Index */
	int IJK_index, XYZ_index, index;

	npixels = dims[0] * dims[1] * dims[2];
        // dims[0] -- nz; dims[1] -- nx; dims[2] -- ny;  // -wdd
	//cout << "dims[0] = " << dims[0] << endl;
	//cout << "dims[1] = " << dims[1] << endl;
	//cout << "dims[2] = " << dims[2] << endl;       // -wdd

	// compute average velocity in the small box around the soruce - Dongzhuo
	double aveVel = averageVelocityAroundSource3D(bz, bx, by, izsrcn, ixsrcn, iysrcn, dims, F);
	double aveAtt = averageAttenuationAroundSource3D(bz, bx, by, izsrcn, ixsrcn, iysrcn, dims, Q);
	//tp cout << "aveVel = " << aveVel << endl;
	// cout << "aveVel = " << aveVel << endl;  //-wdd

	/* Pixels which are processed and have a final distance are frozen */
	Frozen = (bool *)malloc( npixels * sizeof(bool));
	for (q = 0; q < npixels; q++) {
		Frozen[q] = 0; // All nodes are not Frozen - '0' !!!
		T[q] = -1;     // 'T' is finally used to save the traveltime values in each nodes. //-wdd
	        TS[q] = -1;
	}
	if (Ed) {
		for (q = 0; q < npixels; q++) {
			Y[q] = -1;
	        //cout << "Y[0] = " << Y[0] << endl;	// -wdd
		}
	}

	/*Free memory to store neighbours of the (segmented) region */
	//neg_free = 100000;
	neg_free = 1000000;
	neg_pos = 0;

	neg_listx = (double *)malloc( neg_free * sizeof(double) );
	neg_listy = (double *)malloc( neg_free * sizeof(double) );
	neg_listz = (double *)malloc( neg_free * sizeof(double) );

	if (Ed) {
		neg_listo = (double *)malloc( neg_free * sizeof(double) );
		for (q = 0; q < neg_free; q++) {
			neg_listo[q] = 0;
		}
	}

	/* List parameters array */
	listproptt = (int *)malloc(3 * sizeof(int));
	listpropts = (int *)malloc(3 * sizeof(int));
	// Allocates a contiguous block of memory large enough to store 3 integers. -wdd
	// cout << "listprop[3] = " << listprop[3] << endl;  -wdd
	// listprop[0] = 0; listprop[1] = 1; listprop[2] = 2; //-wdd
        // std::cout << sizeof(listprop) / sizeof(int) << std::endl;
	
	/* Make jagged list to store a maximum of 2^64 values */
	listval = (double **)malloc( 64 * sizeof(double *) );
	listtstar = (double **)malloc( 64 * sizeof(double *) );  // -wdd
	// listval[i] - i is between 0, 1, 2.  ???? -wdd
	// listval[i][j] is a 2D array. i from 0 to 1 and 2. ??? j from what to what ? -wdd

	// cout << "listval[0] = " << listval[0] << endl;

	/* Initialize parameter list */
	initialize_list(listval, listproptt);
	initialize_list(listtstar, listpropts);
	// I guess that here 'listprop' is used to set 'listval' to 3 rows. -wdd
	// Based on the follows, here listval[i][j]; i is 0, 1, 2. - wdd
	neg_listv = listval[listproptt[1] - 1];
	neg_listq = listtstar[listpropts[1] - 1];
	// listprop[1] = 1, thus listprop[1] - 1 = 0;
	// neg_listv = listval[0];  -wdd
	//cout << "listprop[1] - 1 = " << listprop[1] - 1 << endl;

	// cout << "neg_listv = " << neg_listv << endl;   // -wdd

	/*(There are 3 pixel classes: */
	/*  - frozen (processed) */
	/*  - narrow band (boundary) (in list to check for the next pixel with smallest distance) */
	/*  - far (not yet used) */
	/* set all starting points to distance zero and frozen */
	/* and add all neighbours of the starting points to narrow list */

	// initial the source with Frozen[] = 1, indicating the source is frozened. -wdd
	for (s = 0; s < dims_sp[1]; s++) {
		// dims_sp[0] = 3; dims_sp[1] = 1; Thus, s = 0;
		/*starting point */
		x = izsrcn;  // z - 125 [151]
		y = ixsrcn;  // x - 75  [151]
		z = iysrcn;  // y - 1   [3]

		XYZ_index = mindex3(x, y, z, dims[0], dims[1]); // dims[0] - nz; dims[1] - nx; -wdd
		// int mindex3(int x, int y, int z, int dim_x, int dim_y) {
		// return z * (dim_x * dim_y) + y * dim_x + x;}   // -wdd
		
		// The order is for 3D: 1st, from up to down, 2nd, from left to right, 3rd, from forth to back; -wdd

		// cout << "XYZ_index = " << XYZ_index << endl;  // 34251.
		// z * (dims[0] * dims[1]) + y * dims[0] + x
		// 1 * 151 * 151 + 75 *151 + 125 = 34251

		Frozen[XYZ_index] = 1; // Frozen the source. -wdd
		T[XYZ_index] = 0;      // At the source location, T = 0. -wdd
	        TS[XYZ_index] = 0;
		if (Ed) {
			Y[XYZ_index] = 0;
		}
	}

	for (s = 0; s < dims_sp[1]; s++) { 
		// dims_sp[1] = 1; thus, s = 0; Just loop one node. -wdd
		/*starting point */
		x = izsrcn;
		y = ixsrcn;
		z = iysrcn;

		XYZ_index = mindex3(x, y, z, dims[0], dims[1]);

		for (w = 0; w < 6; w++) {
			// 0, 1, 2, 3, 4, 5
			/*Location of neighbour */
			// ne[18] = { -1, 0, 0, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, -1, 0, 0, 1};  // -wdd
			
			i = x + ne[w];     // { -1, 0, 0, 1, 0, 0}  // '-1' - up, '1' -down;    -wdd
			j = y + ne[w + 6]; // { 0, -1, 0, 0, 1, 0}  // '-1' - left, '1' -right; -wdd
			k = z + ne[w + 12];// { 0, 0, -1, 0, 0, 1}  // '-1' - forth, '1' -back; -wdd

			IJK_index = mindex3(i, j, k, dims[0], dims[1]);
			// k * (dims[0] * dims[1]) + j * dims[0] + i  // -wdd

			/*Check if current neighbour is not yet frozen and inside the */

			/*picture */
			if (isntfrozen3d(i, j, k, dims, Frozen)) {
				// Here 'isntfrozen3d'- should be a subroutine to determine if this node is frozen. -wdd
				//Tt = (1 / (max(F[IJK_index], eps))); // 'Tt' is a temporary time (T) or not. -wdd
				Tt = (dzxy / (max(F[IJK_index], eps))); // 'Tt' is a temporary time (T) or not. -wdd
				//TSt = (1 / (max(F[IJK_index] * Q[IJK_index], eps)));
				TSt = (dzxy / (max(F[IJK_index] * Q[IJK_index], eps)));
				/*Update distance in neigbour list or add to neigbour list */
				if (T[IJK_index] > 0) {  // At beginning, T is '-1'. -wdd
							 // This section cannot be executed, because it located at the source location. -wdd
							 // Each node can only be calculated once. IJK_index can only once. -wdd
					if (neg_listv[(int)T[IJK_index]] > Tt) {
						// T[IJK_index]=neg_pos, thus here 'T' is integer. -wdd
						// T - the number(index) of 'IJK_index' in the select matrix. -wdd
						// neg_listv[] invoke the previous Tt value. -wdd
						listupdate(listval, listproptt, (int)T[IJK_index], Tt);
						listupdate(listtstar, listpropts, (int)TS[IJK_index], TSt);
						// cout << "IJK_index = " << IJK_index << endl; // -wdd
					}
				} else {
					/*If running out of memory at a new block */
					if (neg_pos >= neg_free) {  // indicates that the total number of node within the vector, which then can be used to select min values. -wdd 
						neg_free += 1000000; // At beginning, neg_free = 100000; // -wdd
								    // If it exceeds 100000, then again allocate one more '100000'. -wdd
						neg_listx = (double *)realloc(neg_listx, neg_free * sizeof(double) );
						neg_listy = (double *)realloc(neg_listy, neg_free * sizeof(double) );
						neg_listz = (double *)realloc(neg_listz, neg_free * sizeof(double) );
						// Why the index in 3 directions (dimensions) are specified as double 'neg_listx', 'neg_listy' and 'neg_listz'. -wdd
						if (Ed) {
							neg_listo = (double *)realloc(neg_listo, neg_free * sizeof(double) );
						}
					}
					//cout << "Tt = " << Tt << endl;
					//cout << "listprop[2] = " << listprop[2] << endl;
					// listprop[0] = 0, 1, 2, 3, 4, 5. -wdd
					// listprop[1] = 1, 1, 1, 2, 2, 3. -wdd
					// listprop[2] = 2, 2, 2, 4, 4, 8. -wdd  Just see the follows. This section is no-use !!! -wdd
					list_add(listval, listproptt, Tt);
					list_add(listtstar, listpropts, TSt);  // -wdd
				        //cout << "listval[1][1] = " << listval[1][1] << endl;
					//cout << "listprop[2] = " << listprop[2] << endl;
                                        // listprop[0] = 1, 2, 3, 4, 5, 6. -wdd
					// listprop[1] = 1, 1, 2, 2, 3, 3. -wdd
					// listprop[2] = 2, 2, 4, 4, 8, 8. -wdd
					// Thus, we can understand that 'listprop[0]' is the total number of elements in the vector 'listval'. -wdd
					// Because 6 neighboring nodes around the sources, thus finally listprop[0] = 0. -wdd
					//Here 'list_add'- should be a subroutine to add the 'Tt' for each node to this matrix. -wdd
					neg_listv = listval[listproptt[1] - 1];
					neg_listq = listtstar[listpropts[1] - 1];  // -wdd
					// listprop[1] - 1 = 0, 0, 1, 1, 2, 2. -wdd
					// Pick the first row[0] of listval. -wdd
                                        //cout << "listprop[1] = " << listprop[1] << endl;  //  -wdd
					//cout << "neg_listv[2] = " << neg_listv[2] << endl;
					neg_listx[neg_pos] = i;
					neg_listy[neg_pos] = j;
					neg_listz[neg_pos] = k;  // i, j, and k are index in 3 directions. -wdd
					T[IJK_index] = neg_pos;  // At beginning, 'neg_pos = 0'.  -wdd
					TS[IJK_index] = neg_pos;
					neg_pos++;
				}
			}
		}
	}
	// stack - neg_pos - number of all neigbouring nodes. -wdd
	// The above 'subroutine' is used to calculate the Tt in the neighbouring 6 nodes, and add them into the matrix. -wdd
	// But did not do any other further processing.  -wdd

	//int n = 6; // Assume the size of the array is known
	//for (int i = 0; i < n; i++) {
    	//cout << "listval[2][" << i << "] = " << listval[2][i] << endl;
	//} // -wdd
	// listval[2][:] = 0.166667.
	// we only have: listval[0][:],listval[1][:],listval[2][:]     -wdd
        // We guess that the real 'Tt' is saved ib the 'listval[2][:]' - wdd

	//cout << "listval[0][:] = " << listval[0][:] << endl;

	// cout << "listval = " << listval << endl;
        // cout << "listprop[1] - 1 = " << listprop[1] - 1 << endl;
	// listprop[1] - 1 = 2. -wdd
	/*Loop through all pixels of the image */
	//npixels = 40;
	for (itt = 0; itt < (npixels); itt++) { /* */
		/*Get the pixel from narrow list (boundary list) with smallest */
		/*distance value and set it to current pixel location */
		index = list_minimum(listval, listproptt);
		// 'index' indicates which column in the 2D matrix 'listval'. with the dimension: '3*N'.
		//  'N' indicates these nodes added. -wdd
		//cout << "index = " << index << endl;
		//cout << "listprop[1] - 1 = " << listprop[1] - 1 << endl;
		neg_listv = listval[listproptt[1] - 1];
                neg_listq = listtstar[listpropts[1] - 1];
		// 'neg_listv' is used to extract the neighboring Tt values. -wdd
		// neg_listv - all the neighboring Tt values calculated in previous round !!! -wdd
		/* Stop if pixel distance is infinite (all pixels are processed) */
		if (IsInf(neg_listv[index])) {
			break;
		}
		//if (IsInf(neg_listq[index])) {
                //        break;
                //}       // -wdd
		/*index=minarray(neg_listv, neg_pos); */
		x = (int)neg_listx[index]; y = (int)neg_listy[index]; z = (int)neg_listz[index];

		XYZ_index = mindex3(x, y, z, dims[0], dims[1]);
		Frozen[XYZ_index] = 1;   // Frozen this node.
		T[XYZ_index] = neg_listv[index];  // neg_listv[index] is the Tt value of the Frozen node. -wdd
						  // Thus, T[XYZ_index] is the real Tt of the Frozen node 'XYZ_index'. -wdd
		TS[XYZ_index] = neg_listq[index];
		//cout << "neg_listv[index] = " << neg_listv[index] << endl; // -wdd
		//cout << "neg_listv[index] = " << neg_listv[index] << endl;
		//cout << "neg_listq[index] = " << neg_listq[index] << endl;  // -wdd
		//cout << "neg_listv[index]/neg_listq[index] = " << neg_listv[index]/neg_listq[index] << endl;  // -wdd
		//neg_listv[index] = 0.166667  - Tt. -- T.  -wdd
		// listval[0/1/2][index]  -wdd
		//
		if (Ed) {
			Y[XYZ_index] = neg_listo[index];
		}

		/*Remove min value by replacing it with the last value in the array */
		list_remove_replace(listval, listproptt, index) ;
		list_remove_replace(listtstar, listpropts, index) ;  // -wdd
		// remove the specific comuln indicated by 'index'. -wdd
		// remove the index column, and move the last column to the index column. -wdd
		// Thus, we need to replace their corresponding location: neg_listx,neg_listy,neg_listz. -wdd
		// After executing this 'list_remove_replace', then both the structure 'listval' and 'listprop' will be changed. -wdd
		neg_listv = listval[listproptt[1] - 1];
		neg_listq = listtstar[listpropts[1] - 1];   // -wdd
		// After removing the index column, thus we need to generate a new neg_listv including all other Tt values. -wdd
	        //cout << "index = " << index << ", neg_pos - 1 = " << neg_pos - 1 << endl;
		if (index < (neg_pos - 1)) {  // The min Tt is not the last column. -wdd
                        //cout << "neg_listx[neg_pos - 1] = " << neg_listx[neg_pos - 1] << endl;
                        //cout << "neg_listx[index] = " << neg_listx[index] << ", neg_listx[neg_pos - 1] = " << neg_listx[neg_pos - 1] << endl;
	 		neg_listx[index] = neg_listx[neg_pos - 1];
			neg_listy[index] = neg_listy[neg_pos - 1];
			neg_listz[index] = neg_listz[neg_pos - 1];
			// Replace the location, because it has been replaced with the last column. -wdd
			//cout << "neg_listx[index] = " << neg_listx[index] << ", neg_listx[neg_pos - 1] = " << neg_listx[neg_pos - 1] << endl;
                        // It is replaced by the last value in the array. -wdd
			if (Ed) {
				neg_listo[index] = neg_listo[neg_pos - 1];
			}
			T[(int)mindex3((int)neg_listx[index], (int)neg_listy[index], (int)neg_listz[index], dims[0], dims[1])] = index;
			TS[(int)mindex3((int)neg_listx[index], (int)neg_listy[index], (int)neg_listz[index], dims[0], dims[1])] = index;
			// Because it has been replaced. Thus, it is necessary to also change its (original in the last column) corresponding 'index'.
		}

		// Here we make a brief summary:
		// T - can represent two parameters:
		// If it is not forzen, T represents the location or index with a vector, where then find the min value. -wdd
		// If it is Frozen, T represents the real Tt in this node, and keep unchanged in the follows. -wdd 
		neg_pos = neg_pos - 1;
		// Because 'index' column has been removed and Frozen, thus we need to minus 1 from original vector. -wdd

		/*Loop through all 6 neighbours of current pixel */
		for (w = 0; w < 6; w++) { // 0, 1, 2, 3, 4, 5. -wdd
			/*Location of neighbour */
			i = x + ne[w]; j = y + ne[w + 6]; k = z + ne[w + 12]; // Here we know that 'x', 'y' and 'z' is the index in 3D for the frozen node. -wdd
			IJK_index = mindex3(i, j, k, dims[0], dims[1]);       // 'i', 'j' and 'k' are the index in the 3D for its 6 neighboring nodes. -wdd

			/*Check if current neighbour is not yet frozen and inside the */
			/*picture */
			if (isntfrozen3d(i, j, k, dims, Frozen)) {
				// analyticaly compute the time table near the source - Dongzhuo
				// originally
				// Tt = CalculateDistance3D(T, F[IJK_index], dims, i, j, k, usesecond, usecross, Frozen);
				if (abs(izsrcn - i) <= bz && abs(ixsrcn - j) <= bx && abs(iysrcn - i) <= by) {
					// indicates that the node is located in the velocity average regions. -wdd
					//Tt = CalculateTConstant3D(aveVel, izsrcn, ixsrcn, iysrcn, i, j, k);
					Tt = CalculateTConstant3D(aveVel, izsrcn, ixsrcn, iysrcn, dzxy, i, j, k);
					//TSt = CalculateTSConstant3D(aveVel, aveAtt, izsrcn, ixsrcn, iysrcn, i, j, k);
					TSt = CalculateTSConstant3D(aveVel, aveAtt, izsrcn, ixsrcn, iysrcn, dzxy, i, j, k);
					//TSt = CalculateTSConstant3D(aveVel, aveAtt, izsrcn, ixsrcn, iysrcn, i, j, k);  // -wdd
				} else {
					// Tt = CalculateDistance3D(T, F[IJK_index], dims, i, j, k, usesecond, usecross, Frozen);
					// should multiply the dimension of the grid - DongzhuoCalculateTConstant3D
					//Tt = CalculateDistance3D(T, F[IJK_index], dims, i, j, k, usesecond, usecross, Frozen);
					CalculateDis3DResult result = CalculateDistance3D(T, TS, F[IJK_index], Q[IJK_index], dims, dzxy, i, j, k, usesecond, usecross, Frozen);
					// If it is not located in the source regions, then use another method to calculate the 'Tt'. -wdd
					Tt = result.Tt; TSt = result.TSt;
					//cout << "Tt/TSt = " << Tt/TSt << endl;
				}

				if (Ed) {
					// should multiply the dimension of the grid - Dongzhuo
					//Ty = CalculateDistance3D(Y, 1, dims, i, j, k, usesecond, usecross, Frozen);
					CalculateDis3DResult result = CalculateDistance3D(Y, Y, 1, 1, dims, dzxy, i, j, k, usesecond, usecross, Frozen);
					Tt = result.Tt; TSt = result.TSt;
					//cout << "Tt/TSt = " << Tt/TSt << endl;
				}

				/*Update distance in neigbour list or add to neigbour list */
				IJK_index = mindex3(i, j, k, dims[0], dims[1]);

				if ((T[IJK_index] > -1) && T[IJK_index] <= listproptt[0]) {  // T[IJK_index] > -1 - indicates that this node has an existing Tt (T is not -1).
											   // T[IJK_index] <= listprop[0] ensure that no exceed the max nodes number. -wdd
                                        //cout << "listprop[0] = " << listprop[0] << endl;
					if (neg_listv[(int)T[IJK_index]] > Tt) {
						//cout << "Tt = " << Tt << endl;
						// T[IJK_index] is used to output the index of this node, which has a known or calculated Tt. -wdd
                                                //cout << "neg_listv[(int)T[IJK_index]] = " << neg_listv[(int)T[IJK_index]] << ", Tt = " << Tt << endl;
						listupdate(listval, listproptt, (int)T[IJK_index], Tt);
						listupdate(listtstar,  listpropts,    (int)TS[IJK_index], TSt);  // -wdd
					}
				} else {
					/*If running out of memory at a new block */
					if (neg_pos >= neg_free) {
						neg_free += 1000000; // If it exceeds 100000, then we allocate another 100000. -wdd
						neg_listx = (double *)realloc(neg_listx, neg_free * sizeof(double) );
						neg_listy = (double *)realloc(neg_listy, neg_free * sizeof(double) );
						neg_listz = (double *)realloc(neg_listz, neg_free * sizeof(double) );
						if (Ed) {
							neg_listo = (double *)realloc(neg_listo, neg_free * sizeof(double) );
						}
					}
					list_add(listval, listproptt, Tt);
					list_add(listtstar, listpropts, TSt);  // -wdd

					neg_listv = listval[listproptt[1] - 1];
					neg_listq = listtstar[listpropts[1] - 1];  // -wdd

					// neg_listv - to select the previous calculated 'Tt' values and save them into a vector 'neg_listv'. -wdd
					// Then, we select the minimum value from the vector 'neg_listv'. -wdd
					neg_listx[neg_pos] = i; neg_listy[neg_pos] = j; neg_listz[neg_pos] = k;
					// Save the index in 3D (three dimensions) to 'neg_listx', 'neg_listy', and 'neg_listz'. -wdd
					if (Ed) {
						neg_listo[neg_pos] = Ty;
					}

					T[IJK_index] = neg_pos; // Also allocate a new number of this node into the vector. -wdd
					TS[IJK_index] = neg_pos; // -wdd
					// While add a new neighbouring nodes into the vector, then the index should also add '1'. -wdd
					neg_pos++;
				}
			}
		}

	}
	/* Free memory */
	/* Destroy parameter list */
	destroy_list(listval, listproptt);
	destroy_list(listtstar, listpropts);

	free(neg_listx);
	free(neg_listy);
	free(neg_listz);
	free(Frozen);
}

//struct TConstant3DResult {
//struct CalculateDis2DResult {
//    double Tt;
//    double Tst;
//};

//double CalculateDistance2D(double *T, double *TS, double Fij, double Qij, int *dims, int i, int j, bool usesecond, bool usecross, bool *Frozen) {
//CalculateDis2DResult CalculateDistance2D(double *T, double *TS, double Fij, double Qij, int *dims, int i, int j, bool usesecond, bool usecross, bool *Frozen) {
CalculateDis2DResult CalculateDistance2D(double *T, double *TS, double Fij, double Qij, int *dims, double dzxy, int i, int j, bool usesecond, bool usecross, bool *Frozen) {
	/* Derivatives */
	double Tm[4] = {0, 0, 0, 0};  // initializes an array Tm of type double with four elements, each set to 0. -wdd
	double TSm[4] = {0, 0, 0, 0}; // -wdd
	double Tm2[4] = {0, 0, 0, 0}; // initializes an array Tm2 of type double with four elements, each set to 0. -wdd
	double TSm2[4] = {0, 0, 0, 0};// -wdd
	double Coeff[3]; // declares an array named Coeff with three elements of type double. -wdd
        double Coeff_ts[2];

	/* local derivatives in distance image */
	double Tpatch_2_3, Txm2, Tpatch_4_3, Txp2;
	double Tpatch_3_2, Tym2, Tpatch_3_4, Typ2;
	double Tpatch_2_2, Tr1m2, Tpatch_4_4, Tr1p2;
	double Tpatch_2_4, Tr2m2, Tpatch_4_2, Tr2p2;

	double TSpatch_2_3, TSxm2, TSpatch_4_3, TSxp2;
	double TSpatch_3_2, TSym2, TSpatch_3_4, TSyp2;
	double TSpatch_2_2, TSr1m2, TSpatch_4_4, TSr1p2;
	double TSpatch_2_4, TSr2m2, TSpatch_4_2, TSr2p2;  // -wdd

	/* Return values root of polynomial */
	double ansroot[2] = {0, 0};

	/* Loop variables  */
	int q, t;

	/* Derivative checks */
	bool ch1, ch2;

	/* Order derivatives in a certain direction */
	int Order[4] = {0, 0, 0, 0};

	/* Current location */
	int in, jn;

	/* Constant cross term */
	const double c1 = 0.5;

	double Tt, Tt2;
        double TSt, TSt2;
	/*Get First order derivatives (only use frozen pixel)  */
	in = i - 1; jn = j + 0; if (isfrozen2d(in, jn, dims, Frozen)) {
		// Here 'isfrozen2d' is different with 'isntfrozen2d' used previously. -wdd
		Tpatch_2_3 = T[mindex2(in, jn, dims[0])]; TSpatch_2_3 = TS[mindex2(in, jn, dims[0])];
        // Tpatch_2_3 - left -wdd		
	// suggests you’re accessing an element in a 2D array (flattened into a 1D array) called T using a 
	// function mindex2 to calculate the index.             -wdd
	// T is a 1D array that represents a 2D grid of values. -wdd
	// mindex2 is a helper function that computes the index in this 1D array T based on two coordinates (in, jn) 
	// and the row size dims[0] (often the number of columns in the 2D representation).  -wdd
	} else {
		Tpatch_2_3 = INF; TSpatch_2_3 = INF;	
	}
	in = i + 0; jn = j - 1; if (isfrozen2d(in, jn, dims, Frozen)) {
		Tpatch_3_2 = T[mindex2(in, jn, dims[0])];  TSpatch_3_2 = TS[mindex2(in, jn, dims[0])];
        // Tpatch_2_3 - up -wdd
	} else {
		Tpatch_3_2 = INF; TSpatch_3_2 = INF;
	}
	in = i + 0; jn = j + 1; if (isfrozen2d(in, jn, dims, Frozen)) {
		Tpatch_3_4 = T[mindex2(in, jn, dims[0])];  TSpatch_3_4 = TS[mindex2(in, jn, dims[0])];
        // Tpatch_2_3 - down -wdd
	} else {
		Tpatch_3_4 = INF; TSpatch_3_4 = INF;
	}
	in = i + 1; jn = j + 0; if (isfrozen2d(in, jn, dims, Frozen)) {
		Tpatch_4_3 = T[mindex2(in, jn, dims[0])];  TSpatch_4_3 = TS[mindex2(in, jn, dims[0])];
        // Tpatch_4_3 - right -wdd
	} else {
		Tpatch_4_3 = INF; TSpatch_4_3 = INF;
	}

	if (usecross) {
		in = i - 1; jn = j - 1; if (isfrozen2d(in, jn, dims, Frozen)) {
			Tpatch_2_2 = T[mindex2(in, jn, dims[0])];  TSpatch_2_2 = TS[mindex2(in, jn, dims[0])];
	        // usecoss - means that the angle of rotation is 45 and 135 to use the cross points. -wdd 		
		} else {
			Tpatch_2_2 = INF;  TSpatch_2_2 = INF;
		}
		in = i - 1; jn = j + 1; if (isfrozen2d(in, jn, dims, Frozen)) {
			Tpatch_2_4 = T[mindex2(in, jn, dims[0])];  TSpatch_2_4 = TS[mindex2(in, jn, dims[0])];
		} else {
			Tpatch_2_4 = INF;  TSpatch_2_4 = INF;
		}
		in = i + 1; jn = j - 1; if (isfrozen2d(in, jn, dims, Frozen)) {
			Tpatch_4_2 = T[mindex2(in, jn, dims[0])];  TSpatch_4_2 = TS[mindex2(in, jn, dims[0])];
		} else {
			Tpatch_4_2 = INF;  TSpatch_4_2 = INF;
		}
		in = i + 1; jn = j + 1; if (isfrozen2d(in, jn, dims, Frozen)) {
			Tpatch_4_4 = T[mindex2(in, jn, dims[0])];  TSpatch_4_4 = TS[mindex2(in, jn, dims[0])];
		} else {
			Tpatch_4_4 = INF;  TSpatch_4_4 = INF;
		// Tpatch_2_2, Tpatch_2_4, Tpatch_4_2 and Tpatch_4_4 are the times used in these points. -wdd
		// They are similar with that adopted in the rotated finite difference methods.          -wdd
		}
	}

	// to be more visual, these order should be like this:
	//             Tpatch_2_2    Tpatch_3_2     Tpatch_4_2
	//             Tpatch_2_3   (Tpatch_3_3)    Tpatch_4_3
	//             Tpatch_2_4    Tpatch_3_4     Tpatch_4_4
             
	//                  Tpatch_3_3 --- Tpatch_i_j

	/*The values in order is 0 if no neighbours in that direction  */
	/*1 if 1e order derivatives is used and 2 if second order  */
	/*derivatives are used  */
	Order[0] = 0; Order[1] = 0; Order[2] = 0; Order[3] = 0;
	// Order[0] - 0 - indicates the horizontal direction.
	// Order[1] - 1 - indicates the vertical direction.
	// Order[2] - 2 - indicates the inclined direction with angle 45.
	// Order[3] - 3 - indicates the inclined direction with angle 135.
	// The values of 'Order[0], Order[1], Order[2] and Order[3]' indicates that if and how this terms show in the calculations. -wdd
	/*Make 1e order derivatives in x and y direction  */
	Tm[0] = min( Tpatch_2_3 , Tpatch_4_3);
        if (Tm[0] == Tpatch_2_3){
		TSm[0] = TSpatch_2_3;
	} else {
		TSm[0] = TSpatch_4_3;
	}
	if (IsFinite(Tm[0])) {
		Order[0] = 1;   // in the horizontal direction. -wdd
	} // If both 'Tpatch_2_3' and 'Tpatch_4_3' are INF, then Order[]=0, indicating that this them is 0. -wdd
	Tm[1] = min( Tpatch_3_2 , Tpatch_3_4);
        if (Tm[1] == Tpatch_3_2){
		TSm[1] = TSpatch_3_2;
	} else {
	        TSm[1] = TSpatch_3_4;
	}
	if (IsFinite(Tm[1])) {
		Order[1] = 1;   // in the vertical direction. -wdd
	} // If both 'Tpatch_3_2' and 'Tpatch_3_4' are INF, then Order[]=0, indicating that this them is 0. -wdd
	/*Make 1e order derivatives in cross directions  */
	if (usecross) {
		Tm[2] = min( Tpatch_2_2 , Tpatch_4_4);
		if (Tm[2] == Tpatch_2_2){
			TSm[2] = TSpatch_2_2;
		} else {
			TSm[2] = TSpatch_4_4;
		}
	       	if (IsFinite(Tm[2])) {
			Order[2] = 1; // in the -45' direction
		} // If both 'Tpatch_2_2' and 'Tpatch_4_4' are INF, then Order[]=0, indicating that this them is 0. -wdd
		Tm[3] = min( Tpatch_2_4 , Tpatch_4_2); 
		if (Tm[3] == Tpatch_2_4){
			TSm[3] = TSpatch_2_4;
		} else {
			TSm[3] = TSpatch_4_2;
		}  // -wdd
		if (IsFinite(Tm[3])) {
			Order[3] = 1; // in the +45' direction
		} // If both 'Tpatch_2_4' and 'Tpatch_4_2' are INF, then Order[]=0, indicating that this them is 0. -wdd
	}
	// Finally, we can obtain Tm[0], Tm[1], Tm[2] and Tm[3].  -wdd
        // Order[0], Order[1], Order[2] and Order[3] equals to 1 indicates only the one-order difference is used. -wdd 
	 
        // The above is used to select minimum T value in horizontal, vertical, and 2 inclined directions. -wdd
	// In each direction, only one T value in each direction are selected. -wdd
        // It is define as Tm[0], Tm=[1], Tm[2] and Tm[3], respectively in horizontal, vertical, and two inclined directions. -wdd
	
	/*Make 2e order derivatives  */
	if (usesecond) {
		/*Get Second order derivatives (only use frozen pixel) */
		in = i - 2; jn = j + 0; if (isfrozen2d(in, jn, dims, Frozen)) {
			Txm2 = T[mindex2(in, jn, dims[0])];  TSxm2 = TS[mindex2(in, jn, dims[0])];
		} else {
			Txm2 = INF;  TSxm2 = INF;
		}
		in = i + 2; jn = j + 0; if (isfrozen2d(in, jn, dims, Frozen)) {
			Txp2 = T[mindex2(in, jn, dims[0])];  TSxp2 = TS[mindex2(in, jn, dims[0])];
		} else {
			Txp2 = INF;  TSxp2 = INF;
		}
		in = i + 0; jn = j - 2; if (isfrozen2d(in, jn, dims, Frozen)) {
			Tym2 = T[mindex2(in, jn, dims[0])];  TSym2 = TS[mindex2(in, jn, dims[0])];
		} else {
			Tym2 = INF;  TSym2 = INF;
		}
		in = i + 0; jn = j + 2; if (isfrozen2d(in, jn, dims, Frozen)) {
			Typ2 = T[mindex2(in, jn, dims[0])];  TSyp2 = TS[mindex2(in, jn, dims[0])];
		} else {
			Typ2 = INF;  TSyp2 = INF;
		} // Allocate the specific values to the nodes involved in the second order differences. -wdd
		  // If it is frozen, use its values; If not, specify it with 'INF'. -wdd
		  // These are for the horizontal and vertical nodes (2-order). -wdd

		if (usecross) {
			in = i - 2; jn = j - 2; if (isfrozen2d(in, jn, dims, Frozen)) {
				Tr1m2 = T[mindex2(in, jn, dims[0])];  TSr1m2 = TS[mindex2(in, jn, dims[0])];
			} else {
				Tr1m2 = INF;  TSr1m2 = INF;
			}
			in = i - 2; jn = j + 2; if (isfrozen2d(in, jn, dims, Frozen)) {
				Tr2m2 = T[mindex2(in, jn, dims[0])];  TSr2m2 = TS[mindex2(in, jn, dims[0])];
			} else {
				Tr2m2 = INF;  TSr2m2 = INF;
			}
			in = i + 2; jn = j - 2; if (isfrozen2d(in, jn, dims, Frozen)) {
				Tr2p2 = T[mindex2(in, jn, dims[0])];  TSr2p2 = TS[mindex2(in, jn, dims[0])];
			} else {
				Tr2p2 = INF;  TSr2p2 = INF;
			}
			in = i + 2; jn = j + 2; if (isfrozen2d(in, jn, dims, Frozen)) {
				Tr1p2 = T[mindex2(in, jn, dims[0])];  TSr1p2 = TS[mindex2(in, jn, dims[0])];
			} else {
				Tr1p2 = INF;  TSr1p2 = INF;
			}
		} // Allocate the specific values to the nodes involved in the second order differences. -wdd
		  // If it is frozen, use its values; If not, specify it with 'INF'. -wdd
		  // These are for these nodes in the 2 inclined directions (2-order). -wdd

		// Similarly, here we can summarize it like follows:
		// to be more visual, these order should be like this:
		//      Tr1m2                   Tym2                      Tr2p2
                //             Tpatch_2_2    Tpatch_3_2     Tpatch_4_2
                //      Txm2   Tpatch_2_3   (Tpatch_3_3)    Tpatch_4_3    Txp2
                //             Tpatch_2_4    Tpatch_3_4     Tpatch_4_4
		//      Tr2m2                   Typ2                      Tr1p2

                //              Tpatch_3_3 --- Tpatch_i_j   -wdd

		Tm2[0] = 0; Tm2[1] = 0; Tm2[2] = 0; Tm2[3] = 0;
		TSm2[0] = 0; TSm2[1] = 0; TSm2[2] = 0; TSm2[3] = 0;   // -wdd
		/*pixels with a pixeldistance 2 from the center must be */
		/*lower in value otherwise use other side or first order */
		ch1 = (Txm2 < Tpatch_2_3) && IsFinite(Tpatch_2_3); ch2 = (Txp2 < Tpatch_4_3) && IsFinite(Tpatch_4_3);
		// Here 'Txm2 < Tpatch_2_3' and 'Txp2 < Tpatch_4_3' make sense, because the middle T should be determined by the neighboring T at the source sides.
		// Thus, this condition should be guaranteed at the two sizes of the sources. But only one side should be used (see follows). -wdd 
                // This case is for the horizontal direction. -wdd

		// performs a logical operation involving two conditions and assigns the result to ch1.  -wdd
		// Txm2 < Tpatch_2_3: This part compares the value of Txm2 with Tpatch_2_3. It returns true if Txm2 is less than Tpatch_2_3 and false otherwise. -wdd
		// For ch2, it is the same.  -wdd
		
		// Order[0] - indicates the horizontal direction.
                // Order[1] - indicates the vertical direction.
                // Order[2] - indicates the inclined direction with angle 45.
                // Order[3] - indicates the inclined direction with angle 135.

		if (ch1 && ch2) {
			Tm2[0] = min( (4.0 * Tpatch_2_3 - Txm2) / 3.0 , (4.0 * Tpatch_4_3 - Txp2) / 3.0);  Order[0] = 2;
			TSm2[0] = min( (4.0 * TSpatch_2_3 - TSxm2) / 3.0 , (4.0 * TSpatch_4_3 - TSxp2) / 3.0);  Order[0] = 2;  // -wdd
		} else if (ch1) {
			Tm2[0] = (4.0 * Tpatch_2_3 - Txm2) / 3.0; Order[0] = 2;
			TSm2[0] = (4.0 * TSpatch_2_3 - TSxm2) / 3.0; Order[0] = 2;
		} else if (ch2) {
			Tm2[0] = (4.0 * Tpatch_4_3 - Txp2) / 3.0; Order[0] = 2;
			TSm2[0] = (4.0 * TSpatch_4_3 - TSxp2) / 3.0; Order[0] = 2;
		}

		ch1 = (Tym2 < Tpatch_3_2) && IsFinite(Tpatch_3_2); ch2 = (Typ2 < Tpatch_3_4) && IsFinite(Tpatch_3_4);
		// Here 'Tym2 < Tpatch_3_2' and 'Typ2 < Tpatch_3_4' make sense, because the middle T should be determined by the neighboring T at the source sides. -wdd
		// This case is for the vertical direction. -wdd

		if (ch1 && ch2) {
			Tm2[1] = min( (4.0 * Tpatch_3_2 - Tym2) / 3.0 , (4.0 * Tpatch_3_4 - Typ2) / 3.0); Order[1] = 2;
			TSm2[1] = min( (4.0 * TSpatch_3_2 - TSym2) / 3.0 , (4.0 * TSpatch_3_4 - TSyp2) / 3.0); Order[1] = 2;
		} else if (ch1) {
			Tm2[1] = (4.0 * Tpatch_3_2 - Tym2) / 3.0; Order[1] = 2;
			TSm2[1] = (4.0 * TSpatch_3_2 - TSym2) / 3.0; Order[1] = 2;
		} else if (ch2) {
			Tm2[1] = (4.0 * Tpatch_3_4 - Typ2) / 3.0; Order[1] = 2;
			TSm2[1] = (4.0 * TSpatch_3_4 - TSyp2) / 3.0; Order[1] = 2;
		}
		if (usecross) {
			ch1 = (Tr1m2 < Tpatch_2_2) && IsFinite(Tpatch_2_2); ch2 = (Tr1p2 < Tpatch_4_4) && IsFinite(Tpatch_4_4);
			// This case is for the inclined direction 1 (maybe 45 direction). -wdd
			if (ch1 && ch2) {
				Tm2[2] = min( (4.0 * Tpatch_2_2 - Tr1m2) / 3.0 , (4.0 * Tpatch_4_4 - Tr1p2) / 3.0); Order[2] = 2;
				TSm2[2] = min( (4.0 * TSpatch_2_2 - TSr1m2) / 3.0 , (4.0 * TSpatch_4_4 - TSr1p2) / 3.0); Order[2] = 2;
			} else if (ch1) {
				Tm2[2] = (4.0 * Tpatch_2_2 - Tr1m2) / 3.0; Order[2] = 2;
				TSm2[2] = (4.0 * TSpatch_2_2 - TSr1m2) / 3.0; Order[2] = 2;
			} else if (ch2) {
				Tm2[2] = (4.0 * Tpatch_4_4 - Tr1p2) / 3.0; Order[2] = 2;
				TSm2[2] = (4.0 * TSpatch_4_4 - TSr1p2) / 3.0; Order[2] = 2;
			}

			ch1 = (Tr2m2 < Tpatch_2_4) && IsFinite(Tpatch_2_4); ch2 = (Tr2p2 < Tpatch_4_2) && IsFinite(Tpatch_4_2);
			// This case is for the inclined direction 2 (maybe 135 direction). -wdd
			if (ch1 && ch2) {
				Tm2[3] = min( (4.0 * Tpatch_2_4 - Tr2m2) / 3.0 , (4.0 * Tpatch_4_2 - Tr2p2) / 3.0); Order[3] = 2;
				TSm2[3] = min( (4.0 * TSpatch_2_4 - TSr2m2) / 3.0 , (4.0 * TSpatch_4_2 - TSr2p2) / 3.0); Order[3] = 2;
			} else if (ch1) {
				Tm2[3] = (4.0 * Tpatch_2_4 - Tr2m2) / 3.0; Order[3] = 2;
				TSm2[3] = (4.0 * TSpatch_2_4 - TSr2m2) / 3.0; Order[3] = 2;
			} else if (ch2) {
				Tm2[3] = (4.0 * Tpatch_4_2 - Tr2p2) / 3.0; Order[3] = 2;
				TSm2[3] = (4.0 * TSpatch_4_2 - TSr2p2) / 3.0; Order[3] = 2;
			}
		}
	} // Order[0], Order[1], Order[2] and Order[3] equals to 2 indicates only the second-order difference is used. -wdd
	/*Calculate the distance using x and y direction */
	Coeff[0] = 0; Coeff[1] = 0; Coeff[2] = -1 / (max(pow2(Fij), eps));
        // Without using 'usecross' - ---- Order[0],Order[1]
	for (t = 0; t < 2; t++) {    // t = 0, 1. indicates that we only use the horizontal and vertical directions.   -wdd
		switch (Order[t]) {  // Order[0] = 1 <For the 1 order derivative>, Order[0] = 2 <For the 2 order derivative>. -wdd
				     // Because when we use second order, Order[0], Order[1] will equal to 2.
		case 1: // 1 order  -wdd
			Coeff[0] += 1; Coeff[1] += -2 * Tm[t]; Coeff[2] += pow2(Tm[t]);
			break;
		case 2: // 2 order  -wdd
			Coeff[0] += (2.2500); Coeff[1] += -2.0 * Tm2[t] * (2.2500); Coeff[2] += pow2(Tm2[t]) * (2.2500);
			break;
		}
	}
	roots(Coeff, ansroot);
	// calling a function named roots, which is likely intended to compute the roots of a polynomial defined by the 
	// coefficients in the Coeff array. The results (roots) are stored in the ansroot variable.  -wdd
	Tt = max(ansroot[0], ansroot[1]);
	/*Calculate the distance using the cross directions */
	if (usecross) {  // ---- Order[2],Order[3]
        // With using 'usecross' - ---- Order[2],Order[3] 
		/* Original Equation */
		/*    Coeff[0]=0; Coeff[1]=0; Coeff[2]=-1/(max(pow2(Fij),eps)) */
		Coeff[0] += 0; Coeff[1] += 0; Coeff[2] += -1 / (max(pow2(Fij), eps));
		for (t = 2; t < 4; t++) {    // t = 0, 1 - indicates not used cross and only use horizontal and vertical directions. -wdd
					     // t = 2, 3 - indicates use the horizontal, vertical, and the two inclined directions;  -wdd
					     // While use cross, it implicit denotes that the horizontal and vertical are used. -wdd
					     // But while do not use cross, it explicit denotes that the inclined are not used. -wdd
			switch (Order[t]) {  // Order[t] - Order[2] - 1 <For the 1 order derivative>, Order[t] - Order[3] - 2 <For the 1 order derivative>.  -wdd
			case 1:  // 1 order  -wdd // The value of Order[t] can only be 1 or 2. It means that during the calculation, we only choose one kind. -wdd
				Coeff[0] += c1; Coeff[1] += -2.0 * c1 * Tm[t]; Coeff[2] += c1 * pow2(Tm[t]);
				break;
			case 2: // 2 order  -wdd
				Coeff[0] += c1 * 2.25; Coeff[1] += -2 * c1 * Tm2[t] * (2.25); Coeff[2] += pow2(Tm2[t]) * c1 * 2.25;
				break;
			}
		}
		if (Coeff[0] > 0) {
			roots(Coeff, ansroot);
			Tt2 = max(ansroot[0], ansroot[1]);
			/*Select minimum distance value of both stensils */
			Tt = min(Tt, Tt2);
		}
	}
	/*Upwind condition check, current distance must be larger */
	/*then direct neighbours used in solution */
	/*(Will this ever happen?) */
	if (usecross) {
		for (q = 0; q < 4; q++) {  // q = 0, 1, 2, ... 3;  -wdd
			if (IsFinite(Tm[q]) && (Tt < Tm[q])) {
				Tt = Tm[minarray(Tm, 4)] + (1 / (max(Fij, eps)));
			}
		}
	} else {
		for (q = 0; q < 2; q++) { // q = 0, 1. -wdd
			if (IsFinite(Tm[q]) && (Tt < Tm[q])) {
				Tt = Tm[minarray(Tm, 2)] + (1 / (max(Fij, eps)));
			}
		}
	}

	// to calculate t* data. // -wdd
	Coeff_ts[0] = 0;  Coeff_ts[1] = -1 / (max(pow2(Fij) * Qij, eps));
	//Tt - Tm[t]

	for (t = 0; t < 2; t++) {    // t = 0, 1. indicates that we only use the horizontal and vertical directions.   -wdd
                switch (Order[t]) {  // Order[0] = 1 <For the 1 order derivative>, Order[0] = 2 <For the 2 order derivative>. -wdd
                                     // Because when we use second order, Order[0], Order[1] will equal to 2.
                case 1: // 1 order  -wdd
                        //Coeff[0] += 1; Coeff[1] += -2 * Tm[t]; Coeff[2] += pow2(Tm[t]);
			Coeff_ts[0] += Tt - Tm[t];  Coeff_ts[1] += -1 * (Tt - Tm[t]) * TSm[t];
                        break;
                case 2: // 2 order  -wdd
                        // Coeff[0] += (2.2500); Coeff[1] += -2.0 * Tm2[t] * (2.2500); Coeff[2] += pow2(Tm2[t]) * (2.2500);
			Coeff_ts[0] += (2.2500)*(Tt - Tm2[t]);  Coeff_ts[1] += -1 * (2.2500) * (Tt - Tm2[t]) * TSm[t];
                        break;
                }
        }

	TSt = -1 * Coeff_ts[1] / Coeff_ts[0];

	//cout << "Tt/TSt = " << Tt/TSt << endl;  // -wdd

        return {Tt, TSt};
	//return Tt;
}

//double CalculateDistance3D(double *T, double Fijk, int *dims, int i, int j, int k, bool usesecond, bool usecross, bool *Frozen) {
//CalculateDis3DResult CalculateDistance3D(double *T, double *TS, double Fijk, double Qijk, int *dims, int i, int j, int k, bool usesecond, bool usecross, bool *Frozen) {
CalculateDis3DResult CalculateDistance3D(double *T, double *TS, double Fijk, double Qijk, int *dims, double dzxy, int i, int j, int k, bool usesecond, bool usecross, bool *Frozen) {
	/* Loop variables */
	int q, t;

	/* Current location */
	int in, jn, kn;

	/* Derivatives */
	double Tm[18] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	double TSm[18] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};  // -wdd
	double Tm2[18] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	double TSm2[18] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};  // -wdd
	double Coeff[3];
	double Coeff_ts[2];

	/* local derivatives in distance image */
	double Txm1, Txm2, Txp1, Txp2;
	double Tym1, Tym2, Typ1, Typ2;
	double Tzm1, Tzm2, Tzp1, Tzp2;

	double TSxm1, TSxm2, TSxp1, TSxp2;
        double TSym1, TSym2, TSyp1, TSyp2;
        double TSzm1, TSzm2, TSzp1, TSzp2;  // -wdd

	/* local cross derivatives in distance image */
	double Tr2t1m1, Tr2t1m2, Tr2t1p1, Tr2t1p2;
	double Tr2t2m1, Tr2t2m2, Tr2t2p1, Tr2t2p2;
	double Tr2t3m1, Tr2t3m2, Tr2t3p1, Tr2t3p2;
	double Tr3t1m1, Tr3t1m2, Tr3t1p1, Tr3t1p2;
	double Tr3t2m1, Tr3t2m2, Tr3t2p1, Tr3t2p2;
	double Tr3t3m1, Tr3t3m2, Tr3t3p1, Tr3t3p2;
	double Tr4t1m1, Tr4t1m2, Tr4t1p1, Tr4t1p2;
	double Tr4t2m1, Tr4t2m2, Tr4t2p1, Tr4t2p2;
	double Tr4t3m1, Tr4t3m2, Tr4t3p1, Tr4t3p2;
	double Tr5t1m1, Tr5t1m2, Tr5t1p1, Tr5t1p2;
	double Tr5t2m1, Tr5t2m2, Tr5t2p1, Tr5t2p2;
	double Tr5t3m1, Tr5t3m2, Tr5t3p1, Tr5t3p2;
	double Tr6t1m1, Tr6t1m2, Tr6t1p1, Tr6t1p2;
	double Tr6t2m1, Tr6t2m2, Tr6t2p1, Tr6t2p2;
	double Tr6t3m1, Tr6t3m2, Tr6t3p1, Tr6t3p2;

	double TSr2t1m1, TSr2t1m2, TSr2t1p1, TSr2t1p2;
        double TSr2t2m1, TSr2t2m2, TSr2t2p1, TSr2t2p2;
        double TSr2t3m1, TSr2t3m2, TSr2t3p1, TSr2t3p2;
        double TSr3t1m1, TSr3t1m2, TSr3t1p1, TSr3t1p2;
        double TSr3t2m1, TSr3t2m2, TSr3t2p1, TSr3t2p2;
        double TSr3t3m1, TSr3t3m2, TSr3t3p1, TSr3t3p2;
        double TSr4t1m1, TSr4t1m2, TSr4t1p1, TSr4t1p2;
        double TSr4t2m1, TSr4t2m2, TSr4t2p1, TSr4t2p2;
        double TSr4t3m1, TSr4t3m2, TSr4t3p1, TSr4t3p2;
        double TSr5t1m1, TSr5t1m2, TSr5t1p1, TSr5t1p2;
        double TSr5t2m1, TSr5t2m2, TSr5t2p1, TSr5t2p2;
        double TSr5t3m1, TSr5t3m2, TSr5t3p1, TSr5t3p2;
        double TSr6t1m1, TSr6t1m2, TSr6t1p1, TSr6t1p2;
        double TSr6t2m1, TSr6t2m2, TSr6t2p1, TSr6t2p2;
        double TSr6t3m1, TSr6t3m2, TSr6t3p1, TSr6t3p2;

	double Tt, Tt2;
	double TSt, TSt2;

	/* Return values root of polynomial */
	double ansroot[2] = {0, 0};

	/* Order derivatives in a certain direction */
	int Order[18] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	/* Neighbours 4x2 */
	int ne[18] = { -1,  0,  0, 1, 0, 0, 0, -1,  0, 0, 1, 0, 0,  0, -1, 0, 0, 1};

	/* Stencil constants */
	double G1[18] = {1, 1, 1, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 0.3333333333333, 0.3333333333333, 0.5, 0.3333333333333, 0.3333333333333};
	double G2[18] = {2.250, 2.250, 2.250, 2.250, 1.125, 1.125, 2.250, 1.125, 1.125, 2.250, 1.125, 1.125, 1.125, 0.750, 0.750, 1.125, 0.750, 0.750};


	/*Get First order derivatives (only use frozen pixel) */
	in = i - 1; jn = j + 0; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
		Txm1 = T[mindex3(in, jn, kn, dims[0], dims[1])];   TSxm1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
	} else {
		Txm1 = INF;  TSxm1 = INF;
	}
	in = i + 1; jn = j + 0; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
		Txp1 = T[mindex3(in, jn, kn, dims[0], dims[1])];   TSxp1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
	} else {
		Txp1 = INF;  TSxp1 = INF;
	}
	in = i + 0; jn = j - 1; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
		Tym1 = T[mindex3(in, jn, kn, dims[0], dims[1])];   TSym1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
	} else {
		Tym1 = INF;  TSym1 = INF;
	}
	in = i + 0; jn = j + 1; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
		Typ1 = T[mindex3(in, jn, kn, dims[0], dims[1])];   TSyp1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
	} else {
		Typ1 = INF;  TSyp1 = INF;
	}
	in = i + 0; jn = j + 0; kn = k - 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
		Tzm1 = T[mindex3(in, jn, kn, dims[0], dims[1])];   TSzm1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
	} else {
		Tzm1 = INF;  TSzm1 = INF;
	}
	in = i + 0; jn = j + 0; kn = k + 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
		Tzp1 = T[mindex3(in, jn, kn, dims[0], dims[1])];   TSzp1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
	} else {
		Tzp1 = INF;  TSzp1 = INF;
	}

	if (usecross) {
		Tr2t1m1 = Txm1;  TSr2t1m1 = TSxm1;
		Tr2t1p1 = Txp1;  TSr2t1p1 = TSxp1;
		in = i - 0; jn = j - 1; kn = k - 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr2t2m1 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr2t2m1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr2t2m1 = INF;  TSr2t2m1 = INF;
		}
		in = i + 0; jn = j + 1; kn = k + 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr2t2p1 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr2t2p1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr2t2p1 = INF;  TSr2t2p1 = INF;
		}
		in = i - 0; jn = j - 1; kn = k + 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr2t3m1 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr2t3m1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr2t3m1 = INF;  TSr2t3m1 = INF;
		}
		in = i + 0; jn = j + 1; kn = k - 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr2t3p1 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr2t3p1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr2t3p1 = INF;  TSr2t3p1 = INF;
		}
		Tr3t1m1 = Tym1;  TSr3t1m1 = TSym1;
		Tr3t1p1 = Typ1;  TSr3t1p1 = TSyp1;
		in = i - 1; jn = j + 0; kn = k + 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr3t2m1 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr3t2m1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr3t2m1 = INF;  TSr3t2m1 = INF;
		}
		in = i + 1; jn = j + 0; kn = k - 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr3t2p1 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr3t2p1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr3t2p1 = INF;  TSr3t2p1 = INF;
		}
		in = i - 1; jn = j - 0; kn = k - 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr3t3m1 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr3t3m1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr3t3m1 = INF;  TSr3t3m1 = INF;
		}
		in = i + 1; jn = j + 0; kn = k + 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr3t3p1 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr3t3p1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr3t3p1 = INF;  TSr3t3p1 = INF;
		}
		Tr4t1m1 = Tzm1;  TSr4t1m1 = TSzm1;
		Tr4t1p1 = Tzp1;  TSr4t1p1 = TSzp1;
		in = i - 1; jn = j - 1; kn = k - 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr4t2m1 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr4t2m1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr4t2m1 = INF;  TSr4t2m1 = INF;
		}
		in = i + 1; jn = j + 1; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr4t2p1 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr4t2p1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr4t2p1 = INF;  TSr4t2p1 = INF;
		}
		in = i - 1; jn = j + 1; kn = k - 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr4t3m1 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr4t3m1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr4t3m1 = INF;  TSr4t3m1 = INF;
		}
		in = i + 1; jn = j - 1; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr4t3p1 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr4t3p1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr4t3p1 = INF;  TSr4t3p1 = INF;
		}
		Tr5t1m1 = Tr3t3m1;  TSr5t1m1 = TSr3t3m1;
		Tr5t1p1 = Tr3t3p1;  TSr5t1p1 = TSr3t3p1;
		in = i - 1; jn = j - 1; kn = k + 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr5t2m1 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr5t2m1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr5t2m1 = INF;  TSr5t2m1 = INF;
		}
		in = i + 1; jn = j + 1; kn = k - 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr5t2p1 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr5t2p1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr5t2p1 = INF;  TSr5t2p1 = INF;
		}
		in = i - 1; jn = j + 1; kn = k + 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr5t3m1 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr5t3m1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr5t3m1 = INF;  TSr5t3m1 = INF;
		}
		in = i + 1; jn = j - 1; kn = k - 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr5t3p1 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr5t3p1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr5t3p1 = INF;  TSr5t3p1 = INF;
		}
		Tr6t1m1 = Tr3t2p1;  TSr6t1m1 = TSr3t2p1;
		Tr6t1p1 = Tr3t2m1;  TSr6t1p1 = TSr3t2m1;
		in = i - 1; jn = j - 1; kn = k - 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr6t2m1 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr6t2m1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr6t2m1 = INF;  TSr6t2m1 = INF;
		}
		in = i + 1; jn = j + 1; kn = k + 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr6t2p1 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr6t2p1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr6t2p1 = INF;  TSr6t2p1 = INF;
		}
		in = i - 1; jn = j + 1; kn = k - 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr6t3m1 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr6t3m1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr6t3m1 = INF;  TSr6t3m1 = INF;
		}
		in = i + 1; jn = j - 1; kn = k + 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr6t3p1 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr6t3p1 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr6t3p1 = INF;  TSr6t3p1 = INF;
		}
	}

	/*The values in order is 0 if no neighbours in that direction */
	/*1 if 1e order derivatives is used and 2 if second order */
	/*derivatives are used */


	/*Make 1e order derivatives in x and y direction */
	Tm[0] = min(Txm1 , Txp1); 
	if (Tm[0] == Txm1){
		TSm[0] = TSxm1;
	} else {
		TSm[0] = TSxp1;
	}
	if (IsFinite(Tm[0])) {
		Order[0] = 1;
	} else {
		Order[0] = 0;
	}

	Tm[1] = min(Tym1 , Typ1);
        if (Tm[1] == Tym1){
                TSm[1] = TSym1;
        } else {
                TSm[1] = TSyp1;
        }
	if (IsFinite(Tm[1])) {
		Order[1] = 1;
	} else {
		Order[1] = 0;
	}

	Tm[2] = min(Tzm1 , Tzp1); 
	if (Tm[2] == Tzm1){
                TSm[2] = TSzm1;
        } else {
                TSm[2] = TSzp1;
        }
	if (IsFinite(Tm[2])) {
		Order[2] = 1;
	} else {
		Order[2] = 0;
	}

	/*Make 1e order derivatives in cross directions */
	if (usecross) {
		Tm[3] = Tm[0];   Order[3] = Order[0];
                TSm[3] = TSm[0]; Order[3] = Order[0];  // -wdd

		Tm[4] = min(Tr2t2m1 , Tr2t2p1);
	        if (Tm[4] == Tr2t2m1){
                TSm[4] = TSr2t2m1;
                } else {
                TSm[4] = TSr2t2p1;
                }	
		if (IsFinite(Tm[4])) {
			Order[4] = 1;
		} else {
			Order[4] = 0;
		}

		Tm[5] = min(Tr2t3m1 , Tr2t3p1);
	        if (Tm[5] == Tr2t3m1){
                TSm[5] = TSr2t3m1;
                } else {
                TSm[5] = TSr2t3p1;
                }	
		if (IsFinite(Tm[5])) {
			Order[5] = 1;
		} else {
			Order[5] = 0;
		}

		Tm[6] = Tm[1]; Order[6] = Order[1];
		TSm[6] = TSm[1]; Order[6] = Order[1];  // -wdd

		Tm[7] = min( Tr3t2m1 , Tr3t2p1);
	        if (Tm[7] == Tr3t2m1){
                TSm[7] = TSr3t2m1;
                } else {
                TSm[7] = TSr3t2p1;
                }	
		if (IsFinite(Tm[7])) {
			Order[7] = 1;
		} else {
			Order[7] = 0;
		}

		Tm[8] = min( Tr3t3m1 , Tr3t3p1);
	        if (Tm[8] == Tr3t3m1){
                TSm[8] = TSr3t3m1;
                } else {
                TSm[8] = TSr3t3p1;
                }	
		if (IsFinite(Tm[8])) {
			Order[8] = 1;
		} else {
			Order[8] = 0;
		}

		Tm[9] = Tm[2]; Order[9] = Order[2];
		TSm[9] = TSm[2]; Order[9] = Order[2];  // -wdd

		Tm[10] = min( Tr4t2m1 , Tr4t2p1);
	        if (Tm[10] == Tr4t2m1){
                TSm[10] = TSr4t2m1;
                } else {
                TSm[10] = TSr4t2p1;
                }	
		if (IsFinite(Tm[10])) {
			Order[10] = 1;
		} else {
			Order[10] = 0;
		}

		Tm[11] = min( Tr4t3m1 , Tr4t3p1);
	        if (Tm[11] == Tr4t3m1){
                TSm[11] = TSr4t3m1;
                } else {
                TSm[11] = TSr4t3p1;
                }	
		if (IsFinite(Tm[11])) {
			Order[11] = 1;
		} else {
			Order[11] = 0;
		}

		Tm[12] = Tm[8]; Order[12] = Order[8];
		TSm[12] = TSm[8]; Order[12] = Order[8]; // -wdd

		Tm[13] = min( Tr5t2m1 , Tr5t2p1);
	        if (Tm[13] == Tr5t2m1){
                TSm[13] = TSr5t2m1;
                } else {
                TSm[13] = TSr5t2p1;
                }	
		if (IsFinite(Tm[13])) {
			Order[13] = 1;
		} else {
			Order[13] = 0;
		}

		Tm[14] = min( Tr5t3m1 , Tr5t3p1);
	        if (Tm[14] == Tr5t3m1){
                TSm[14] = TSr5t3m1;
                } else {
                TSm[14] = TSr5t3p1;
                }	
		if (IsFinite(Tm[14])) {
			Order[14] = 1;
		} else {
			Order[14] = 0;
		}

		Tm[15] = Tm[7]; Order[15] = Order[7];
		TSm[15] = TSm[7]; Order[15] = Order[7];  // -wdd

		Tm[16] = min( Tr6t2m1 , Tr6t2p1);
	        if (Tm[16] == Tr6t2m1){
                TSm[16] = TSr6t2m1;
                } else {
                TSm[16] = TSr6t2p1;
                }	
		if (IsFinite(Tm[16])) {
			Order[16] = 1;
		} else {
			Order[16] = 0;
		}

		Tm[17] = min( Tr6t3m1 , Tr6t3p1);
	        if (Tm[17] == Tr6t3m1){
                TSm[17] = TSr6t3m1;
                } else {
                TSm[17] = TSr6t3p1;
                }	
		if (IsFinite(Tm[17])) {
			Order[17] = 1;
		} else {
			Order[17] = 0;
		}
	}

	// here here here // -wdd 2024.12.4 20:54
	
	/*Make 2e order derivatives */
	if (usesecond) {
		/*Get Second order derivatives (only use frozen pixel) */
		/*Get First order derivatives (only use frozen pixel) */
		in = i - 2; jn = j + 0; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Txm2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSxm2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Txm2 = INF;  TSxm2 = INF;
		}
		in = i + 2; jn = j + 0; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Txp2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSxp2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Txp2 = INF;  TSxp2 = INF;
		}
		in = i + 0; jn = j - 2; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tym2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSym2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tym2 = INF;  TSym2 = INF;
		}
		in = i + 0; jn = j + 2; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Typ2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSyp2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Typ2 = INF;  TSyp2 = INF;
		}
		in = i + 0; jn = j + 0; kn = k - 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tzm2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSzm2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tzm2 = INF;  TSzm2 = INF;
		}
		in = i + 0; jn = j + 0; kn = k + 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tzp2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSzp2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tzp2 = INF;  TSzp2 = INF;
		}

		if (usecross) {
			Tr2t1m2 = Txm2;  TSr2t1m2 = TSxm2;
			Tr2t1p2 = Txp2;  TSr2t1p2 = TSxp2;
			in = i - 0; jn = j - 2; kn = k - 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr2t2m2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr2t2m2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr2t2m2 = INF;  TSr2t2m2 = INF;
			}
			in = i + 0; jn = j + 2; kn = k + 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr2t2p2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr2t2p2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr2t2p2 = INF;  TSr2t2p2 = INF;
			}
			in = i - 0; jn = j - 2; kn = k + 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr2t3m2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr2t3m2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr2t3m2 = INF;  TSr2t3m2 = INF;
			}
			in = i + 0; jn = j + 2; kn = k - 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr2t3p2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr2t3p2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr2t3p2 = INF;  TSr2t3p2 = INF;
			}
			Tr3t1m2 = Tym2;  TSr3t1m2 = TSym2;
			Tr3t1p2 = Typ2;  TSr3t1p2 = TSyp2;
			in = i - 2; jn = j + 0; kn = k + 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr3t2m2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr3t2m2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr3t2m2 = INF;  TSr3t2m2 = INF;
			}
			in = i + 2; jn = j + 0; kn = k - 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr3t2p2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr3t2p2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr3t2p2 = INF;  TSr3t2p2 = INF;
			}
			in = i - 2; jn = j - 0; kn = k - 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr3t3m2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr3t3m2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr3t3m2 = INF;  TSr3t3m2 = INF;
			}
			in = i + 2; jn = j + 0; kn = k + 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr3t3p2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr3t3p2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr3t3p2 = INF;  TSr3t3p2 = INF;
			}
			Tr4t1m2 = Tzm2;  TSr4t1m2 = TSzm2;
			Tr4t1p2 = Tzp2;  TSr4t1p2 = TSzp2;
			in = i - 2; jn = j - 2; kn = k - 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr4t2m2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr4t2m2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr4t2m2 = INF;  TSr4t2m2 = INF;
			}
			in = i + 2; jn = j + 2; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr4t2p2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr4t2p2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr4t2p2 = INF;  TSr4t2p2 = INF;
			}
			in = i - 2; jn = j + 2; kn = k - 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr4t3m2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr4t3m2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr4t3m2 = INF;  TSr4t3m2 = INF;
			}
			in = i + 2; jn = j - 2; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr4t3p2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr4t3p2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr4t3p2 = INF;  TSr4t3p2 = INF;
			}
			Tr5t1m2 = Tr3t3m2;  TSr5t1m2 = TSr3t3m2;
			Tr5t1p2 = Tr3t3p2;  TSr5t1p2 = TSr3t3p2;
			in = i - 2; jn = j - 2; kn = k + 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr5t2m2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr5t2m2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr5t2m2 = INF;  TSr5t2m2 = INF;
			}
			in = i + 2; jn = j + 2; kn = k - 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr5t2p2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr5t2p2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr5t2p2 = INF;  TSr5t2p2 = INF;
			}
			in = i - 2; jn = j + 2; kn = k + 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr5t3m2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr5t3m2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr5t3m2 = INF;  TSr5t3m2 = INF;
			}
			in = i + 2; jn = j - 2; kn = k - 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr5t3p2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr5t3p2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr5t3p2 = INF;  TSr5t3p2 = INF;
			}
			Tr6t1m2 = Tr3t2p2;  TSr6t1m2 = TSr3t2p2;
			Tr6t1p2 = Tr3t2m2;  TSr6t1p2 = TSr3t2m2;
			in = i - 2; jn = j - 2; kn = k - 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr6t2m2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr6t2m2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr6t2m2 = INF;  TSr6t2m2 = INF;
			}
			in = i + 2; jn = j + 2; kn = k + 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr6t2p2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr6t2p2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr6t2p2 = INF;  TSr6t2p2 = INF;
			}
			in = i - 2; jn = j + 2; kn = k - 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr6t3m2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr6t3m2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr6t3m2 = INF;  TSr6t3m2 = INF;
			}
			in = i + 2; jn = j - 2; kn = k + 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr6t3p2 = T[mindex3(in, jn, kn, dims[0], dims[1])];  TSr6t3p2 = TS[mindex3(in, jn, kn, dims[0], dims[1])];
			} else { 
				Tr6t3p2 = INF;  TSr6t3p2 = INF;
			}
		}


		/*pixels with a pixeldistance 2 from the center must be */
		/*lower in value otherwise use other side or first order */

		Tm2[0] = second_derivative(Txm1, Txm2, Txp1, Txp2);  
		TSm2[0] = second_derivative_tt_ts(Txm1, Txm2, Txp1, Txp2, TSxm1, TSxm2, TSxp1, TSxp2);
		if (IsInf(Tm2[0])) {
			Tm2[0] = 0;   TSm2[0] = 0;
		} else {
			Order[0] = 2; // Order[0] = 2;
		}

		Tm2[1] = second_derivative(Tym1, Tym2, Typ1, Typ2);  
		TSm2[1] = second_derivative_tt_ts(Tym1, Tym2, Typ1, Typ2, TSym1, TSym2, TSyp1, TSyp2);
		if (IsInf(Tm2[1])) {
			Tm2[1] = 0;   TSm2[1] = 0;
		} else {
			Order[1] = 2; // Order[1] = 2; 
		}

		Tm2[2] = second_derivative(Tzm1, Tzm2, Tzp1, Tzp2);  
		TSm2[2] = second_derivative_tt_ts(Tzm1, Tzm2, Tzp1, Tzp2, TSzm1, TSzm2, TSzp1, TSzp2);
		if (IsInf(Tm2[2])) {
			Tm2[2] = 0;   TSm2[2] = 0;
		} else {
			Order[2] = 2;  // Order[2] = 2;
		}

		if (usecross) {
			Tm2[3] = Tm2[0];   Order[3] = Order[0];
			TSm2[3] = TSm2[0]; Order[3] = Order[0]; // -wdd
			Tm2[4] = second_derivative(Tr2t2m1, Tr2t2m2, Tr2t2p1, Tr2t2p2);  TSm2[4] = second_derivative(TSr2t2m1, TSr2t2m2, TSr2t2p1, TSr2t2p2);
			if (IsInf(Tm2[4])) {
				Tm2[4] = 0;   TSm2[4] = 0;
			} else {
				Order[4] = 2;
			}
			Tm2[5] = second_derivative(Tr2t3m1, Tr2t3m2, Tr2t3p1, Tr2t3p2);  TSm2[5] = second_derivative(TSr2t3m1, TSr2t3m2, TSr2t3p1, TSr2t3p2);
			if (IsInf(Tm2[5])) {
				Tm2[5] = 0;   TSm2[5] = 0;
			} else {
				Order[5] = 2;
			}

			Tm2[6] = Tm2[1];   Order[6] = Order[1];
			TSm2[6] = TSm2[1]; Order[6] = Order[1]; // -wdd
			Tm2[7] = second_derivative(Tr3t2m1, Tr3t2m2, Tr3t2p1, Tr3t2p2);  TSm2[7] = second_derivative(TSr3t2m1, TSr3t2m2, TSr3t2p1, TSr3t2p2);
			if (IsInf(Tm2[7])) {
				Tm2[7] = 0;   TSm2[7] = 0;
			} else {
				Order[7] = 2;
			}
			Tm2[8] = second_derivative(Tr3t3m1, Tr3t3m2, Tr3t3p1, Tr3t3p2);  TSm2[8] = second_derivative(TSr3t3m1, TSr3t3m2, TSr3t3p1, TSr3t3p2);
			if (IsInf(Tm2[8])) {
				Tm2[8] = 0;   TSm2[8] = 0;
			} else {
				Order[8] = 2;
			}

			Tm2[9] = Tm2[2];   Order[9] = Order[2];
			TSm2[9] = TSm2[2]; Order[9] = Order[2];  // -wdd
			Tm2[10] = second_derivative(Tr4t2m1, Tr4t2m2, Tr4t2p1, Tr4t2p2);  TSm2[10] = second_derivative(TSr4t2m1, TSr4t2m2, TSr4t2p1, TSr4t2p2);
			if (IsInf(Tm2[10])) {
				Tm2[10] = 0;   TSm2[10] = 0;
			} else {
				Order[10] = 2;
			}
			Tm2[11] = second_derivative(Tr4t3m1, Tr4t3m2, Tr4t3p1, Tr4t3p2);  TSm2[11] = second_derivative(TSr4t3m1, TSr4t3m2, TSr4t3p1, TSr4t3p2);
			if (IsInf(Tm2[11])) {
				Tm2[11] = 0;   TSm2[11] = 0;
			} else {
				Order[11] = 2;
			}

			Tm2[12] = Tm2[8];   Order[12] = Order[8];
			TSm2[12] = TSm2[8]; Order[12] = Order[8];  // -wdd
			Tm2[13] = second_derivative(Tr5t2m1, Tr5t2m2, Tr5t2p1, Tr5t2p2);   TSm2[13] = second_derivative(TSr5t2m1, TSr5t2m2, TSr5t2p1, TSr5t2p2);
			if (IsInf(Tm2[13])) {
				Tm2[13] = 0;   TSm2[13] = 0;
			} else {
				Order[13] = 2;
			}
			Tm2[14] = second_derivative(Tr5t3m1, Tr5t3m2, Tr5t3p1, Tr5t3p2);   TSm2[14] = second_derivative(TSr5t3m1, TSr5t3m2, TSr5t3p1, TSr5t3p2);
			if (IsInf(Tm2[14])) {
				Tm2[14] = 0;   TSm2[14] = 0;
			} else {
				Order[14] = 2;
			}

			Tm2[15] = Tm2[7]; Order[15] = Order[7];
			TSm2[15] = TSm2[7]; Order[15] = Order[7];  // -wdd
			Tm2[16] = second_derivative(Tr6t2m1, Tr6t2m2, Tr6t2p1, Tr6t2p2);   TSm2[16] = second_derivative(TSr6t2m1, TSr6t2m2, TSr6t2p1, TSr6t2p2);
			if (IsInf(Tm2[16])) {
				Tm2[16] = 0;   TSm2[16] = 0;
			} else {
				Order[16] = 2;
			}
			Tm2[17] = second_derivative(Tr6t3m1, Tr6t3m2, Tr6t3p1, Tr6t3p2);   TSm2[17] = second_derivative(TSr6t3m1, TSr6t3m2, TSr6t3p1, TSr6t3p2);
			if (IsInf(Tm2[17])) {
				Tm2[17] = 0;   TSm2[17] = 0;
			} else {
				Order[17] = 2;
			}
		}

	}

	/*Calculate the distance using x and y direction */
	Coeff[0] = 0; Coeff[1] = 0; Coeff[2] = -1 / (max(pow2(Fijk), eps));

	for (t = 0; t < 3; t++) {  // t = 0, 1, 2
		switch (Order[t]) {
		case 1:
			//Coeff[0] += G1[t]; Coeff[1] += -2.0 * Tm[t] * G1[t]; Coeff[2] += pow2(Tm[t]) * G1[t];
			Coeff[0] += G1[t] * (1 / (dzxy * dzxy)); Coeff[1] += -2.0 * Tm[t] * (1 / (dzxy * dzxy)) * G1[t]; Coeff[2] += pow2(Tm[t] / dzxy) * G1[t];
			break;
		case 2:
			//Coeff[0] += G2[t]; Coeff[1] += -2.0 * Tm2[t] * G2[t]; Coeff[2] += pow2(Tm2[t]) * G2[t];
			Coeff[0] += G2[t] * (1 / (dzxy * dzxy)); Coeff[1] += -2.0 * Tm2[t] * (1 / (dzxy * dzxy)) * G2[t]; Coeff[2] += pow2(Tm2[t] / dzxy) * G2[t];
			break;
		}
	}

	// double G1[18] = {1, 1, 1, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 0.3333333333333, 0.3333333333333, 0.5, 0.3333333333333, 0.3333333333333};
        // double G2[18] = {2.250, 2.250, 2.250, 2.250, 1.125, 1.125, 2.250, 1.125, 1.125, 2.250, 1.125, 1.125, 1.125, 0.750, 0.750, 1.125, 0.750, 0.750};

	roots(Coeff, ansroot);
	Tt = max(ansroot[0], ansroot[1]);

	/*Calculate the distance using the cross directions */
	if (usecross) {

		for (q = 1; q < 6; q++) {
			/* Original Equation */
			/*    Coeff[0]=0; Coeff[1]=0; Coeff[2]=-1/(max(pow2(Fijk),eps)) */
			Coeff[0] += 0; Coeff[1] += 0; Coeff[2] += -1 / (max(pow2(Fijk), eps));

			for (t = q * 3; t < ((q + 1) * 3); t++) {
				switch (Order[t]) {
				case 1:
					Coeff[0] += G1[t]; Coeff[1] += -2.0 * Tm[t] * G1[t]; Coeff[2] += pow2(Tm[t]) * G1[t];
					//Coeff[0] += G1[t] * (1 / (dzxy * dzxy)); Coeff[1] += -2.0 * Tm[t] * (1 / (dzxy * dzxy)) * G1[t]; Coeff[2] += pow2(Tm[t] / dzxy) * G1[t];
					break;
				case 2:
					//Coeff[0] += G2[t]; Coeff[1] += -2.0 * Tm2[t] * G2[t]; Coeff[2] += pow2(Tm2[t]) * G2[t];
					Coeff[0] += G2[t]; Coeff[1] += -2.0 * Tm2[t] * G2[t]; Coeff[2] += pow2(Tm2[t]) * G2[t];
					break;
				}
			}
			/*Select maximum root solution and minimum distance value of both stensils */
			if (Coeff[0] > 0) {
				roots(Coeff, ansroot);
				Tt2 = max(ansroot[0], ansroot[1]);
				Tt = min(Tt, Tt2);
			}
		}
	}

	/*Upwind condition check, current distance must be larger */
	/*then direct neighbours used in solution */
	/*(Will this ever happen?) */
	if (usecross) {
		for (q = 0; q < 18; q++) {
			if (IsFinite(Tm[q]) && (Tt < Tm[q])) {
				Tt = Tm[minarray(Tm, 18)] + (1 / (max(Fijk, eps)));
			}
		}
	} else {
		for (q = 0; q < 3; q++) {  // q = 0, 1, 2
			if (IsFinite(Tm[q]) && (Tt < Tm[q])) {
				Tt = Tm[minarray(Tm, 3)] + (1 / (max(Fijk, eps)));
				//cout << "Tt = " << Tt << endl; // -wdd
			}
		}
	}

	// to calculate t* data. // -wdd
	Coeff_ts[0] = 0;  Coeff_ts[1] = -1 / (max(pow2(Fijk) * Qijk, eps));

	for (t = 0; t < 3; t++) {    // t = 0, 1, 2. indicates that we only use the horizontal and vertical directions.   -wdd
                switch (Order[t]) {  // Order[0] = 1 <For the 1 order derivative>, Order[0] = 2 <For the 2 order derivative>. -wdd
                                     // Because when we use second order, Order[0], Order[1] will equal to 2.
                case 1: // 1 order  -wdd
                        //Coeff[0] += 1; Coeff[1] += -2 * Tm[t]; Coeff[2] += pow2(Tm[t]);
                        //Coeff_ts[0] += Tt - Tm[t];  Coeff_ts[1] += -1 * (Tt - Tm[t]) * TSm[t];
			//Coeff_ts[0] += Tt - Tm[t]; Coeff_ts[1] += -1 * (Tt - Tm[t]) * TSm[t];
			Coeff_ts[0] += (Tt - Tm[t]) * (1 / (dzxy * dzxy)); Coeff_ts[1] += -1 * (Tt - Tm[t]) * (1 / (dzxy * dzxy)) * TSm[t];
                        break;
                case 2: // 2 order  -wdd
                        // Coeff[0] += (2.2500); Coeff[1] += -2.0 * Tm2[t] * (2.2500); Coeff[2] += pow2(Tm2[t]) * (2.2500);
                        // Coeff_ts[0] += (2.2500)*(Tt - Tm2[t]);  Coeff_ts[1] += -1 * (2.2500) * (Tt - Tm2[t]) * TSm[t];
			Coeff_ts[0] += (2.2500)*(Tt - Tm2[t]) * (1 / (dzxy * dzxy));  Coeff_ts[1] += -1 * (2.2500) * (Tt - Tm2[t]) * (1 / (dzxy * dzxy)) * TSm2[t];
                        break;
                }
        }

	TSt = -1 * Coeff_ts[1] / Coeff_ts[0];
        //if (usecross) {
        //        for (q = 0; q < 18; q++) {
        //                if (IsFinite(Tm[q]) && (Tt < Tm[q])) {
        //                        Tt = Tm[minarray(Tm, 18)] + (1 / (max(Fijk, eps)));
        //                }
        //        }
        //} else {
        //        for (q = 0; q < 3; q++) {
        //                if (IsFinite(Tm[q]) && (Tt < Tm[q])) {
        //                        Tt = Tm[minarray(Tm, 3)] + (1 / (max(Fijk, eps)));
        //                }
        //        }
        //}
	
	//cout << "Tt = " << Tt << endl;  // -wdd

	//TSt = -1 * Coeff_ts[1] / Coeff_ts[0];
	//TSt = Tt;
	//return Tt;
	return {Tt, TSt};
}

double second_derivative(double Txm1, double Txm2, double Txp1, double Txp2) {
	bool ch1, ch2;
	double Tm;
	Tm = INF;
	ch1 = (Txm2 < Txm1) && IsFinite(Txm1); ch2 = (Txp2 < Txp1) && IsFinite(Txp1);
	if (ch1 && ch2) {
		Tm = min( (4.0 * Txm1 - Txm2) / 3.0 , (4.0 * Txp1 - Txp2) / 3.0);
	} else if (ch1) {
		Tm = (4.0 * Txm1 - Txm2) / 3.0;
	} else if (ch2) {
		Tm = (4.0 * Txp1 - Txp2) / 3.0;
	}
	return Tm;
}

double second_derivative_tt_ts(double Txm1, double Txm2, double Txp1, double Txp2, double TSxm1, double TSxm2, double TSxp1, double TSxp2) {
        bool ch1, ch2;
        double Tm;
	double TSm;
	double Tml, Tmr, TSml, TSmr;
        Tm = INF; TSm = INF;
        ch1 = (Txm2 < Txm1) && IsFinite(Txm1); ch2 = (Txp2 < Txp1) && IsFinite(Txp1);
        if (ch1 && ch2) {
		Tml =  (4.0 * Txm1  - Txm2)  / 3.0; Tmr =  (4.0 * Txp1  - Txp2)  / 3.0;
		TSml = (4.0 * TSxm1 - TSxm2) / 3.0; TSmr = (4.0 * TSxp1 - TSxp2) / 3.0;
                Tm = min(Tml, Tmr);
		if (Tm == Tml) {
			TSm = TSml;
		} else {
			TSm = TSmr;
		}
        } else if (ch1) {
                Tm =  (4.0 * Txm1  - Txm2)  / 3.0;
		TSm = (4.0 * TSxm1 - TSxm2) / 3.0;  // -wdd
        } else if (ch2) {
                Tm =  (4.0 * Txp1  - Txp2)  / 3.0;
		TSm = (4.0 * TSxp1 - TSxp2) / 3.0;     // -wdd
        }
        return TSm;
}

// calculate the average velocity near the source, note that the indexs are different from the original code. Here
// we have the z, x, y system
// Dongzhuo June 5 2015
double averageVelocityAroundSource2D(int bz, int bx, int izsrcn, int ixsrcn, int *dims, double *F) {  
	// bz = 2; bx = 2.  -wdd
	// izsrcn - the index of source in x direction; ixsrcn - the index of source in x direction. -wdd
	double count = 0.0;
	double velSum = 0.0;
	for (int i = 0; i <= 2 * bx; i++ ) {
		// i = 0, 1, 2, 3, 4. -wdd
		for (int k = 0; k <= 2 * bz; k++) {
			// k = 0, 1, 2, 3, 4. -wdd
			if (k >= 0 && k <= dims[0] && i >= 0 && i <= dims[1]) {
				// It is absolute that this condition is naturaly satisfied due to that it usually has a larger element number. -wdd
				// Good. Now we know that 'dims[0]' and 'dims[1]' are respectively the total number of grid
				// in z and x directions. -wdd 'dims[0]' - z; 'dims[1]' - x;
				velSum = velSum + F[mindex2(izsrcn + (k - bz), ixsrcn + (i - bx), dims[0])];
				// 'k - bz' - '-2', '-1', '0', '1', '2'; 'i - bx' - '-2', '-1', '0', '1', '2';
				// here 'mindex2' - int mindex2(int row, int col, int num_cols) {return row * num_cols + col;}
				// we should determine here it is 'dims[0]'-total row; or 'dims[1]'-total column. -wdd
				// And how to define 'izsrcn'-row number of source; and 'ixsrcn'-column number of source. -wdd
				// What is 'F' ? Is it the velocity model, which rewrite the 2D velocity model into the vector format. -wdd
				// Here we guess that 'F' is velocity model, which is written as a vector format. Row1: 1eft to right,
				// Then Row2: left to right; then, row3: left to right.  -wdd
				//
				std::cout << "dims[0] = " << dims[0] << std::endl;

				count ++;
			}
		}
	}
	double aveVel = velSum / count;
	return aveVel;
}

double averageAttenuationAroundSource2D(int bz, int bx, int izsrcn, int ixsrcn, int *dims, double *Q) {
	double count = 0.0;
	double attSum = 0.0;
	for (int i = 0; i <= 2 * bx; i++ ) {
		for (int k = 0; k <= 2 * bz; k++) {
			if (k >= 0 && k <= dims[0] && i >= 0 && i <= dims[1]) {
				attSum = attSum + Q[mindex2(izsrcn + (k - bz), ixsrcn + (i - bx), dims[0])];
				std::cout << "dims[0] = " << dims[0] << std::endl;
				count ++;
			}
		}
	}
	double aveAtt = attSum / count;
	return aveAtt;
}

// calculate time table near source based upon the average veocity in the box
// Dongzhuo June 5 2015
double CalculateTConstant2D(double aveVel, int izsrcn, int ixsrcn, int zloc, int xloc) {
	return sqrt(pow((zloc - izsrcn), 2) + pow((xloc - ixsrcn), 2)) / aveVel;
	// pow((zloc - izsrcn), 2): Computes the square of the difference in the z-coordinates (zloc and izsrcn). -wdd
	// pow((xloc - ixsrcn), 2): Computes the square of the difference in the x-coordinates (xloc and ixsrcn). -wdd
	// sqrt(...): Takes the square root of the sum of the squared differences.
	// The function 'pow' is a mathematical function in C++ that calculates the power of a number. -wdd
	// Here why not consider the size of the element, or the size of the grid in x and z directions. -wdd
}

double CalculateTSConstant2D(double aveVel, double aveAtt, int izsrcn, int ixsrcn, int zloc, int xloc) {
	//return sqrt(pow((zloc - izsrcn), 2) + pow((xloc - ixsrcn), 2)) / aveVel;
	return sqrt(pow((zloc - izsrcn), 2) + pow((xloc - ixsrcn), 2)) / (aveVel * aveAtt);
}

// calculate the average velocity near the source, note that the indexs are different from the original code. Here
// we have the z, x, y system
// Dongzhuo June 5 2015
double averageVelocityAroundSource3D(int bz, int bx, int by, int izsrcn, int ixsrcn, int iysrcn, int *dims, double *F) {
	double count = 0.0;
	double velSum = 0.0;
	for (int j = 0; j <= 2 * by; j++) {
		for (int i = 0; i <= 2 * bx; i++ ) {
			// 2 * by (bx, bz) indicates the two sides of the source. // -wdd
			for (int k = 0; k <= 2 * bz; k++) {
				if (k >= 0 && k <= dims[0] && i >= 0 && i <= dims[1] && j >= 0 && j <= dims[2]) {
					// dims[0] - nz; dims[1] - nx; dims[2] - ny; -wdd
					velSum = velSum + F[mindex3(izsrcn + (k - bz), ixsrcn + (i - bx), iysrcn + (j - by), dims[0], dims[1])];
					// dims[0] -- nz;  dims[1] -- nx;
					// cout << "izsrcn = " << izsrcn << endl; // izsrcn - 125
					// cout << "ixsrcn = " << ixsrcn << endl; // ixsrcn - 75
					// cout << "iysrcn = " << iysrcn << endl; // iysrcn - 1
					// Thus, 'F' should be the velocity vector (modified by the 3D velocity model)  // -wdd				
					count ++;
				}
			}
		}
	}
	double aveVel = velSum / count;
	return aveVel;
}

double averageAttenuationAroundSource3D(int bz, int bx, int by, int izsrcn, int ixsrcn, int iysrcn, int *dims, double *Q) {
	double count = 0.0;
	double attSum = 0.0;
	for (int j = 0; j <= 2 * by; j++) {
		for (int i = 0; i <= 2 * bx; i++ ) {
			for (int k = 0; k <= 2 * bz; k++) {
				if (k >= 0 && k <= dims[0] && i >= 0 && i <= dims[1] && j >= 0 && j <= dims[2]) {
					attSum = attSum + Q[mindex3(izsrcn + (k - bz), ixsrcn + (i - bx), iysrcn + (j - by), dims[0], dims[1])];
					count ++;
				}
			}
		}
	}
	double aveAtt = attSum / count;
	return aveAtt;
}

// calculate time table near source based upon the average veocity in the box
// Dongzhuo June 5 2015
//double CalculateTConstant3D(double aveVel, int izsrcn, int ixsrcn, int iysrcn, int zloc, int xloc, int yloc) {
double CalculateTConstant3D(double aveVel, int izsrcn, int ixsrcn, int iysrcn, double dzxy, int zloc, int xloc, int yloc) {
	//return sqrt(pow((zloc - izsrcn), 2) + pow((xloc - ixsrcn), 2) + pow((yloc - iysrcn), 2)) / aveVel;
	return sqrt(pow(((zloc - izsrcn) * dzxy), 2) + pow(((xloc - ixsrcn) * dzxy), 2) + pow(((yloc - iysrcn) * dzxy), 2)) / aveVel;
}

//double CalculateTSConstant3D(double aveVel, double aveAtt, int izsrcn, int ixsrcn, int iysrcn, int zloc, int xloc, int yloc) {
double CalculateTSConstant3D(double aveVel, double aveAtt, int izsrcn, int ixsrcn, int iysrcn, double dzxy, int zloc, int xloc, int yloc) {
	//return sqrt(pow((zloc - izsrcn), 2) + pow((xloc - ixsrcn), 2) + pow((yloc - iysrcn), 2)) / aveVel;
        //return sqrt(pow((zloc - izsrcn), 2) + pow((xloc - ixsrcn), 2) + pow((yloc - iysrcn), 2)) / (aveVel * aveAtt);
        return sqrt(pow(((zloc - izsrcn) * dzxy), 2) + pow(((xloc - ixsrcn) * dzxy), 2) + pow(((yloc - iysrcn) * dzxy), 2)) / (aveVel * aveAtt);
}
