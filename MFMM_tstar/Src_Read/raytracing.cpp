// ---------------
//   Dongzhuo Li
//   March, 2015
// ---------------


// This program is to trace the rays
// and compute the ray length in each
// cell giving a time table.

// by Shaoyu Lu
// Oct, 2013 at Stanford University

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <stdlib.h>
#include "raytracing_aux.h"

using namespace std;

/* Main program */

void raytracing(int dimension, int nx, int ny, int nz, double dx, double dy, double dz, double xmin, double ymin, double zmin, \
                int invnx, int invny, int invnz, double invdx, double invdy, double invdz, double invxmin, double invymin, double invzmin, \
                double xsrc, double ysrc, double zsrc, double xstart, double ystart, double zstart, double *timeTable, \
                vector< vector<double> > &rayPoints) {

	char DimensionIndex;


	if (dimension == 2) {
		DimensionIndex = 'a';
	}
	if (dimension == 3) {
		DimensionIndex = 'b';
	}

	double x[nx], y[ny], z[nz];
	ConstructGrids(x, y, z, nx, ny, nz, dx, dy, dz, xmin, ymin, zmin);
	double xmax = x[nx - 1];
	double ymax = y[ny - 1];
	double zmax = z[nz - 1];
	double dist = dx / 5.0;
	double totalRayLength = 0.0;
	TraceRayConstDist(rayPoints, totalRayLength, timeTable, dist, xsrc, ysrc, zsrc,
	                  xstart, ystart, zstart, x, y, z, nx, ny, nz);

}


/*
 * Function: TraceRayConstDist
 * Usage: TraceRayConstDist(<matrix> rayPoints,totalRayLength,
 * <cube> timeTable,xsrc, ysrc, zsrc, xstart, ystart, zstart,
 * <vec> x, <vec> y, <vec> z);
 * Outputs: <matrix> rayPoints: (rayx, rayy, rayz)
 *          totalRayLength: raylength from start point to source
 *          rayx, rayy, rayz: coordinates of intercept of the
 *          ray with cell boundaries
 * --------------------------
 * Trace the ray from starting point to source point
 * with timetable gradients and compute the raylenth
 * in each cell
 */
void TraceRayConstDist(
  vector< vector<double> > &rayPoints,
  double &totalRayLength,
  double *timeTable,
  double dist,
  double xsrc,
  double ysrc,
  double zsrc,
  double xstart,
  double ystart,
  double zstart,
  double *x,
  double *y,
  double *z,
  int nx,
  int ny,
  int nz
) {
	vector<double> onePoint(3);
	double xmin, ymin, zmin, xmax, ymax, zmax, dx, dy, dz;
	ExtractGridInfor(dx, dy, dz, xmin, ymin, zmin, xmax, ymax, zmax, x, y, z, nx, ny, nz);
	// Locate the source box
	int ixsrcn, iysrcn, izsrcn;
	LocateNode(ixsrcn, iysrcn, izsrcn, dx, dy, dz, xmin, ymin, zmin, xsrc, ysrc, zsrc);
	vector<int> sourceBoxIndexs = SourceBox(nx, ny, nz, ixsrcn, iysrcn, izsrcn);
	// Initialize from start point
	double xp = xstart, yp = ystart, zp = zstart;
	onePoint[0] = xp;
	onePoint[1] = yp;
	onePoint[2] = zp;
	rayPoints.push_back(onePoint);

	// Find next point while not in the sourcebox
	while (!TouchSourceBox(sourceBoxIndexs, x, y, z, xp, yp, zp)) {
		if (!IsInBound(xmin, xmax, ymin, ymax, zmin, zmax, xp, yp, zp)) {
			cout << xp <<" "<< yp <<" "<< zp <<" "<<"out of the velocity domain!!!" << endl;
			break;
		}
		double gx, gy, gz, xnex, ynex, znex;
		ComputeGradient(gx, gy, gz, timeTable, xp, yp, zp, dx, dy, dz, \
		                xmin, xmax, ymin, ymax, zmin, zmax, nx, ny, nz);
		if (!FindNextPointConstDist(xnex, ynex, znex, xp, yp, zp, gx, gy, gz, dist)) {
			cout << "Stop!" << endl; //Arrived the source point
			break;
		} else {
			xp = xnex;
			yp = ynex;
			zp = znex;
		}
		onePoint[0] = xp;
		onePoint[1] = yp;
		onePoint[2] = zp;
		rayPoints.push_back(onePoint);
	}
	// go straightly to the source point
	vector<double> srcPoint(3);
	srcPoint[0] = xsrc;
	srcPoint[1] = ysrc;
	srcPoint[2] = zsrc;

	rayPoints.push_back(srcPoint);

	totalRayLength = (rayPoints.size() - 2) * dist +
	                 OneLength(onePoint, srcPoint);

}



/*
 * Function: ComputeRayLength
 * Usage: ComputeRayLength(<matrix> CrossCellIndex,
 * <matrix> rayLength, <matrix> newRayPoints,
 * <matrix> interceptPoints, <matrix> rayPoints,
 * <char> DimensionIndex, <vec> invx, <vec> invy,
 * <vec> invz);
 * Outputs: CrossCellIndex: ixc, iyc, izc
 *          rayLength: in each crossed cell
 *          newRayPoints: rayPoints with interceptions
 *          interceptPoints: points on cell edges
 * --------------------------
 * Computes raylength in each cell and gives
 * cell indexs
 */
void ComputeRayLength(
  double *G1,
  vector< vector<double> > &newRayPoints,
  vector< vector<double> > &rayPoints,
  char DimensionIndex,
  double *invx,
  double *invy,
  double *invz,
  int invnx,
  int invny,
  int invnz
) {
#define G1CUBE(icz,icx,icy) G1[icy*(invnz-1)*(invnx-1)+icx*(invnz-1)+icz]
#define G1MAT(icz,icx) G1[icx*(invnz-1)+icz]
	vector< vector<double> > interceptPoints;
	double invdx, invdy, invdz, invxmin, invymin, invzmin,
	       invxmax, invymax, invzmax;
	ExtractGridInfor(invdx, invdy, invdz, invxmin, invymin, invzmin, invxmax,
	                 invymax, invzmax, invx, invy, invz, invnx, invny, invnz);
	// compute interceptions of the ray with cell boundaries
	newRayPoints = InterceptInvGrids(rayPoints, DimensionIndex, invx, invy, invz, invnx, invny, invnz);
	int i = 0, istart = 0;
	double secLength = 0.0;
	while (i < newRayPoints.size() - 2) {
		vector<double> onePoint1(newRayPoints[i]);
		vector<double> onePoint2(newRayPoints[i + 1]);
		vector<double> onePoint3(newRayPoints[i + 2]);
		secLength += OneLength(onePoint1, onePoint2);
		i++;
		//cout << "i = " << i << endl;
		// whether to stop accumulate length
		double xp1 = onePoint1[0];
		double yp1 = onePoint1[1];
		double zp1 = onePoint1[2];
		double xp = onePoint2[0];
		double yp = onePoint2[1];
		double zp = onePoint2[2];
		if (!IsInBound(invxmin, invxmax, invymin, invymax,
		               invzmin, invzmax, xp1, yp1, zp1)) {
			// don't accumulate length outside of inverse boundary
			secLength = 0.0;
			istart ++;
			continue;
		}
		int ixp = LocateGrid(xp, invdx, invxmin, invxmax);
		int iyp = LocateGrid(yp, invdy, invymin, invymax);
		int izp = LocateGrid(zp, invdz, invzmin, invzmax);
		double xintercept = invx[ixp];
		double yintercept = invy[iyp];
		double zintercept = invz[izp];
		bool jumpJudge = trueIntercept(onePoint1, onePoint2, onePoint3, DimensionIndex,
		                               xintercept, yintercept, zintercept);
		if (jumpJudge) {
			interceptPoints.push_back(onePoint2);
			int ixc, iyc, izc;
			LocateCell(ixc, iyc, izc, invdx, invdy, invdz,
			           invxmin, invymin, invzmin,
			           (onePoint1[0] + onePoint2[0]) / 2.0,
			           (onePoint1[1] + onePoint2[1]) / 2.0,
			           (onePoint1[2] + onePoint2[2]) / 2.0);
			if (invny == 1) {
				G1MAT(izc, ixc) = secLength;
			} else {
				G1CUBE(izc, ixc, iyc) = secLength;
			}
			secLength = 0.0;
			/*
			mexPrintf("\nInterception point (x = %f, y = %f, z = %f)\n",
			        onePoint2[0],onePoint2[1],onePoint2[2]);
			mexPrintf("\nG1(izc = %d, ixc = %d, iyc = %d) = %f\n",
			        izc+1,ixc+1,iyc+1,secLength);
			 */
			/*
			                cout << "Interception point (x = " << onePoint2[0] << ", y = " <<
			                         onePoint2[1] << ", z = " << onePoint2[2] << " )" << endl;
			                cout << "G1(izc=" << izc+1 << ",ixc=" << ixc+1 <<
			                       ",iyc=" << iyc+1 << ") = " << secLength << endl;
			 */

		}
	}
	// record the length for last point if it is in the boundary
	vector<double> endPoint1(newRayPoints[newRayPoints.size() - 1]);
	if (IsInBound(invxmin, invxmax, invymin, invymax,
	              invzmin, invzmax, endPoint1[0], endPoint1[1], endPoint1[2])) {
		vector<double> nearEndPoint1(newRayPoints[newRayPoints.size() - 2]);
		vector<double> onePoint2(interceptPoints[interceptPoints.size() - 1]);
		secLength += OneLength(nearEndPoint1, endPoint1);
		int ixc, iyc, izc;
		LocateCell(ixc, iyc, izc, invdx, invdy, invdz,
		           invxmin, invymin, invzmin,
		           (endPoint1[0] + onePoint2[0]) / 2.0,
		           (endPoint1[1] + onePoint2[1]) / 2.0,
		           (endPoint1[2] + onePoint2[2]) / 2.0);
		if (invny == 1) {
			G1MAT(izc, ixc) = secLength;
		} else {
			G1CUBE(izc, ixc, iyc) = secLength;
		}
		/*
		mexPrintf("\nLast interception point (x = %f, y = %f, z = %f)\n",
		                onePoint2[0],onePoint2[1],onePoint2[2]);
		mexPrintf("\nEnd point (x = %f, y = %f, z = %f)\n",
		                endPoint1[0],endPoint1[1],endPoint1[2]);
		mexPrintf("\nG1(izc = %d, ixc = %d, iyc = %d) = %f\n",
		                izc+1,ixc+1,iyc+1,secLength);
		*/
		/*
		         cout << "Last interception point (x = " << onePoint2[0] << ", y = " <<
		                onePoint2[1] << ", z = " << onePoint2[2] << " )" << endl;
		        cout << "End point (x = " << endPoint1[0] << ", y = " <<
		                endPoint1[1] << ", z = " << endPoint1[2] << " )" << endl;
		        cout << "G1(izc=" << izc+1 << ",ixc=" << ixc+1 <<
		               ",iyc=" << iyc+1 << ") = " << secLength << endl;
		*/
		interceptPoints.insert(interceptPoints.begin(), newRayPoints[istart]);
	}
}



/*
 * Function: trueIntercept
 * Usage: trueIntercept(newRayPoints,currentIndex,endPointOrNot,
 * DimensionIndex,xintercept,yintercept,zintercept);
 * ------------------------
 * Judge whether the current point is the true
 * intercept point or not
 */
bool trueIntercept(
  vector<double> prevPoint,
  vector<double> currentPoint,
  vector<double> nextPoint,
  char DimensionIndex,
  double xintercept,
  double yintercept,
  double zintercept
) {
	// current point on x plane?
	bool curOnXplane = (fabs(currentPoint[0] - xintercept) < TOL);
	// current point on y plane? (2D case always on y plane)
	bool curOnYplane = (fabs(currentPoint[1] - yintercept) < TOL);
	// current point on z plane?
	bool curOnZplane = (fabs(currentPoint[2] - zintercept) < TOL);

	// see if current point is a true interception
	if (DimensionIndex == 'a') { //2D case
		// judge the simplist cases
		if (!curOnXplane && !curOnZplane) {
			return false;// if not on a mesh, return false
		} else if (curOnXplane && curOnZplane) {
			return true;// if on a node, return true
		}
		// goes to more complex case
		if (curOnXplane) { //on x plane but not on z plane
			bool prevOnXplane = (fabs(prevPoint[0] - xintercept) < TOL);
			bool nextOnXplane = (fabs(nextPoint[0] - xintercept) < TOL);
			if (prevOnXplane && nextOnXplane) {//both prev and next on x plane
				return false;
			} else {
				return true;
			}
		} else { //on z plane but not on x plane
			bool prevOnZplane = (fabs(prevPoint[2] - zintercept) < TOL);
			bool nextOnZplane = (fabs(nextPoint[2] - zintercept) < TOL);
			if (prevOnZplane && nextOnZplane) {//both prev and next on z plane
				return false;
			} else {
				return true;
			}
		}

	} else { //3D case
		// if not on a plane, then return false
		if (!curOnXplane && !curOnYplane && !curOnZplane) return false;
		if (curOnXplane) {// if on X plane
			if (curOnYplane || curOnZplane) {
				return true; //on the edge of cube
			} else {
				bool prevOnXplane = (fabs(prevPoint[0] - xintercept) < TOL);
				bool nextOnXplane = (fabs(nextPoint[0] - xintercept) < TOL);
				if (prevOnXplane && nextOnXplane) {//both prev and next on x plane
					return false;
				} else {
					return true;
				}
			}
		} else if (curOnYplane) {// if on Y plane
			if (curOnXplane || curOnZplane) {
				return true; //on the edge of cube
			} else {
				bool prevOnYplane = (fabs(prevPoint[1] - yintercept) < TOL);
				bool nextOnYplane = (fabs(nextPoint[1] - yintercept) < TOL);
				if (prevOnYplane && nextOnYplane) {//both prev and next on y plane
					return false;
				} else {
					return true;
				}
			}
		} else {// if on Z plane
			if (curOnXplane || curOnYplane) {
				return true; //on the edge of cube
			} else {
				bool prevOnZplane = (fabs(prevPoint[2] - zintercept) < TOL);
				bool nextOnZplane = (fabs(nextPoint[2] - zintercept) < TOL);
				if (prevOnZplane && nextOnZplane) {//both prev and next on z plane
					return false;
				} else {
					return true;
				}
			}
		}
	}
}



/*
 * Function: PrintArray
 * Usage: PrintArray(arraypointer);
 * -------------------------
 * Print out the elements in array
 */
void PrintArray(int *ptr, int length) {
	for (int i = 0; i < length; i++) {
		cout << *ptr << "   " ;
		ptr++;
	}
	cout << endl;
}

void PrintArray(double *ptr, int length) {
	for (int i = 0; i < length; i++) {
		cout << *ptr << "   " ;
		ptr++;
	}
	cout << endl;
}


/*
 * Function: PrintArray2D
 * Usage: PrintArray2D(arraypointer,nrow,ncol)
 * ------------------------
 * Print out Array as a Matlab matrix
 */
void PrintArray2D(int *ptr, int nrow, int ncol) {
#define ARRAY2D(irow,icol) ptr[icol*nrow+irow]
	for (int irow = 0; irow < nrow; irow++) {
		for (int icol = 0; icol < ncol; icol++) {
			cout << ARRAY2D(irow, icol) << "   ";
		}
		cout << endl;
	}
	cout << endl;
}


void PrintArray2D(double *ptr, int nrow, int ncol) {
#define ARRAY2D(irow,icol) ptr[icol*nrow+irow]
	for (int irow = 0; irow < nrow; irow++) {
		for (int icol = 0; icol < ncol; icol++) {
			cout << ARRAY2D(irow, icol) << "   ";
		}
		cout << endl;
	}
	cout << endl;
}


/*
 * Function: PrintArray3D
 * Usage: PrintArray3D(arraypointer,nrow,ncol,nslide)
 * ------------------------
 * Print out Array as a Matlab cube
 */
void PrintArray3D(int *ptr, int nrow, int ncol, int nslide) {
#define ARRAY3D(irow,icol,islide) ptr[islide*nrow*ncol+icol*nrow+irow]
	for (int islide = 0; islide < nslide; islide++) {
		cout << "slide #" << islide + 1 << endl;
		for (int icol = 0; icol < ncol; icol++) {
			for (int irow = 0; irow < nrow; irow++) {
				cout << ARRAY3D(irow, icol, islide) << "   ";
			}
			cout << endl;
		}
		cout << endl;
	}
	cout << endl;
}


void PrintArray3D(double *ptr, int nrow, int ncol, int nslide) {
#define ARRAY3D(irow,icol,islide) ptr[islide*nrow*ncol+icol*nrow+irow]
	for (int islide = 0; islide < nslide; islide++) {
		cout << "slide #" << islide + 1 << endl;
		for (int icol = 0; icol < ncol; icol++) {
			for (int irow = 0; irow < nrow; irow++) {
				cout << ARRAY3D(irow, icol, islide) << "   ";
			}
			cout << endl;
		}
		cout << endl;
	}
	cout << endl;
}


/*
 * Function Vector2Array
 * Usage: Vector2Array(ptr,vector);
 * -------------------------
 * Copy vector to Array
 */
void Vector2Array(double *ptr, vector<double> vec) {
	int length = vec.size();
	for (int i = 0; i < length; i++) {
		ptr[i] = vec[i];
	}
}


void Vector2Array(int *ptr, vector<int> vec) {
	int length = vec.size();
	for (int i = 0; i < length; i++) {
		ptr[i] = vec[i];
	}
}


/*
 * Function Matrix2Array2D
 * Usage: Matrix2Array2D(ptr,matrix);
 * -------------------------
 * Copy matrix to 2D array
 */
void Matrix2Array2D(double *ptr, vector< vector <double> > matrix) {
	int length1 = matrix.size();
	int length2 = matrix[0].size();
#define ARRAY2Dmat(i,j) ptr[(i)*length2+j]
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			ARRAY2Dmat(i, j) = matrix[i][j];
		}
	}
}


void Matrix2Array2D(int *ptr, vector< vector <int> > matrix) {
	int length1 = matrix.size();
	int length2 = matrix[0].size();
#define ARRAY2Dmat(i,j) ptr[(i)*length2+j]
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			ARRAY2Dmat(i, j) = matrix[i][j];
		}
	}
}


/*
 * Function: PrintVector
 * Usage: PrintVector(vector);
 * -------------------------
 * Print out the vector values one by one
 */
void PrintVector(vector<double> vec) {
	int length = vec.size();
	for (int i = 0; i < length; i++) {
		cout << vec[i] << ", ";
	}
	cout << endl;
}

void PrintVector(vector<int> vec) {
	int length = vec.size();
	for (int i = 0; i < length; i++) {
		cout << vec[i] << ", ";
	}
	cout << endl;
}


/*
 * Function: PrintMatrix
 * Usage: PrintMatrix(vector<vector>);
 * -------------------------
 * Print out the matrix values line by line
 */
void PrintMatrix(vector< vector<double> > &matrix) {
	int length1 = matrix.size();
	for (int i = 0; i < length1; i++) {
		PrintVector(matrix[i]);
	}
	cout << endl;
}

void PrintMatrix(vector< vector<int> > &matrix) {
	int length1 = matrix.size();
	for (int i = 0; i < length1; i++) {
		PrintVector(matrix[i]);
	}
	cout << endl;
}


/*
 * Function: PrintCube
 * Usage: PrintCube(vector<vector<vector>>);
 * -------------------------
 * Print out the cube values slide by slide
 */
void PrintCube(vector< vector< vector<double> > > &cube) {
	int length2 = cube.size();
	for (int i = 0; i < length2; i++) {
		cout << "Slide #" << i + 1 << endl;
		PrintMatrix(cube[i]);
	}
	cout << endl;
}

void PrintCube(vector< vector< vector<int> > > &cube) {
	int length2 = cube.size();
	for (int i = 0; i < length2; i++) {
		cout << "Slide #" << i + 1 << endl;
		PrintMatrix(cube[i]);
	}
	cout << endl;
}


/*
 * Function: ConstructGrids
 * Usage: ConstructGrids(<vec>& x, <vec>& y, <vec>& z, nx,...
 * ny, nz, dx, dy, dz, xmin, ymin, zmin);
 * Outputs: x, y, z
 * --------------------------
 * Construct Model dimensions along x, y, z direction
 * vector x, y, z are the grids along x, y, z direction
 * for 2D model, ny=1 and dy=ymin=0.0
 */
void ConstructGrids(
  double *x,
  double *y,
  double *z,
  int nx,
  int ny,
  int nz,
  double dx,
  double dy,
  double dz,
  double xmin,
  double ymin,
  double zmin
) {
	for (int i = 0; i < nx; i++) {
		x[i] = xmin + i * dx;
	}
	for (int i = 0; i < ny; i++) {
		y[i] = ymin + i * dy;
	}
	for (int i = 0; i < nz; i++) {
		z[i] = zmin + i * dz;
	}
}


/*
 * Function: ExtractGridInfor
 * Usage: ExtractGridInfor(nx,ny,nz,dx,dy,dz,
 * xmin,ymin,zmin,xmax,ymax,zmax,x,y,z);
 * Outputs: nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,
 * xmax,ymax,zmax
 * --------------------------
 * Extract the grid information from
 * x, y, z vectors
 */
void ExtractGridInfor(
  double &dx,
  double &dy,
  double &dz,
  double &xmin,
  double &ymin,
  double &zmin,
  double &xmax,
  double &ymax,
  double &zmax,
  double *x,
  double *y,
  double *z,
  int nx,
  int ny,
  int nz
) {
	xmin = x[0], ymin = y[0], zmin = z[0];
	xmax = x[nx - 1], ymax = y[ny - 1], zmax = z[nz - 1];
	dx = x[1] - x[0], dz = z[1] - z[0], dy = dx;
	if (ny > 1) {
		dy = y[1] - y[0];
	}
}


/*
 * Function: TimeTableConstVel
 * Usage: TimeTableConstVel(constVelValue, <vec> x,...
 * <vec> y, <vec> z, xsrc, ysrc, zsrc);
 * Output: timeTable[iy][ix][iz]
 * --------------------------
 * Compute the analytical traveltime from source point
 * to every point on the velocity grid and return a
 * timeTable[iy][ix][iz]
 */
void TimeTableConstVel(
  double *timeTable,
  double vel,
  double *x,
  double *y,
  double *z,
  int nx,
  int ny,
  int nz,
  double xsrc,
  double ysrc,
  double zsrc
) {

#define TIMETABLE3D(iz,ix,iy) timeTable[iy*nz*nx+ix*nz+iz] // need to be modified

	for (int iy = 0; iy < ny; iy++) {
		for (int ix = 0; ix < nx; ix++) {
			for (int iz = 0; iz < nz; iz++) {
				TIMETABLE3D(iz, ix, iy) = sqrt(pow(x[ix] - xsrc, 2.0) +
				                               pow(y[iy] - ysrc, 2.0) +
				                               pow(z[iz] - zsrc, 2.0)) / vel;
			}
		}
	}

}


/*
 * Function: TakeInStartPoint
 * Usage: TakeInStartPoint(xstart,ystart,zstart,xmin,
 * xmax,ymin,ymax,zmin,zmax,DimensionIndex);
 * Output: xstart, ystart, zstart
 * --------------------------
 * Ask for start point from user
 */
void TakeInStartPoint(
  double &xstart,
  double &ystart,
  double &zstart,
  double xmin,
  double xmax,
  double ymin,
  double ymax,
  double zmin,
  double zmax,
  char DimensionIndex
) {
	cout << "Input the coordinates of starting point inside boundary of:" << endl;
	if (DimensionIndex == 'b') {
		cout << "x: [" << xmin << "~" << xmax << "], " <<
		     "z: [" << zmin << "~" << zmax << "], " <<
		     "y: [" << ymin << "~" << ymax << "]. " << endl;
	} else {
		cout << "x: [" << xmin << "~" << xmax << "], " <<
		     "z: [" << zmin << "~" << zmax << "]. " << endl;
	}

	while (true) {
		cout << "x = "; cin >> xstart; cout << endl;
		if (xstart >= xmin && xstart <= xmax) break;
		cout << "Inputted xstart exceeds grid boundary [" <<
		     xmin << "~" << xmax << "]" << endl;
		cout << "Please input again." << endl;
	}
	while (true) {
		cout << "z = "; cin >> zstart; cout << endl;
		if (zstart >= zmin && zstart <= zmax) break;
		cout << "Inputted zstart exceeds grid boundary[" <<
		     zmin << "~" << zmax << "]" << endl;
		cout << "Please input again." << endl;
	}
	if (DimensionIndex == 'b') {
		while (true) {
			cout << "y = "; cin >> ystart; cout << endl;
			if (ystart >= ymin && ystart <= ymax) break;
			cout << "Inputted ystart exceeds grid boundary[" <<
			     ymin << "~" << ymax << "]" << endl;
			cout << "Please input again." << endl;
		}
	} else {
		ystart = 0.0;
	}
	cout << "The inputted starting point is (" <<
	     xstart << "," << ystart << "," << zstart << ")" << endl;
	cout << endl;
}



/*
 * Function: LocateNode
 * Usage: LocateNode(ixn,iyn,izn,dx,dy,dz,xmin,ymin,zmin,xp,yp,zp);
 * Outputs: ixn, iyn, izn: indexes of closest node
 * --------------------------
 * Returns the index of closest node
 */
void LocateNode(
  int &ixn,
  int &iyn,
  int &izn,
  double dx,
  double dy,
  double dz,
  double xmin,
  double ymin,
  double zmin,
  double xp,
  double yp,
  double zp
) {
	ixn = int(round((xp - xmin) / dx));
	izn = int(round((zp - zmin) / dz));
	iyn = int(round((yp - ymin) / dy));
}


/*
 * Function: SourceBox
 * Usage: SourceBox(nx,ny,nz,xsrc,ysrc,zsrc);
 * Outputs: x,y,z indexes of sourcebox
 * --------------------------
 * Returns the sourcebox edges
 */
vector<int> SourceBox(
  int nx,
  int ny,
  int nz,
  int ixsrcn,
  int iysrcn,
  int izsrcn
) {
	vector<int> sourceBoxIndexs(6, 0);
	sourceBoxIndexs[0] = max(0, ixsrcn - BOXWIDTH);
	sourceBoxIndexs[1] = min(nx - 1, ixsrcn + BOXWIDTH);
	sourceBoxIndexs[2] = max(0, iysrcn - BOXWIDTH);
	sourceBoxIndexs[3] = min(ny - 1, iysrcn + BOXWIDTH);
	sourceBoxIndexs[4] = max(0, izsrcn - BOXWIDTH);
	sourceBoxIndexs[5] = min(nz - 1, izsrcn + BOXWIDTH);
	return sourceBoxIndexs;
}


/*
 * Function: IsInBound
 * Usage: IsInBound(xmin,xmax,ymin,ymax,zmin,zmax,xp,yp,zp);
 * Output: true if in the boundary
 * --------------------------
 * Judge if current point is in the boundary
 */
bool IsInBound(
  double xmin,
  double xmax,
  double ymin,
  double ymax,
  double zmin,
  double zmax,
  double xp,
  double yp,
  double zp
) {
	return (xp >= xmin && xp <= xmax && yp >= ymin && yp <= ymax && zp >= zmin && zp <= zmax);
}


/*
 * Function: TouchSourceBox
 * Usage: TouchSourceBox(<vec> sourceBoxIndexs,
 * <vec> x, <vec> y, <vec> z, xp, yp, zp);
 * Output: true if the current point has touched
 * the source box boundary
 * --------------------------
 * Judge if current point is on the source
 * box boundary
 */
bool TouchSourceBox(
  vector<int> sourceBoxIndexs,
  double *x,
  double *y,
  double *z,
  double xp,
  double yp,
  double zp
) {
	int ixmax = sourceBoxIndexs[0];
	int ix2 = sourceBoxIndexs[1];
	int iymax = sourceBoxIndexs[2];
	int iy2 = sourceBoxIndexs[3];
	int izmax = sourceBoxIndexs[4];
	int iz2 = sourceBoxIndexs[5];

	// cout << "ixmax = " << x[ixmax] << endl;
	// cout << "ix2 = " << x[ix2] << endl;
	// cout << "iymax = " << y[iymax] << endl;
	// cout << " iy2 = " << y[iy2] << endl;
	// cout << "izmax = " << z[izmax] << endl;
	// cout << "iz2 = " << z[iz2] << endl;

	// cout << "xp = " << xp << endl;
	// cout << "yp = " << yp << endl;
	// cout << "zp = " << zp << endl;



	return xp >= x[ixmax] && xp <= x[ix2] &&
	       yp >= y[iymax] && yp <= y[iy2] &&
	       zp >= z[izmax] && zp <= z[iz2];
}


/*
 * Function: LocateCell
 * Usage: LocateCell(ixc,iyc,izc,dx,dy,dz,xmin,ymin,zmin,xp,yp,zp);
 * Outputs: ixc, iyc, izc: indexes of closest node
 * --------------------------
 * Returns the index of current cell
 */
void LocateCell(
  int &ixc,
  int &iyc,
  int &izc,
  double dx,
  double dy,
  double dz,
  double xmin,
  double ymin,
  double zmin,
  double xp,
  double yp,
  double zp
) {
	ixc = int(floor((xp - xmin) / dx));
	izc = int(floor((zp - zmin) / dz));
	iyc = int(floor((yp - ymin) / dy));
}


/*
 * Function: InterpTableValue
 * Usage: InterpTableValue(TableValues[iy][ix][iz],
 * xp,yp,zp,dx,dy,dz,xmin,ymin,zmin);
 * Output: interped value on (xp,yp,zp)
 * --------------------------
 * Gives interped value on (xp,yp,zp)
 */
double InterpTableValue(
  double *TableValues,
  double xp,
  double yp,
  double zp,
  double dx,
  double dy,
  double dz,
  double xmin,
  double ymin,
  double zmin,
  int nx,
  int ny,
  int nz
) {
#define TABLE3D(iz,ix,iy) TableValues[iy*nz*nx+ix*nz+iz] // original
	int ixc, iyc, izc, ixc2, iyc2, izc2;
	LocateCell(ixc, iyc, izc, dx, dy, dz, xmin, ymin, zmin, xp, yp, zp);
	double t1, t2, t3, t4, t5, t6, t7, t8;
	double xfrac = xp - ixc * dx - xmin;
	double yfrac = yp - iyc * dy - ymin;
	double zfrac = zp - izc * dz - zmin;
	double yzfrac = yfrac * zfrac;
	double tleft, tright, interpValue;
	// Keep inside the boundary
	if (fabs(xfrac) < TOL) {
		ixc2 = ixc; //on x grid edge
	} else {
		ixc2 = ixc + 1;
	}
	if (fabs(yfrac) < TOL) {
		iyc2 = iyc; //on y grid edge
	} else {
		iyc2 = iyc + 1;
	}
	if (fabs(zfrac) < TOL) {
		izc2 = izc; //on z grid edge
	} else {
		izc2 = izc + 1;
	}

	t1 = TABLE3D(izc, ixc, iyc); //[iyc][ixc][izc];
	t2 = TABLE3D(izc, ixc2, iyc); //[iyc][ixc2][izc];
	t5 = TABLE3D(izc2, ixc, iyc); //[iyc][ixc][izc2];
	t6 = TABLE3D(izc2, ixc2, iyc); //[iyc][ixc2][izc2];
	t3 = TABLE3D(izc, ixc, iyc2); //[iyc2][ixc][izc];
	t4 = TABLE3D(izc, ixc2, iyc2); //[iyc2][ixc2][izc];
	t7 = TABLE3D(izc2, ixc, iyc2); //[iyc2][ixc][izc2];
	t8 = TABLE3D(izc2, ixc2, iyc2); //[iyc2][ixc2][izc2];
	tleft = ((t3 - t1) * yfrac * dz + (t5 - t1) * zfrac * dy +
	         (t1 - t3 + t7 - t5) * yzfrac) / (dy * dz) + t1;
	tright = ((t4 - t2) * yfrac * dz + (t6 - t2) * zfrac * dy +
	          (t2 - t4 + t8 - t6) * yzfrac) / (dy * dz) + t2;
	interpValue = (dx - xfrac) / dx * tleft + xfrac / dx * tright;
	return interpValue;
}


/*
 * Function: ComputeGradient
 * Usage: ComputeGradient(gx,gy,gz,TabelValues[iy][ix][iz],
 * xp,yp,zp,xmin,xmax,ymin,ymax,zmin,zmax);
 * Output: gx, gy, gz: partial derivative along x, y, z
 * --------------------------
 * Gives partial derivative along each direction
 */
void ComputeGradient(
  double &gx,
  double &gy,
  double &gz,
  double *timeTable,
  double xp,
  double yp,
  double zp,
  double dx,
  double dy,
  double dz,
  double xmin,
  double xmax,
  double ymin,
  double ymax,
  double zmin,
  double zmax,
  int nx,
  int ny,
  int nz
) {
	double x1 = max(xmin, xp - dx / 2.0);
	double x2 = min(xmax, xp + dx / 2.0);
	double ddx = x2 - x1;
	double Tx1 = InterpTableValue(timeTable, x1, yp, zp, dx, dy, dz, xmin, ymin, zmin, nx, ny, nz);
	double Tx2 = InterpTableValue(timeTable, x2, yp, zp, dx, dy, dz, xmin, ymin, zmin, nx, ny, nz);
	gx = (Tx2 - Tx1) / ddx;
	double z1 = max(zmin, zp - dz / 2.0);
	double z2 = min(zmax, zp + dz / 2.0);
	double ddz = z2 - z1;
	double Tz1 = InterpTableValue(timeTable, xp, yp, z1, dx, dy, dz, xmin, ymin, zmin, nx, ny, nz);
	double Tz2 = InterpTableValue(timeTable, xp, yp, z2, dx, dy, dz, xmin, ymin, zmin, nx, ny, nz);
	gz = (Tz2 - Tz1) / ddz;
	if (ny == 1) {
		gy = 0.0;
	} else {
		double y1 = max(ymin, yp - dy / 2.0);
		double y2 = min(ymax, yp + dy / 2.0);
		double ddy = y2 - y1;
		double Ty1 = InterpTableValue(timeTable, xp, y1, zp, dx, dy, dz, xmin, ymin, zmin, nx, ny, nz);
		double Ty2 = InterpTableValue(timeTable, xp, y2, zp, dx, dy, dz, xmin, ymin, zmin, nx, ny, nz);
		gy = (Ty2 - Ty1) / ddy;
	}
}


/*
 * Function: FindNextPointConstDist
 * Usage: FindNextPointConstDist(xnex,ynex,znex,
 * xp,yp,zp,gx,gy,gz,dist);
 * Output: xnex, ynex, znex: next point in ray
 * --------------------------
 * Gives the next point with constant distance - dist
 */
int FindNextPointConstDist(
  double &xnex,
  double &ynex,
  double &znex,
  double xp,
  double yp,
  double zp,
  double gx,
  double gy,
  double gz,
  double dist
) {
	double gall = sqrt(pow(gx, 2.0) + pow(gy, 2.0) + pow(gz, 2.0));
	if (gall < TOL) {
		cout << "Arrived the source point!" << endl;
		xnex = xp;
		ynex = yp;
		znex = zp;
		return 0;
	} else { //go to the opposite direction of gradient
		xnex = xp - gx / gall * dist;
		ynex = yp - gy / gall * dist;
		znex = zp - gz / gall * dist;
		return 1;
	}
}


/*
 * Function: InterceptInvGrids
 * Usage: <matrix> newRayPoints = InterceptInvGrids
 * (<matrix> rayPoints, <vec> invx, <vec> invy,
 * <vec> invz);
 * Output: <matrix> newRayPoints: including
 * interception points
 * --------------------------
 * Inserts interception points into rayPoints
 */
vector< vector<double> > InterceptInvGrids(
  vector< vector<double> > &rayPoints,
  char DimensionIndex,
  double *invx,
  double *invy,
  double *invz,
  int invnx,
  int invny,
  int invnz
) {
	vector< vector<double> > newRayPoints(rayPoints);
	double invxmin, invymin, invzmin, invxmax, invymax, invzmax, invdx, invdy, invdz;
	ExtractGridInfor(invdx, invdy, invdz, invxmin, invymin, invzmin, invxmax, invymax,
	                 invzmax, invx, invy, invz, invnx, invny, invnz);
	// x loop
	int length = newRayPoints.size();
	int i = 0;
	while (i < length - 1) {
		int j = 0;
		vector<double> onePoint1(newRayPoints[i]);
		vector<double> onePoint2(newRayPoints[i + 1]);
		double x1 = onePoint1[0];
		double x2 = onePoint2[0];
		int ix1 = LocateIndex(x1, invdx, invxmin, invxmax);
		int ix2 = LocateIndex(x2, invdx, invxmin, invxmax);
		// insert interceptions
		while (ix1 != ix2) {
			vector<double> newPoint(onePoint1);
			double z1 = onePoint1[2];
			double z2 = onePoint2[2];
			double xgrid;
			//cout << "ix1 = " << ix1 << ", ix2 = " << ix2;
			if (ix1 > ix2) {
				xgrid = invx[ix1];
				ix1 --;
			} else {
				xgrid = invx[ix1 + 1];
				ix1 ++;
			}
			//cout << ", new ix1 = " << ix1 << ", ix2 = " << ix2 << endl;
			//cout << "x1 = " << x1 << ", x2 = " << x2 << ", xgrid = " << xgrid << endl;
			newPoint[0] = xgrid;
			if (DimensionIndex == 'b') {
				double y1 = onePoint1[1];
				double y2 = onePoint2[1];
				double yinter = InterceptLines(x1, x2, xgrid, y1, y2);
				int iyinter = LocateGrid(yinter, invdy, invymin, invymax);
				if (fabs(yinter - invy[iyinter]) < TOL) {
					yinter = invy[iyinter];
				}
				newPoint[1] = yinter;
			}
			double zinter = InterceptLines(x1, x2, xgrid, z1, z2);
			int izinter = LocateGrid(zinter, invdz, invzmin, invzmax);
			if (fabs(zinter - invz[izinter]) < TOL) {
				zinter = invz[izinter];
			}
			newPoint[2] = zinter;
			if ((fabs(newPoint[0] - onePoint1[0]) < TOL &&
			     fabs(newPoint[1] - onePoint1[1]) < TOL &&
			     fabs(newPoint[2] - onePoint1[2]) < TOL) ||
			    (fabs(newPoint[0] - onePoint2[0]) < TOL &&
			     fabs(newPoint[1] - onePoint2[1]) < TOL &&
			     fabs(newPoint[2] - onePoint2[2]) < TOL)) {
				continue;
			}
			j++;
			newRayPoints.insert(newRayPoints.begin() + i + j, newPoint);
		}
		i = i + j + 1;
		length = length + j;
	}

	// z loop
	length = newRayPoints.size();
	i = 0;
	while (i < length - 1) {
		int j = 0;
		vector<double> onePoint1(newRayPoints[i]);
		vector<double> onePoint2(newRayPoints[i + 1]);
		double z1 = onePoint1[2];
		double z2 = onePoint2[2];
		int iz1 = LocateIndex(z1, invdz, invzmin, invzmax);
		int iz2 = LocateIndex(z2, invdz, invzmin, invzmax);
		// insert interceptions
		while (iz1 != iz2) {
			vector<double> newPoint(onePoint1);
			double x1 = onePoint1[0];
			double x2 = onePoint2[0];
			//cout << "iz1 = " << iz1 << ", iz2 = " << iz2;
			double zgrid;
			if (iz1 > iz2) {
				zgrid = invz[iz1];
				iz1--;
			} else {
				zgrid = invz[iz1 + 1];
				iz1++;
			}
			//cout << ", new iz1 = " << iz1 << ", iz2 = " << iz2 << endl;
			//cout << "z1 = " << z1 << ", z2 = " << z2 << ", zgrid = " << zgrid << endl;
			newPoint[2] = zgrid;
			if (DimensionIndex == 'b') {
				double y1 = onePoint1[1];
				double y2 = onePoint2[1];
				double yinter = InterceptLines(z1, z2, zgrid, y1, y2);
				int iyinter = LocateGrid(yinter, invdy, invymin, invymax);
				if (fabs(yinter - invy[iyinter]) < TOL) {
					yinter = invy[iyinter];
				}
				newPoint[1] = yinter;
			}
			double xinter = InterceptLines(z1, z2, zgrid, x1, x2);
			int ixinter = LocateGrid(xinter, invdx, invxmin, invxmax);
			if (fabs(xinter - invx[ixinter]) < TOL) {
				xinter = invx[ixinter];
			}
			newPoint[0] = xinter;
			if ((fabs(newPoint[0] - onePoint1[0]) < TOL &&
			     fabs(newPoint[1] - onePoint1[1]) < TOL &&
			     fabs(newPoint[2] - onePoint1[2]) < TOL) ||
			    (fabs(newPoint[0] - onePoint2[0]) < TOL &&
			     fabs(newPoint[1] - onePoint2[1]) < TOL &&
			     fabs(newPoint[2] - onePoint2[2]) < TOL)) {
				continue;
			}
			j++;
			newRayPoints.insert(newRayPoints.begin() + i + j, newPoint);
		}
		i = i + j + 1;
		length = length + j;
	}

	// y loop
	if (DimensionIndex == 'b') {
		length = newRayPoints.size();
		i = 0;
		while (i < length - 1) {
			int j = 0;
			vector<double> onePoint1(newRayPoints[i]);
			vector<double> onePoint2(newRayPoints[i + 1]);
			double y1 = onePoint1[1];
			double y2 = onePoint2[1];
			int iy1 = LocateIndex(y1, invdy, invymin, invymax);
			int iy2 = LocateIndex(y2, invdy, invymin, invymax);
			// insert interceptions
			while (iy1 != iy2) {
				vector<double> newPoint(onePoint1);
				double ygrid;
				if (iy1 > iy2) {
					ygrid = invy[iy1];
					iy1--;
				} else {
					ygrid = invy[iy1 + 1];
					iy1++;
				}
				newPoint[1] = ygrid;
				double x1 = onePoint1[0];
				double x2 = onePoint2[0];
				double xinter = InterceptLines(y1, y2, ygrid, x1, x2);
				int ixinter = LocateGrid(xinter, invdx, invxmin, invxmax);
				if (fabs(xinter - invx[ixinter]) < TOL) {
					xinter = invx[ixinter];
				}
				newPoint[0] = xinter;
				double z1 = onePoint1[2];
				double z2 = onePoint2[2];
				double zinter = InterceptLines(y1, y2, ygrid, z1, z2);
				int izinter = LocateGrid(zinter, invdz, invzmin, invzmax);
				if (fabs(zinter - invz[izinter]) < TOL) {
					zinter = invz[izinter];
				}
				newPoint[2] = zinter;
				if ((fabs(newPoint[0] - onePoint1[0]) < TOL &&
				     fabs(newPoint[1] - onePoint1[1]) < TOL &&
				     fabs(newPoint[2] - onePoint1[2]) < TOL) ||
				    (fabs(newPoint[0] - onePoint2[0]) < TOL &&
				     fabs(newPoint[1] - onePoint2[1]) < TOL &&
				     fabs(newPoint[2] - onePoint2[2]) < TOL)) {
					continue;
				}
				j++;
				newRayPoints.insert(newRayPoints.begin() + i + j, newPoint);
			}
			i = i + j + 1;
			length = length + j;
		}
	}
	return newRayPoints;
}




/*
 * Function: InterceptLines
 * Usage: yintercept = InterceptLines(x1,x2,xgrid,y1,y2);
 * Output: intercept location along y (assuming
 * x direction is gridded)
 * --------------------------
 * Compute the Intercept
 */
double InterceptLines(
  double x1,
  double x2,
  double xgrid,
  double y1,
  double y2
) {
	if (abs(x1 - x2) < TOL) {
		return (y1 + y2) / 2.0; // actually we need to avoid this
	} else {
		return ((x1 - xgrid) * y2 - (x2 - xgrid) * y1) / (x1 - x2);
	}
}


/*
 * Function: LocateIndex
 * Usage: ixp = LocateIndex(xp,dx,xmin,xmax);
 * Output: ixp: index of cell xp in <vec> x
 * --------------------------
 * Find out the ixp that x[ixp] <= xp < x[ixp+1]
 */
int LocateIndex(
  double xp,
  double dx,
  double xmin,
  double xmax
) {
	if (xp > xmax) {
		xp = xmax;
	} else if (xp < xmin) {
		xp = xmin;
	}
	return int(floor((xp - xmin) / dx));
}


/*
 * Function: LocateGrid
 * Usage: ixp = LocateGrid(xp,dx,xmin,xmax);
 * Output: ixp: index of closest xp in <vec> x
 * --------------------------
 * Find out the ixp that x[ixp] <= xp < x[ixp+1]
 */
int LocateGrid(
  double xp,
  double dx,
  double xmin,
  double xmax
) {
	if (xp < xmin) {
		xp = xmin;
	} else if (xp > xmax) {
		xp = xmax;
	}
	return int(round((xp - xmin) / dx));
}


/*
 * Function: OneLength
 * Usage: sectLength = OneLength(<vec> onePoint1,
 * <vec> onePoint2);
 * Output: sectLength: dist(Point1, Point2)
 * --------------------------
 * Computes distance between two adjent points
 */
double OneLength(
  vector<double> onePoint1,
  vector<double> onePoint2
) {
	return sqrt(pow(onePoint1[0] - onePoint2[0], 2.0) +
	            pow(onePoint1[1] - onePoint2[1], 2.0) +
	            pow(onePoint1[2] - onePoint2[2], 2.0));
}


