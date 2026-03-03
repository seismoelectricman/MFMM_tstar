// ---------------
//   Dongzhuo Li
//   March, 2015
// ---------------


#ifndef RAYTRACINGHEADERS_H
#define RAYTRACINGHEADERS_H


/* Constants */

#define BOXWIDTH 2
#define TOL 0.000001


//const int BOXWIDTH = 1; // source box width
//const double TOL = 0.000000001; // compare difference tolerance

/* Function prototypes */


/*
 * Function: PrintArray
 * Usage: PrintArray(arraypointer);
 * -------------------------
 * Print out the elements in array
 */
void PrintArray(int *ptr, int length);

void PrintArray(double *ptr, int length);

/*
 * Function: PrintArray2D
 * Usage: PrintArray2D(arraypointer,nrow,ncol)
 * ------------------------
 * Print out Array as a Matlab matrix
 */
void PrintArray2D(int *ptr, int nrow, int ncol);

void PrintArray2D(double *ptr, int nrow, int ncol);


/*
 * Function: PrintArray3D
 * Usage: PrintArray3D(arraypointer,nrow,ncol,nslide)
 * ------------------------
 * Print out Array as a Matlab cube
 */
void PrintArray3D(int *ptr, int nrow, int ncol, int nslide);

void PrintArray3D(double *ptr, int nrow, int ncol, int nslide);


/*
 * Function Vector2Array
 * Usage: Vector2Array(ptr,vector);
 * -------------------------
 * Copy vector to Array
 */
void Vector2Array(double *ptr, std::vector<double> vec);

void Vector2Array(int *ptr, std::vector<int> vec);


/*
 * Function Matrix2Array2D
 * Usage: Matrix2Array2D(ptr,matrix);
 * -------------------------
 * Copy matrix to 2D array
 */
void Matrix2Array2D(double *ptr, std::vector< std::vector <double> > matrix);

void Matrix2Array2D(int *ptr, std::vector< std::vector <int> > matrix);


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
  std::vector< std::vector<double> > &rayPoints,
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
);



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
  std::vector< std::vector<double> > &newRayPoints,
  std::vector< std::vector<double> > &rayPoints,
  char DimensionIndex,
  double *invx,
  double *invy,
  double *invz,
  int invnx,
  int invny,
  int invnz
);

/*
 * Function: trueIntercept
 * Usage: trueIntercept(newRayPoints,currentIndex,endPointOrNot,
 * DimensionIndex,xintercept,yintercept,zintercept);
 * ------------------------
 * Judge whether the current point is the true
 * intercept point or not
 */
bool trueIntercept(
  std::vector<double> prevPoint,
  std::vector<double> currentPoint,
  std::vector<double> nextPoint,
  char DimensionIndex,
  double xintercept,
  double yintercept,
  double zintercept
);


/*
 * Function: Printstd::vector
 * Usage: Printstd::vector(std::vector);
 * -------------------------
 * Print out the std::vector values one by one
 */
void PrintVector(std::vector<double> vec);

void PrintVector(std::vector<int> vec);


/*
 * Function: PrintMatrix
 * Usage: PrintMatrix(std::vector<std::vector>);
 * -------------------------
 * Print out the matrix values line by line
 */
void PrintMatrix(std::vector< std::vector<double> > &matrix);

void PrintMatrix(std::vector<std::vector<int> > &matrix);


/*
 * Function: PrintCube
 * Usage: PrintCube(std::vector<std::vector<std::vector>>);
 * -------------------------
 * Print out the cube values slide by slide
 */
void PrintCube(std::vector< std::vector< std::vector<double> > > &cube);

void PrintCube(std::vector< std::vector< std::vector<int> > > &cube);


/*
 * Function: ConstructGrids
 * Usage: ConstructGrids(<vec>& x, <vec>& y, <vec>& z, nx,...
 * ny, nz, dx, dy, dz, xmin, ymin, zmin);
 * Outputs: x, y, z
 * --------------------------
 * Construct Model dimensions along x, y, z direction
 * std::vector x, y, z are the grids along x, y, z direction
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
);


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
);


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
);


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
);


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
);


/*
 * Function: SourceBox
 * Usage: SourceBox(nx,ny,nz,xsrc,ysrc,zsrc);
 * Outputs: x,y,z indexes of sourcebox
 * --------------------------
 * Returns the sourcebox edges
 */
std::vector<int> SourceBox(
  int nx,
  int ny,
  int nz,
  int ixsrcn,
  int iysrcn,
  int izsrcn
);


/*
 * Function: TouchSourceBox
 * Usage: TouchSourceBox(<vec> sourceBoxIndexs,
 * <vec> x, <vec> y, <vec> z, xp, yp, zp);
 * Outputs: true if the current point has touched
 * the source box boundary
 * --------------------------
 * Judge if current point is on the source
 * box boundary
 */
bool TouchSourceBox(
  std::vector<int> sourceBoxIndexs,
  double *x,
  double *y,
  double *z,
  double xp,
  double yp,
  double zp
);


/*
 * Function: IsInBound
 * Usage: IsInBound(xmin,xmax,ymin,ymax,zmin,zmax,xp,yp,zp);
 * Outputs: true if in the boundary
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
);


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
);


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
);


/*
 * Function: ComputeGradient
 * Usage: ComputeGradient(gx,gy,gz,TabelValues[iy][ix][iz],xp,yp,zp);
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
);


/*
 * Function: FindNextPointConstDist
 * Usage: FindNextPointConstDist(xnex,ynex,znex,
 * xp,yp,zp,gx,gy,gz,dist);
 * Output: xnex, ynex, znex: next point in ray
 * --------------------------
 * Gives the next point with constant distance
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
);


/*
 * Function: InterceptInvGrids
 * Usage: <matrix> newRayPoints = InterceptInvGrids
 * (<matrix> rayPoints,<char> DimensionIndex,
 * <vec> invx, <vec> invy,
 * <vec> invz);
 * Output: <matrix> newRayPoints: including
 * interception points
 * --------------------------
 * Inserts interception points into rayPoints
 */
std::vector< std::vector<double> > InterceptInvGrids(
  std::vector< std::vector<double> > &rayPoints,
  char DimensionIndex,
  double *invx,
  double *invy,
  double *invz,
  int invnx,
  int invny,
  int invnz
);


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
);


/*
 * Function: LocateIndex
 * Usage: ixp = LocateIndex(xp,dx,xmin,xmax);
 * Output: ixp: index of xp in <vec> x
 * --------------------------
 * Find out the ixp that x[ixp] <= xp < x[ixp+1]
 */
int LocateIndex(
  double xp,
  double dx,
  double xmin,
  double xmax
);


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
);


/*
 * Function: OneLength
 * Usage: sectionLength = OneLength(<vec> onePoint1,
 * <vec> onePoint2);
 * Output: sectionLength = dist(Point1, Point2)
 * --------------------------
 * Compute distance between two adjent points
 */
double OneLength(
  std::vector<double> onePoint1,
  std::vector<double> onePoint2
);



#endif // RAYTRACINGHEADERS_H
