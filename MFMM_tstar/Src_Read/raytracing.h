// ---------------
//   Dongzhuo Li
//   March, 2015
// ---------------

#ifndef RAYTRACING_H
#define RAYTRACING_H


using namespace std;

void raytracing(int dimension, int nx, int ny, int nz, double dx, double dy, double dz, double xmin, double ymin, double zmin, \
                int invnx, int invny, int invnz, double invdx, double invdy, double invdz, double invxmin, double invymin, double invzmin, \
                double xsrc, double ysrc, double zsrc, double xstart, double ystart, double zstart, double *timeTable, \
                vector< vector<double> > &rayPoints);

void localDimension(double xsrc, double ysrc, double zsrc, double xrec, double yrec, double zrec, double dxyz_loc, int &dim_x, int &dim_y, int &dim_z);

void coordinateBack(vector< vector<double> > &rayPoints, double theta, double xsrc, double ysrc, double zsrc);

void coordinateTransform(int nx, int ny, int nz, double xmin, double ymin, double zmin, \
                         double dxyz_glb, double xsrc, double ysrc, double zsrc, double xrec, double yrec, \
                         double zrec, int dim_x, int dim_y, int dim_z, double dxyz_loc, double *velModel, \
                         double &loc_xsrc, double &loc_ysrc, double &loc_zsrc, double &loc_xrec, double &loc_yrec, \
                         double &loc_zrec, int & ixsrcn, int &iysrcn, int &izsrcn, \
                         double  &xmin_loc, double &ymin_loc, double &zmin_loc, \
                         double &theta, double *locVelModel);

#endif
