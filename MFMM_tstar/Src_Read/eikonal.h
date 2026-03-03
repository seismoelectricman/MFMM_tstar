// ---------------
//   Dongzhuo Li
//   March, 2015
//
//   Modified from the fast marching eikonal solver
//   written by D.Kroon University of Twente (Oct 2010)
// ---------------

#ifndef EIKONAL_H
#define EIKONAL_H

// #include <iostream>
using namespace std;

//void msfm2dCpp(double *T, double *F, int izsrcn, int ixsrcn, int *dims, bool usesecond, bool usecross);
//void msfm2dCpp(double *T, double *TS, double *F, double *Q, int izsrcn, int ixsrcn, int *dims, bool usesecond, bool usecross);
void msfm2dCpp(double *T, double *TS, double *F, double *Q, int izsrcn, int ixsrcn, double dzxy, int *dims, bool usesecond, bool usecross);

//void msfm3dCpp(double *T, double *F, int izsrcn, int ixsrcn, int iysrcn, int *dims, bool usesecond, bool usecross);
//void msfm3dCpp(double *T, double *TS, double *F, double *Q, int izsrcn, int ixsrcn, int iysrcn, int *dims, bool usesecond, bool usecross);
void msfm3dCpp(double *T, double *TS, double *F, double *Q, int izsrcn, int ixsrcn, int iysrcn, double dzxy, int *dims, bool usesecond, bool usecross);

#endif

