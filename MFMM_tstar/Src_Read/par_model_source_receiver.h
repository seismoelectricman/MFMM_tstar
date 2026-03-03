// ---------------
//   Dongzhuo Li
//   March, 2015
// ---------------

#ifndef PAR_MODEL_SOURCE_RECEIVER_H
#define PAR_MODEL_SOURCE_RECEIVER_H

using namespace std;

void read_parafile(int &dimension, int &nx, int &ny, int &nz, double &dx, double &dy, double &dz, double &xmin, double &ymin, double &zmin, int &invnx, int &invny, int &invnz, double &invdx, double &invdy, double &invdz, double &invxmin, double &invymin, double &invzmin, string &velModel_file, string &attModel_file, string &source_file, string &receiver_file);

void velModelLoad(int Dimension, string &velModel_file, double *velModel);

void attModelLoad(int Dimension, string &attModel_file, double *attModel);

void sourceLoad(int dimension, string &source_file, vector <double *> &SourcePoints);

void receiverLoad(int dimension, string &receiver_file, map < int, vector< vector <double> > > &ReceiverPoints);

#endif
