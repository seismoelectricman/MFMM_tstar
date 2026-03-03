// ---------------
//   Dongzhuo Li
//   March, 2015
// ---------------

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>

using namespace std;

void receiverLoad(int dimension, string &receiver_file, map <int, vector<vector<double> > > &ReceiverPoints) {
	// map <int, vector<vector<double> > > &ReceiverPoints
	// The purpose of this declaration is to allow a function to receive a reference to a map where each integer key maps to a vector 
	// of vectors of doubles. This structure can be useful for various applications, such as handling multi-dimensional data points
	//  associated with specific keys (like indices or identifiers).  -wdd
	string sourceIndex, receiverIndex, receiverTrueIndex, x , y , z, picktime;
	vector <double> newReceiverPoint;

	ifstream rfile;
	rfile.open(receiver_file);
	string line;
	stringstream ss;

	if (dimension == 2) {
		while (getline(rfile, line)) {
			ss.clear();
			ss << line;
			ss >> sourceIndex >> receiverIndex >> z >> x;
			newReceiverPoint.push_back(stod(receiverIndex));
			newReceiverPoint.push_back(stod(z));
			newReceiverPoint.push_back(stod(x));
			ReceiverPoints[stoi(sourceIndex)].push_back(newReceiverPoint);
			newReceiverPoint.clear();
		}
	}

	if (dimension == 3) {
		while (getline(rfile, line)) {
			ss.clear();
			ss << line; //is used in conjunction with std::stringstream to insert the contents of the variable line into the string stream object ss.
			ss >> sourceIndex >> receiverIndex >> receiverTrueIndex >> z >> x >> y >> picktime;
			// The statement ss >> sourceIndex >> receiverIndex >> receiverTrueIndex >> z >> x >> y >> picktime; 
			// is used in C++ to extract multiple values from a std::stringstream into respective variables.
			newReceiverPoint.push_back(stod(receiverIndex));
			// The purpose of this statement is to take a string variable receiverIndex, convert it to a double, 
			// and then add this double value to the newReceiverPoint vector. -wdd
			newReceiverPoint.push_back(stod(z));
			// The purpose of this statement is to convert the string z (which is expected to represent a numeric value) into a double and 
			// then add that double to the newReceiverPoint vector. -wdd
			newReceiverPoint.push_back(stod(x));
			newReceiverPoint.push_back(stod(y));
			ReceiverPoints[stoi(sourceIndex)].push_back(newReceiverPoint);
                        // The purpose of this statement is to store a new receiver point (represented by newReceiverPoint, which is expected to be a vector of 
			// doubles) into the ReceiverPoints data structure under the key that corresponds to sourceIndex. This is often used in applications where 
			// you want to group or categorize data points based on some index. -wdd
			newReceiverPoint.clear();
		}
	}

	rfile.close();
}
