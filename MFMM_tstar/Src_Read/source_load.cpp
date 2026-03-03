// ---------------
//   Dongzhuo Li
//   March, 2015
// ---------------

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

void sourceLoad(int dimension, string &source_file, vector<double * > &SourcePoints) {
	// vector<double *> is a vector (dynamic array) of double * pointers, meaning each element in this 
	// vector is a pointer to a double value.   -wdd
	int i = 0;
	string index, x , y , z; // This means each of these variables (index, x, y, z) will be able to store a string value. -wdd

	ifstream sfile; // ifstream: This stands for input file stream; 
			// sfile: This is the variable name for the file input stream object. -wdd
	sfile.open(source_file);
	string line;
	stringstream ss; // in C++ creates an object named ss of the type std::stringstream. -wdd

	if (dimension == 2) {
		while (getline(sfile, line)) {
			ss.clear();
			ss << line;
			ss >> index >> z >> x;
			SourcePoints.push_back(new double[3]);
			SourcePoints.at(i)[0] = stod(index);
			SourcePoints.at(i)[1] = stod(z);
			SourcePoints.at(i)[2] = stod(x);
			i++;
		}
	}

	if (dimension == 3) {
		while (getline(sfile, line)) {
			ss.clear(); // ss.clear() clears all of these state flags, allowing the stream to be used again without errors. -wdd
			ss << line; // in C++ is used to insert the contents of the variable line into a std::stringstream object ss.
			ss >> index >> z >> x >> y;
			// in C++ is used to extract values from a std::stringstream object (ss) into multiple variables (index, z, x, and y). -wdd
			SourcePoints.push_back(new double[4]); // This part of the statement dynamically allocates an array of 4 double values on the heap.
			SourcePoints.at(i)[0] = stod(index); // The statement assigns the converted double value from the string 
							     // index to the 'first' element of the array at index i in the SourcePoints vector. -wdd
			SourcePoints.at(i)[1] = stod(z);  // The statement assigns the converted double value from the string z to the 'second' element 
							  // of the array at index i in the SourcePoints vector.  -wdd
			SourcePoints.at(i)[2] = stod(x); // This statement assigns the converted double value from the string x to the 'third' element 
							 // of the array at index i in the SourcePoints vector.
			SourcePoints.at(i)[3] = stod(y);
			i++;
		}
	}

	sfile.close();
}
