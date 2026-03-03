// ---------------
//   Dongzhuo Li
//   March, 2015
// ---------------

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>

using namespace std; // is a directive in C++ that tells the compiler to use the standard namespace (std). Here's a detailed explanation. -wdd

void velModelLoad(int dimension, string &velModel_file, double *velModel) {
	int i = 0;
	string a, b , c , d;

	ifstream vfile;
	vfile.open(velModel_file);
	string line,dimxyz;
	stringstream ss;

	if (dimension == 2) {
                getline(vfile, dimxyz);
		while (getline(vfile, line)) {
			ss.clear();
			ss << line;
			ss >> a >> b >> c;
			velModel[i] = stod(c);
			i++;
		}
	}

	if (dimension == 3) {
		while (getline(vfile, line)) {
			ss.clear();
			ss << line; // is used to insert the contents of the string variable line into a std::stringstream object ss.
			ss >> d; // in C++ is used to extract a value from a std::stringstream object ss and store it in a variable d.
			velModel[i] = stod(d); // The purpose of this statement is to take a string representation of a numeric value (d),
					       // convert it to a double, and store that double in the velModel at the specified index i.
			i++;
		}
	}

	vfile.close();
}
