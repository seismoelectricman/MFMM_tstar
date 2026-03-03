#include <iostream>  // stands for standard input-output stream for C++. -wdd
#include <fstream>  // stream class to both read and write from/to files. -wdd
#include <string> 

using namespace std; 
// The code 'using namespace std;' tells your compiler that if you find a method being called that 
// is in 'std' class, then immediately use method in that class.

void read_parafile(int &dimension, int &nx, int &ny, int &nz, double &dx, double &dy, double &dz, double &xmin, \
                   double &ymin, double &zmin, int &invnx, int &invny, int &invnz, double &invdx, double &invdy, double &invdz, \
                   double &invxmin, double &invymin, double &invzmin, string &velModel_file, string &attModel_file, string &source_file, string &receiver_file) {

	string line;
	ifstream parafile;  // use an ifstream when you only want to perform input, an ofstream for output only,
			    // and an fstream for a stream on which you want to perform both input and output. -wdd

	parafile.open("./Par/Par_file");

	while (!parafile.eof()) {
		// eof() is a member function of the file stream class that returns true if the end of the file has been reached; otherwise, it returns false.
		// !parafile.eof() checks if the file is not at the end, meaning there is still data to read. -wdd
		getline(parafile, line);
		// reads an entire line from a file stream (parafile) and stores it in the string variable line. -wdd
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			// line != "": Checks if line is not an empty string. This ensures the line contains some data.
			// line.at(0) != '#': Checks if the first character of line is not the # symbol. This is commonly used to skip 
			// comments in files where lines beginning with # are comments.
			// !isspace(line.at(0)): Checks if the first character of line is not whitespace. The function isspace returns true if the character 
			// is a whitespace character (such as a space, tab, or newline); the ! negates it, so this condition passes if the first character is not whitespace. -wdd
			dimension = stoi(line);
			// stoi, is part of the C++ Standard Library and stands for 'string to integer'.
			cout << "dimension = " << dimension << endl;
			break;
		}
	} 

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			nx = stoi(line);
			cout << "nx = " << nx << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			ny = stoi(line);
			cout << "ny = " << ny << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			nz = stoi(line);
			cout << "nz = " << nz << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			dx = stod(line);
			// in C++ is used to convert a string (line) to a double-precision floating-point number and assign it to the variable dx.
			cout << "dx = " << dx << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			dy = stod(line);
			cout << "dy = " << dy << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			dz = stod(line);
			cout << "dz = " << dz << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			xmin = stod(line);
			cout << "xmin = " << xmin << endl;
			break;
		}
	}


	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			ymin = stod(line);
			cout << "ymin = " << ymin << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			zmin = stod(line);
			cout << "zmin = " << zmin << endl;
			break;
                }
        }

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			velModel_file = line;
			// velModel_file -- ./Par/velocity3d
			cout << "velocity model =" << velModel_file << endl;
			break;
		}
	}

	while (!parafile.eof()) {
                getline(parafile, line);
                if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
                        attModel_file = line;
                        // velModel_file -- ./Par/velocity3d
                        cout << "attenuation model =" << attModel_file << endl;
                        break;
                }
        }   // -wdd

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			source_file = line;
			// source_file -- ./Par/sources
			cout << "source file = " << source_file << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			receiver_file = line;
			// receiver_file - ./Par/receivers
			cout << "receiver file = " << receiver_file << endl;
			break;
		}
	}


	parafile.close();
}
