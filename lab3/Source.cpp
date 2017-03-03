#include <iostream>
#include <iomanip>
#include <fstream>
#include "TSLE.h"

using namespace std;

int main(int argc, const char* argv[])
{
	char * inFile = new char[20];
	char * outFile = new char[20];
	if (argc > 1) strcpy(inFile, argv[1]);
	else strcpy(inFile, "input.txt");
	if (argc > 2) strcpy(outFile, argv[2]);
	else strcpy(outFile, "output.txt");

	ifstream in(inFile);
	ofstream out(outFile);

	int order;
	in >> order;

	double ** A;
	A = new double*[order];

	for (int i = 0; i < order; i++)
	{
		A[i] = new double[order];
		for (int j = 0; j < order; j++)
			in >> A[i][j];
	}

	double * b = new double[order];

	for (int i = 0; i < order; i++)
		in >> b[i];

	double eps;
	in >> eps;

	TSLE matr(A, order);
	
	char ** iterations = new char*[1];

	double * x = matr.zeidel(b,eps,iterations);

	out << "Solution:\nx=(";
	out << setprecision(-(int)log10(eps / 10));
	for (int i = 0; i < order; i++)
	{
		out << x[i];
		if (i != order - 1)
			out << " , ";
	}
	out << ")\n---------------------------------------\n";
	out << *iterations << endl;



	delete[] inFile;
	delete[] outFile;

	for (int i = 0; i < order; i++)
		delete[] A[i];
	delete[] A;

	delete[] b;

	in.close();
	out.close();

	return 0;
}