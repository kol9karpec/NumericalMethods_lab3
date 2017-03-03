#pragma once

#include <iostream>
#include <cstring>
#include <fstream>
#include <cmath>

using namespace std;

#define DEFAULT_SIZE 5

class TSLE
{
private:
	double ** A;
	int order;
	//int * b;
public:
	TSLE();
	TSLE(TSLE&);
	TSLE(double ** A, int order);

	double * residualsVect(double * b, double * x);
	double * zeidel(double * b, double eps,char ** iterations);

	char * print();
	~TSLE();
};



