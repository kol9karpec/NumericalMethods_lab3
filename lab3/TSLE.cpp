#include "TSLE.h"

TSLE::TSLE()
{
	order = DEFAULT_SIZE;
	A = new double*[order];
	for (int i = 0; i < order; i++)
	{
		A[i] = new double[order];
		for (int j = 0; j < order; j++)
			A[i][j] = 0;
	}
	/*b = new int[order];
	for (int i = 0; i < order; i++)
		b[i] = 0;*/
}

TSLE::TSLE(TSLE & sle)
{
	order = sle.order;
	A = new double*[order];
	for (int i = 0; i < order; i++)
	{
		A[i] = new double[order];
		for (int j = 0; j < order; j++)
			A[i][j] = sle.A[i][j];
	}
	/*b = new int[order];
	for (int i = 0; i < order; i++)
		b[i] = sle.b[i];*/
}

TSLE::TSLE(double ** A, int order)
{
	this->order = order;
	this->A = new double*[this->order];
	for (int i = 0; i < this->order; i++)
	{
		this->A[i] = new double[this->order];
		for (int j = 0; j < this->order; j++)
			this->A[i][j] = A[i][j];
	}
	/*this->b = new int[this->order];
	for (int i = 0; i < this->order; i++)
		this->b[i] = b[i];*/
}

double * TSLE::residualsVect(double * b, double * x)
{
	double * res = new double[order];
	for (int i = 0; i < order; i++)
	{
		double rightPart=0;
		for (int j = 0; j < order; j++)
			rightPart += A[i][j] * x[j];
		res[i] = b[i] - rightPart;
	}
	return res;
}

double * TSLE::zeidel(double * b, double eps, char ** iterations)
{
	double * p = new double[order];
	double * x = new double[order];

	//start approximation
	for (int i = 0; i < order; i++)
		x[i] = 0;
	
	*iterations = new char[50000];
	strcpy(*iterations, "");

	int iter = 0;

	while(true)
	{
		//saving previous x
		for (int i = 0; i < order; i++)
			p[i] = x[i];

		for (int i = 0; i < order; i++)
		{
			double var = 0;
			for (int j = 0; j < i; j++) //using founded values
				var += (A[i][j] * x[j]);

			//skipping i==j (a[i][j] == 0)

			for (int j = i + 1; j < order; j++) // using old values
				var += (A[i][j] * p[j]);

			x[i] = (b[i] - var) / A[i][i];
		}

		char * buff = new char[256];
		strcat(*iterations, "Iteration #");
		strcat(*iterations, _itoa(iter,buff,10));
		strcat(*iterations, "\nx=(");
		for (int i = 0; i < order; i++)
		{
			strcpy(buff, "");
			sprintf(buff, "%2.6f", x[i]);
			strcat(*iterations, buff);
			if (i != order - 1)
				strcat(*iterations, " , ");
		}

		strcat(*iterations, ")\neps=(");
		double * res = residualsVect(b, x);

		for (int i = 0; i < order; i++)
		{
			strcpy(buff, "");
			sprintf(buff, "%2.6f", res[i]);
			strcat(*iterations, buff);
			if (i != order - 1)
				strcat(*iterations, " , ");
		}

		strcat(*iterations, ")\nNorma eps = ");

		double norm = 0;

		for (int i = 0; i < order; i++)
			norm += pow(res[i],2);

		norm = sqrt(norm);

		strcpy(buff, "");
		sprintf(buff, "%2.6f", norm);
		strcat(*iterations, buff);
		strcat(*iterations, "\n\n");

		if (norm < eps)
			break;

		delete[] buff;

		iter++;
	} 

	delete[] p;
	p = NULL;

	return x;
}

char * TSLE::print()
{
	int size = (this->order * 6)*this->order;
	char * result = new char[size];
	strcpy(result, "");
	for (int i = 0; i < this->order; i++)
	{
		for (int j = 0; j < this->order; j++)
		{
			char * buff = new char[256];
			sprintf(buff, "%2.2f", A[i][j]);
			strcat(result,buff);
			delete(buff);
			if(j < this->order-1) strcat(result,"\t");
		}
		if (i < this->order - 1) strcat(result,"\n");
	}

	return result;
}

TSLE::~TSLE()
{
	for (int i = 0; i < order; i++)
	{
		delete[](A[i]);
		A[i] = nullptr;
	}
	delete[](A);
	A = nullptr;
	/*delete(b);
	b = nullptr;*/
}
