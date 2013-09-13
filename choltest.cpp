/* 
*************************************************************
Copyright © 2013 Igor Kohanovsky e-mail: Igor.Kohanovsky@gmail.com
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*************************************************************
*/

// Test.cpp - Test of the Cholesky decomposition updating algorithms 
#include <iostream>
#include <iomanip>

#include <memory>
#include <random>

#include "./spd/spdchol.h"
#include "./helper/helper1.h"

#include "./service/stopwatch.h"

using std::cin;
using std::cout;
using std::endl;

using namespace mns;

typedef Defs<double>::VectorT VectorT;
typedef Defs<double>::SpdMatrixT SpdMatrixT;
typedef SpdChol<double> SpdCholT;
typedef Helper1<double> Helper1T;

typedef std::vector<std::vector<double>> Matrix2D;

Status GetRandomMatrix(int n, int ix, SpdMatrixT& orig, SpdMatrixT& reduced);
void Compress(int n, int ix, SpdMatrixT& orig, SpdMatrixT& reduced);
void SetPrintParams(int width, int precision, std::ios::fmtflags fmt=std::ios::fixed);

std::ostream& operator << (std::ostream& os, const mns::Status& obj)
{
   os << static_cast<std::underlying_type<mns::Status>::type>(obj);
   return os;
}

bool Test(int n, int ix);

int main(int argc, char* argv[])
{
	int n, ix;
	cout << "Matrix dimension (2 - 10000): ";
	cin >> n;
    if (cin.fail())
    {
       cin.clear();
	   cout << n << " is Not a number";
	   cin.get();
	   exit(1);
    }
	cout << "Row/column being deleted(1 - 10000): ";
	cin >> ix;
    if (cin.fail())
    {
       cin.clear();
	   cout << ix << " is Not a number";
	   cin.get();
	   exit(1);
    }
	if ( !Test(n, ix) ) 
		cout << "Test did not pass..." << endl;
	else
		cout << "Test passed..." << endl;

	cout << endl << "Hit <Return> key to exit..." << endl;
	cin.clear();
	cin.ignore(cin.rdbuf()->in_avail());
	cin.get();
	return 0;
}

bool Test(int n, int ix)
{
	// Test of the updating Cholesky factorization of the random SPD test matrix after the symmetric column/row deletion
	// Get a Test matrix
	SetPrintParams(2, 0);
	cout << endl << n << " x " << n << " Test matrix" << endl;
	cout << "Deleting row/column: " << ix << endl;

	SetPrintParams(7, 3);

	SpdMatrixT orig(n * (n + 1) / 2);
	SpdMatrixT reduced(n * (n + 1) / 2);

	StopWatch sw;
	if ( GetRandomMatrix(n, ix - 1, orig, reduced) != Status::Success )
	{
		cout << "Cannot get a test matrix, incorrect input parameter (matrix dimension)." << endl;
		return false;
	}
	cout <<"Time of the matrices generating: "  << sw.Elapsed() << " (s)\n";

	SpdCholT cholReduced(std::move(reduced), n - 1); 
	sw.Restart();
	Status rc = cholReduced.Factorize();
	cout <<"Time of the matrix factorizing: "  << sw.Elapsed() << " (s)\n";
	if ( rc != Status::Success )
	{
		cout << "Cannot factorize the Reduced test matrix, it is not positive-definite one. rc = " << rc << endl;
		return false;
	}
	cout << "The Reduced test matrix is factorized, rc = " << rc << endl; 

	SpdCholT cholOrig(std::move(orig), n); 
	rc = cholOrig.Factorize();
	if ( rc != Status::Success )
	{
		cout << "Cannot factorize the original test matrix, it is not positive-definite one. rc = " << rc << endl;
		return false;
	}
	cout << "The Original test matrix is factorized, rc = " << rc << endl;

	// Get assessment of the test matrix condition number
	cout.unsetf(std::ios::fixed);
	SetPrintParams(7, 1, std::ios::scientific);
	cout << "The Original matrix reciprocal condition number estimate = " << cholOrig.GetRCond() << endl;

	sw.Restart();
	rc = cholOrig.UpdateDel(ix - 1);
	cout <<"Time of the matrix factor recalculating: "  << sw.Elapsed() << " (s)\n";
	if ( rc != Status::Success )
	{
		cout << "Cannot calculate the updated Cholesky factor of the Test matrix. rc = " << rc << endl;
		return false;
	}
	cout << "The updated Cholesky factor of the Original matrix is calculated. rc = " << rc << endl;

    SpdMatrixT diff(cholOrig.GetMatrix());
	SpdMatrixT fac = cholReduced.GetMatrix();

	int n1 = n - 1;
	Size_T sz = n1 * (n1 + 1) / 2;
	for ( Size_T i = 0; i < sz; ++i ) 
	{
		diff[i] -= fac[i];
	}

	Helper1T hlp;
	auto rn = hlp.GetVectorNorm2(sz, diff);
	cout << "Norm of diference of the Orginal(recalculated) and Reduced factors: " << rn << endl;

	VectorT b1(n1, 1);
	rc = cholReduced.Solve(b1);

	VectorT b2(n1, 1);
	rc = cholOrig.Solve(b2);
	for ( Size_T i = 0; i < n1; ++i ) 
	{
		b1[i] -= b2[i];
	}
	rn = hlp.GetVectorNorm2(n1, b1);
	cout << "Norm of diference of the solutions: " << rn << endl;

	return true;
}

Status GetRandomMatrix(int n, int ix, SpdMatrixT& orig, SpdMatrixT& reduced)
{
	typedef std::mt19937 RNG;   // the Mersenne Twister with a popular choice of parameters 
								//	http://stackoverflow.com/questions/7114043/random-number-generation-in-c11-how-to-generate-how-do-they-work
	std::uniform_int_distribution<int> dist(-5, 5);
	Matrix2D tmp(n);
	RNG rng;
	rng.seed();
	for( int i = 0; i < n; ++i )
	{
		tmp[i].resize(2*n);
		for( int j = 0; j < 2*n; ++j )
		{
			tmp[i][j] = dist(rng);
		}
	}

	Size_T ij = 0;
	for( int i = 0; i < n; ++i )
	{
		for( int j = 0; j <= i; ++j )
		{
			double s = 0.0;
			for( int k = 0; k < 2*n; ++k )
			{
				s += tmp[i][k] * tmp[j][k]; 
			}
			orig[ij++] = s;
		}
	}

	Compress(n, ix, orig, reduced);
	return Status::Success;
}

void Compress(int n, int ix, SpdMatrixT& orig, SpdMatrixT& reduced)
{
	Size_T ij = 0;
	for( int i = 0; i < n; ++i )
	{
		if ( i == ix )
		{
			continue;
		}
		for( int j = 0; j <= i; ++j )
		{
			if ( j == ix )
			{
				continue;
			}
			reduced[ij++] = orig[j + ((Size_T)i) * (i + 1) / 2];
		}
	}
}

void SetPrintParams(int width, int precision, std::ios::fmtflags fmt)
{
	cout.setf(fmt);
	cout.width(width);
	cout.precision(precision);
}
