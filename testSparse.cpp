#include<iostream>
#include<eigen3/Eigen/Sparse>
#include<vector>

using namespace std;
using namespace Eigen;

typedef Triplet<double> Trip;
typedef SparseMatrix<double> SpMat;

int main()
{
	vector<Trip> T;
	SpMat X (3,3);

	T.push_back(Trip(0,0,11.0));
	T.push_back(Trip(0,1,12.0));
	T.push_back(Trip(1,0,21.0));
	T.push_back(Trip(1,1,22.0));

	X.setFromTriplets(T.begin(),T.end());

	cout << X << endl;
}

	
