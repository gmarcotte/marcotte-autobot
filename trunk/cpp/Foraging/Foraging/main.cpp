#include <iostream>

#include <utils.hpp>

using namespace std;


int main(int argc, char* argv[])
{
	// Define the cone
	double coneOrg[3] = {0.0, 0.0, 2.0};
	double coneDir[3] = {0.0, 0.0, -1.0};
	double coneAngle = 0.2;
	
	// Define the plane
	double planeOrg[3] = {0.0, 0.0, 0.0};
	double planeDir1[3] = {1.0, 0.0, 0.0};
	double planeDir2[3] = {0.0, 1.0, 0.0};
	
	double* ellipse = conePlaneIntersection(coneOrg, coneDir, coneAngle, planeOrg, planeDir1, planeDir2);

	cout << "Ellipse translated by (" << ellipse[0] << ", " << ellipse[1] << ")" << endl;
	cout << "Covariances are: (" << ellipse[2] << ", " << ellipse[3] << ", " << ellipse[5] << ")" << endl;
	delete [] ellipse;
	return 0;
}