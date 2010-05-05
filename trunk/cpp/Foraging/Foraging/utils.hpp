// Definitions of utility functions for various common purposes


#ifndef UTILS_HPP
#define UTILS_HPP
#include <gsl/gsl_rng.h>

#include <windows.h>

const double PI = 3.1415926;

// Random number sampling
void init_rng();
gsl_rng* get_rng();
void free_rng();
bool rand_bool(double p = 0.5);
double rand_double_uniform(double min, double max);



// Custom base exception class for error handling
class Exception
{
private:
	const char* _msg;
public:
	Exception(const char* msg = "Error!") : _msg(msg) {}
	const char* what() {return this->_msg;}
};


/** Normalize a vector
Inputs:
	N: An integer, the vector length.
	src: An N-length vector of doubles to be normalized.
	dest: An N-length vector of doubles in which to store the normalized vector (can be the same as src).
*/
void normalize(int N, double* src, double* dest);


/**Compute the intersection between a cone and a plane.
Inputs:
	coneOrg: A 3-length double vector with the 3-D cartesian coordinates of the cone vertex
	coneDir: A 3-length double vector (non-zero magnitude) indicating the direction of the cone
	coneAngle: A double indicating the width of the cone (in radians)
	planeOrg: A 3-length double vector with the 3-D cartesian coordinates of an origin point for the plane
	planeDir1 and planeDir2: 3-length double vectors (non-zero magnitude) defining the plane
	Note that the normal vector of the plane is (planeDir1 x planeDir2)
	ellipse: A 5-length double vector into which the output will be written.

Output:
	A 6-length double vector.  The first two elements indicate the translation from the plane origin to
	the center of the ellipse, with respect to the plane basis defined by planeDir1 and planeDir2.
	The next two elements are the axes of the ellipse.
	The final element is the rotation of the ellipse from the plane basis.
*/
void conePlaneIntersection(double* coneOrg, double* coneDir, double coneAngle,
							 double* planeOrg, double* planeDir1, double* planeDir2,
							 double* ellipse);


// Create points to simulate ellipse using beziers
void EllipseToBezier(RECT& r, POINT* cCtlPt);


// Rotate an array of points about a given point by a given angle
void Rotate(double radians, const POINT& c, POINT* vCtlPt, UINT Cnt);


// Compute the cross product of two vectors
void CrossProduct(double* v1, double* v2, double* dest);



#endif //UTILS_HPP