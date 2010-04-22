// Definitions of utility functions for various common purposes


#ifndef UTILS_HPP
#define UTILS_HPP

// Custom base exception class for error handling
class Exception
{
private:
	const char* _msg;
public:
	Exception(const char* msg = "Error!") : _msg(msg) {}
	const char* what() {return this->_msg;}
};


/**Compute the intersection between a cone and a plane.
Inputs:
	coneOrg: A 3-length double vector with the 3-D cartesian coordinates of the cone vertex
	coneDir: A 3-length double vector (non-zero magnitude) indicating the direction of the cone
	coneAngle: A double indicating the width of the cone (in radians)
	planeOrg: A 3-length double vector with the 3-D cartesian coordinates of an origin point for the plane
	planeDir1 and planeDir2: 3-length double vectors (non-zero magnitude) defining the plane
	Note that the normal vector of the plane is (planeDir1 x planeDir2)
	ellipse: A 6-length double vector into which the output will be written.

Output:
	A 6-length double vector.  The first two elements indicate the translation from the plane origin to
	the center of the ellipse, with respect to the plane basis defined by planeDir1 and planeDir2.
	The final four elements define the rotation and axes of the ellipse as the covariance matrix of
	a 2D Gaussian distribution w.r.t the plane basis.
*/
void conePlaneIntersection(double* coneOrg, double* coneDir, double coneAngle,
							 double* planeOrg, double* planeDir1, double* planeDir2,
							 double* ellipse);



#endif //UTILS_HPP