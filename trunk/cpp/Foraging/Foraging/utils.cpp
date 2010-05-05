// Implementations of utility functions
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_rng.h>
#include <cmath>
#include <ctime>
#include <stdio.h>
#include <windows.h>

#include <utils.hpp>

using namespace std;



// Implement the random number generation as a singleton
gsl_rng* randGen = NULL;

void init_rng()
{
	randGen = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(randGen, (unsigned long)time(NULL));
}

gsl_rng* get_rng()
{
	if (randGen == NULL)
		throw Exception("Error: Random number generator has not been initialized.");

	return randGen;
}

void free_rng()
{
	gsl_rng_free(randGen);
	randGen = NULL;
}

bool rand_bool(double p)
{
	if (p < 0.0 || p > 1.0)
		throw Exception("Error: Probabilities must be in [0, 1]");

	double val = gsl_rng_uniform(get_rng());
	return (val <= p);
}

double rand_double_uniform(double min, double max)
{
	double val = gsl_rng_uniform(get_rng());
	return (min + val*(max-min));
}

/********************************************************************/



void normalize(int N, double* src, double* dest)
{
	double sum = 0.0;
	for (int i=0; i<N; i++)
		sum += src[i]*src[i];
	double abs = sqrt(sum);

	for (int i=0; i<N; i++)
		dest[i] = src[i] / abs;
}

// Compute the determinant of N x N matrix
double _determinant(int N, gsl_matrix* mat)
{
	int sign;
	double det;
	gsl_matrix* ludecomp = gsl_matrix_alloc(N, N);
	gsl_matrix_memcpy(ludecomp, mat);
	gsl_permutation* perm = gsl_permutation_alloc(N);
	gsl_linalg_LU_decomp(ludecomp, perm, &sign);
	det = gsl_linalg_LU_det(ludecomp, sign);
	gsl_matrix_free(ludecomp);
	gsl_permutation_free(perm);
	return det;
}

// Compute the inverse of N x N matrix
void _inverse(int N, gsl_matrix* mat, gsl_matrix* inv)
{
	int sign;
	gsl_matrix* ludecomp = gsl_matrix_alloc(N, N);
	gsl_matrix_memcpy(ludecomp, mat);
	gsl_permutation* perm = gsl_permutation_alloc(N);
	gsl_linalg_LU_decomp(ludecomp, perm, &sign);
	gsl_linalg_LU_invert(ludecomp, perm, inv);
	gsl_matrix_free(ludecomp);
	gsl_permutation_free(perm);
}


// Compute u'*A*v for vectors u, v and matrix A
double _vecTMatVecMultiply(int M, int N, gsl_matrix* vec1, gsl_matrix* mat, gsl_matrix* vec2)
{
	gsl_matrix* temp = gsl_matrix_alloc(1, N);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, vec1, mat, 0.0, temp);  // temp = vec1' * mat
	
	gsl_matrix* temp2 = gsl_matrix_alloc(1,1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp, vec2, 0.0, temp2); // temp2 = temp * vec2

	double ret = gsl_matrix_get(temp2, 0, 0);
	gsl_matrix_free(temp);
	gsl_matrix_free(temp2);
	return ret;
}


// Print a visual representation of a matrix
void _printMatrix(int M, int N, gsl_matrix* mat, char* name)
{
	printf("%s:\n", name);
	for (int i=0; i<M; i++)
	{
		printf("[ ");
		for (int j=0; j<N; j++)
		{
			double val = gsl_matrix_get(mat, i, j);
			printf("%.4f ", val);
		}
		printf(" ]\n");
	}
	printf("\n\n");
}


// Compute the elliptical intersection between a cone and a plane
void conePlaneIntersection(double* coneOrg, double* coneDir, double coneAngle,
							 double* planeOrg, double* planeDir1, double* planeDir2,
							 double* ellipse)
{
	// Normalize directional vectors
	double* unitConeDir = new double[3];
	double* unitPlaneDir1 = new double[3];
	double* unitPlaneDir2 = new double[3];
	normalize(3, coneDir, unitConeDir);
	normalize(3, planeDir1, unitPlaneDir1);
	normalize(3, planeDir2, unitPlaneDir2);

	// Setup matrix views for input arrays
	gsl_matrix_view mvConeOrg = gsl_matrix_view_array(coneOrg, 3, 1);
	gsl_matrix_view mvConeDir = gsl_matrix_view_array(unitConeDir, 3, 1);
	gsl_matrix_view mvPlaneOrg = gsl_matrix_view_array(planeOrg, 3, 1);
	gsl_matrix_view mvPlaneDir1 = gsl_matrix_view_array(unitPlaneDir1, 3, 1);
	gsl_matrix_view mvPlaneDir2 = gsl_matrix_view_array(unitPlaneDir2, 3, 1);


	// Declare and allocate non-temporary matrices and other memory structures
	gsl_matrix* matM = gsl_matrix_alloc(3, 3);
	gsl_matrix* matConePlaneDir = gsl_matrix_alloc(3, 1);
	gsl_matrix* matC = gsl_matrix_alloc(3, 3);
	gsl_eigen_symmv_workspace* eigWorkspace = gsl_eigen_symmv_alloc(2);
	gsl_vector* eigVals = gsl_vector_alloc(2);
	gsl_matrix* matEigVals = gsl_matrix_calloc(2, 2);
	gsl_matrix* matEigVecs = gsl_matrix_alloc(2, 2);
	gsl_matrix* matEigValsInv = gsl_matrix_alloc(2, 2);
	gsl_matrix* matT = gsl_matrix_alloc(2, 1);
	gsl_matrix* matH = gsl_matrix_calloc(3, 3);
	gsl_matrix* matC_c = gsl_matrix_alloc(3,3);
	gsl_matrix* matSigma_c = gsl_matrix_calloc(2, 2);
	gsl_matrix* matSigma = gsl_matrix_alloc(2, 2);


	// M = coneDir * coneDir' - cos(coneAngle)^2*eye(3)
	double cosAngSq = cos(coneAngle)*cos(coneAngle);
	gsl_matrix_set_identity(matM);

	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 
				   1.0, &mvConeDir.matrix, &mvConeDir.matrix, 
				   -cosAngSq, matM);
	//_printMatrix(3, 3, matM, "M");

	// conePlaneDir = planeOrg - coneOrg;
	gsl_matrix_memcpy(matConePlaneDir, &mvPlaneOrg.matrix); 
	gsl_matrix_sub(matConePlaneDir, &mvConeOrg.matrix);
	//_printMatrix(3, 1, matConePlaneDir, "conePlaneDir");
	
	// Build the conic matrix defined by:
	//c1 = planeDir1' * M * planeDir1;
	//c2 = planeDir1' * M * planeDir2;
	//c3 = planeDir2' * M * planeDir2;
	//c4 = conePlaneDir' * M * planeDir1;
	//c5 = conePlaneDir' * M * planeDir2;
	//c6 = conePlaneDir' * M * conePlaneDir;
	// and C = [c1 c2 c4; c2 c3 c5; c4 c5 c6]
	double c1 = _vecTMatVecMultiply(3, 3, &mvPlaneDir1.matrix, matM, &mvPlaneDir1.matrix);
	double c2 = _vecTMatVecMultiply(3, 3, &mvPlaneDir1.matrix, matM, &mvPlaneDir2.matrix);
	double c3 = _vecTMatVecMultiply(3, 3, &mvPlaneDir2.matrix, matM, &mvPlaneDir2.matrix);
	double c4 = _vecTMatVecMultiply(3, 3, matConePlaneDir, matM, &mvPlaneDir1.matrix);
	double c5 = _vecTMatVecMultiply(3, 3, matConePlaneDir, matM, &mvPlaneDir2.matrix);
	double c6 = _vecTMatVecMultiply(3, 3, matConePlaneDir, matM, matConePlaneDir);

	gsl_matrix_set(matC, 0, 0, c1);
	gsl_matrix_set(matC, 0, 1, c2);
	gsl_matrix_set(matC, 0, 2, c4);
	gsl_matrix_set(matC, 1, 0, c2);
	gsl_matrix_set(matC, 1, 1, c3);
	gsl_matrix_set(matC, 1, 2, c5);
	gsl_matrix_set(matC, 2, 0, c4);
	gsl_matrix_set(matC, 2, 1, c5);
	gsl_matrix_set(matC, 2, 2, c6);

	//_printMatrix(3, 3, matC, "C");

	// Useful submatrices of C
	gsl_matrix_view mvC_R = gsl_matrix_submatrix(matC, 0, 0, 2, 2);
	gsl_matrix_view mvC_t = gsl_matrix_submatrix(matC, 0, 2, 2, 1);
	gsl_matrix_view mvC_delta = gsl_matrix_submatrix(matC, 2, 2, 1, 1);

	
	// Check that this is actually an ellipse, not a parabola or hyperbola:
	// if (det(C)~=0 && det(C_R)>0 && det(C)/(c1+c3)<0)==false)
	double detC = _determinant(3, matC);
	double detC_R = _determinant(2, &mvC_R.matrix);
    if (!((detC != 0.0) && (detC_R > 0) && (detC/(c1+c3) < 0))) 
		throw Exception("Error: plane-cone intersection is not an ellipse");

	// [R,Lambda] = eig(C_R);
	gsl_matrix* matC_Rtemp = gsl_matrix_alloc(2, 2);
	gsl_matrix_memcpy(matC_Rtemp, &mvC_R.matrix);
	gsl_eigen_symmv(matC_Rtemp, eigVals, matEigVecs, eigWorkspace);
	gsl_matrix_free(matC_Rtemp);

	gsl_matrix_set(matEigVals, 0, 0, gsl_vector_get(eigVals, 0));
	gsl_matrix_set(matEigVals, 1, 1, gsl_vector_get(eigVals, 1));

	//_printMatrix(2, 2, matEigVals, "Lambda");
	//_printMatrix(2, 2, matEigVecs, "R");

	// t = -R * inv(Lambda) * R' * C_t;
	_inverse(2, matEigVals, matEigValsInv);
	//_printMatrix(2, 2, matEigValsInv, "inv(R)");
	gsl_matrix* matTemp1 = gsl_matrix_alloc(2, 2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, matEigVecs, matEigValsInv, 0.0, matTemp1);
	gsl_matrix* matTemp2 = gsl_matrix_alloc(2, 2);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, matTemp1, matEigVecs, 0.0, matTemp2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, matTemp2, &mvC_t.matrix, 0.0, matT);
	gsl_matrix_free(matTemp1);
	gsl_matrix_free(matTemp2);
	//_printMatrix(2, 1, matT, "t");

	// H = [R t; 0 0 1];
	gsl_matrix_memcpy(&gsl_matrix_submatrix(matH, 0, 0, 2, 2).matrix, matEigVecs);
	gsl_matrix_memcpy(&gsl_matrix_submatrix(matH, 0, 2, 2, 1).matrix, matT);
	gsl_matrix_set(matH, 2, 2, 1.0);
	//_printMatrix(3, 3, matH, "H");

	// C_c = H' * C * H;
	matTemp1 = gsl_matrix_alloc(3,3);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, matH, matC, 0.0, matTemp1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, matTemp1, matH, 0.0, matC_c);
	gsl_matrix_free(matTemp1);
	//_printMatrix(3, 3, matC_c, "C_c");

	// a = sqrt(-C_c(3,3)/C_c(1,1));
	// b = sqrt(-C_c(3,3)/C_c(2,2));
	double a = sqrt(-gsl_matrix_get(matC_c, 2, 2) / gsl_matrix_get(matC_c, 0, 0));
	double b = sqrt(-gsl_matrix_get(matC_c, 2, 2) / gsl_matrix_get(matC_c, 1, 1));

	// Sigma_c = [a^2 0; 0 b^2];
	gsl_matrix_set(matSigma_c, 0, 0, a*a);
	gsl_matrix_set(matSigma_c, 1, 1, b*b);
	//_printMatrix(2, 2, matSigma_c, "Sigma_c");

	// Sigma = R * Sigma_c * R';
	matTemp1 = gsl_matrix_alloc(2, 2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, matEigVecs, matSigma_c, 0.0, matTemp1);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, matTemp1, matEigVecs, 0.0, matSigma);
	gsl_matrix_free(matTemp1);
	//_printMatrix(2, 2, matSigma, "Sigma");

	double theta = atan2(gsl_matrix_get(matEigVecs, 1, 0), gsl_matrix_get(matEigVecs, 0, 0));

	// Fill in ellipse with values in matT and matSigma
	ellipse[0] = gsl_matrix_get(matT, 0, 0);
	ellipse[1] = gsl_matrix_get(matT, 1, 0);
	ellipse[2] = a;
	ellipse[3] = b;
	ellipse[4] = theta;


	// Clean up
	delete [] unitConeDir;
	delete [] unitPlaneDir1;
	delete [] unitPlaneDir2;
	gsl_matrix_free(matM);
	gsl_matrix_free(matConePlaneDir);
	gsl_matrix_free(matC);
	gsl_eigen_symmv_free(eigWorkspace);
	gsl_vector_free(eigVals);
	gsl_matrix_free(matEigVals);
	gsl_matrix_free(matEigVecs);
	gsl_matrix_free(matEigValsInv);
	gsl_matrix_free(matT);
	gsl_matrix_free(matH);
	gsl_matrix_free(matC_c);
	gsl_matrix_free(matSigma_c);
	gsl_matrix_free(matSigma);
}


// Create points to simulate ellipse using beziers
void EllipseToBezier(RECT& r, POINT* cCtlPt)
{
    // MAGICAL CONSTANT to map ellipse to beziers
    //  			2/3*(sqrt(2)-1) 
    const double EToBConst =	0.2761423749154; 

    SIZE offset;
	offset.cx = (int)((r.right-r.left) * EToBConst);
	offset.cy = (int)((r.bottom-r.top) * EToBConst);
//  Use the following line instead for mapping systems where +ve Y is upwards
//  CSize offset((int)(r.Width() * EToBConst), -(int)(r.Height() * EToBConst));

    POINT centre;
	centre.x = (r.left + r.right) / 2;
	centre.y = (r.top + r.bottom) / 2;

    cCtlPt[0].x  =                            //------------------------/
    cCtlPt[1].x  =                            //                        /
    cCtlPt[11].x =                            //        2___3___4       /
    cCtlPt[12].x = r.left;                    //     1             5    /
    cCtlPt[5].x  =                            //     |             |    /
    cCtlPt[6].x  =                            //     |             |    /
    cCtlPt[7].x  = r.right;                   //     0,12          6    /
    cCtlPt[2].x  =                            //     |             |    /
    cCtlPt[10].x = centre.x - offset.cx;      //     |             |    /
    cCtlPt[4].x  =                            //    11             7    /
    cCtlPt[8].x  = centre.x + offset.cx;      //       10___9___8       /
    cCtlPt[3].x  =                            //                        /
    cCtlPt[9].x  = centre.x;                  //------------------------*

    cCtlPt[2].y  =
    cCtlPt[3].y  =
    cCtlPt[4].y  = r.top;
    cCtlPt[8].y  =
    cCtlPt[9].y  =
    cCtlPt[10].y = r.bottom;
    cCtlPt[7].y  =
    cCtlPt[11].y = centre.y + offset.cy;
    cCtlPt[1].y  =
    cCtlPt[5].y  = centre.y - offset.cy;
    cCtlPt[0].y  =
    cCtlPt[12].y =
    cCtlPt[6].y  = centre.y;
}


void Rotate(double radians, const POINT& c, POINT* vCtlPt, UINT Cnt)
{    
    for (int i = Cnt-1; i>=0; --i)
    {
		double x = vCtlPt[i].x - c.x;
		double y = vCtlPt[i].y - c.y;
		double r = sqrt(x*x + y*y);
		double theta = atan2(y, x);
		theta += radians;
		vCtlPt[i].x = (int)(r*cos(theta) + c.x);
		vCtlPt[i].y = (int)(r*sin(theta) + c.y);
    }
}


void CrossProduct(double* v1, double* v2, double* dest)
{
	dest[0] = v1[1]*v2[2] - v1[2]*v2[1];
	dest[1] = v1[2]*v2[0] - v1[0]*v2[2];
	dest[2] = v1[0]*v2[1] - v1[1]*v2[0];
}
