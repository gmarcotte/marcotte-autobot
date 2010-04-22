#include <iostream>

#include <utils.hpp>
#include <environment.hpp>

using namespace std;



void testConePlaneIntersect()
{
	// Define the cone
	double coneOrg[3] = {0.762, -3.219, 6.212};
	double coneDir[3] = {-.213, .459, -.71};
	double coneAngle = 0.1467;
	
	// Define the plane
	double planeOrg[3] = {-1.271, 3.929, .789};
	double planeDir1[3] = {1.2, 2.4, -1.9};
	double planeDir2[3] = {-3.2, 1.69, 0.828};
	
	double* ellipse = new double[6];
	conePlaneIntersection(coneOrg, coneDir, coneAngle, planeOrg, planeDir1, planeDir2, ellipse);

	cout << "Ellipse translated by (" << ellipse[0] << ", " << ellipse[1] << ")" << endl;
	cout << "Covariances are: (" << ellipse[2] << ", " << ellipse[3] << ", " << ellipse[4] << ", " << ellipse[5] << ")" << endl;
	delete [] ellipse;
}

void testDiscreteTarget()
{
	double rewards[2] = {0.0, 8.0};
	double probs[2] = {0.3, 0.7};
	DiscreteTarget* dt = new DiscreteTarget(2, probs, rewards);
	int count[2] = {0, 0};
	int num_runs = 1000;
	for (int i=0; i<num_runs; i++)
	{
		double r = dt->Payout();
		//printf("Payout: %f\n", r);
		if (r == 0.0)
			count[0]++;
		else if (r == 8.0)
			count[1]++;
		else
			printf("Error: Got invalid reward\n");
	}
	printf("After %d runs, got %d 0.0 rewards and %d 8.0 rewards\n", num_runs, count[0], count[1]);
	delete dt;
}

int main(int argc, char* argv[])
{
	testDiscreteTarget();
	return 0;
}