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
	
	double* ellipse = new double[5];
	conePlaneIntersection(coneOrg, coneDir, coneAngle, planeOrg, planeDir1, planeDir2, ellipse);

	cout << "Ellipse translated by (" << ellipse[0] << ", " << ellipse[1] << ")" << endl;
	cout << "Axes are: (" << ellipse[2] << ", " << ellipse[3] << endl;
	cout << "Rotation is: " << ellipse[4] << " radians." << endl;
	delete [] ellipse;
}

void testDiscreteTarget()
{
	double rewards[2] = {0.0, 8.0};
	double probs[2] = {0.3, 0.7};
	DiscreteTarget* dt = new DiscreteTarget(2, probs, rewards, RED);
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
	printf("Target Type: %c\n", dt->AsChar());
	printf("After %d runs, got %d 0.0 rewards and %d 8.0 rewards\n", num_runs, count[0], count[1]);
	delete dt;
}

void testField()
{
	double reward = 2.0;
	Field* f = new Field(10, 20, 10.0);
	int clr;
	for (int i=0; i<10; i++)
	{
		for (int j=0; j<20; j++)
		{
			clr = f->GetIndex(i, j) % 3;
			f->AddTarget(i, j, new ConstantTarget(reward, (COLOR)clr));
		}
	}
	printf(f->AsString());
	delete f;
}


void testBasicBinaryField()
{
	int N = 50;
	double nPct = .1;
	double rPct = .6;
	double rRwd = 14.0;
	double rPrb = 0.2;
	double bPct = 0.3;
	double bRwd = 8.0;
	double bPrb = 0.9;
	BasicBinaryField* bbf = new BasicBinaryField(N, nPct, rPct, rRwd, rPrb, bPct, bRwd, bPrb);
	//printf(bbf->AsString());
	printf("Red: %d, Blue: %d, Gray: %d\n", bbf->NumRed(), bbf->NumBlue(), bbf->NumNeutral());
	printf("Expected Payout: %f\n", rPct*rRwd*rPrb + bPct*bRwd*bPrb);
	double payout = 0.0;
	gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
	int trials = 1000000;
	for (int i=0; i<trials; i++)
	{
		int row = gsl_rng_uniform_int(rng, N);
		int col = gsl_rng_uniform_int(rng, N);
		payout += bbf->Sample(row, col);
	}
	double mean = payout / (double)trials;
	printf("Actual Mean Payout: %f\n", mean);
	gsl_rng_free(rng);
	delete bbf;
}
