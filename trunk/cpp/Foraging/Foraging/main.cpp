#include <config.h>
#include <environment.hpp>
#include <tasks.hpp>
#include <controller.hpp>
#include <utils.hpp>

#include <iostream>

using namespace std;


int main(int argc, char* argv)
{

	BasicBinaryField* bbf;
	Forager* fgr;
	ForagingTask* ft;
	NivController* nvc;

	init_rng();
	bbf = new BasicBinaryField(fieldSide, neutralPct, redPct, redRwd, redPrb, bluePct, blueRwd, bluePrb, 1.0);
	fgr = new Forager(30.0, 30.0, 8.0, -1.0, 0.0, -0.2, foragerSpeed, visualConeAngle, projectionRadius, projection_dR, projection_dTheta);
	ft = new ForagingTask(bbf, fgr, (double)fieldSide, 8.0, 9.0, 100);
	nvc = new NivController();
	fgr->SetController(nvc);

	cout << "Starting Run...." << endl;
	ft->Run(sim_dT);
	cout << "Run Complete...." << endl;

	delete ft;
	delete bbf;
	delete fgr;
	delete nvc;
}